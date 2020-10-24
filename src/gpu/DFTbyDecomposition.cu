/*
 * DFTbyDecomposition.cpp
 *
 *  Created on: Oct 21, 2020
 *      Author: himesb
 */

#include "gpu_core_headers.h"
#include "/groups/himesb/cufftdx/include/cufftdx.hpp"
using namespace cufftdx;

__global__ void DFT_R2C_WithPaddingKernel(cufftReal* input_values, cufftComplex* output_values, int4 dims_in, int4 dims_out, float C);
__global__ void DFT_C2C_WithPaddingKernel_strided(cufftComplex* input_values, int4 dims_in, int4 dims_out, float C);
__global__ void DFT_R2C_WithPaddingKernel_strided(cufftReal* input_values, cufftComplex* output_values, int4 dims_in, int4 dims_out, float C);
__global__ void DFT_C2C_WithPaddingKernel(cufftComplex* input_values, int4 dims_in, int4 dims_out, float C);
__global__ void DFT_C2C_WithPaddingKernel_rdx2(cufftComplex* input_values, int4 dims_in, int4 dims_out, float C);

template<class FFT>
__global__ void block_fft_kernel_R2C_strided(typename FFT::input_type* input_values, typename FFT::output_type* output_values, int4 dims_in, int4 dims_out, float twid_constant, int n_sectors);

using FFT_256          = decltype(Block() + Size<256>() + Type<fft_type::c2c>() +
                     Precision<float>() + ElementsPerThread<2>() + FFTsPerBlock<1>() + SM<700>());
using FFT_16          = decltype(Block() + Size<16>() + Type<fft_type::c2c>() +
                     Precision<float>() + ElementsPerThread<2>() + FFTsPerBlock<1>() + SM<700>());

DFTbyDecomposition::DFTbyDecomposition() // @suppress("Class members should be properly initialized")
{
	is_set_gpu_images = false;
	is_set_twiddles = false;
//	is_set_outputs = false;
}

DFTbyDecomposition::~DFTbyDecomposition()
{
	if (is_set_twiddles)
	{
		cudaErr(cudaFree(twiddles));
	}
//	if (is_set_outputs)
//	{
//		cudaErr(cudaFree(output_real));
//		cudaErr(cudaFree(output_imag));
//	}
}

DFTbyDecomposition::DFTbyDecomposition(const DFTbyDecomposition &other)
{
	// TODO Auto-generated constructor stub

}

DFTbyDecomposition& DFTbyDecomposition::operator=(
		const DFTbyDecomposition &other) {
	// TODO Auto-generated method stub

}

void DFTbyDecomposition::InitTestCase(int wanted_input_size_x, int wanted_input_size_y, int wanted_output_size_x, int wanted_output_size_y)
{
	dims_input = make_int2(wanted_input_size_x, wanted_input_size_y);
	dims_output = make_int2(wanted_output_size_x, wanted_output_size_y);

	// In practice we'll give a pointer to the arrays in some GpuImages
}

void DFTbyDecomposition::SetGpuImages(Image& cpu_input, Image& cpu_output)
{

	// Should be in real space, TODO add check
	input_image.CopyFromCpuImage(cpu_input);
	input_image.CopyHostToDevice();


	// Initialize to Fourier space
	output_image.CopyFromCpuImage(cpu_output);
	output_image.Allocate((int)dims_output.x, (int)dims_output.y, 1, false);
	output_image.Zeros();



	is_set_gpu_images = true;

}


void DFTbyDecomposition::DFT_R2C_WithPadding()
{

	// FIXME when adding real space complex images
	MyAssertTrue( input_image.is_in_memory_gpu, "Input image is in not on the GPU!");
	MyAssertTrue( output_image.is_in_memory_gpu, "Output image is in not on the GPU!");



	int threadsPerBlock = input_image.dims.x; // FIXME make sure its a multiple of 32
	int gridDims = input_image.dims.y;
//	dim3 gridDims = dim3((output_image.dims.w/2 + threadsPerBlock - 1) / threadsPerBlock,
//					  	1, 1);
//  output_image.dims.y
	int shared_mem = sizeof(float)*input_image.dims.x;
	float C = -2*PIf/output_image.dims.x;
	DFT_R2C_WithPaddingKernel<< <gridDims, threadsPerBlock, shared_mem, cudaStreamPerThread>> > ( input_image.real_values_gpu,  output_image.complex_values_gpu, input_image.dims, output_image.dims, C);
	cudaStreamSynchronize(cudaStreamPerThread);



}

__global__ void DFT_R2C_WithPaddingKernel(cufftReal* input_values, cufftComplex* output_values, int4 dims_in, int4 dims_out, float C)
{

//	// Initialize the shared memory, assuming everying matches the input data X size in
	extern __shared__ float s[];
	// Avoid N*k type conversion and multiplication
	float* data = s;
//	float* coeff= (float*)&data[dims_in.x];


	int x = threadIdx.x;
	int pixel_out = (dims_out.w/2)*blockIdx.x;


	data[x] = __ldg((const float *)&input_values[dims_in.w*blockIdx.x + x]);
	__syncthreads();
//
//	 Loop over N updating the actual twiddle value along the way. This might lead to accuracy problems.
	float sum_real;
	float sum_imag;
	float twi_r;
	float twi_i;
	float coeff;

	for (int k = threadIdx.x; k < dims_out.w/2; k+=blockDim.x)
	{
		coeff = C*(float)k;
		sum_real = 0.0f;
		sum_imag = 0.0f;
		for (int n = 0; n < dims_in.x; n++)
		{
			__sincosf(coeff*n,&twi_i,&twi_r);
			sum_real = __fmaf_rn(data[n],twi_r,sum_real);
			sum_imag = __fmaf_rn(data[n],twi_i,sum_imag);
		}

		// Not sure if an async write, or storage to a shared mem temp would be faster.
		output_values[pixel_out + k].x = sum_real;
		output_values[pixel_out + k].y = sum_imag;
	}


	return;

}


void DFTbyDecomposition::DFT_C2C_WithPadding()
{

	// FIXME when adding real space complex images
	MyAssertTrue( input_image.is_in_memory_gpu, "Input image is in not on the GPU!");
	MyAssertTrue( output_image.is_in_memory_gpu, "Output image is in not on the GPU!");


	int threadsPerBlock = input_image.dims.x; // FIXME make sure its a multiple of 32
	int gridDims = output_image.dims.w/2;

	int shared_mem = sizeof(cufftComplex)*input_image.dims.x;

	float C = -2*PIf/output_image.dims.x;
	DFT_C2C_WithPaddingKernel<< <gridDims, threadsPerBlock, shared_mem, cudaStreamPerThread>> > ( output_image.complex_values_gpu, input_image.dims, output_image.dims, C);
	cudaStreamSynchronize(cudaStreamPerThread);



}

__global__ void DFT_C2C_WithPaddingKernel(cufftComplex* inplace_image, int4 dims_in, int4 dims_out, float C)
{

	// Initialize the shared memory, assuming everying matches the input data X size in
	// Check that setting cudaFuncSetSharedMemConfig  to 8byte makes any diff for complex reads
	extern __shared__ cufftComplex c[];
	cufftComplex* data = c;


	int x = threadIdx.x;
	int pixel_out = (dims_out.w/2)*blockIdx.x;

	data[x] = __ldg((const cufftComplex *)&inplace_image[pixel_out + x]);
	__syncthreads();
//
//	 Loop over N updating the actual twiddle value along the way. This might lead to accuracy problems.
	cufftComplex sum;
	float twi_r;
	float twi_i;
	float coeff;
	float tmp;

	for (int k = threadIdx.x; k < dims_out.w/2; k+=blockDim.x)
	{
		coeff = C*(float)k;
		sum.x = 0.0f;
		sum.y = 0.0f;
		for (int n = 0; n < dims_in.y; n++)
		{
			__sincosf(coeff*n,&twi_i,&twi_r);
			tmp = data[n].x * twi_i;
			sum.x += __fmaf_rn(data[n].x, twi_r, -twi_i * data[n].y);
			sum.y += __fmaf_rn(data[n].y, twi_r, tmp);
		}

		// Not sure if an async write, or storage to a shared mem temp would be faster.
//		inplace_image[pixel_out + k].x = sum_real;
//		inplace_image[pixel_out + k].y = sum_imag;
		inplace_image[pixel_out + k] = sum;
	}



	return;

}


void DFTbyDecomposition::DFT_R2C_WithPadding_strided()
{

	// FIXME when adding real space complex images
	MyAssertTrue( input_image.is_in_memory_gpu, "Input image is in not on the GPU!");
	MyAssertTrue( output_image.is_in_memory_gpu, "Output image is in not on the GPU!");



	int threadsPerBlock = input_image.dims.y; // FIXME make sure its a multiple of 32
	int gridDims = input_image.dims.x;
//	dim3 gridDims = dim3((output_image.dims.w/2 + threadsPerBlock - 1) / threadsPerBlock,
//					  	1, 1);
//  output_image.dims.y
	int shared_mem = sizeof(float)*input_image.dims.y;
	float C = -2*PIf/output_image.dims.y;
	DFT_R2C_WithPaddingKernel_strided<< <gridDims, threadsPerBlock, shared_mem, cudaStreamPerThread>> > ( input_image.real_values_gpu,  output_image.complex_values_gpu, input_image.dims, output_image.dims, C);
	cudaStreamSynchronize(cudaStreamPerThread);



}

__global__ void DFT_R2C_WithPaddingKernel_strided(cufftReal* input_values, cufftComplex* output_values, int4 dims_in, int4 dims_out, float C)
{

//	// Initialize the shared memory, assuming everying matches the input data X size in
	extern __shared__ float s[];
	// Avoid N*k type conversion and multiplication
	float* data = s;
//	float* coeff= (float*)&data[dims_in.x];


	int y = threadIdx.x;
	int pixel_in = blockIdx.x + y * (dims_in.w);

	data[y] = __ldg((const cufftReal *)&input_values[pixel_in]);
	__syncthreads();
//

//
//	 Loop over N updating the actual twiddle value along the way. This might lead to accuracy problems.
	float sum_real;
	float sum_imag;
	float twi_r;
	float twi_i;
	float coeff;

	for (int k = threadIdx.x; k < dims_out.y; k+=blockDim.x)
	{
		coeff = C*(float)k;
		sum_real = 0.0f;
		sum_imag = 0.0f;
		for (int n = 0; n < dims_in.x; n++)
		{
			__sincosf(coeff*n,&twi_i,&twi_r);
			sum_real = __fmaf_rn(data[n],twi_r,sum_real);
			sum_imag = __fmaf_rn(data[n],twi_i,sum_imag);
		}

		// Not sure if an async write, or storage to a shared mem temp would be faster.
		// Not sure if an async write, or storage to a shared mem temp would be faster.
		output_values[blockIdx.x + k * (dims_out.w/2)].x = sum_real;
		output_values[blockIdx.x + k * (dims_out.w/2)].y = sum_imag;
	}


	return;

}


void DFTbyDecomposition::DFT_C2C_WithPadding_strided()
{

	// FIXME when adding real space complex images
	MyAssertTrue( input_image.is_in_memory_gpu, "Input image is in not on the GPU!");
	MyAssertTrue( output_image.is_in_memory_gpu, "Output image is in not on the GPU!");


	int threadsPerBlock = input_image.dims.y; // FIXME make sure its a multiple of 32
	int gridDims = output_image.dims.w/2;

	int shared_mem = sizeof(cufftComplex)*input_image.dims.y;

	float C = -2*PIf/output_image.dims.y;
	DFT_C2C_WithPaddingKernel_strided<< <gridDims, threadsPerBlock, shared_mem, cudaStreamPerThread>> > ( output_image.complex_values_gpu, input_image.dims, output_image.dims, C);
	cudaStreamSynchronize(cudaStreamPerThread);



}

__global__ void DFT_C2C_WithPaddingKernel_strided(cufftComplex* inplace_image, int4 dims_in, int4 dims_out, float C)
{

	// Initialize the shared memory, assuming everying matches the input data X size in
	// Check that setting cudaFuncSetSharedMemConfig  to 8byte makes any diff for complex reads
	extern __shared__ cufftComplex c[];
	cufftComplex* data = c;


	int y = threadIdx.x;
	int pixel_in = blockIdx.x + y * (dims_out.w/2);


	data[y] = __ldg((const cufftComplex *)&inplace_image[pixel_in]);
	__syncthreads();
//
//	 Loop over N updating the actual twiddle value along the way. This might lead to accuracy problems.
	float sum_real;
	float sum_imag;
	float twi_r;
	float twi_i;
	float coeff;
	float tmp;

	for (int k = threadIdx.x; k < dims_out.y; k+=blockDim.x)
	{
		coeff = C*(float)k;
		sum_real = 0.0f;
		sum_imag = 0.0f;
		for (int n = 0; n < dims_in.y; n++)
		{
			__sincosf(coeff*n,&twi_i,&twi_r);
			tmp = data[n].x * twi_i;
			sum_real += __fmaf_rn(data[n].x, twi_r, -twi_i * data[n].y);
			sum_imag += __fmaf_rn(data[n].y, twi_r, tmp);
		}

		// Not sure if an async write, or storage to a shared mem temp would be faster.
		inplace_image[blockIdx.x + k * (dims_out.w/2)].x = sum_real;
		inplace_image[blockIdx.x + k * (dims_out.w/2)].y = sum_imag;
	}


	return;

}

void DFTbyDecomposition::DFT_C2C_WithPadding_rdx2()
{

	// FIXME when adding real space complex images
	MyAssertTrue( input_image.is_in_memory_gpu, "Input image is in not on the GPU!");
	MyAssertTrue( output_image.is_in_memory_gpu, "Output image is in not on the GPU!");


	int threadsPerBlock = input_image.dims.x; // FIXME make sure its a multiple of 32
	int gridDims = output_image.dims.w/2;

	int shared_mem = sizeof(cufftComplex)*input_image.dims.x;

	float C = -2*PIf/output_image.dims.x*2;
	DFT_C2C_WithPaddingKernel_rdx2<< <gridDims, threadsPerBlock, shared_mem, cudaStreamPerThread>> > ( output_image.complex_values_gpu, input_image.dims, output_image.dims, C);
	cudaStreamSynchronize(cudaStreamPerThread);



}

__global__ void DFT_C2C_WithPaddingKernel_rdx2(cufftComplex* inplace_image, int4 dims_in, int4 dims_out, float C)
{

	// Initialize the shared memory, assuming everying matches the input data X size in
	// Check that setting cudaFuncSetSharedMemConfig  to 8byte makes any diff for complex reads
	extern __shared__ cufftComplex c[];
	cufftComplex* data = c;


	int x = threadIdx.x;
	int pixel_out = (dims_out.w/2)*blockIdx.x;

	data[x] = __ldg((const cufftComplex *)&inplace_image[pixel_out + x]);
	__syncthreads();
//
//	 Loop over N updating the actual twiddle value along the way. This might lead to accuracy problems.
	cufftComplex sum;
	cufftComplex eve;
	float twi_r;
	float twi_i;
	float coeff;
	float tmp;

	for (int k = threadIdx.x; k < dims_out.w/4; k+=blockDim.x)
	{
		// get the even DFT
		coeff = C*(float)k;
		sum.x = 0.0f;
		sum.y = 0.0f;
		for (int n = 0; n < dims_in.y; n+=2)
		{
			__sincosf(coeff*n,&twi_i,&twi_r);
			tmp = data[n].x * twi_i;
			sum.x += __fmaf_rn(data[n].x, twi_r, -twi_i * data[n].y);
			sum.y += __fmaf_rn(data[n].y, twi_r, tmp);
		}

		eve = sum;

		// get the odd DFT
		sum.x = 0.0f;
		sum.y = 0.0f;
		for (int n = 1; n < dims_in.y; n+=2)
		{
			__sincosf(coeff*n,&twi_i,&twi_r);
			tmp = data[n].x * twi_i;
			sum.x += __fmaf_rn(data[n].x, twi_r, -twi_i * data[n].y);
			sum.y += __fmaf_rn(data[n].y, twi_r, tmp);
		}

		// Get the twiddle for the combined radix
		__sincosf(coeff/2.0f,&twi_i,&twi_r);
		// Multiply the odd
		tmp = sum.x * twi_i;
		sum.x = __fmaf_rn(sum.x, twi_r, -twi_i * sum.y);
		sum.y = __fmaf_rn(sum.y, twi_r, tmp);

		inplace_image[pixel_out + k].x = eve.x + sum.x;
		inplace_image[pixel_out + k].y = eve.y + sum.y;

		inplace_image[pixel_out + k + dims_out.w/4].x = eve.x - sum.x;
		inplace_image[pixel_out + k + dims_out.w/4].y = eve.y - sum.y;

	}



	return;

}


void DFTbyDecomposition::FFT_R2C_WithPadding_strided()
{

	// This is the first set of 1d ffts when the input data are real valued, accessing the strided dimension. Since we need the full length, it will actually run a C2C xform

	// FIXME when adding real space complex images
	MyAssertTrue( input_image.is_in_memory_gpu, "Input image is in not on the GPU!");
	MyAssertTrue( output_image.is_in_memory_gpu, "Output image is in not on the GPU!");

	// Elements per thread must be [2,32]
    const int ept = 2;

    // FFts per block. Might be able to re-use twiddles but prob more mem intensive. TODO test me and also evaluate memory size
    const int ffts_per_block = 1; // 1 is the default.

    // For now consider the simplest launch params, where one input element is handled per thread.
    MyAssertFalse(input_image.dims.y % ept, "The elements per thread is not a divisor of the input y-dimension.");
	int threadsPerBlock = input_image.dims.y / ept; // FIXME make sure its a multiple of 32
	int gridDims = input_image.dims.x;

	// For the twiddle factors ahead of the P size ffts
	float CN = -2*PIf/output_image.dims.y;
	int   IQ = output_image.dims.y / input_image.dims.y; // FIXME assuming for now this is already divisible

    // FFT is defined, its: size, type, direction, precision. Block() operator informs that FFT
    // will be executed on block level. Shared memory is required for co-operation between threads.


	if (input_image.dims.y == 256)
	{
	    using FFT = decltype(FFT_256() + Direction<fft_direction::forward>() );
		int shared_mem = sizeof(FFT::value_type)*(input_image.dims.y) + FFT::shared_memory_size;
		block_fft_kernel_R2C_strided<FFT><< <gridDims, threadsPerBlock, shared_mem, cudaStreamPerThread>> > ( (typename FFT::input_type*)input_image.real_values_gpu,  (typename FFT::output_type*)output_image.complex_values_gpu, input_image.dims, output_image.dims, CN,IQ);



	}
	else if (input_image.dims.y == 16)
	{
	    using FFT = decltype(FFT_16() + Direction<fft_direction::forward>() );
		int shared_mem = sizeof(float)*(2+input_image.dims.y) + FFT::shared_memory_size;
		block_fft_kernel_R2C_strided<FFT><< <gridDims, threadsPerBlock, shared_mem, cudaStreamPerThread>> > ( (typename FFT::input_type*)input_image.real_values_gpu,  (typename FFT::output_type*)output_image.complex_values_gpu, input_image.dims, output_image.dims, CN,IQ);



	}
	else
	{
		exit(-1);
	}




}

template<class FFT>
//__launch_bounds__(FFT::max_threads_per_block) __global__
__global__ void block_fft_kernel_R2C_strided(typename FFT::input_type* input_values, typename FFT::output_type* output_values, int4 dims_in, int4 dims_out, float twid_constant, int n_sectors)
{

//	// Initialize the shared memory, assuming everying matches the input data X size in
    using complex_type = typename FFT::value_type;

	extern __shared__  complex_type real_data[];
	complex_type* shared_mem_work = (complex_type*)&real_data[dims_in.y];


	// Memory used by FFT
	complex_type twiddle;
    complex_type thread_data[FFT::elements_per_thread];

    // To re-map the thread index to the data
    int input_MAP[FFT::elements_per_thread];
    // To re-map the decomposed frequency to the full output frequency
    int output_MAP[FFT::elements_per_thread];
    // For a given decomposed fragment
    float twiddle_factors_args[FFT::elements_per_thread];

    // This way reads are
    int i;

    for (i = 0; i < FFT::elements_per_thread; i++)
    {
    	// index into the input data
    	input_MAP[i] = threadIdx.x + i * (size_of<FFT>::value / FFT::elements_per_thread);
		output_MAP[i] = n_sectors * input_MAP[i];
		twiddle_factors_args[i] = twid_constant * input_MAP[i];

		// Unpack the floats and move from shared mem into the register space.in R2C this would happen anyway as a preprocessing step.
		real_data[input_MAP[i]].x = __ldg((const cufftReal *)&input_values[blockIdx.x + input_MAP[i] * (dims_in.w)]);
		real_data[input_MAP[i]].y = 0.0f;
    }
	__syncthreads();


	// this data will be re-used for each n_sectors FFTs
    for (i = 0; i < FFT::elements_per_thread; i++)
    {
		thread_data[i] = real_data[input_MAP[i]];
    }


    // For loop zero the twiddles don't need to be computed
    FFT().execute(thread_data, shared_mem_work);

    // The memory access is strided anyway so just send to global
    for (i = 0; i < FFT::elements_per_thread; i++)
    {
        output_values[blockIdx.x + output_MAP[i] * (dims_out.w/2)].x = (float)thread_data[i].x;
        output_values[blockIdx.x + output_MAP[i] * (dims_out.w/2)].y = (float)thread_data[i].y;
    }

    // For the other fragments we need the initial twiddle
	for (int fft_fragment = 0; fft_fragment < n_sectors; fft_fragment++)
	{
	    for (i = 0; i < FFT::elements_per_thread; i++)
	    {
			// Pre shift with twiddle
			__sincosf(twiddle_factors_args[0]*fft_fragment,&twiddle.x,&twiddle.y);
			twiddle *= real_data[input_MAP[i]]; // Only the inplace operators are included in cufftdx::types TODO expand
		    thread_data[i] = twiddle;
	    }

	      FFT().execute(thread_data, shared_mem_work);

		for (i = 0; i < FFT::elements_per_thread; i++)
		{
		      output_values[blockIdx.x + (fft_fragment + output_MAP[i]) * (dims_out.w/2)].x = (float)thread_data[i].x;
		      output_values[blockIdx.x + (fft_fragment + output_MAP[i]) * (dims_out.w/2)].y = (float)thread_data[i].y;
		}

	}



	return;

}
