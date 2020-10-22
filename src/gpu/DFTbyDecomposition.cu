/*
 * DFTbyDecomposition.cpp
 *
 *  Created on: Oct 21, 2020
 *      Author: himesb
 */

#include "gpu_core_headers.h"

__global__ void DFT_R2C_WithPaddingKernel(cufftReal* input_values, cufftComplex* output_values, int4 dims_in, int4 dims_out, float C);
__global__ void DFT_C2C_WithPaddingKernel(cufftComplex* input_values, int4 dims_in, int4 dims_out, float C);



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

	pre_checkErrorsAndTimingWithSynchronization(cudaStreamPerThread);

	int div = 1;
	int threadsPerBlock = input_image.dims.x; // FIXME make sure its a multiple of 32
	int gridDims = input_image.dims.y;
//	dim3 gridDims = dim3((output_image.dims.w/2 + threadsPerBlock - 1) / threadsPerBlock,
//					  	1, 1);
//  output_image.dims.y
	int shared_mem = sizeof(float)*input_image.dims.x;
wxPrintf("DIMS DIMS %d\n",output_image.dims.w/2);
	float C = -2*PIf/output_image.dims.x;
	DFT_R2C_WithPaddingKernel<< <gridDims, threadsPerBlock, shared_mem, cudaStreamPerThread>> > ( input_image.real_values_gpu,  output_image.complex_values_gpu, input_image.dims, output_image.dims, C);
	cudaStreamSynchronize(cudaStreamPerThread);


	checkErrorsAndTimingWithSynchronization(cudaStreamPerThread);

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

	pre_checkErrorsAndTimingWithSynchronization(cudaStreamPerThread);

	int threadsPerBlock = input_image.dims.y; // FIXME make sure its a multiple of 32
	int gridDims = output_image.dims.w/2;

	int shared_mem = sizeof(cufftComplex)*input_image.dims.y;

	float C = -2*PIf/output_image.dims.y;
	DFT_C2C_WithPaddingKernel<< <gridDims, threadsPerBlock, shared_mem, cudaStreamPerThread>> > ( output_image.complex_values_gpu, input_image.dims, output_image.dims, C);
	cudaStreamSynchronize(cudaStreamPerThread);


	checkErrorsAndTimingWithSynchronization(cudaStreamPerThread);

}

__global__ void DFT_C2C_WithPaddingKernel(cufftComplex* inplace_image, int4 dims_in, int4 dims_out, float C)
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
