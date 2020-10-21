/*
 * DFTbyDecomposition.cpp
 *
 *  Created on: Oct 21, 2020
 *      Author: himesb
 */

#include "gpu_core_headers.h"

__global__ void DFT_R2C_WithPaddingKernel(cufftReal* input_values, cufftComplex* output_values, float* input_twiddles, int4 dims_in, int4 dims_out);


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

void DFTbyDecomposition::SetTwiddleAndOutputs()
{
	MyAssertTrue(is_set_gpu_images, "You did not SetGpuImages which must happen before setting Twiddles");


	int n_twiddles = output_image.dims.w;
	float* tmp_twiddles = new float[n_twiddles];
	float C = {2 * PIf / output_image.dims.x};
	for (int i = 0; i < n_twiddles/2 ; i+=2)
	{
		sincosf(C * (float)i, &tmp_twiddles[i+1], &tmp_twiddles[i]);
	}

	// If twiddles or chunks of output are changed, then nT,nI,nO should be set and tracked. then shared mem = each of these summed, and should be passed to the kernel to know where in shared mem things liv.
	cudaErr(cudaMalloc((void **) &twiddles, sizeof(float)*n_twiddles));

	cudaErr(cudaMemcpyAsync(twiddles,tmp_twiddles, n_twiddles, cudaMemcpyHostToDevice, cudaStreamPerThread));
	is_set_twiddles = true;

	delete [] tmp_twiddles;

	cudaErr(cudaStreamSynchronize(cudaStreamPerThread));
//	cudaErr(cudaMallocManaged((void **)&output_real, sizeof(float)*input_image.dims.x));
//	cudaErr(cudaMallocManaged((void **)&output_imag, sizeof(float)*input_image.dims.x));
//	is_set_outputs = true;

	shared_mem = sizeof(float)*input_image.dims.x + sizeof(float)*input_image.dims.w*2;


}

void DFTbyDecomposition::DFT_R2C_WithPadding()
{

	// FIXME when adding real space complex images
	MyAssertTrue( input_image.is_in_memory_gpu, "Input image is in not on the GPU!");
	MyAssertTrue( output_image.is_in_memory_gpu, "Output image is in not on the GPU!");

	pre_checkErrorsAndTimingWithSynchronization(cudaStreamPerThread);

	int div = 1;
	int threadsPerBlock = input_image.dims.x;
	dim3 gridDims = dim3((output_image.dims.w/2 + threadsPerBlock - 1) / threadsPerBlock,
					  	1, 1);
//  output_image.dims.y
	wxPrintf("Threads n grids %d %d %d %d shared mem %d\n",threadsPerBlock,gridDims.x,gridDims.y,gridDims.z,shared_mem);

	output_image.printVal("Before",0);
	output_image.printVal("Before",10);

	DFT_R2C_WithPaddingKernel<< <gridDims, threadsPerBlock, shared_mem, cudaStreamPerThread>> > ( input_image.real_values_gpu,  output_image.complex_values_gpu, twiddles, input_image.dims, output_image.dims);
	cudaStreamSynchronize(cudaStreamPerThread);
	output_image.printVal("After",0);
	output_image.printVal("After",10);


	checkErrorsAndTimingWithSynchronization(cudaStreamPerThread);

}

__global__ void DFT_R2C_WithPaddingKernel(cufftReal* input_values, cufftComplex* output_values, float* input_twiddles, int4 dims_in, int4 dims_out)
{

//	// Initialize the shared memory, assuming everying matches the input data X size in
	extern __shared__ float s[];
//    __shared__ float s[2056];

	float *data   = s;
	float *twiddles = (float*)&data[dims_in.x];
//	float *output   = (float*)&twiddles[dims_out.w];


	int x = threadIdx.x;
	int r = 2*x;
	int i = 2*x+1;
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	if (x > dims_in.x)
	if (k > dims_out.w/2) return;

	// Every thread works on one spatial freq (k) looping over the N real inputs
	// There are N threads in a block
	// Each real vector is divided into NOUT/N
	// There are input dims y blocks in the Y direction

	data[x] = __ldg((const float *)&input_values[dims_in.w*blockIdx.y + x]);
	twiddles[r] = __ldg((const float *)&input_twiddles[k]);
	twiddles[i] = __ldg((const float *)&input_twiddles[k]);

	__syncthreads();
//
//	 Loop over N updating the actual twiddle value along the way. This might lead to accuracy problems.
	float sum_real = 0.0f;
	float sum_imag = 0.0f;
	float twi_real = twiddles[r];
	float twi_imag = twiddles[i];
	float tmp;


	for (int n = 0; n < dims_in.x; n++)
	{
//		sum_real += data[n] * twi_real;
//		sum_imag += data[n] * twi_imag;
//		tmp = twi_real;
//		twi_real = twi_real*twiddles[r] - twi_imag*twiddles[i];
//		twi_imag = twi_imag*twiddles[r] + tmp * twiddles[i];

//		twi_real = __fmaf_rn(twi_real,twiddles[r],-twi_imag*twiddles[i]);
//		twi_imag = __fmaf_rn(twi_imag,twiddles[r], tmp*twiddles[i]);
		sum_real += data[n] * cosf(-2*PIf*n*k/dims_out.x);
		sum_imag += data[n] * sinf(-2*PIf*n*k/dims_out.x);
	}

//	output[r] = sum_real;
//	output[i] = sum_imag;
//
	__syncthreads();

//		if (k==3) {output_values[0].x=(float)dims_out.w*blockIdx.y;}
	//	if (k==3) {output_values[0].y=(float)k;}

		output_values[dims_out.w*blockIdx.y + k].x = sum_real;//output[r];
		output_values[dims_out.w*blockIdx.y + k].y = sum_imag;//output[r];

//		output_values[dims_out.w*blockIdx.y + x].y =  -3;//sum_imag;//output[i];




	return;


}
