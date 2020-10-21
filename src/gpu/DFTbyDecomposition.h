/*
 * DFTbyDecomposition.h
 *
 *  Created on: Oct 21, 2020
 *      Author: himesb
 */


#ifndef SRC_GPU_DFTBYDECOMPOSITION_H_
#define SRC_GPU_DFTBYDECOMPOSITION_H_

class DFTbyDecomposition
{

public:

	DFTbyDecomposition();
	virtual ~DFTbyDecomposition();
	DFTbyDecomposition(const DFTbyDecomposition &other);
	DFTbyDecomposition& operator=(const DFTbyDecomposition &other);

	// These properties are for testing and should not be in the final class
	StopWatch timer;
	GpuImage baseline_image;
	GpuImage input_image;
	GpuImage output_image;
	int2 dims_input;
	int2 dims_output;

	int numel;
	// These will point to the memory in the GpuImages, but won't have memory allocated for them.
	cufftReal* input_image_real;
	cufftComplex* output_image_complex;


	// To be used in development where a fake image is created.
	void InitTestCase(int wanted_input_size_x, int wanted_input_size_y, int wanted_output_size_x, int wanted_output_size_y);
	void SetGpuImages(Image& cpu_input, Image& cpu_output);
};

#endif /* SRC_GPU_DFTBYDECOMPOSITION_H_ */
