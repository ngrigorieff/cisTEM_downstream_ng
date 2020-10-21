/*
 * DFTbyDecomposition.cpp
 *
 *  Created on: Oct 21, 2020
 *      Author: himesb
 */

#include "gpu_core_headers.h"

DFTbyDecomposition::DFTbyDecomposition() // @suppress("Class members should be properly initialized")
{

}

DFTbyDecomposition::~DFTbyDecomposition()
{

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

}
