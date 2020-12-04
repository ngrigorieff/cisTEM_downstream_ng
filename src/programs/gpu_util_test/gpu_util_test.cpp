#include "../../core/core_headers.h"


class
GpuUtilTest : public MyApp
{

	public:

	bool DoCalculation();
	void DoInteractiveUserInput();
	void TemplateMatchingStandalone(int nGPUs, int nThreads);
	void createImageAddOne();
	void DFTbyDecomp();
	void FFTwithRotation();

	private:
};


IMPLEMENT_APP(GpuUtilTest)

// override the DoInteractiveUserInput

void GpuUtilTest::DoInteractiveUserInput()
{

}

// override the do calculation method which will be what is actually run..

bool GpuUtilTest::DoCalculation()
{

  wxPrintf("GpuUtilTest is running!\n");

//  this->createImageAddOne();
  int nThreads = 1;
  int nGPUs = 1;
//  this->TemplateMatchingStandalone(nThreads, nGPUs);
 // this->DFTbyDecomp();
 this->FFTwithRotation();
  int gpuID = 0;
  wxPrintf("I made it here\n");

  

}	

void GpuUtilTest::TemplateMatchingStandalone(int nThreads, int nGPUs)
{

	int number_of_jobs_per_image_in_gui = 1;
	nThreads = 2;
	nGPUs = 1;
	int minPos = 0;
	int maxPos = 60;
	int incPos = 60 / (nThreads*nGPUs); // FIXME

	ProgressBar *my_progress;
	long total_correlation_positions = 1000; // Not calculated properly: FIXME
	long current_correlation_position = 0;

	if (is_running_locally == true)
	{
		my_progress = new ProgressBar(total_correlation_positions);
	}

//	DeviceManager gpuDev(nGPUs);
//    omp_set_num_threads(nThreads * gpuDev.nGPUs);  // create as many CPU threads as there are CUDA devices
//	#pragma omp parallel
//    {

		TemplateMatchingCore GPU[nThreads];
//    	TemplateMatchingCore GPU(number_of_jobs_per_image_in_gui);


		// Check the number of available gpus
		DeviceManager gpuDev;
		gpuDev.Init(nGPUs);


		#pragma omp parallel num_threads(nThreads)
		{

			int tIDX = ReturnThreadNumberOfCurrentThread();

	    	Image template_reconstruction;
	    	Image projection_filter;
	    	Image input_image;
	    	Image current_projection;
	    	Image padded_reference;
	    	Image max_intensity_projection;
	    	Image best_psi;
	    	Image best_phi;
	    	Image best_theta;
	    	Image best_defocus;
	    	Image best_pixel_size;

	    	ImageFile template_reconstruction_file;

	    	template_reconstruction_file.OpenFile("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/template_reconstruction.mrc", false);
	    	template_reconstruction.ReadSlices(&template_reconstruction_file, 1, template_reconstruction_file.ReturnNumberOfSlices());

			projection_filter.QuickAndDirtyReadSlice("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/projection_filter.mrc",1);
			input_image.QuickAndDirtyReadSlice("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/input_image.mrc",1);
			current_projection.QuickAndDirtyReadSlice("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/current_projection.mrc",1);
			padded_reference.QuickAndDirtyReadSlice("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/padded_reference.mrc",1);


			input_image.Resize(4096,4096,1,0.0f);
			padded_reference.CopyFrom(&input_image);
			// These are all blank to start
			max_intensity_projection.CopyFrom(&input_image);
			best_psi.CopyFrom(&input_image);
			best_phi.CopyFrom(&input_image);
			best_theta.CopyFrom(&input_image);
			best_pixel_size.CopyFrom(&input_image);
			best_defocus.CopyFrom(&input_image);

	//		// These should be in Fourier space, but were ifft to save
			template_reconstruction.ForwardFFT();
			template_reconstruction.SwapRealSpaceQuadrants();
			input_image.ForwardFFT();
			input_image.SwapRealSpaceQuadrants();
			projection_filter.ForwardFFT();



			// These also were set up prior to entering the GPU loop
			EulerSearch	global_euler_search;
			AnglesAndShifts angles;

			float angular_step = 2.5f;
			float psi_step = 1.5f;
			float pixel_size = 1.5;
			float pixel_size_search_range = 0.0f;
			float pixel_size_step = 0.001f;
			float defocus_search_range = 0.0f;
			float defocus_step = 200.0f;
			float defocus1 = 19880.0f;
			float defocus2 = 18910.0f;
			long first_search_position = 0 + (tIDX*incPos);
			long last_search_position = incPos + (tIDX*incPos);

			if (tIDX == (nThreads*nGPUs - 1)) last_search_position = maxPos;

			float high_resolution_limit_search = 2.0f * pixel_size;
			int best_parameters_to_keep = 1;
			float psi_start = 0.0f;
			float psi_max = 360.0f;
			ParameterMap parameter_map; // needed for euler search init
			parameter_map.SetAllTrue();
			global_euler_search.InitGrid("O", angular_step, 0.0, 0.0, psi_max, psi_step, psi_start, pixel_size / high_resolution_limit_search, parameter_map, best_parameters_to_keep);
			global_euler_search.CalculateGridSearchPositions(false);

			wxDateTime 	overall_start;
			wxDateTime 	overall_finish;
			overall_start = wxDateTime::Now();
			gpuDev.SetGpu(tIDX);

			int max_padding = 0;
			const float histogram_min = -20.0f;
			const float histogram_max = 50.0f;
			const int histogram_number_of_points = 1024;
			float histogram_step;
			float histogram_min_scaled, histogram_step_scaled;
			histogram_step = (histogram_max - histogram_min) / float(histogram_number_of_points);

			histogram_min_scaled = histogram_min / double(sqrt(input_image.logical_x_dimension * input_image.logical_y_dimension));
			histogram_step_scaled = histogram_step / double(sqrt(input_image.logical_x_dimension * input_image.logical_y_dimension));

			GPU[tIDX].Init(this, template_reconstruction, input_image, current_projection,
					pixel_size_search_range, pixel_size_step, pixel_size,
					defocus_search_range, defocus_step, defocus1, defocus2,
					psi_max, psi_start, psi_step,
					angles, global_euler_search,
					histogram_min_scaled, histogram_step_scaled, histogram_number_of_points,
					max_padding, first_search_position, last_search_position,
					my_progress, total_correlation_positions, is_running_locally);

			int size_i = 0;
			int defocus_i = 0;


			GPU[tIDX].RunInnerLoop(projection_filter, size_i, defocus_i, tIDX, current_correlation_position);

			long* histogram_data = new long[GPU[tIDX].histogram.histogram_n_bins];
			for (int iBin = 0; iBin < GPU[tIDX].histogram.histogram_n_bins; iBin++)
			{
				histogram_data[iBin] = 0;
			}
			GPU[tIDX].histogram.CopyToHostAndAdd(histogram_data);
			std::string fileNameOUT4 = "/tmp/tmpMip" + std::to_string(tIDX) + ".mrc";
			max_intensity_projection.QuickAndDirtyWriteSlice(fileNameOUT4,1,true,1.5);

			wxPrintf("\n\n\tTimings: Overall: %s\n",(wxDateTime::Now()-overall_start).Format());


    } // end of omp block
}

void GpuUtilTest::createImageAddOne()
{


	bool do_all_tests = false;
	bool do_shift = false;
	bool do_fft = true;
	bool do_scale = false;
	bool do_swap = false;
	bool do_pad = false;


	int wanted_number_of_gpus = 1;
	int wanted_number_threads_per_gpu = 1;

	DeviceManager gpuDev(wanted_number_of_gpus);



//
//
//	wxPrintf("Found %d gpus to use\n",gDev.nGPUs);
//    ContextManager CM[gDev.nGPUs];
//
//
	#pragma omp parallel num_threads(wanted_number_threads_per_gpu * gpuDev.nGPUs)
    {

		int threadIDX = ReturnThreadNumberOfCurrentThread();
		gpuDev.SetGpu(threadIDX);

		Image cpu_image_half;
		Image cpu_work_half;
		GpuImage d_image_half;

		Image cpu_image_full;
		Image cpu_work_full;
		GpuImage d_image_full;

		EmpiricalDistribution eDist;

		float t_rmsd;


//		GpuImage mask;
//		cpu_image_full.QuickAndDirtyReadSlice("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/oval_full.mrc",1);
//		d_image_full.Init(cpu_image_full);
//		d_image_full.CopyHostToDevice();
//		d_image_full.ForwardFFT(false);
//		float wm =  d_image_full.ReturnSumSquareModulusComplexValuesMask();
//		float nm = d_image_full.ReturnSumSquareModulusComplexValues();
//		wxPrintf("Compare the slow loop and the cublas %3.3e, %3.3e\n",wm,nm);
//		mask.Allocate(512,512,1,true);
//		mask.Wait();
//		mask.ReturnSumSquareModulusComplexValuesMask();
//		mask.mask_CSOS->printVal("after copy val 0 is",0);
//		mask.mask_CSOS->printVal("after copy val 100 is",100);
//		mask.mask_CSOS->QuickAndDirtyWriteSlices("/tmp/mask.mrc",1,1);


		if (do_shift || do_all_tests)
		{

			// full size image

			// Read in the CPU image
			cpu_image_full.QuickAndDirtyReadSlice("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/oval_full.mrc",1);

			// Copy of the image to operate on using cpu method
			cpu_work_full.CopyFrom(&cpu_image_full);

			// Initialize gpu image, and copy to device
			d_image_full.Init(cpu_image_full);
			d_image_full.CopyHostToDevice();

			cpu_work_full.PhaseShift(10, -200,0);
			d_image_full.PhaseShift(10,-200,0);

			d_image_full.Wait();
//			d_image_full.CopyDeviceToHost(true, true);
			d_image_full.CopyDeviceToHost(true, true);

			d_image_full.Wait();
			d_image_full.QuickAndDirtyWriteSlices("/tmp/oval_full_shift_fromGPU.mrc",1,1);

			cpu_image_full.QuickAndDirtyWriteSlice("/tmp/oval_full_shift.mrc",1);

			cpu_work_full.SubtractImage(&cpu_image_full);

			t_rmsd = sqrtf(cpu_work_full.ReturnSumOfSquares(0,0,0,0,false));
			wxPrintf("RMSD between cpu phase shift and gpu phase shift for 512 x 512 is %3.3e\n", t_rmsd);


			// Half size image

			cpu_image_half.QuickAndDirtyReadSlice("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/oval_half.mrc",1);
			// Copy of the image to operate on using cpu method
			cpu_work_half.CopyFrom(&cpu_image_half);

			// Initialize gpu image, and copy to device
			d_image_half.Init(cpu_image_half);
			d_image_half.CopyHostToDevice();

			cpu_work_half.PhaseShift(10, -200,0);
			d_image_half.PhaseShift(10,-200,0);

			d_image_half.Wait();
//			d_image_half.CopyDeviceToHost(true, true);
			d_image_half.CopyDeviceToHost(false, true);

			d_image_half.Wait();
			d_image_half.QuickAndDirtyWriteSlices("/tmp/oval_half_shift_fromGPU.mrc",1,1);

			cpu_image_half.QuickAndDirtyWriteSlice("/tmp/oval_half_shift.mrc",1);

			cpu_work_half.SubtractImage(&cpu_image_half);

			t_rmsd = sqrtf(cpu_work_half.ReturnSumOfSquares(0,0,0,0,false));
			wxPrintf("RMSD between cpu phase shift and gpu phase shift for 256 x 512 is %3.3e\n", t_rmsd);


		}


		if (do_fft || do_all_tests)
		{

			wxDateTime start;
			start = wxDateTime::Now();

			Image c;
			c.Allocate(514,514,1,true);


			// full size image

			// Read in the CPU image
			cpu_image_full.QuickAndDirtyReadSlice("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/oval_full.mrc",1);
			cpu_image_full.AddConstant(-cpu_image_full.ReturnAverageOfRealValues(0.0f,false));

			// Copy of the image to operate on using cpu method
			cpu_work_full.CopyFrom(&cpu_image_full);

			// Re-use the image and library contexts - this needs to have debug asserts added'			'
			d_image_full.Init(cpu_image_full);

			d_image_full.CopyHostToDevice();

//			d_image_full.Mean();
//			d_image_full.AddConstant(-d_image_full.img_mean);
			bool doNorm = false;
			cpu_work_full.AddConstant(-cpu_work_full.ReturnAverageOfRealValues(0.0f,false));

			cpu_work_full.MultiplyByConstant(1.0f/sqrtf(cpu_work_full.ReturnSumOfSquares()));
			wxPrintf("cpu var before fft is %3.3e\n", (cpu_work_full.ReturnSumOfSquares()));

			cpu_work_full.ForwardFFT(doNorm);
			wxPrintf("cpu var after fft no norm is %3.3e or *= / n^2 %3.3e\n", (cpu_work_full.ReturnSumOfSquares()), cpu_work_full.ReturnSumOfSquares()/cpu_work_full.number_of_real_space_pixels/cpu_work_full.number_of_real_space_pixels);
			cpu_work_full.MultiplyByConstant(1.0f/sqrtf(cpu_work_full.ReturnSumOfSquares()));
			cpu_work_full.BackwardFFT();
			wxPrintf("cpu var after ifft with norm is %3.3e or *= / n^2 %3.3e\n", (cpu_work_full.ReturnSumOfSquares()), cpu_work_full.ReturnSumOfSquares()/(cpu_work_full.number_of_real_space_pixels)/cpu_work_full.number_of_real_space_pixels);


			d_image_full.MultiplyByConstant(1.0f/sqrtf(d_image_full.ReturnSumOfSquares()));
			wxPrintf("gpu var before fft is %3.3e\n", (d_image_full.ReturnSumOfSquares()));

			d_image_full.ForwardFFT(doNorm);
			wxPrintf("gpu var after fft no norm is %3.3e or *= / n %3.3e\n", (d_image_full.ReturnSumSquareModulusComplexValues()), d_image_full.ReturnSumSquareModulusComplexValues()/(d_image_full.number_of_real_space_pixels));
			d_image_full.MultiplyByConstant(1.0f/sqrtf(d_image_full.ReturnSumSquareModulusComplexValues()));
			d_image_full.BackwardFFT();
			wxPrintf("gpu var after ifft with norm is %3.3e or *= / n %3.3e\n", (d_image_full.ReturnSumOfSquares()), d_image_full.ReturnSumOfSquares()/(d_image_full.number_of_real_space_pixels));

			exit(0);



			cpu_work_full.ForwardFFT(true);
			cpu_work_full.BackwardFFT();

			d_image_full.ForwardFFT(true);
			d_image_full.BackwardFFT();

			d_image_full.Wait();

			d_image_full.CopyDeviceToHost(true, true);
			d_image_full.Wait();

			cpu_image_full.QuickAndDirtyWriteSlice("/tmp/oval_full_fft_ifft.mrc",1);
			cpu_work_full.QuickAndDirtyWriteSlice("/tmp/oval_full_work_fft_ifft.mrc",1);

			cpu_work_full.SubtractImage(&cpu_image_full);
//			cpu_work_full.DivideByConstant(cpu_work_full.number_of_real_space_pixels);

			t_rmsd = sqrtf(cpu_work_full.ReturnSumOfSquares(0,0,0,0,false));
			wxPrintf("RMSD between cpu fft/ifft and gpu  fft/ifft  for 512 x 512 is %3.3e\n", t_rmsd);


			// Half size image

			cpu_image_half.QuickAndDirtyReadSlice("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/oval_half.mrc",1);

			// Copy of the image to operate on using cpu method
			cpu_work_half.CopyFrom(&cpu_image_half);

			// Re-use the image and library contexts - this needs to have debug asserts added
			d_image_half.CopyHostToDevice();

			cpu_work_half.ForwardFFT(true);
			cpu_work_half.BackwardFFT();

			d_image_half.ForwardFFT(true);
			d_image_half.BackwardFFT();

			d_image_half.Wait();

			d_image_half.CopyDeviceToHost(true, true);
			d_image_half.Wait();

			cpu_image_half.QuickAndDirtyWriteSlice("/tmp/oval_half_fft_ifft.mrc",1);
			cpu_work_half.QuickAndDirtyWriteSlice("/tmp/oval_half_work_fft_ifft.mrc",1);


			cpu_work_half.SubtractImage(&cpu_image_half);

			t_rmsd = sqrtf(cpu_work_half.ReturnSumOfSquares(0,0,0,0,false));
			wxPrintf("RMSD between cpu fft/ifft and gpu  fft/ifft  for 256 x 512 is %3.3e\n", t_rmsd);


		}

		if (do_scale || do_all_tests)
		{

			// full size image

			// Read in the CPU image
			cpu_image_full.QuickAndDirtyReadSlice("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/oval_full.mrc",1);

			// Copy of the image to operate on using cpu method
			cpu_work_full.CopyFrom(&cpu_image_full);
			// Re-use the image and library contexts - this needs to have debug asserts added
			d_image_full.CopyHostToDevice();

			cpu_work_full.MultiplyByConstant(PIf);
			d_image_full.MultiplyByConstant(PIf);
			d_image_full.Wait();




			wxPrintf("Real sums are cpu: %f gpu: %f\n",cpu_work_full.ReturnSumOfRealValues(), d_image_full.ReturnSumOfRealValues());
			d_image_full.Wait();

			cpu_work_full.ForwardFFT(true);
			d_image_full.ForwardFFT(true);

			wxPrintf("Complex sums are cpu: %4.4e gpu: %4.4e\n",cpu_work_full.ReturnSumOfSquares(), d_image_full.ReturnSumSquareModulusComplexValues());
			d_image_full.Wait();

			cpu_work_full.BackwardFFT();
			d_image_full.BackwardFFT();

			d_image_full.Wait();



		}

		if (do_pad || do_all_tests)
		{

			// full size image

			// Read in the CPU image
			cpu_image_full.QuickAndDirtyReadSlice("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/oval_full.mrc",1);

			// Copy of the image to operate on using cpu method
			cpu_work_full.CopyFrom(&cpu_image_full);
			// Re-use the image and library contexts - this needs to have debug asserts added
			d_image_full.CopyHostToDevice();

			Image padded;
			Image padded_work;
			padded.Allocate(1024,768,1,true);
			padded.SetToConstant(0.0f);
			padded_work.CopyFrom(&padded);
			GpuImage d_padded(padded);
			d_padded.CopyHostToDevice();

			cpu_work_full.ClipInto(&padded_work,1,false,0.0f,10,-30,0);

			d_image_full.ClipInto(&d_padded,1,false,0.0f,10,-30,0);
			d_padded.Wait();
			d_padded.CopyDeviceToHost(true, true);

			//			d_fft.ClipInto(&d_padded,-1,0,128,0);



			d_padded.Wait();


			padded_work.SubtractImage(&padded);
//			cpu_work_full.DivideByConstant(cpu_work_full.number_of_real_space_pixels);

			t_rmsd = sqrtf(padded_work.ReturnSumOfSquares(0,0,0,0,false));
			wxPrintf("RMSD between padded images 512/512-->1024x768 is %3.3e\n", t_rmsd);



		}
//

		if (do_swap || do_all_tests)
		{

			// full size image

			// Read in the CPU image
			cpu_image_full.QuickAndDirtyReadSlice("/groups/grigorieff/home/himesb/cisTEM_2/cisTEM/trunk/gpu/include/oval_full.mrc",1);

			// Copy of the image to operate on using cpu method
			cpu_work_full.CopyFrom(&cpu_image_full);
			// Re-use the image and library contexts - this needs to have debug asserts added
			d_image_full.CopyHostToDevice();

			cpu_work_full.SwapRealSpaceQuadrants();
			d_image_full.SwapRealSpaceQuadrants();
			d_image_full.Wait();




			wxPrintf("Real sums after swapping are cpu: %f gpu: %f\n",cpu_work_full.ReturnSumOfRealValues(), d_image_full.ReturnSumOfRealValues());
			d_image_full.Wait();




		}

//	wxPrintf("Making an image of zeros in host memory\n");
//	Image zeros;
//	zeros.Allocate(512,512,512,true);
//	zeros.SetToConstant(1.f);
//
//	Image zeros2;
//	zeros2.real_values = zeros.real_values;
//
//	  wxPrintf("\n\nhost ZEROS2 %p with value %f\n\n", &zeros2.real_values,  zeros2.real_values[0]);
//
//
//	wxPrintf("Initializing a GpuImage\n");
//
//	GpuImage d_zeros;
//
//	wxPrintf("Set up a GpuImage with size %d\n", d_zeros.logical_x_dimension);
//	d_zeros.CopyFromCpuImage(zeros);
//
//
//	wxPrintf("copied the Image into the GpuImage with size %d\n", d_zeros.logical_x_dimension);
//
////	d_zeros.CopyVolumeHostToDevice();
//	d_zeros.CopyHostToDevice();
//
//	wxPrintf("I copied to the device\n");
//
//	wxPrintf("Now I'll try to multipy by 5\n");
//	d_zeros.MultiplyByScalar(5.0f);
//
//	wxPrintf("Now I'll try to copy to the host\n");
////	d_zeros.CopyVolumeDeviceToHost();
//	d_zeros.CopyDeviceToHost();
//
//
//	wxPrintf("Checking the values in the copied array %f\n",d_zeros.real_values[10]);
//
//	zeros.real_values = d_zeros.real_values;
//
//	zeros.QuickAndDirtyWriteSlices("TestGpuOutx5.mrc",1,512);

    } // end of parallel omp block

}

void GpuUtilTest::DFTbyDecomp()
{

	bool complex_strided = false;
	bool do_rotate = true;

	DFTbyDecomposition DFT;
	int wanted_input_size_x = 512;
	int wanted_input_size_y = wanted_input_size_x;
	int wanted_output_size_x = 4096;
	int wanted_output_size_y = 4096;
	int wanted_number_of_iterations = 1000;

	DFT.InitTestCase(wanted_input_size_x,wanted_input_size_y,wanted_output_size_x,wanted_output_size_y);

	Image cpu_image_in, cpu_image_out;
	Image gpu_image_in, gpu_image_out;
	Image baseline_image, baseline_gpu_image;
	GpuImage d_gpu_in, d_gpu_image_out;
	cpu_image_in.Allocate(wanted_input_size_x,wanted_input_size_y,true);
	cpu_image_out.Allocate(wanted_output_size_x,wanted_output_size_y,true);

	RandomNumberGenerator rg(PIf);
	for (int current_pixel=0; current_pixel < cpu_image_in.real_memory_allocated; current_pixel++)
	{
		cpu_image_in.real_values[current_pixel] = rg.GetNormalRandomSTD(0.0,1.0);
	}

	// Copy of gpu images (on host)
	gpu_image_in.CopyFrom(&cpu_image_in);
	gpu_image_out.CopyFrom(&cpu_image_out);


	cpu_image_in.QuickAndDirtyWriteSlice("checkRandom.mrc", 1, true, 1.0f);

	// Check the accuracy relative for MKL and cuFFT
	baseline_image.CopyFrom(&cpu_image_in);
	baseline_gpu_image.CopyFrom(&gpu_image_in);


	cpu_image_in.ForwardFFT(true);
	cpu_image_in.BackwardFFT();
	baseline_image.SubtractImage(&cpu_image_in);
	float t_rmsd = sqrtf(baseline_image.ReturnSumOfSquares(0,0,0,0,false));
	wxPrintf("RMSD between cpu image before/after FFT pair is %3.3e\n", t_rmsd);

	// Refresh the cpu input
	for (int current_pixel=0; current_pixel < cpu_image_in.real_memory_allocated; current_pixel++)
	{
		cpu_image_in.real_values[current_pixel] = baseline_gpu_image.real_values[current_pixel];
	}

	d_gpu_in.CopyFromCpuImage(gpu_image_in);
	d_gpu_in.CopyHostToDevice();
	d_gpu_in.ForwardFFT(true);
	d_gpu_in.BackwardFFT();
	d_gpu_in.CopyDeviceToHost(true, true);
	d_gpu_in.Wait();
	baseline_gpu_image.SubtractImage(&gpu_image_in);
	t_rmsd = sqrtf(baseline_gpu_image.ReturnSumOfSquares(0,0,0,0,false));
	wxPrintf("RMSD between gpu image before/after FFT pair is %3.3e\n", t_rmsd);

	// Refresh the gpu input
	for (int current_pixel=0; current_pixel < cpu_image_in.real_memory_allocated; current_pixel++)
	{
		gpu_image_in.real_values[current_pixel] = cpu_image_in.real_values[current_pixel];
	}


	// Associate the test images in the DFT object. The input is copied device to host, the output is allocated directly on the GPU
	DFT.SetGpuImages(gpu_image_in, gpu_image_out);
	DFT.AllocateRotatedBuffer();
	cudaDeviceSynchronize();


	if (complex_strided || do_rotate)
	{
		MyPrintWithDetails("In ROTATE");
		DFT.FFT_R2C_WithPadding(do_rotate);
//		DFT.DFT_R2C_WithPadding();

		MyPrintWithDetails("");

	}
	else
	{
		MyPrintWithDetails("");
//		DFT.DFT_R2C_WithPadding_strided();

		DFT.FFT_R2C_WithPadding_strided(do_rotate);
	}
	MyPrintWithDetails("");
//
//	// Check the first dimension
//	Image first_dim, last_dim;
//	first_dim.Allocate(cpu_image_out.logical_x_dimension,1,1,true);
//	last_dim.Allocate(cpu_image_out.logical_x_dimension,1,1,true);
//	int last_index = (cpu_image_in.logical_x_dimension + cpu_image_in.padding_jump_value)*(cpu_image_in.logical_y_dimension-2);
//	int actual_pixel;
//	for (int current_pixel=0; current_pixel < first_dim.real_memory_allocated; current_pixel++)
//	{
//		if (complex_strided)
//		{
//			actual_pixel = current_pixel;
//		}
//		else
//		{
//			actual_pixel = current_pixel * (first_dim.logical_x_dimension + first_dim.padding_jump_value);
//		}
//		if (current_pixel < cpu_image_in.logical_x_dimension)
//		{
//			first_dim.real_values[current_pixel] = cpu_image_in.real_values[actual_pixel];
//			last_dim.real_values[current_pixel] = 0.0f;//cpu_image_in.real_values[last_index+current_pixel];
//
//		}
//		else
//		{
//			first_dim.real_values[current_pixel] = 0.0f;
//			last_dim.real_values[current_pixel] = 0.0f;
//
//		}
//	}
//
//	first_dim.ForwardFFT(false);
//	last_dim.ForwardFFT(false);
//	wxPrintf("Half dim %ld\n",first_dim.real_memory_allocated/2);
//	for (int current_pixel=0; current_pixel < 16; current_pixel++)
//	{
//		wxPrintf("Transformed Vals %d, %3.3e,%3.3e  %3.3e :LAST: %d %3.3e,%3.3e  %3.3e\n",
//				current_pixel, first_dim.real_values[2*current_pixel],first_dim.real_values[2*current_pixel+1],
//				sqrtf(powf(first_dim.real_values[2*current_pixel],2)+powf(first_dim.real_values[2*current_pixel+1],2)),
//				current_pixel + last_index,last_dim.real_values[2*current_pixel],last_dim.real_values[2*current_pixel+1],
//				sqrtf(powf(last_dim.real_values[2*current_pixel],2)+powf(last_dim.real_values[2*current_pixel+1],2)));
//
//	}
//
//	wxPrintf("\n\n");


//	DFT.output_image.CopyDeviceToHost(false,false);
//	wxPrintf("\n\n");
//
//	if (complex_strided)
//	{
//		last_index = (gpu_image_out.logical_upper_bound_complex_x + 1)*(cpu_image_in.logical_y_dimension-2);
//
//		for (int current_pixel=0; current_pixel < 16; current_pixel++)
//		{
//			wxPrintf("Transformed Vals %d, %3.3e,%3.3e  %3.3e :LAST: %d %3.3e,%3.3e  %3.3e\n",
//					current_pixel, gpu_image_out.real_values[2*current_pixel],gpu_image_out.real_values[2*current_pixel+1],
//					sqrtf(powf(gpu_image_out.real_values[2*current_pixel],2)+powf(gpu_image_out.real_values[2*current_pixel+1],2)),
//					current_pixel + last_index,gpu_image_out.real_values[2*(current_pixel + last_index)],gpu_image_out.real_values[2*(current_pixel + last_index)+1],
//					sqrtf(powf(gpu_image_out.real_values[2*(current_pixel + last_index)],2)+powf(gpu_image_out.real_values[2*(current_pixel + last_index)+1],2)));
//
//		}
//	}
//	else
//	{
//		last_index = (gpu_image_out.logical_upper_bound_complex_x + 1)*(cpu_image_in.logical_y_dimension-2);
//		int actual_pixel;
//		for (int current_pixel=0; current_pixel < 16; current_pixel++)
//		{
//			actual_pixel = current_pixel*(DFT.output_image.dims.w);
//			wxPrintf("Transformed Vals %d, %3.3e,%3.3e  %3.3e :LAST: %d %3.3e,%3.3e  %3.3e\n",
//					current_pixel, gpu_image_out.real_values[actual_pixel],gpu_image_out.real_values[actual_pixel+1],
//					sqrtf(powf(gpu_image_out.real_values[actual_pixel],2)+powf(gpu_image_out.real_values[actual_pixel+1],2)),
//					current_pixel + last_index,gpu_image_out.real_values[2*(current_pixel + last_index)],gpu_image_out.real_values[2*(current_pixel + last_index)+1],
//					sqrtf(powf(gpu_image_out.real_values[2*(current_pixel + last_index)],2)+powf(gpu_image_out.real_values[2*(current_pixel + last_index)+1],2)));
//
//		}
//	}


MyPrintWithDetails("");
	// Complete the second dimension and calc cpu 2d xform to compare
	if (complex_strided && ! do_rotate)
	{
		DFT.FFT_C2C_WithPadding_strided(do_rotate);
//		DFT.DFT_C2C_WithPadding_strided();
	}
	else
	{
//		DFT.DFT_C2C_WithPadding();
//		DFT.FFT_C2C_rotate(true, true);

		DFT.FFT_C2C_WithPadding(do_rotate);
	}

	if (do_rotate)
	{
		MyPrintWithDetails("Doing r2c inverse");

		DFT.FFT_C2C_rotate(true,false);
		DFT.FFT_C2R_rotate(true);
	}

	cpu_image_in.Resize(wanted_output_size_x, wanted_output_size_y, 1, 0.0f);
	cpu_image_in.ForwardFFT(false);
	cpu_image_in.PhaseShift(-(wanted_output_size_x/2-wanted_input_size_x/2), -(wanted_output_size_y/2-wanted_input_size_y/2), 0);
	cpu_image_in.QuickAndDirtyWriteSlice("cpu_shift.mrc", 1, true, 1.0);
	DFT.output_image.CopyDeviceToHost(false,false);
//	gpu_image_out.is_in_real_space = false; //
    gpu_image_out.QuickAndDirtyWriteSlice("DFT_xformed.mrc", 1, false, 1.0);
	wxPrintf("2d cpu\n\n");

//	for (int current_pixel=0; current_pixel < 16; current_pixel++)
//	{
//		wxPrintf("Transformed Vals %d, %3.3e,%3.3e  %3.3e :LAST: %d %3.3e,%3.3e  %3.3e\n",
//				current_pixel, cpu_image_in.real_values[2*current_pixel],cpu_image_in.real_values[2*current_pixel+1],
//				sqrtf(powf(cpu_image_in.real_values[2*current_pixel],2)+powf(cpu_image_in.real_values[2*current_pixel+1],2)),
//				current_pixel + last_index,cpu_image_in.real_values[2*(current_pixel + last_index)],cpu_image_in.real_values[2*(current_pixel + last_index)+1],
//				sqrtf(powf(cpu_image_in.real_values[2*(current_pixel + last_index)],2)+powf(cpu_image_in.real_values[2*(current_pixel + last_index)+1],2)));
//	}
//	wxPrintf("2d GPU \n\n");
//
//	for (int current_pixel=0; current_pixel < 16; current_pixel++)
//	{
//		wxPrintf("Transformed Vals %d, %3.3e,%3.3e  %3.3e :LAST: %d %3.3e,%3.3e  %3.3e\n",
//				current_pixel, gpu_image_out.real_values[2*current_pixel],gpu_image_out.real_values[2*current_pixel+1],
//				sqrtf(powf(gpu_image_out.real_values[2*current_pixel],2)+powf(gpu_image_out.real_values[2*current_pixel+1],2)),
//				current_pixel + last_index,gpu_image_out.real_values[2*(current_pixel + last_index)],gpu_image_out.real_values[2*(current_pixel + last_index)+1],
//				sqrtf(powf(gpu_image_out.real_values[2*(current_pixel + last_index)],2)+powf(gpu_image_out.real_values[2*(current_pixel + last_index)+1],2)));
//	}
//	wxPrintf("\n\n");

	// Timing
	StopWatch timer;
	cpu_image_in.BackwardFFT();
    GpuImage paddedGpu;
    paddedGpu.CopyFromCpuImage(cpu_image_in);
    paddedGpu.CopyHostToDevice();
    int nLoops = 60000;
//	for (int iLoop = 0; iLoop < nLoops; iLoop++)
//	{
//		timer.start("padded");
//		paddedGpu.ForwardFFT(false);
//		paddedGpu.Wait();
//		timer.lap("padded");
//		if (iLoop == 0)     paddedGpu.QuickAndDirtyWriteSlices("FFT_forward.mrc",1,1);
//
//		paddedGpu.BackwardFFT();
//		paddedGpu.AddConstant(rg.GetNormalRandom());
//
//
//	}

	for (int iLoop = 0; iLoop < nLoops; iLoop++)
	{
		timer.start("GPU");
//		if (complex_strided)
//		{
//			DFT.FFT_R2C_WithPadding(do_rotate);
//			DFT.FFT_C2C_WithPadding_strided(do_rotate);
//		}
//		else
//		{
//			DFT.FFT_R2C_WithPadding_strided(do_rotate);
//			DFT.FFT_C2C_WithPadding(do_rotate);
//		}

			DFT.FFT_R2C_WithPadding(do_rotate);
			DFT.FFT_C2C_WithPadding(do_rotate);



		timer.lap("GPU");
	}
	DFT.output_image.Wait();


	timer.print_times();


}

void GpuUtilTest::FFTwithRotation()
{


	unsigned int n_terms = 16;
	unsigned int n_bits_in_final_term = (unsigned int)(std::log2(float(n_terms)));
	const unsigned int n_bits = sizeof(int)*8;
	const unsigned int n_butterflies = 2;
	unsigned int r,x,i,idx;
	// Each (nth) radix 2 butterfly results in the least significant bit going to the
	// (nth) most significant in log2(n_terms) bits. A full FFT will result in a full
	// bit reversal. Intermediate results are calculated as follows:
	for (idx=0; idx< n_terms; idx++)
	{
		r = 0;
		x = idx;

		// Fill in r from the least significant, marching to the most significant. We only need the
		// n_butterflies least significant to be flipped.
	    for(i = 0; i < n_butterflies; i++)
	    {
	    	// I don't think the () are necessary in c++
	        r = r << 1 | (x & 1);
	        x >>= 1;
	    }
	    x=idx;
	    r <<= n_bits_in_final_term - n_butterflies;
	    r += (x >> n_butterflies);
	    wxPrintf("bit reversed from r a %d %d\n", idx, r);
	}




	exit(-1);
	DFTbyDecomposition DFT;
	StopWatch timer;
	bool do_rotation = true;

	int nLoops = 1e4;
	float t_rmsd;

	int wanted_input_size_x = 4096;
	int wanted_input_size_y = 4096;
	int wanted_output_size_x = 1*wanted_input_size_x;
	int wanted_output_size_y = 1*wanted_input_size_y;
	int wanted_number_of_iterations = 1;

	DFT.InitTestCase(wanted_input_size_x,wanted_input_size_y,wanted_output_size_x,wanted_output_size_y);


	Image regular_fft, rotated_fft, rotated_fft_inv, buffer;
	GpuImage d_regular_fft, d_rotated_fft;
	regular_fft.Allocate(wanted_input_size_x,wanted_input_size_y,true);
	rotated_fft.Allocate(wanted_output_size_x,wanted_output_size_y,true);
	rotated_fft_inv.Allocate(wanted_output_size_x,wanted_output_size_y,true);


	RandomNumberGenerator rg(PIf);
	for (int current_pixel=0; current_pixel < regular_fft.real_memory_allocated; current_pixel++)
	{
		regular_fft.real_values[current_pixel] = rg.GetNormalRandomSTD(0.0,1.0);
	}
	rotated_fft.CopyFrom(&regular_fft);
//	rotated_fft_inv.CopyFrom(&regular_fft);
//
//	int cp = 0;
//	for (int iY=0; iY < regular_fft.logical_y_dimension; iY++)
//	{
//		for (int iX=0; iX < regular_fft.logical_x_dimension; iX++)
//		{
//			regular_fft.real_values[cp] = iY* regular_fft.logical_x_dimension + iX;
//			cp++;
//		}
//		regular_fft.real_values[cp] = 0.0;
//		cp++;
//		regular_fft.real_values[cp] = 0.0;
//		cp++;
//
//
//	}
//	rotated_fft.CopyFrom(&regular_fft);


	// Copy of gpu images (on host)
	d_regular_fft.CopyFromCpuImage(regular_fft);
	d_regular_fft.CopyHostToDevice();
//	d_rotated_fft.CopyFromCpuImage(rotated_fft);
//	d_rotated_fft.CopyHostToDevice();

	// Get baseline cpu
	buffer.CopyFrom(&regular_fft);
	buffer.QuickAndDirtyWriteSlice("fft_in.mrc", 1, false, 1.0);
	buffer.ForwardFFT(false);
	buffer.BackwardFFT();
	buffer.MultiplyByConstant(1.0f/buffer.real_memory_allocated);
	buffer.QuickAndDirtyWriteSlice("fft_out_cpu.mrc", 1, false, 1.0);


	// Warm up for regular fft
	d_regular_fft.ForwardFFT(false);
	d_regular_fft.BackwardFFT();
	d_regular_fft.MultiplyByConstant(1.0/d_regular_fft.real_memory_allocated);
	d_regular_fft.CopyDeviceToHost(false, false);
	d_regular_fft.Wait();
	regular_fft.QuickAndDirtyWriteSlice("fft_out_cufft.mrc", 1, false, 1.0);



	buffer.SubtractImage(&regular_fft);
	t_rmsd = sqrtf(buffer.ReturnSumOfSquares(0, 0, 0, 0, false));
	wxPrintf("RMSD between regular gpu fft/ifft pair and cpu fft/ifft pair is %3.3e\n", t_rmsd);
	// reset for comparison to rotated FFT
	buffer.AddImage(&regular_fft);
	MyPrintWithDetails("");

	// Record FFTs timing
	timer.start("cufft");
	for (int iLoop = 0; iLoop < nLoops; iLoop++)
	{
		d_regular_fft.ForwardFFT(false);
		d_regular_fft.BackwardFFT();
	}
	d_regular_fft.Wait();
	timer.lap("cufft");

	DFT.SetGpuImages(rotated_fft,rotated_fft_inv);
	DFT.input_image.QuickAndDirtyWriteSlices("input_out.mrc", 1, 1);
	DFT.output_image.QuickAndDirtyWriteSlices("output_out.mrc", 1, 1);
	DFT.AllocateRotatedBuffer();
	cudaDeviceSynchronize();

//
//	DFT.test_main();
//	exit(-1);
//	// Warm up
	DFT.FFT_R2C_rotate(do_rotation);
	cudaDeviceSynchronize();
	DFT.FFT_C2C_rotate(do_rotation,true);
	cudaDeviceSynchronize();
	DFT.FFT_C2C_rotate(do_rotation,false);
	cudaDeviceSynchronize();
	DFT.FFT_C2R_rotate(do_rotation);
	cudaDeviceSynchronize();
//
	DFT.output_image.MultiplyByConstant(1.0/DFT.output_image.real_memory_allocated);
//	DFT.output_image.MultiplyByConstant(1.0/DFT.output_image.dims.x);


	DFT.output_image.CopyDeviceToHost(false, false);
	DFT.output_image.Wait();
	cudaDeviceSynchronize();

	rotated_fft_inv.QuickAndDirtyWriteSlice("fft_out_rotfft.mrc", 1, false, 1.0);



	buffer.SubtractImage(&rotated_fft_inv);
	t_rmsd = sqrtf(buffer.ReturnSumOfSquares(0, 0, 0, 0, false));
	wxPrintf("RMSD between regular gpu fft/ifft pair and cpu fft/ifft pair is %3.3e\n", t_rmsd);

	// Record FFTs timing
	timer.start("rot_fft");
	for (int iLoop = 0; iLoop < nLoops; iLoop++)
	{
		DFT.FFT_R2C_rotate(do_rotation);
		DFT.FFT_C2C_rotate(do_rotation,true);
		DFT.FFT_C2C_rotate(do_rotation,false);
		DFT.FFT_C2R_rotate(do_rotation);

	}
	DFT.output_image.Wait();
	timer.lap("rot_fft");
	timer.print_times();


}
