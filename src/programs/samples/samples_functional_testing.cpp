#include <wx/wx.h>
#include <wx/app.h>
#include <wx/cmdline.h>
#include <cstdio>
#include "wx/socket.h"

#include "shared/shared_functions.cpp"
#include "samples_functional_testing.hpp"
#include "../../core/core_headers.h"
#include "0_Simple/disk_io_image.cpp"


// embedded images..

// #include "../console_test/hiv_image_80x80x1.cpp"
// #include "../console_test/hiv_images_shift_noise_80x80x10.cpp"
// #include "../console_test/sine_128x128x1.cpp"

// #define PrintResult(result)	PrintResultSlave(result, __LINE__);
// #define FailTest {if (test_has_passed == true) PrintResultSlave(false, __LINE__); test_has_passed = false;}//#include "samples_functional_testing.hpp"



// TODO //
// TEST 3D's
// ODD Sized IMAGES
// NON SQUARE IMAGES






IMPLEMENT_APP(SamplesTestingApp);

bool SamplesTestingApp::OnInit()
{

	wxPrintf("Starting samples testing.\n");



	DoDiskIOImageTests();


	wxPrintf("Samples testing done.\n");
	return false;
}


void SamplesTestingApp::DoInteractiveUserInput()
{
	 UserInput *my_input = new UserInput("Simulator", 0.25);
	 test_has_passed = false;
	 delete my_input;
}


bool SamplesTestingApp::DoCalculation()
{

	//wxPrintf("")

	wxPrintf("\n");

	// Do tests..
	wxPrintf("Hello World\n");
#ifdef ENABLEGPU
	wxPrintf("The functional samples have been compiled with GPU support!\n");
#else
	wxPrintf("The functional samples NOT have been compiled with GPU support!\n");
#endif



	wxPrintf("\n\n\n");




	return false;
}




