#include "../../core/core_headers.h"

// embedded images..

#include "../console_test/hiv_image_80x80x1.cpp"
#include "../console_test/hiv_images_shift_noise_80x80x10.cpp"
#include "../console_test/sine_128x128x1.cpp"

#define PrintResult(result)	PrintResultSlave(result, __LINE__);
#define FailTest {if (test_has_passed == true) PrintResultSlave(false, __LINE__); test_has_passed = false;}//#include "samples_functional_testing.hpp"



// TODO //
// TEST 3D's
// ODD Sized IMAGES
// NON SQUARE IMAGES




class SamplesApp : public MyApp
{
	wxString hiv_image_80x80x1_filename;
	wxString hiv_images_80x80x10_filename;
	wxString sine_wave_128x128x1_filename;
	wxString numeric_text_filename;
	wxString temp_directory;

	public:
		bool DoCalculation();
		void DoInteractiveUserInput();

		bool test_has_passed;


};


IMPLEMENT_APP(SamplesApp);


void SamplesApp::DoInteractiveUserInput()
{
	 UserInput *my_input = new UserInput("Simulator", 0.25);
	 test_has_passed = false;
	 delete my_input;
}

bool SamplesApp::DoCalculation()
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


