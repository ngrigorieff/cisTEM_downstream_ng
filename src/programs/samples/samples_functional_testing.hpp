#include "../../core/core_headers.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This block is just from console_test.cpp - Evaluate what is actually needed after image selections have been made. TODO
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// embedded images..

#include "../console_test/hiv_image_80x80x1.cpp"
#include "../console_test/hiv_images_shift_noise_80x80x10.cpp"
#include "../console_test/sine_128x128x1.cpp"

#define PrintResult(result)	PrintResultSlave(result, __LINE__);
#define FailTest {if (test_has_passed == true) PrintResultSlave(false, __LINE__); test_has_passed = false;}//#include "samples_functional_testing.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "0_Simple/disk_io_image.hpp"

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
		// This block taken verbatim from console_test except the use of HomeDir rather than TmpDir to avoid conflicts if different users leave the temp images. TODO add an auto remove for the temp files
		void WriteEmbeddedFiles();
		void WriteEmbeddedArray(const char *filename, const unsigned char *array, long length);
		void WriteNumericTextFile(const char *filename);
		// end verbatim


		bool test_has_passed;


};


IMPLEMENT_APP(SamplesApp);
