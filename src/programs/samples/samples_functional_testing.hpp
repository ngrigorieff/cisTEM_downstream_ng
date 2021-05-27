#ifndef SAMPLES_FUNCTION_TESTING_HPP_
#define SAMPLES_FUNCTION_TESTING_HPP_

#include "../../core/core_headers.h"
#include <cstdio>
#include <list>



#include "classes/TestFile.cpp"
#include "classes/EmbeddedTestFile.cpp"
#include "classes/NumericTestFile.cpp"
#include "classes/TestResult.hpp"

#include "0_Simple/disk_io_image.cpp"


#include "../console_test/hiv_image_80x80x1.cpp"
#include "../console_test/hiv_images_shift_noise_80x80x10.cpp"
#include "../console_test/sine_128x128x1.cpp"


class SamplesTestingApp : public wxAppConsole
{

	// constructor: set file names and temp folder, write embbeded files to harddrive.
  wxString temp_directory = wxFileName::GetHomeDir();
  wxString hiv_image_80x80x1_filename = 	temp_directory + "/hiv_image_80x80x1.mrc";
  wxString hiv_images_80x80x10_filename = 	temp_directory + "/hiv_images_shift_noise_80x80x10.mrc";  
  wxString sine_wave_128x128x1_filename = 	temp_directory + "/sine_wave_128x128x1.mrc";  
  wxString numeric_text_filename = 			temp_directory + "/numbers.num";

	std::list<TestFile*> testFiles;
	
	public:
		SamplesTestingApp();
		~SamplesTestingApp();
    	bool OnInit();
		void WriteFiles();
		//void WriteEmbeddedArray(const char *filename, const unsigned char *array, long length);
		//void WriteNumericTextFile(const char *filename);
   
		bool DoCalculation();
		void DoInteractiveUserInput();

};


#endif
