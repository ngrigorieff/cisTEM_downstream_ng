#include "samples_functional_testing.hpp"

// TODO //
// TEST 3D's
// ODD Sized IMAGES
// NON SQUARE IMAGES







void SamplesApp::DoInteractiveUserInput()
{
	 UserInput *my_input = new UserInput("Samples", 0.1);
	 test_has_passed = false;
	 delete my_input;
}

bool SamplesApp::DoCalculation()
{

	wxPrintf("\n\n\n     **   ");
	if (OutputIsAtTerminal() == true) wxPrintf(ANSI_UNDERLINE "ProjectX Library Tester" ANSI_UNDERLINE_OFF);
	else wxPrintf("ProjectX Library Tester");
	wxPrintf("   **\n");

	//wxPrintf("")

	WriteEmbeddedFiles();
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This block is taken verbatim from console_test except the use of HomeDir rather than TmpDir to avoid conflicts if different users leave the temp images. TODO add an auto remove for the temp files
void SamplesApp::WriteEmbeddedFiles()
{
	temp_directory = wxFileName::GetHomeDir();
	wxPrintf("\nWriting out embedded test files to '%s'...", temp_directory);
	fflush(stdout);

	hiv_image_80x80x1_filename = temp_directory + "/hiv_image_80x80x1.mrc";
	hiv_images_80x80x10_filename = temp_directory + "/hiv_images_shift_noise_80x80x10.mrc";
	sine_wave_128x128x1_filename = temp_directory + "/sine_wave_128x128x1.mrc";

	WriteEmbeddedArray(hiv_image_80x80x1_filename, hiv_image_80x80x1_array, sizeof(hiv_image_80x80x1_array));
	WriteEmbeddedArray(hiv_images_80x80x10_filename, hiv_images_shift_noise_80x80x10_array, sizeof(hiv_images_shift_noise_80x80x10_array));
	WriteEmbeddedArray(hiv_images_80x80x10_filename, hiv_images_shift_noise_80x80x10_array, sizeof(hiv_images_shift_noise_80x80x10_array));
	WriteEmbeddedArray(sine_wave_128x128x1_filename, sine_128x128x1_array, sizeof(sine_128x128x1_array));

	numeric_text_filename = temp_directory + "/numbers.num";
	WriteNumericTextFile(numeric_text_filename);

	wxPrintf("done!\n");


}

void SamplesApp::WriteEmbeddedArray(const char *filename, const unsigned char *array, long length)
{

	FILE *output_file = NULL;
	output_file = fopen(filename, "wb+");

	if (output_file == NULL)
	{
		wxPrintf(ANSI_COLOR_RED "\n\nError: Can't open output file %s.\n", filename);
		wxPrintf(ANSI_COLOR_RESET "\n\nError: Can't open output file %s.\n", filename);
		DEBUG_ABORT;

	}

	 fwrite (array , sizeof(unsigned char), length, output_file);

	 fclose(output_file);
}

void SamplesApp::WriteNumericTextFile(const char *filename)
{

	FILE *output_file = NULL;
	output_file = fopen(filename, "wb+");

	if (output_file == NULL)
	{
		wxPrintf(ANSI_COLOR_RED "\n\nError: Can't open output file %s.\n", filename);
		wxPrintf(ANSI_COLOR_RESET "\n\nError: Can't open output file %s.\n", filename);
		DEBUG_ABORT;

	}

	fprintf(output_file, "# This is comment, starting with #\n");
	fprintf(output_file, "C This is comment, starting with C\n");
	fprintf(output_file, "%f %f %f %f %f\n%f %f %f %f %f\n", 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.1, 8.3, 9.4, 10.5);
	fprintf(output_file, "# The next line will be blank, but contain 5 spaces\n     \n");
	fprintf(output_file, "%f %f %f %f %f\n", 11.2, 12.7, 13.2, 14.1, 15.8);
	fprintf(output_file, "   # This comment line starts with #, but not at the first character\n");
	fprintf(output_file, "   C This comment line starts with C, but not at the first character\n");
	fprintf(output_file, "C The next line will have varying spaces between the datapoints\n");
	fprintf(output_file, "   %f %f   %f       %f          %f\n", 16.1245, 17.81003, 18.5467, 19.7621, 20.11111);

	fclose(output_file);
}
// end console test blocks
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



