#include "samples_functional_testing.hpp"



IMPLEMENT_APP(SamplesTestingApp);

bool SamplesTestingApp::OnInit() {

  wxPrintf("Starting samples testing.\n\n");

  DoDiskIOImageTests(hiv_images_80x80x10_filename, temp_directory);

  wxPrintf("\n\nSamples testing done.\n");
  return false;
}



void SamplesTestingApp::WriteFiles() {

  wxPrintf("\nWriting out embedded test files to '%s'...\n", temp_directory);
  fflush(stdout);

  testFiles.push_back(new EmbeddedTestFile(hiv_image_80x80x1_filename, hiv_image_80x80x1_array, sizeof(hiv_image_80x80x1_array)));
  testFiles.push_back(new EmbeddedTestFile(hiv_images_80x80x10_filename, hiv_images_shift_noise_80x80x10_array, sizeof(hiv_images_shift_noise_80x80x10_array)));
  testFiles.push_back(new EmbeddedTestFile(sine_wave_128x128x1_filename, sine_128x128x1_array, sizeof(sine_128x128x1_array)));
  testFiles.push_back(new NumericTestFile(numeric_text_filename));

  wxPrintf("done writing files!\n");
}


SamplesTestingApp::SamplesTestingApp() {

  WriteFiles();
}

SamplesTestingApp::~SamplesTestingApp() {
	// destructor: remove all files written to harddrive.
	wxPrintf("Removing test files from '%s'... \n", temp_directory);

  for(auto &it:testFiles) delete it;
  testFiles.clear();
  // removeFile(hiv_image_80x80x1_filename.mb_str());
  // removeFile(hiv_images_80x80x10_filename.mb_str());
  // removeFile(sine_wave_128x128x1_filename.mb_str());
  // removeFile(numeric_text_filename.mb_str());
  wxPrintf("done!\n");
}