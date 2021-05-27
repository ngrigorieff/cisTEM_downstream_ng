/*
 * disk_io_image.hpp
 *
 *  Created on: May 27, 2021
 *      Author: B.A. Himes, Shiran Dror
 *
 *      Goal:
 *      	The purpose of this recipe, is to demonstrate the different
 * means of reading/writing gpu images to disk using cisTEM.
 * 
 *      Prerequisites: 
 *          As of the writing, gpu io must go through the pcie bus via the cpu, so disk_io_image.cpp is a prereq.
 *
 *      Scope:
 *      	At the time of writing, all disk i/o must go through the cpu,
 * but hopefully we can do direct reads to the gpu in the future. A pointer to
 * this recipe will be <here>.
 *
 *      Background:
 *      	The GpuImage class in cisTEM attempts to parallel the Image class closely, sharing many of the same methods.
 * Perhaps the most significant difference is the handling of the separate memory spaces on the GPU and the CPU.
 * This recipe will demonstrate the different methods in the GpuImage class available to simplify memory management on the GPU.
 * 
 * 
 * 
 *
 *
 */


bool DoDiskIOGpuImageTests(wxString hiv_image_80x80x10_filename, wxString temp_directory) {

  TestResult tr;
  bool passed = true, allPassed = true;
  wxPrintf("  Starting GPU I/O image tests.\n\n");

  // Assuming the disk_io_image tests have passed, get the hiv stack to work with
  Image test_image;
  MRCFile input_file;
  input_file.OpenFile(hiv_image_80x80x10_filename.ToStdString(), false); 
  test_image.ReadSlice(&input_file, 1);

/* Basic test (delete on completeion of draft code.)
  testName = "TestName"; // If read slice fails, obviously this fails
    // Init variable
    try {
      // Try things

    }
    catch (...) {passed = false;} 
    tr.PrintResults(testName, passed);
*/



  testName = "  GPU I/O image all tests";
    allPassed = tr.ReturnAllPassed();
    passed = tr.ReturnAllPassed();
    tr.PrintResults(testName, passed); // This will set tr.ReutrnAllPassed = true, so we need to record the value in this scope in allPassed for return
    return allPassed;
}