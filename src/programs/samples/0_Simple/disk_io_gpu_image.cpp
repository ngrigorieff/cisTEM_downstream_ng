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
#include "disk_io_gpu_image.hpp"

bool DoDiskIOGpuImageTests(wxString hiv_image_80x80x10_filename, wxString temp_directory) {

  TestResult tr;
  bool passed = true, allPassed = true;
  wxString testName;
  wxPrintf("  Starting GPU I/O image tests.\n\n");

  // Assuming the disk_io_image tests have passed, get the hiv stack to work with
  Image test_image;
  MRCFile input_file;
  input_file.OpenFile(hiv_image_80x80x10_filename.ToStdString(), false); 
  test_image.ReadSlice(&input_file, 1);

  testName = "Allocate GpuImage on GPU"; // If read slice fails, obviously this fails
    GpuImage my_gpu_image;
    int wanted_x_size = 128;
    int wanted_y_size = 128;
    int wanted_z_size = 1;  
    try {
      // Allocate memory for am (int) x,y,z image and specify it will be in realspace to start. 
      my_gpu_image.Allocate(wanted_x_size,wanted_y_size,wanted_z_size,true);
    }
    catch (...) {passed = false;} 
    tr.PrintResults(testName, passed);

  testName = "Set GpuImage to constant";
    float wanted_constant = 1.0f;
    try {
      my_gpu_image.SetToConstant(wanted_constant);
    }
    catch (...) {passed = false;} 
    tr.PrintResults(testName, passed);

  testName = "Return sum";
    float expected_sum;
    try {
      expected_sum = float(wanted_x_size * wanted_y_size * wanted_z_size) * wanted_constant;
      float sum = my_gpu_image.ReturnSumOfRealValues();
      if ( ! FloatsAreAlmostTheSame(expected_sum, sum) ) passed = false;
    }
    catch (...) {passed = false;} 
    tr.PrintResults(testName, passed);

  testName = "Save GPU image";
    std::string quick_and_dirty_filename = temp_directory.ToStdString() + "/gpu_quick.mrc";
    try {
      my_gpu_image.QuickAndDirtyWriteSlices(quick_and_dirty_filename, 1, 1);
    }
    catch (...) {passed = false;} 
    tr.PrintResults(testName, passed); 

  testName = "Read back in GPU image";
    Image my_cpu_image;
    MRCFile re_open = MRCFile(quick_and_dirty_filename,false);
    try {
      my_cpu_image.ReadSlice(&re_open, 1);
      float sum = my_cpu_image.ReturnSumOfRealValues();
      if ( ! FloatsAreAlmostTheSame(expected_sum, sum) ) passed = false;
    }
    catch (...) {passed = false;} 
    tr.PrintResults(testName, passed);

    // A gpu image may be associated directly with a cpu image. In this case, the cpu image is "pinned" to improve data transfer to the gpu. This page-locked state degrades caching performance to a small degree.
    testName = "Create GPU image from CPU image"; // If read slice fails, obviously this fails
    GpuImage my_cpu_associated_gpu_image;
    try {
      // CopyFromCpuImage automatically allocates the memory on the device, pins the memory on the host. It DOES NOT transfer the actual image data itself
      // Alternatively, you could use the constructor
      // GpuImage my_cpu_associated_gpu_image = GPUImage(my_cpu_image);
      my_cpu_associated_gpu_image.CopyFromCpuImage(my_cpu_image);
    }
    catch (...) {passed = false;} 
    tr.PrintResults(testName, passed);  

    testName = "Copy real values to GPU image from CPU image"; // If read slice fails, obviously this fails
    try {

      my_cpu_associated_gpu_image.CopyHostToDevice();
    }
    catch (...) {passed = false;} 
    tr.PrintResults(testName, passed);  

    testName = "Copy back device to host"; // If read slice fails, obviously this fails
    try {
      // Both boolean options default to true, freeing up resources when done.
      bool free_gpu_memory = true;
      bool unpin_host_memory = true;
      my_cpu_associated_gpu_image.CopyDeviceToHost(free_gpu_memory, unpin_host_memory);
    }
    catch (...) {passed = false;} 
    tr.PrintResults(testName, passed);  

    // NOTE!
    // If you create an image on the gpu as we did with my_gpu_image, and want to send it to a host side image (which was implicitly done with my_gpu_image.QuickAndDirtyWriteSLices)
    // You can use an alternate method
    // 	void CopyDeviceToHost(Image &cpu_image, bool should_block_until_complete = false, bool free_gpu_memory = true);

    testName = "Delete MRC file";
      try {
        re_open.CloseFile();
        passed = remove(re_open.filename.mb_str()) == 0;
      } catch (...) {
        passed = false;
      }
      tr.PrintResults(testName, passed);

/* Basic test (delete on completeion of draft code.)
  testName = "TestName"; // If read slice fails, obviously this fails
    // Init variable
    try {
      // Try things

    }
    catch (...) {passed = false;} 
    tr.PrintResults(testName, passed);
*/



  testName = "GPU I/O image all tests";
    allPassed = tr.ReturnAllPassed();
    passed = tr.ReturnAllPassed();
    tr.PrintResults(testName, passed); // This will set tr.ReutrnAllPassed = true, so we need to record the value in this scope in allPassed for return
    return allPassed;
}