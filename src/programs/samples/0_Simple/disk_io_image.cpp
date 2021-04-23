/*
 * disk_io_image.hpp
 *
 *  Created on: Apr 1, 2021
 *      Author: B.A. Himes, Shiran Dror
 *
 *      Goal:
 *      	The purpose of this recipe, is to demonstrate the different means of reading/writing images to disk using cisTEM.
 *
 *      Scope:
 *      	At the time of writing, all disk i/o must go through the cpu, but hopefully we can do direct reads to the gpu in the future. A pointer to this recipe will be <here>.
 *
 *      Background:
 *      	The Image class in cisTEM can be used to represent 1,2 or 3D images, in memory, with extensions to higher dimensions only available via arrays of Images (or pointers to them.)
 *      	To read an image in from disk, we need to know about its representation, the properties of which can be accessed via the Image_file class. This provides an interface with the
 *      	primary types supported in cisTEM: MRC, TIF(F), DM, EER. We predominantly use the MRC format, which is very similar to CCP4. There are class specializations of mrc_file and mrc_header.
 *
 *
 *
 *
 *
 */

// #ifndef SRC_PROGRAMS_SAMPLES_0_SIMPLE_DISK_IO_IMAGE_HPP_
// #define SRC_PROGRAMS_SAMPLES_0_SIMPLE_DISK_IO_IMAGE_HPP_
// #endif
/* SRC_PROGRAMS_SAMPLES_0_SIMPLE_DISK_IO_IMAGE_HPP_ */

#include "quick_and_dirty_test.cpp"
#include "dimensions_test.cpp"
#include "mrc_file_test.cpp"
#include "pixel_size_test.cpp"
#include "write_mrc_file_to_disk.cpp"
#include "delete_mrc_file.cpp"

#include "do_test.cpp"

bool DoDiskIOImageTests()
{
    bool passed = true;
    // Assuming this is called from samples_functional_testing.cpp, you will have a file written to disk "${HOME}/hiv_image_80x80x1.mrc"
    wxString temp_directory = wxFileName::GetHomeDir();
    wxString hiv_image_80x80x1_filename = temp_directory + "/hiv_image_80x80x1.mrc";
    MRCFile input_file(std::string(hiv_image_80x80x1_filename.mb_str()), false);
    wxString temp_filename = temp_directory + "/tmp1.mrc";
    MRCFile output_file(std::string(temp_filename.mb_str()), false);
    Image test_image;
    test_image.ReadSlice(&input_file, 1);

    wxPrintf(" Starting disk I/O image tests.\n");

    // The easiest way to read or write and image is simply to use the "quick and dirty methods" in Image class
    // TODO instantiate an image object, read it in using quick and dirt read slice.
  wxString testName = "Quick and dirty instantiation";
    passed = DoTest(testName, hiv_image_80x80x1_filename, QuickAndDirtyTests) && passed;
    


    // At times, we may want to have information about the image, prior to reading it into memory. The most general way to do this is to create an Image_file object, which handles all supported types (see Background)
    // TODO create image file type, refer to console_test  MyTestApp::TestMRCFunctions()
    // TODO run size checks as in above
     testName = "Dimensions";
    passed = DoTest(testName, input_file, CheckDimensions) && passed;

    // We can do a little bit more with MRC files, particularly modifying the header information. Here, we'll set the pixel size to 2
    // TODO Same as image file, create, read, add checks as in MyTestApp::TestMRCFunctions()
    testName = "MRC file";
    passed = DoTest(testName, test_image, TestMRCFile) && passed;

    // TODO Additionally add check on pixel size (should == 1) but use the approximate check that is in  MyTestApp::TestMRCFunctions() for pixel values
    // TODO Modify pixel size, set to 2
    testName = "Pixel";
    passed = DoTest(testName, input_file, TestPixelSizeAndChange) && passed;

    // We'll skip the quick and dirty write slices, and now write a temporary file with our modified pixel size.
    // TODO use mrc_file method write slice to disk, to write tmp1.mrc (in the temp directory above, you'll need a new wxString too

    // Alternatively, you can pass a pointer to your mrc_file object to the image object method WriteSlices
    // TODO write out tmp2.mrc using the Image class method
    testName = "Write temp file";
    passed = DoTest(testName, input_file, output_file, TestWriteMRCFile) && passed;

    // TODO read in both tmp images, and check that the pixel size has been changed appropriately (with FAIL etc.)
    testName = "Change and compare pixel sizes";
    passed = DoTest(testName, input_file, output_file, TestComparePixelSize) && passed;

    // TODO remove the tmp images from disk
    testName = "Delete MRC file";
    passed = DoTest(testName, output_file, DeleteMRCFile) && passed;

    // TODO call end test and ensure the printout indicates this test (disk_io_image) has pass/failed.
    testName = "Disk I/O image";
    passed ? PrintSuccess(testName.mb_str()) : PrintFailure(testName.mb_str());
    return passed;
}
