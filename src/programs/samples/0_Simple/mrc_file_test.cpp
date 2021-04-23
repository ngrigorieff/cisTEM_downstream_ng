bool TestMRCFile(Image &test_image)
{
	bool expectedIsInRealSpace = false;
	int expectedLogicalXDimension = 80;
	int expectedLogicalYDimension = 80;
	int expectedLogicalZDimension = 1;

	bool passed = true;




	// check dimensions and type

	passed = Assert("test_image.is_in_real_space", test_image.is_in_real_space, expectedIsInRealSpace, false) && passed;
	passed = Assert("test_image.logical_x_dimension", test_image.logical_x_dimension, expectedLogicalXDimension, true) && passed;
	passed = Assert("test_image.logical_y_dimension", test_image.logical_y_dimension, expectedLogicalYDimension, true) && passed;
	passed = Assert("test_image.logical_z_dimension", test_image.logical_z_dimension, expectedLogicalZDimension, true) && passed;

	// check first and last pixel...

	//wxPrintf(std::to_string(test_image.real_values[0]).c_str());
	passed = Assert("DoublesAreAlmostTheSame(test_image.real_values[0], -0.340068)",
					DoublesAreAlmostTheSame(test_image.real_values[0], -0.340068),
					false, false) &&
			 passed;

	passed = Assert("DoublesAreAlmostTheSame(test_image.real_values[test_image.real_memory_allocated - 3], 0.637069)",
					DoublesAreAlmostTheSame(test_image.real_values[test_image.real_memory_allocated - 3], 0.637069),
					false, false) &&
			 passed;

	return passed;
}
