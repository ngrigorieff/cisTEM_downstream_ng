bool TestPixelSizeAndChange(MRCFile &input_file)
{
	input_file.SetPixelSize(1.0);

	return Assert("DoublesAreAlmostTheSame(input_file.ReturnPixelSize(), 1.0)", DoublesAreAlmostTheSame(input_file.ReturnPixelSize(), 1.0), true, true);
}

bool TestComparePixelSize(MRCFile &input_file, MRCFile &output_file) {

	output_file.SetPixelSize(2.0);

	return Assert("input_file.ReturnPixelSize()", input_file.ReturnPixelSize(), 1.0, true) &&
    		Assert("output_file.ReturnPixelSize()", output_file.ReturnPixelSize(), 2.0, true);
}