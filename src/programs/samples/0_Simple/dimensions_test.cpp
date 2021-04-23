

bool CheckDimensions(MRCFile &input_file)
{
	int expectedReturnNumberOfSlices = 1, expectedReturnXSize = 80, expectedReturnYSize = 80, expectedReturnZSize = 1;
	bool passed = true;

	passed = Assert("input_file.ReturnNumberOfSlices()", input_file.ReturnNumberOfSlices(), expectedReturnNumberOfSlices, true) && passed;
	passed = Assert("input_file.ReturnXSize()", input_file.ReturnXSize(), expectedReturnXSize, true) && passed;
	passed = Assert("input_file.ReturnYSize()", input_file.ReturnYSize(), expectedReturnYSize, true) && passed;
	passed = Assert("input_file.ReturnZSize()", input_file.ReturnZSize(), expectedReturnZSize, true) && passed;

	return passed;
}