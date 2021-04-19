
bool DoTest(wxString testName, wxString str, std::function<bool(wxString &)> testFunction)
{
    wxString indendation = "  ";
    testName = indendation + testName + " test";
    bool result = testFunction(str);
    result ? PrintSuccess(testName.mb_str()) : PrintFailure(testName.mb_str());
    return result;
}
bool DoTest(wxString testName, Image &test_image, std::function<bool(Image &)> testFunction)
{
    wxString indendation = "  ";
    testName = indendation + testName + " test";
    bool result = testFunction(test_image);
    result ? PrintSuccess(testName.mb_str()) : PrintFailure(testName.mb_str());
    return result;
}
bool DoTest(wxString testName, MRCFile &input_file, std::function<bool(MRCFile &)> testFunction)
{
    wxString indendation = "  ";
    testName = indendation + testName + " test";
    bool result = testFunction(input_file);
    result ? PrintSuccess(testName.mb_str()) : PrintFailure(testName.mb_str());
    return result;
}
bool DoTest(wxString testName, MRCFile &input_file, MRCFile &output_file, std::function<bool(MRCFile &, MRCFile &)> testFunction)
{
    wxString indendation = "  ";
    testName = indendation + testName + " test";
    bool result = testFunction(input_file, output_file);
    result ? PrintSuccess(testName.mb_str()) : PrintFailure(testName.mb_str());
    return result;
}