bool TestWriteMRCFile(MRCFile &input_file, MRCFile &output_file)
{
    Image test_image;

    try
    {
        test_image.ReadSlice(&input_file, 1);
        test_image.WriteSlice(&output_file, 1);
        return true;
    }
    
    catch (...)
    {
        return false;
    }
}