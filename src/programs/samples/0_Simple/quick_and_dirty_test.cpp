bool QuickAndDirtyTests(wxString path)
{
    Image quick_image;

    try
    {
        quick_image.QuickAndDirtyReadSlice(std::string(path.mb_str()), 1);
        return true;
    }
    catch (...)
    {
        return false;
    }
}