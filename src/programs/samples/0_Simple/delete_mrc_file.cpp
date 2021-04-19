bool DeleteMRCFile(MRCFile &file) {
    file.CloseFile();
     
    return Assert("remove(temp_filename.mb_str())", remove(file.filename.mb_str()), 0, true);
}