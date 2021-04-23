void PrintSuccess(const char *testName);
void PrintFailure(const char *testName);
std::string BoolToString(bool var);
bool Assert(const char *functionName, int functionOutput, int expectedResult, bool shouldBeEqual);
bool Assert(const char *functionName, float functionOutput, float expectedResult, bool shouldBeEqual);
bool Assert(const char *functionName, bool functionOutput, bool expectedResult, bool shouldBeEqual);
bool AssertInternal(bool result, bool shouldBeEqual, const char *functionName, std::string functionOutput, std::string expectedResult);
void AddErrorMessage(const char *functionName, std::string functionOutput, std::string expectedResult, bool shouldBeEqual);