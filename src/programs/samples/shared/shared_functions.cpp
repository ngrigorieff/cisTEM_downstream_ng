#include "shared_functions.hpp"

std::string BoolToString(bool var)
{
    return var ? "true" : "false";
}


bool Assert(const char *functionName, int functionOutput, int expectedResult, bool shouldBeEqual = false)
{
    return AssertInternal(functionOutput == expectedResult, shouldBeEqual, functionName, std::to_string(functionOutput), std::to_string(expectedResult));
}
bool Assert(const char *functionName, float functionOutput, float expectedResult, bool shouldBeEqual = false)
{
    return AssertInternal(functionOutput == expectedResult, shouldBeEqual, functionName, std::to_string(functionOutput), std::to_string(expectedResult));
}
bool Assert(const char *functionName, bool functionOutput, bool expectedResult, bool shouldBeEqual = false)
{

    return AssertInternal(functionOutput == expectedResult, shouldBeEqual, functionName, BoolToString(functionOutput), BoolToString(expectedResult));
}

bool AssertInternal(bool result, bool shouldBeEqual, const char *functionName, std::string functionOutput, std::string expectedResult)
{
    if (result != shouldBeEqual)
    {
        AddErrorMessage(functionName, functionOutput, expectedResult, shouldBeEqual);
        return false;
    }
    return true;
}

void AddErrorMessage(const char *functionName, std::string functionOutput, std::string expectedResult, bool shouldBeEqual = false)
{
    std::string errorMessage = "";
    std::string indentation = "   ";
    errorMessage += indentation + "Error in " + functionName + ". ";
    errorMessage += "Value is: " + functionOutput;
    errorMessage += " And should have been ";
    errorMessage += shouldBeEqual ? "equal to " : "different than ";
    errorMessage += expectedResult;
    errorMessage += "\n";
    wxPrintf(errorMessage.c_str());
}

void PrintSuccess(const char *testName)
{
    wxPrintf(testName);
    wxPrintf(": [Success]\n");
}
void PrintFailure(const char *testName)
{
    wxPrintf(testName);
    wxPrintf(": [Failed]\n");
}