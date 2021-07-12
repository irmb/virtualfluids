#include <gmock/gmock.h>
#include <filesystem>

#include <basics/utilities/UbFileInputASCII.h>


TEST(UbFileInputASCIITest, readIntegerAfterString)
{
    // assuming that the config files is stored parallel to this file.
    std::filesystem::path filePath = __FILE__;
    filePath.replace_filename("UbFileInputASCIITest.cfg");

    UbFileInputASCII sut {filePath.string()};

    const int actual = sut.readIntegerAfterString("test =");

    EXPECT_THAT(actual, testing::Eq(1));
}
