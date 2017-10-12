#include "gmock/gmock.h"

#include "StringUtil.h"

TEST(StringUtilTest, findAndReplace)
{
	std::string source0 = "2232fh";
	std::string source1 = "find_2232fh";
	std::string source2 = "2q232_find_2232fh";
	std::string source3 = "2q232_find";
	std::string source4 = "find_2q232_find";
	std::string source5 = "find_2q2find2_find";

	std::string find = "find";
	std::string replace = "replace";

    EXPECT_THAT(StringUtil::findAndReplace(source0, find, replace), "2232fh");
    EXPECT_THAT(StringUtil::findAndReplace(source1, find, replace), "replace_2232fh");
    EXPECT_THAT(StringUtil::findAndReplace(source2, find, replace), "2q232_replace_2232fh");
    EXPECT_THAT(StringUtil::findAndReplace(source3, find, replace), "2q232_replace");
    EXPECT_THAT(StringUtil::findAndReplace(source4, find, replace), "replace_2q232_replace");
    EXPECT_THAT(StringUtil::findAndReplace(source5, find, replace), "replace_2q2replace2_replace");
}

TEST(StringUtilTest, makeUpper)
{
	std::string source = "hello123world";
	EXPECT_THAT(StringUtil::makeUpper(source), "HELLO123WORLD");
}

TEST(StringUtilTest, makeLower)
{
	std::string source = "HELLO123WORLD";
	EXPECT_THAT(StringUtil::makeLower(source), "hello123world");
}

TEST(StringUtilTest, contains_shouldReturnTrue)
{
    std::string source = "fsd8998--++search++d";
    const char* find = "search";
    EXPECT_TRUE(StringUtil::contains(source, find));
}

TEST(StringUtilTest, contains_shouldReturnFalse)
{
    std::string source = "fsd8998--++++d";
    const char* find = "search";
    EXPECT_FALSE(StringUtil::contains(source, find));
}

TEST(StringUtilTest, pad)
{
    std::string source = "start";
    char pad = '1';
    int length = 6;
    EXPECT_THAT(StringUtil::pad(source, pad, length), "start1");
}

TEST(StringUtilTest, trimAStringDeletesTrimAtTheStartAndTheEnd)
{
    std::string source = " for after ";
    std::string trim = " ";
    EXPECT_THAT(StringUtil::trim(source, trim), "for after");
}

TEST(StringUtilTest, convertStringToInt)
{
    std::string source = " 1 ";
    EXPECT_THAT(StringUtil::toInt(source), 1);
}

TEST(StringUtilTest, convertStringToFloat)
{
    std::string source = " 1.2 ";
    EXPECT_THAT(StringUtil::toFloat(source), testing::FloatEq(1.2f));
}

TEST(StringUtilTest, convertStringToDouble)
{
    std::string source = " 1.2 ";
    EXPECT_THAT(StringUtil::toDouble(source), testing::DoubleEq(1.2));
}

TEST(StringUtilTest, convertStringToBool)
{
    std::string source0 = "false";
    EXPECT_FALSE(StringUtil::toBool(source0));

    std::string source = "true";
    EXPECT_TRUE(StringUtil::toBool(source));
}

TEST(StringUtilTest, toVector)
{
    std::string source = " 1 2 3 ";
    auto values = StringUtil::toVector(source);

    EXPECT_THAT(values[0], 1);
    EXPECT_THAT(values[1], 2);
    EXPECT_THAT(values[2], 3);

    EXPECT_THAT(values.size(), 3);
}

TEST(StringUtilTest, toString_implementedWithInt)
{
    int i = 123;
    EXPECT_THAT(StringUtil::toString(i), "123");
}

TEST(StringUtilTest, string_ends_with)
{
    std::string myString = "MyString";
    EXPECT_TRUE(StringUtil::endsWith(myString, "String"));
    EXPECT_FALSE(StringUtil::endsWith(myString, "FOO"));
}
