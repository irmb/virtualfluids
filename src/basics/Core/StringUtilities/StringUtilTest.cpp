//
// Created by Soeren Peters on 24.07.20.
//

#include <gmock/gmock.h>

#include "StringUtil.h"

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(StringUtilTest, endsWith_shouldReturnTrue)
{
    const std::string input {"input_string"};
    const std::string ends_with {"string"};

    ASSERT_TRUE(StringUtil::endsWith(input, ends_with));
}

TEST(StringUtilTest, endsWith_shouldReturnFalse)
{
    const std::string input {"input_string"};
    const std::string ends_with {"string_"};

    ASSERT_FALSE(StringUtil::endsWith(input, ends_with));
}

TEST(StringUtilTest, toIntVector)
{
    const std::string input {"1   2\n3 4"};
    std::vector<int> expected_result {1, 2, 3, 4};

    auto result = StringUtil::toIntVector(input);

    ASSERT_THAT(result,testing::Eq(expected_result));
}

TEST(StringUtilTest, splitIntoStringsWithDelimeter)
{
    const std::string input {"1   2\n3 4\t5"};
    const std::string delimeter {" \n\t"};
    std::vector<std::string> expected_result {"1", "2", "3", "4", "5"};

    auto result = StringUtil::split(input, delimeter);

    ASSERT_THAT(result,testing::Eq(expected_result));
}

TEST(StringUtilTest, toStringVector)
{
    const std::string input {"1   2\n3 4\t5"};
    std::vector<std::string> expected_result {"1", "2", "3", "4", "5"};

    auto result = StringUtil::toStringVector(input);

    ASSERT_THAT(result,testing::Eq(expected_result));
}