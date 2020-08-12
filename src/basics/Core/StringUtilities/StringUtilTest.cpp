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