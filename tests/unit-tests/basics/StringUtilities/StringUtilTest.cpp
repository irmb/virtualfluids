//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup StringUtilities_tests StringUtilities
//! \ingroup basics_tests
//! \{
//! \author Soeren Peters
//=======================================================================================
#include <gmock/gmock.h>

#include <basics/StringUtilities/StringUtil.h>

TEST(StringUtilTest, endsWith_shouldReturnTrue)
{
    const std::string input{ "input_string" };
    const std::string ends_with{ "string" };

    ASSERT_TRUE(StringUtil::endsWith(input, ends_with));
}

TEST(StringUtilTest, endsWith_shouldReturnFalse)
{
    const std::string input{ "input_string" };
    const std::string ends_with{ "string_" };

    ASSERT_FALSE(StringUtil::endsWith(input, ends_with));
}

TEST(StringUtilTest, toIntVector)
{
    const std::string input{ "1   2\n3 4" };
    std::vector<int> expected_result{ 1, 2, 3, 4 };

    auto result = StringUtil::toIntVector(input);

    ASSERT_THAT(result, testing::Eq(expected_result));
}

TEST(StringUtilTest, splitIntoStringsWithDelimeter)
{
    const std::string input{ "1   2\n3 4\t5" };
    const std::string delimeter{ " \n\t" };
    std::vector<std::string> expected_result{ "1", "2", "3", "4", "5" };

    auto result = StringUtil::split(input, delimeter);

    ASSERT_THAT(result, testing::Eq(expected_result));
}

TEST(StringUtilTest, toStringVector)
{
    const std::string input{ "1   2\n3 4\t5" };
    std::vector<std::string> expected_result{ "1", "2", "3", "4", "5" };

    auto result = StringUtil::toStringVector(input);

    ASSERT_THAT(result, testing::Eq(expected_result));
}
//! \}
