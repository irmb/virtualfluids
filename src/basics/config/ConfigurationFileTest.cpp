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
//! \addtogroup config_tests config
//! \ingroup basics_tests
//! \{
//! \author Soeren Peters
//=======================================================================================
#include <gmock/gmock.h>

#include "ConfigurationFile.h"

using namespace vf::basics;

TEST(ConfigurationFileTest, ContainsReturnsTrueForExistingKey)
{
    ConfigurationFile config;
    config.data["key1"] = "value1";
    config.data["key2"] = "value2";

    EXPECT_TRUE(config.contains("key1"));
    EXPECT_TRUE(config.contains("key2"));
}

TEST(ConfigurationFileTest, ContainsReturnsFalseForMissingKey)
{
    ConfigurationFile config;
    config.data["key1"] = "value1";

    EXPECT_FALSE(config.contains("key2"));
}

TEST(ConfigurationFileTest, GetValueReturnsCorrectValue)
{
    ConfigurationFile config;
    config.data["key1"] = "value1";
    config.data["key2"] = "1234";

    EXPECT_EQ(config.getValue<std::string>("key1"), "value1");
    EXPECT_EQ(config.getValue<int>("key2"), 1234);
}

TEST(ConfigurationFileTest, GetValueThrowsExceptionForMissingKey)
{
    ConfigurationFile config;
    config.data["key1"] = "value1";

    EXPECT_THROW(config.getValue<std::string>("key2"), UbException);
}

TEST(ConfigurationFileTest, GetVectorReturnsCorrectValues)
{
    ConfigurationFile config;
    config.data["key1"] = "1, 2, 3";
    config.data["key2"] = "4; 5; 6";

    std::vector<int> v1 = config.getVector<int>("key1");
    std::vector<int> v2 = config.getVector<int>("key2");

    EXPECT_EQ(v1.size(), 3);
    EXPECT_EQ(v1[0], 1);
    EXPECT_EQ(v1[1], 2);
    EXPECT_EQ(v1[2], 3);

    EXPECT_EQ(v2.size(), 3);
    EXPECT_EQ(v2[0], 4);
    EXPECT_EQ(v2[1], 5);
    EXPECT_EQ(v2[2], 6);
}

TEST(ConfigurationFileTest, GetVectorThrowsExceptionForMissingKey)
{
    ConfigurationFile config;
    config.data["key1"] = "1, 2, 3";

    EXPECT_THROW(config.getVector<int>("key2"), UbException);
}

TEST(ConfigurationFileTest, GetValueReturnsDefaultValueForMissingKey)
{
    ConfigurationFile config;
    config.data["key1"] = "value1";

    int defaultValue = 42;
    int value = config.getValue<int>("key2", defaultValue);

    EXPECT_EQ(value, defaultValue);
}

TEST(ConfigurationFileTest, GetValueReturnsCorrectValueForExistingKey)
{
    ConfigurationFile config;
    config.data["key1"] = "42";

    int defaultValue = 0;
    int value = config.getValue<int>("key1", defaultValue);

    EXPECT_EQ(value, 42);
}

TEST(ConfigurationFileTest, FromStringConvertsValueToBool)
{
    ConfigurationFile config;
    config.data["key1"] = "true";
    config.data["key2"] = "false";

    bool valueTrue = config.getValue<bool>("key1");
    bool valueFalse = config.getValue<bool>("key2");

    EXPECT_TRUE(valueTrue);
    EXPECT_FALSE(valueFalse);
}

TEST(ConfigurationFileTest, GetValueThrowsExceptionForWrongTypeConversion)
{
    ConfigurationFile config;
    config.data["key1"] = "text";

    EXPECT_THROW(config.getVector<int>("key1"), UbException);
}

TEST(ConfigurationFileTest, ClearRemovesAllData)
{
    ConfigurationFile config;
    config.data["key1"] = "value1";
    config.data["key2"] = "value2";

    config.clear();

    EXPECT_TRUE(config.data.empty());
}

//! \}
