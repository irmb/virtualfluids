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
//! \addtogroup gpu_utilitie_tests utilities
//! \ingroup gpu_GridGenerator_tests GridGenerator
//! \{
#include <gmock/gmock.h>

#include "communication.h"

using namespace CommunicationDirections;

TEST(communicationTest, isNegative)
{
    EXPECT_TRUE(isNegative(CommunicationDirection::MX));
    EXPECT_TRUE(isNegative(CommunicationDirection::MY));
    EXPECT_TRUE(isNegative(CommunicationDirection::MZ));
    EXPECT_FALSE(isNegative(CommunicationDirection::PX));
    EXPECT_FALSE(isNegative(CommunicationDirection::PY));
    EXPECT_FALSE(isNegative(CommunicationDirection::PZ));
}

TEST(communicationTest, isPositive)
{
    EXPECT_TRUE(isPositive(CommunicationDirection::PX));
    EXPECT_TRUE(isPositive(CommunicationDirection::PY));
    EXPECT_TRUE(isPositive(CommunicationDirection::PZ));
    EXPECT_FALSE(isPositive(CommunicationDirection::MX));
    EXPECT_FALSE(isPositive(CommunicationDirection::MY));
    EXPECT_FALSE(isPositive(CommunicationDirection::MZ));
}

TEST(communicationTest, getNegativeDirectionAlongAxis)
{
    EXPECT_THAT(getNegativeDirectionAlongAxis(Axis::x), CommunicationDirection::MX);
    EXPECT_THAT(getNegativeDirectionAlongAxis(Axis::y), CommunicationDirection::MY);
    EXPECT_THAT(getNegativeDirectionAlongAxis(Axis::z), CommunicationDirection::MZ);
}

TEST(communicationTest, getPositiveDirectionAlongAxis)
{
    EXPECT_THAT(getPositiveDirectionAlongAxis(Axis::x), CommunicationDirection::PX);
    EXPECT_THAT(getPositiveDirectionAlongAxis(Axis::y), CommunicationDirection::PY);
    EXPECT_THAT(getPositiveDirectionAlongAxis(Axis::z), CommunicationDirection::PZ);
}
//! \}
