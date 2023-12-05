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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \author Martin Schoenherr
//=======================================================================================
#include "FilePartCalculator.h"
#include "gpu/core/Utilities/testUtilitiesGPU.h"
#include <gmock/gmock.h>

TEST(FilePartCalculatorTest, calculateNumberOfParts)
{
    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(FilePartCalculator::limitOfNodesForVTK * 2), testing::Eq(3));

    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(FilePartCalculator::limitOfNodesForVTK), testing::Eq(2));

    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(FilePartCalculator::limitOfNodesForVTK - 1), testing::Eq(1));

    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(1), testing::Eq(1));
}

TEST(FilePartCalculatorTest, calculateNumberOfNodesInPart)
{
    uint part = 0;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfNodesInPart(FilePartCalculator::limitOfNodesForVTK + 13, part),
                FilePartCalculator::limitOfNodesForVTK);

    part = 1;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfNodesInPart(FilePartCalculator::limitOfNodesForVTK + 13, part), 13);

    part = 2;
    EXPECT_THROW(FilePartCalculator::calculateNumberOfNodesInPart(FilePartCalculator::limitOfNodesForVTK + 13, part),
                 std::runtime_error);
}

TEST(FilePartCalculatorTest, getStartingPostionOfPart)
{
    uint part = 0;
    EXPECT_THAT(FilePartCalculator::calculateStartingPostionOfPart(part), 0);

    part = 1;
    EXPECT_THAT(FilePartCalculator::calculateStartingPostionOfPart(part), FilePartCalculator::limitOfNodesForVTK);
}
