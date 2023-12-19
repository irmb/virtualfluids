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
//! \addtogroup gpu_Output Output
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "FilePartCalculator.h"
#include <stdexcept>

uint FilePartCalculator::calculateNumberOfNodesInPart(uint numberOfNodes, uint indexOfFilePart)
{
    uint indexOfLastFilePart = FilePartCalculator::calculateNumberOfParts(numberOfNodes) - 1;

    if (indexOfFilePart > indexOfLastFilePart)
        throw std::runtime_error("The number of nodes for a non-existing part can not be calculated");
    if (indexOfFilePart == indexOfLastFilePart)
        return numberOfNodes - (indexOfFilePart * FilePartCalculator::limitOfNodesForVTK);
    return FilePartCalculator::limitOfNodesForVTK;
}

uint FilePartCalculator::calculateNumberOfParts(uint numberOfNodes)
{
    return numberOfNodes / FilePartCalculator::limitOfNodesForVTK + 1;
}

uint FilePartCalculator::calculateStartingPostionOfPart(uint indexOfPart)
{
    return indexOfPart * FilePartCalculator::limitOfNodesForVTK;
}

const uint FilePartCalculator::limitOfNodesForVTK;

//! \}
