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
#include "WriterUtilities.h"
#include "Parameter/Parameter.h"
#include <basics/StringUtilities/StringUtil.h>

std::string WriterUtilities::makePartFileNameEnding(uint level, int ID, int part, int timestep)
{
    return "_lev_" + StringUtil::toString<int>(level) + "_ID_" + StringUtil::toString<int>(ID) + "_Part_" +
           StringUtil::toString<int>(part) + "_t_" + StringUtil::toString<int>(timestep) + ".vtk";
}

bool WriterUtilities::isPeriodicCell(const LBMSimulationParameter& parH, unsigned int baseNodeOfCell,
                                     unsigned int otherNodeInCell)
{
    // perform periodicity check by calculating the length of the grid cell's space diagonal
    const real distance = sqrt(
        pow(parH.coordinateX[otherNodeInCell] - parH.coordinateX[baseNodeOfCell], 2.) +
        pow(parH.coordinateY[otherNodeInCell] - parH.coordinateY[baseNodeOfCell], 2.) +
        pow(parH.coordinateZ[otherNodeInCell] - parH.coordinateZ[baseNodeOfCell], 2.));
    return distance > 1.01 * sqrt(3 * pow(parH.gridSpacing, 2.));
}

void WriterUtilities::getIndicesOfAllNodesInOct(std::array<uint, 8>& nodeIndices, uint baseNodeOfOct,
                                                const LBMSimulationParameter& parH)
{
    nodeIndices[0] = baseNodeOfOct;
    nodeIndices[1] = parH.neighborX[nodeIndices[0]];
    nodeIndices[2] = parH.neighborY[nodeIndices[1]];
    nodeIndices[3] = parH.neighborY[nodeIndices[0]];
    nodeIndices[4] = parH.neighborZ[nodeIndices[0]];
    nodeIndices[5] = parH.neighborZ[nodeIndices[1]];
    nodeIndices[6] = parH.neighborZ[nodeIndices[2]];
    nodeIndices[7] = parH.neighborZ[nodeIndices[3]];
}

void WriterUtilities::calculateRelativeNodeIndexInPart(std::array<uint, 8>& relativePositionInPart,
                                                       const std::array<uint, 8>& indicesOfOct, uint startPositionOfPart)
{
    for (size_t i = 0; i < relativePositionInPart.size(); i++) {
        relativePositionInPart[i] = indicesOfOct[i] - startPositionOfPart;
    }
}

bool WriterUtilities::areAllNodesInOctValidForWriting(const std::array<uint, 8>& indicesOfOct,
                                                      const LBMSimulationParameter& parH, uint endPositionOfPart)
{
    const bool neighborsAreFluid = std::all_of(indicesOfOct.begin(), indicesOfOct.end(),
                                               [&](uint index) { return parH.typeOfGridNode[index] == GEO_FLUID; });

    const bool neighborIsOutOfPart =
        (std::any_of(indicesOfOct.begin(), indicesOfOct.end(), [&](uint index) { return index > endPositionOfPart; }));

    return neighborsAreFluid && !neighborIsOutOfPart;
}
