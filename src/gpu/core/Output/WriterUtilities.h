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
#ifndef WRITER_UTILITIES
#define WRITER_UTILITIES

#include <array>

#include <basics/DataTypes.h>

class Parameter;
struct LBMSimulationParameter;

class WriterUtilities
{
public:
    //! \brief check whether a grid cell is part of a periodic boundary condition
    //! \param baseNodeOfCell is the index of one node of the grid cell
    //! \param otherNodeInCell is the index of the node which is not on the same face of the grid cell as the base node (i.e.
    //! it is on the other end of the space diagonal)
    static bool isPeriodicCell(const LBMSimulationParameter& parH, unsigned int baseNodeOfCell,
                               unsigned int otherNodeInCell);

    //! \brief use the neighbor relations to find the indices of all nodes in an oct cell
    static void getIndicesOfAllNodesInOct(std::array<uint, 8>& nodeIndices, uint baseNodeOfOct,
                                          const LBMSimulationParameter& parH);

    //! \brief calculate the node index relative to the start position of the part
    static void calculateRelativeNodeIndexInPart(std::array<uint, 8>& relativePositionInPart,
                                                 const std::array<uint, 8>& indicesOfOct, uint startPositionOfPart);

    //! \brief check if all nodes in an oct are valid to be written into an output file
    //! \details to be valid the nodes need to be: 1. have the type GEO_FLUID, 2. not be outside the current file part
    //! \param endPositionOfPart specifies the index of the last node in the current file part
    static bool areAllNodesInOctValidForWriting(const std::array<uint, 8>& indicesOfOct, const LBMSimulationParameter& parH,
                                                uint endPositionOfPart);

    //! \brief create the ending of the file name for a file part
    static std::string makePartFileNameEnding(uint level, int processID, int part, int timestep);
};

#endif
