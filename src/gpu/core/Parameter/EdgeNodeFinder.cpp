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
//! \author Anna Wellmann
//! \brief Functions for finding edge nodes in the multi-gpu implementation
//! \details Edge nodes are nodes, which are part of the communication in multiple directions
//! \ref master thesis of Anna Wellmann (p. 54-57)
//=======================================================================================
#include <vector>
#include <optional>

#include "EdgeNodeFinder.h"
#include "Parameter.h"

namespace vf::gpu
{
//! \brief Find nodes that are both received in the x-direction and sent in the y-direction
void findEdgeNodesXY(LBMSimulationParameter& parameterLB);
//! \brief Find nodes that are both received in the x-direction and sent in the z-direction
void findEdgeNodesXZ(LBMSimulationParameter& parameterLB);
//! \brief Find nodes that are both received in the y-direction and sent in the z-direction
void findEdgeNodesYZ(LBMSimulationParameter& parameterLB);
void findEdgeNodes(const std::vector<ProcessNeighbor27> &recvProcessNeighbor,
                   const std::vector<ProcessNeighbor27> &sendProcessNeighbor,
                   std::vector<LBMSimulationParameter::EdgeNodePositions> &edgeNodes);
std::optional<std::pair<int, int>> findIndexInSendNodes(const int nodeIndex,
                                                        const std::vector<ProcessNeighbor27> &sendProcessNeighbor);

void findEdgeNodesCommMultiGPU(Parameter& parameter)
{
    for (int level = 0; level <= parameter.getFine(); level++) {
        findEdgeNodesXY(*parameter.getParH(level));
        findEdgeNodesXZ(*parameter.getParH(level));
        findEdgeNodesYZ(*parameter.getParH(level));
    }
}

void findEdgeNodesXY(LBMSimulationParameter& parameterLB)
{
    findEdgeNodes(parameterLB.recvProcessNeighborX, parameterLB.sendProcessNeighborY,
                  parameterLB.edgeNodesXtoY);
}

void findEdgeNodesXZ(LBMSimulationParameter& parameterLB)
{
    findEdgeNodes(parameterLB.recvProcessNeighborX, parameterLB.sendProcessNeighborZ,
                  parameterLB.edgeNodesXtoZ);
}

void findEdgeNodesYZ(LBMSimulationParameter& parameterLB)
{
    findEdgeNodes(parameterLB.recvProcessNeighborY, parameterLB.sendProcessNeighborZ,
                  parameterLB.edgeNodesYtoZ);
}

void findEdgeNodes(const std::vector<ProcessNeighbor27> &recvProcessNeighbor,
                   const std::vector<ProcessNeighbor27> &sendProcessNeighbor,
                   std::vector<LBMSimulationParameter::EdgeNodePositions> &edgeNodes)
{
    for (uint neighbor = 0; neighbor < (unsigned int)(recvProcessNeighbor.size()); neighbor++) {
        for (int index = 0; index < recvProcessNeighbor[neighbor].numberOfNodes; index++) {
            if (auto sendIndices = findIndexInSendNodes(recvProcessNeighbor[neighbor].index[index], sendProcessNeighbor)) {
                edgeNodes.emplace_back(neighbor, index, sendIndices->first, sendIndices->second);
            }
        }
    }
}

std::optional<std::pair<int, int>> findIndexInSendNodes(const int nodeIndex,
                                                        const std::vector<ProcessNeighbor27> &sendProcessNeighbor)
{
    for (uint neighbor = 0; neighbor < (unsigned int)sendProcessNeighbor.size(); neighbor++) {
        for (int node = 0; node < sendProcessNeighbor[neighbor].numberOfNodes; node++) {
            if (sendProcessNeighbor[neighbor].index[node] == nodeIndex) {
                return std::pair<int, int>(neighbor, node);
            }
        }
    }
    return std::nullopt;
}

} // namespace vf::gpu
