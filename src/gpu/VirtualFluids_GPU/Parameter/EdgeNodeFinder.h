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
//! \file EdgeNodeFinder.h
//! \ingroup Parameter
//! \author Anna Wellmann
//! \brief Functions for finding edge nodes in the multi-gpu implementation
//! \details Edge nodes are nodes, which are part of the communication in multiple directions
//! \ref master thesis of Anna Wellmann (p. 54-57)

//=======================================================================================
#ifndef GPU_EDGENODES_H
#define GPU_EDGENODES_H

#include <vector>

#include "Core/DataTypes.h"
#include "basics/PointerDefinitions.h"
#include "gpu/VirtualFluids_GPU/LBM/LB.h"

class Parameter;

namespace vf::gpu
{
//! \brief Find nodes which are part of communication in multiple coordinate directions
void findEdgeNodesCommMultiGPU(SPtr<Parameter> parameter);
} // namespace vf::gpu

// anonymous namespace
namespace
{
static bool findIndexInSendNodes(int nodeIndex, const std::vector<ProcessNeighbor27>& sendProcessNeighbor, int &indexOfProcessNeighborSend, int &indexInSendBuffer);
static void findEdgeNodesXY(int level, SPtr<Parameter> parameter);
static void findEdgeNodesXZ(int level, SPtr<Parameter> parameter);
static void findEdgeNodesYZ(int level, SPtr<Parameter> parameter);
} // namespace

#endif