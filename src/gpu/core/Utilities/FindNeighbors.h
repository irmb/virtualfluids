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
#ifndef FIND_NEIGHBORS_H
#define FIND_NEIGHBORS_H

#include <map>

#include <lbm/constants/D3Q27.h>

#include "Parameter/Parameter.h"

using namespace vf::lbm::dir;

struct countersForPointerChasing
{
    uint counterInverse;
    uint counterX;
    uint counterY;
    uint counterZ;
};

const std::map<const size_t, const countersForPointerChasing> mapForPointerChasing =
{
    {d000, countersForPointerChasing{0, 0, 0, 0}},
    {dP00, countersForPointerChasing{0, 1, 0, 0}},
    {dM00, countersForPointerChasing{1, 0, 1, 1}},
    {d0P0, countersForPointerChasing{0, 0, 1, 0}},
    {d0M0, countersForPointerChasing{1, 1, 0, 1}},
    {d00P, countersForPointerChasing{0, 0, 0, 1}},
    {d00M, countersForPointerChasing{1, 1, 1, 0}},

    {dPP0, countersForPointerChasing{0, 1, 1, 0}},
    {dMM0, countersForPointerChasing{1, 0, 0, 1}},
    {dPM0, countersForPointerChasing{1, 2, 0, 1}},
    {dMP0, countersForPointerChasing{1, 0, 2, 1}},
    {dP0P, countersForPointerChasing{0, 1, 0, 1}},
    {dM0M, countersForPointerChasing{1, 0, 1, 0}},
    {dP0M, countersForPointerChasing{1, 2, 1, 0}},
    {dM0P, countersForPointerChasing{1, 0, 1, 2}},
    {d0PP, countersForPointerChasing{0, 0, 1, 1}},
    {d0MM, countersForPointerChasing{1, 1, 0, 0}},
    {d0PM, countersForPointerChasing{1, 1, 2, 0}},
    {d0MP, countersForPointerChasing{1, 1, 0, 2}},

    {dPPP, countersForPointerChasing{0, 1, 1, 1}},
    {dMPP, countersForPointerChasing{1, 0, 2, 2}},
    {dPMP, countersForPointerChasing{1, 2, 0, 2}},
    {dMMP, countersForPointerChasing{1, 0, 0, 2}},
    {dPPM, countersForPointerChasing{1, 2, 2, 0}},
    {dMPM, countersForPointerChasing{1, 0, 2, 0}},
    {dPMM, countersForPointerChasing{1, 2, 0, 0}},
    {dMMM, countersForPointerChasing{1, 0, 0, 0}}
};

// Only use for fluid nodes!
inline uint getNeighborIndex(LBMSimulationParameter* parH, const uint position, const int direction)
{
    uint nodeIndex = position;

    if (mapForPointerChasing.at(direction).counterInverse != 0) {
        nodeIndex = parH->neighborInverse[nodeIndex];
    }

    for (uint x = 0; x < mapForPointerChasing.at(direction).counterX; x++) {
        nodeIndex = parH->neighborX[nodeIndex];
    }

    for (uint y = 0; y < mapForPointerChasing.at(direction).counterY; y++) {
        nodeIndex = parH->neighborY[nodeIndex];
    }

    for (uint z = 0; z < mapForPointerChasing.at(direction).counterZ; z++) {
        nodeIndex = parH->neighborZ[nodeIndex];
    }

    return nodeIndex;
}

#endif
