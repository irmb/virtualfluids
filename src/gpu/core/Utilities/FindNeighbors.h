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
//! \addtogroup gpu_Utilities Utilities
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#ifndef FIND_NEIGHBORS_H
#define FIND_NEIGHBORS_H

#include <map>

#include <lbm/constants/D3Q27.h>

#include "Parameter/Parameter.h"

namespace vf::gpu {

struct countersForPointerChasing
{
    uint counterInverse;
    uint counterX;
    uint counterY;
    uint counterZ;
};

const std::map<const size_t, const countersForPointerChasing> mapForPointerChasing =
{

    {vf::lbm::dir::d000, countersForPointerChasing{0, 0, 0, 0}},
    {vf::lbm::dir::dP00, countersForPointerChasing{0, 1, 0, 0}},
    {vf::lbm::dir::dM00, countersForPointerChasing{1, 0, 1, 1}},
    {vf::lbm::dir::d0P0, countersForPointerChasing{0, 0, 1, 0}},
    {vf::lbm::dir::d0M0, countersForPointerChasing{1, 1, 0, 1}},
    {vf::lbm::dir::d00P, countersForPointerChasing{0, 0, 0, 1}},
    {vf::lbm::dir::d00M, countersForPointerChasing{1, 1, 1, 0}},

    {vf::lbm::dir::dPP0, countersForPointerChasing{0, 1, 1, 0}},
    {vf::lbm::dir::dMM0, countersForPointerChasing{1, 0, 0, 1}},
    {vf::lbm::dir::dPM0, countersForPointerChasing{1, 2, 0, 1}},
    {vf::lbm::dir::dMP0, countersForPointerChasing{1, 0, 2, 1}},
    {vf::lbm::dir::dP0P, countersForPointerChasing{0, 1, 0, 1}},
    {vf::lbm::dir::dM0M, countersForPointerChasing{1, 0, 1, 0}},
    {vf::lbm::dir::dP0M, countersForPointerChasing{1, 2, 1, 0}},
    {vf::lbm::dir::dM0P, countersForPointerChasing{1, 0, 1, 2}},
    {vf::lbm::dir::d0PP, countersForPointerChasing{0, 0, 1, 1}},
    {vf::lbm::dir::d0MM, countersForPointerChasing{1, 1, 0, 0}},
    {vf::lbm::dir::d0PM, countersForPointerChasing{1, 1, 2, 0}},
    {vf::lbm::dir::d0MP, countersForPointerChasing{1, 1, 0, 2}},

    {vf::lbm::dir::dPPP, countersForPointerChasing{0, 1, 1, 1}},
    {vf::lbm::dir::dMPP, countersForPointerChasing{1, 0, 2, 2}},
    {vf::lbm::dir::dPMP, countersForPointerChasing{1, 2, 0, 2}},
    {vf::lbm::dir::dMMP, countersForPointerChasing{1, 0, 0, 2}},
    {vf::lbm::dir::dPPM, countersForPointerChasing{1, 2, 2, 0}},
    {vf::lbm::dir::dMPM, countersForPointerChasing{1, 0, 2, 0}},
    {vf::lbm::dir::dPMM, countersForPointerChasing{1, 2, 0, 0}},
    {vf::lbm::dir::dMMM, countersForPointerChasing{1, 0, 0, 0}}
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

}

#endif

//! \}
