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
//! \addtogroup lbm
//! \{
//! \author Henry Korb
//=======================================================================================
#ifndef LBM_ADVECTION_DIFFUSION_EQUILIBRIUM_H
#define LBM_ADVECTION_DIFFUSION_EQUILIBRIUM_H
#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <lbm/constants/D3Q27.h>

namespace vf::lbm::advection_diffusion
{

constexpr real equilibrium(real weight, real concentration, real velocity, real cu_sq)
{
    using namespace vf::basics::constant;
    return weight * concentration * (c1o1 + c3o1 * velocity + c9o2 * velocity * velocity - cu_sq);
}

template <size_t direction>
constexpr real computeEquilibrium(real concentration, real velocityX, real velocityY, real velocityZ)
{
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;
    const real weight = getWeight<direction>();
    const real cu_sq = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;
    const real velocity = getVelocity<direction>(velocityX, velocityY, velocityZ);
    return weight * concentration * (c1o1 + c3o1 * velocity + c9o2 * velocity * velocity - c3o2 * cu_sq);
}
} // namespace vf::lbm::advection_diffusion
#endif