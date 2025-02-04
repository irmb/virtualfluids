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
#ifndef LBM_TURBULENTDIFFUSION_H
#define LBM_TURBULENTDIFFUSION_H
#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

#include <algorithm>
#include <cmath>

namespace vf::lbm::advection_diffusion
{

//! \brief An enumeration for selecting a template of the advection-diffusion turbulence model
enum class TurbulenceModel {
    //! - Not using eddy diffusivity
    None,
    //! - Scales the eddy-viscosity with a constant Prandtl number
    Default,
    //! - Computes diffusivity according to \ref <a
    //! href="https://doi.org/10.1175/1520-0469(1984)041<2052:ALESMF>2.0.CO;2"><b><C. Moeng, 1984
    //! \DOI:10.1175/1520-0469(1984)041<2052:ALESMF>2.0.CO;2>
    Moeng,
    //! Computes diffusivity according to \ref <a href="https://doi.org/10.1007/s10546-017-0288-4" ><b><M. Abkar and P. Moin,
    //! 2017 \DOI:10.1007/s10546-017-0288-4>
    AMDStratified
};

constexpr real calcTurbulentDiffusivityDefault(real turbulentViscosity, real turbulentPrandtlNumber)
{
    return turbulentViscosity / turbulentPrandtlNumber;
}

constexpr real calcTurbulentDiffusivityMoeng(real temperatureGradient, real turbulentViscosity, real buoyancyParameter)
{
    using namespace vf::basics::constant;
    real lengthScale = c1o1;
    const real bruntVaisalaFrequencySquared = buoyancyParameter * temperatureGradient;
    const real subgridTurbulentKineticEnergy = turbulentViscosity * turbulentViscosity;
    const real limit = c76o10 * c76o10 * subgridTurbulentKineticEnergy;
    if (bruntVaisalaFrequencySquared > limit) { // for numerical stability & clips to minimum of 1
        lengthScale = sqrtf((float)limit / (float)bruntVaisalaFrequencySquared);
    }
    const real invTurbulentPrandtlNumber = c1o1 + c2o1 * lengthScale;
    return turbulentViscosity * invTurbulentPrandtlNumber;
}

constexpr real calcTurbulentDiffusivityAMD(real SGSConstant, real dvxdx, real dvxdy, real dvxdz, real dvydx, real dvydy,
                                           real dvydz, real dvzdx, real dvzdy, real dvzdz, real dthetadx, real dthetady,
                                           real dthetadz)
{
    using namespace vf::basics::constant;
    const real numerator = -((dvxdx * dthetadx + dvxdy * dthetady + dvxdz * dthetadz) * dthetadx +
                             (dvydx * dthetadx + dvydy * dthetady + dvydz * dthetadz) * dthetady +
                             (dvzdx * dthetadx + dvzdy * dthetady + dvzdz * dthetadz) * dthetadz);
    const real denominator = dthetadx * dthetadx + dthetady * dthetady + dthetadz * dthetadz;
    if (denominator < cSmallSingle)
        return c0o1;
    constexpr real zero = c0o1;
    constexpr real upperLimit = c1o100;
    return std::clamp(SGSConstant * numerator / denominator, zero, upperLimit);
}
} // namespace vf::lbm::advection_diffusion
#endif
//! \}