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
//! \addtogroup collision
//! \ingroup lbm
//! \{
//! \author Henry Korb, Henrik Asmuth
//======================================================================================

#ifndef TURBULENT_VISCOSITY_INLINES_CUH_
#define TURBULENT_VISCOSITY_INLINES_CUH_

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <algorithm>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

namespace vf::lbm
{

//! \brief An enumeration for selecting a turbulence model
enum class TurbulenceModel {
    //! - No turbulence model
    None,
    //! - Smagorinsky
    Smagorinsky,
    //! - QR model by Verstappen
    QR,
    //! - AMD (Anisotropic Minimum Dissipation) model, see e.g. Rozema et al., Phys. Fluids 27, 085107 (2015),
    //! https://doi.org/10.1063/1.4928700
    AMD,
    //! Computes diffusivity according to \ref <a href="https://doi.org/10.1007/s10546-017-0288-4" ><b><M. Abkar and P. Moin,
    //! 2017 \DOI:10.1007/s10546-017-0288-4>
    AMDStratified
};

inline __host__ __device__ real calcTurbulentViscositySmagorinsky(real SGSConstant, real dxux, real dyuy, real dzuz, real Dxy, real Dxz, real Dyz)
{
    using namespace vf::basics::constant;
    return SGSConstant * SGSConstant * sqrtf(c2o1 * (dxux * dxux + dyuy * dyuy + dzuz * dzuz) + Dxy * Dxy + Dxz * Dxz + Dyz * Dyz);
}

inline __host__ __device__ real calcTurbulentViscosityQR(real SGSConstant, real dxux, real dyuy, real dzuz, real Dxy, real Dxz, real Dyz)
{
    using namespace vf::basics::constant;
    // ! Verstappen's QR model
    //! Second invariant of the strain-rate tensor
    const real secondInvariant = c1o2 * (dxux * dxux + dyuy * dyuy + dzuz * dzuz) + c1o4 * (Dxy * Dxy + Dxz * Dxz + Dyz * Dyz);
    //! Third invariant of the strain-rate tensor (determinant)
    const real thirdInvariant = -dxux * dyuy * dzuz + c1o4 * (-Dxy * Dxz * Dyz + dxux * Dyz * Dyz + dyuy * Dxz * Dxz + dzuz * Dxy * Dxy);

    constexpr real zero = c0o1; // I Don't know why this is necessary, but it is apparently to pass it to std::max ...
    return SGSConstant * std::max(thirdInvariant, zero) / secondInvariant;
}

constexpr real calcDenominatorAMD(real dvxdx, real dvxdy, real dvxdz, real dvydx, real dvydy, real dvydz, real dvzdx,
                                  real dvzdy, real dvzdz)
{
    return dvxdx * dvxdx + dvydx * dvydx + dvzdx * dvzdx + dvxdy * dvxdy + dvydy * dvydy + dvzdy * dvzdy + dvxdz * dvxdz +
           dvydz * dvydz + dvzdz * dvzdz;
}

constexpr real calcNumeratorAMD(real dvxdx, real dvxdy, real dvxdz, real dvydx, real dvydy, real dvydz, real dvzdx,
                                real dvzdy, real dvzdz)
{
    return -((dvxdx * dvxdx + dvxdy * dvxdy + dvxdz * dvxdz) * dvxdx +
             (dvydx * dvydx + dvydy * dvydy + dvydz * dvydz) * dvydy +
             (dvzdx * dvzdx + dvzdy * dvzdy + dvzdz * dvzdz) * dvzdz +
             (dvxdx * dvydx + dvxdy * dvydy + dvxdz * dvydz) * (dvxdy + dvydx) +
             (dvxdx * dvzdx + dvxdy * dvzdy + dvxdz * dvzdz) * (dvxdz + dvzdx) +
             (dvydx * dvzdx + dvydy * dvzdy + dvydz * dvzdz) * (dvydz + dvzdy));
}
constexpr real calcTurbulentViscosityAMD(real SGSConstant, real dvxdx, real dvxdy, real dvxdz, real dvydx, real dvydy,
                                         real dvydz, real dvzdx, real dvzdy, real dvzdz)
{
    using namespace vf::basics::constant;
    const real denominator = calcDenominatorAMD(dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz);
    if (denominator == c0o1)
        return c0o1;
    const real numerator = calcNumeratorAMD(dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz);
    constexpr real zero = c0o1;
    return std::max(zero, SGSConstant * numerator / denominator);
}

constexpr real calcTurbulentViscosityAMDStratified(real SGSConstant, real dvxdx, real dvxdy, real dvxdz, real dvydx,
                                                   real dvydy, real dvydz, real dvzdx, real dvzdy, real dvzdz,
                                                   real buoyancyParameter, real dthetadx, real dthetady, real dthetadz)
{
    using namespace vf::basics::constant;
    const real denominator = calcDenominatorAMD(dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz);
    if (denominator <= cSmallSingle)
        return c0o1;
    const real temp = buoyancyParameter * (dvzdx * dthetadx + dvzdy * dthetady + dvzdz * dthetadz);
    const real numerator = calcNumeratorAMD(dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz);
    constexpr real zero = c0o1;
    constexpr real upperLimit = c1o100;
    return std::clamp(SGSConstant * (numerator + temp) / denominator, zero, upperLimit);
}

constexpr real calculateOmegaWithTurbulentViscosity(real omega, real turbulenceViscosity)
{
    using namespace vf::basics::constant;
    return omega / (c1o1 + c3o1 * omega * turbulenceViscosity);
}
} // namespace vf::lbm

#endif // TURBULENT_VISCOSITY_H

//! \}
