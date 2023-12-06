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

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;

namespace vf::lbm
{

//! \brief An enumeration for selecting a turbulence model
enum class TurbulenceModel {
    //! - Smagorinsky
    Smagorinsky,
    //! - AMD (Anisotropic Minimum Dissipation) model, see e.g. Rozema et al., Phys. Fluids 27, 085107 (2015),
    //! https://doi.org/10.1063/1.4928700
    AMD,
    //! - QR model by Verstappen
    QR,
    //! - No turbulence model
    None
};

inline __host__ __device__ real calcTurbulentViscositySmagorinsky(real Cs, real dxux, real dyuy, real dzuz, real Dxy, real Dxz,
                                                           real Dyz)
{
    return Cs * Cs * sqrt(c2o1 * (dxux * dxux + dyuy * dyuy + dzuz * dzuz) + Dxy * Dxy + Dxz * Dxz + Dyz * Dyz);
}

template <typename T>
__host__ __device__ T max( T a, T b )
{
    return ( a > b ) ? a : b;
}

inline __host__ __device__ real calcTurbulentViscosityQR(real C, real dxux, real dyuy, real dzuz, real Dxy, real Dxz, real Dyz)
{
    // ! Verstappen's QR model
    //! Second invariant of the strain-rate tensor
    real Q = c1o2 * (dxux * dxux + dyuy * dyuy + dzuz * dzuz) + c1o4 * (Dxy * Dxy + Dxz * Dxz + Dyz * Dyz);
    //! Third invariant of the strain-rate tensor (determinant)
    // real R = - dxux*dyuy*dzuz - c1o4*( Dxy*Dxz*Dyz + dxux*Dyz*Dyz + dyuy*Dxz*Dxz + dzuz*Dxy*Dxy );
    real R = -dxux * dyuy * dzuz + c1o4 * (-Dxy * Dxz * Dyz + dxux * Dyz * Dyz + dyuy * Dxz * Dxz + dzuz * Dxy * Dxy);
    return C * max(R, c0o1) / Q;
}

inline __host__ __device__ real calculateOmegaWithturbulentViscosity(real omega, real turbulenceViscosity)
{
    return omega / (c1o1 + c3o1 * omega * turbulenceViscosity);
}

} // namespace vf::lbm

#endif //TURBULENT_VISCOSITY_H
