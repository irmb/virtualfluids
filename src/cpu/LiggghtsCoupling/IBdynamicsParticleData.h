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
//! \file DataSet3D.h
//! \ingroup Data
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef IBdynamicsParticleData_h
#define IBdynamicsParticleData_h

#include<array>

constexpr auto SOLFRAC_MIN = 0.001;
constexpr auto SOLFRAC_MAX = 0.999;

struct IBdynamicsParticleData {
public:
    IBdynamicsParticleData()
        : partId(0), solidFraction(0.)
    {
        uPart[0] = 0.;
        uPart[1] = 0.;
        uPart[2] = 0.;

        hydrodynamicForce[0] = 0.;
        hydrodynamicForce[1] = 0.;
        hydrodynamicForce[2] = 0.;
    };
    int partId;
    double solidFraction;
    std::array<double, 3> uPart;
    std::array<double, 3> hydrodynamicForce;
};

#endif
