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
#ifndef Stress_Device_H
#define Stress_Device_H

#include "LBM/LB.h"

__global__ void StressCompressible_Device(
    real* DD,
    int* k_Q,
    int* k_N,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    real* turbViscosity,
    real* vx,
    real* vy,
    real* vz,
    real* normalX,
    real* normalY,
    real* normalZ,
    real* vx_bc,
    real* vy_bc,
    real* vz_bc,
    real* vx1,
    real* vy1,
    real* vz1,
    int* samplingOffset,
    real* z0,
    bool  hasWallModelMonitor,
    real* u_star_monitor,
    real* Fx_monitor,
    real* Fy_monitor,
    real* Fz_monitor,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep);

__global__ void StressBounceBackCompressible_Device(
    real* DD,
    int* k_Q,
    int* k_N,
    real* QQ,
    unsigned int numberOfBCnodes,
    real* vx,
    real* vy,
    real* vz,
    real* normalX,
    real* normalY,
    real* normalZ,
    real* vx_bc,
    real* vy_bc,
    real* vz_bc,
    real* vx1,
    real* vy1,
    real* vz1,
    int* samplingOffset,
    real* z0,
    bool  hasWallModelMonitor,
    real* u_star_monitor,
    real* Fx_monitor,
    real* Fy_monitor,
    real* Fz_monitor,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep);

__global__ void StressBounceBackPressureCompressible_Device(
    real* DD,
    int* k_Q,
    int* k_N,
    real* QQ,
    unsigned int  numberOfBCnodes,
    real* vx,
    real* vy,
    real* vz,
    real* normalX,
    real* normalY,
    real* normalZ,
    real* vx_el,
    real* vy_el,
    real* vz_el,
    real* vx_w_mean,
    real* vy_w_mean,
    real* vz_w_mean,
    int* samplingOffset,
    real* z0,
    bool  hasWallModelMonitor,
    real* u_star_monitor,
    real* Fx_monitor,
    real* Fy_monitor,
    real* Fz_monitor,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep);

#endif
