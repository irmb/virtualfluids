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
#ifndef Precursor_Device_H
#define Precursor_Device_H

#include "Calculation/Calculation.h"

__global__ void PrecursorNonReflectiveCompressible_Device(
    int* subgridDistanceIndices,
    int numberOfBCnodes,
    int numberOfPrecursorNodes,
    int sizeQ,
    real omega,
    real* distributions,
    real* subgridDistances,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    uint* neighborsNT,
    uint* neighborsNB,
    uint* neighborsST,
    uint* neighborsSB,
    real* weights0PP,
    real* weights0PM,
    real* weights0MP,
    real* weights0MM,
    real* vLast,
    real* vCurrent,
    real velocityX,
    real velocityY,
    real velocityZ,
    real timeRatio,
    real velocityRatio,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep);

__global__ void PrecursorDistributions_Device(
    int* subgridDistanceIndices,
    int numberOfBCNodes,
    int numberOfPrecursorNodes,
    real* distributions,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    uint* neighborsNT,
    uint* neighborsNB,
    uint* neighborsST,
    uint* neighborsSB,
    real* weights0PP,
    real* weights0PM,
    real* weights0MP,
    real* weights0MM,
    real* fsLast,
    real* fsNext,
    real timeRatio,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep);

#endif
