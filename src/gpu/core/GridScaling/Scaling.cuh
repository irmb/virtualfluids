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
//! \author Martin Schoenherr, Soeren Peters
//======================================================================================

#ifndef SCALING_CUH
#define SCALING_CUH

#include <basics/DataTypes.h>

#include "Calculation/Calculation.h"

struct LBMSimulationParameter;

template <bool hasTurbulentViscosity>
void scaleCoarseToFineCompressible(
    LBMSimulationParameter* parameterDeviceC,
    LBMSimulationParameter* parameterDeviceF,
    ICells* interpolationCellsCoarseToFine,
    ICellNeigh& neighborCoarseToFine,
    CUstream_st* stream);

template <bool hasTurbulentViscosity>
void scaleFineToCoarseCompressible(
    LBMSimulationParameter* parameterDeviceC,
    LBMSimulationParameter* parameterDeviceF,
    ICells* icellFC,
    ICellNeigh& neighborFineToCoarse,
    CUstream_st* stream);


void scaleCoarseToFineAdvectionDiffusion(
    real* DC,
    real* DF,
    real* DD27C,
    real* DD27F,
    uint* neighborCX,
    uint* neighborCY,
    uint* neighborCZ,
    uint* neighborFX,
    uint* neighborFY,
    uint* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
    bool isEvenTimestep,
    uint* posCSWB,
    uint* posFSWB,
    uint kCF,
    real nu,
    real diffusivity_fine,
    uint numberOfThreads,
    ICellNeigh neighborCoarseToFine);

void scaleFineToCoarseAdvectionDiffusion(
    real* DC,
    real* DF,
    real* DD27C,
    real* DD27F,
    uint* neighborCX,
    uint* neighborCY,
    uint* neighborCZ,
    uint* neighborFX,
    uint* neighborFY,
    uint* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
    bool isEvenTimestep,
    uint* posC,
    uint* posFSWB,
    uint kFC,
    real nu,
    real diffusivity_coarse,
    uint numberOfThreads,
    ICellNeigh neighborFineToCoarse);

#endif
