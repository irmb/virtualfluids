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
//! \addtogroup gpu_BoundaryConditions BoundaryConditions
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

#include "Calculation/Calculation.h"
#include <cuda_helper/CudaGrid.h>

#include "BoundaryConditions/Velocity/Velocity_Device.cuh"
#include "Parameter/Parameter.h"

void VelocityBounceBack(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
    dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
    dim3 threads(parameterDevice->numberofthreads, 1, 1 );

    VelocityBounceBack_Device<<< grid, threads >>> (
        boundaryCondition->Vx,
        boundaryCondition->Vy,
        boundaryCondition->Vz,
        parameterDevice->distributions.f[0],
        boundaryCondition->k,
        boundaryCondition->q27[0],
        boundaryCondition->numberOfBCnodes,
        parameterDevice->neighborX,
        parameterDevice->neighborY,
        parameterDevice->neighborZ,
        parameterDevice->numberOfNodes,
        parameterDevice->isEvenTimestep);
    getLastCudaError("VelocityBounceBack_Device execution failed");
}

void VelocityInterpolatedIncompressible(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
    dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
    dim3 threads(parameterDevice->numberofthreads, 1, 1 );

    VelocityInterpolatedIncompressible_Device<<< grid, threads >>> (
        boundaryCondition->Vx,
        boundaryCondition->Vy,
        boundaryCondition->Vz,
        parameterDevice->distributions.f[0],
        boundaryCondition->k,
        boundaryCondition->q27[0],
        boundaryCondition->numberOfBCnodes,
        parameterDevice->omega,
        parameterDevice->neighborX,
        parameterDevice->neighborY,
        parameterDevice->neighborZ,
        parameterDevice->numberOfNodes,
        parameterDevice->isEvenTimestep);
    getLastCudaError("VelocityInterpolatedIncompressible_Device execution failed");
}

void VelocityInterpolatedCompressible(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
    dim3 threads(parameterDevice->numberofthreads, 1, 1 );

    VelocityInterpolatedCompressible_Device<<< grid, threads >>> (
        boundaryCondition->Vx,
        boundaryCondition->Vy,
        boundaryCondition->Vz,
        parameterDevice->distributions.f[0],
        boundaryCondition->k,
        boundaryCondition->q27[0],
        boundaryCondition->numberOfBCnodes,
        parameterDevice->omega,
        parameterDevice->neighborX,
        parameterDevice->neighborY,
        parameterDevice->neighborZ,
        parameterDevice->numberOfNodes,
        parameterDevice->isEvenTimestep);
   getLastCudaError("VelocityInterpolatedCompressible_Device execution failed");
}

void VelocityWithPressureInterpolatedCompressible(LBMSimulationParameter *parameterDevice, QforBoundaryConditions *boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

    VelocityWithPressureInterpolatedCompressible_Device<<< grid, threads >>> (
        boundaryCondition->Vx,
        boundaryCondition->Vy,
        boundaryCondition->Vz,
        parameterDevice->distributions.f[0],
        boundaryCondition->k,
        boundaryCondition->q27[0],
        boundaryCondition->numberOfBCnodes,
        parameterDevice->omega,
        parameterDevice->neighborX,
        parameterDevice->neighborY,
        parameterDevice->neighborZ,
        parameterDevice->numberOfNodes,
        parameterDevice->isEvenTimestep);
    getLastCudaError("VelocityWithPressureInterpolatedCompressible_Device execution failed");
}

//! \}
