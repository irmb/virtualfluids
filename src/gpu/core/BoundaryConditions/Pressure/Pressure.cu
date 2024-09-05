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
//! \author Martin Schoenherr, Anna Wellmann
//======================================================================================
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Calculation/Calculation.h"
#include <cuda_helper/CudaGrid.h>

#include "PressureNonEquilibrium.cuh"
#include "Parameter/Parameter.h"

#include <lbm/MacroscopicQuantities.h>

void PressureNonEquilibriumIncompressible(LBMSimulationParameter* parameterDevice, QforDirectionalBoundaryCondition* boundaryCondition)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDevice->numberofthreads, boundaryCondition->numberOfBCnodes);
    dim3 threads(parameterDevice->numberofthreads, 1, 1);

    auto macroscopicQuantityFunctor = [] __device__(const real *const &f /*[27]*/, real &rho, real &vx1, real &vx2, real &vx3) {
        vf::lbm::getIncompressibleMacroscopicValues(f, rho, vx1, vx2, vx3);
    };

    PressureNonEquilibrium_Device<decltype(macroscopicQuantityFunctor)><<< grid, threads >>> (
        boundaryCondition->RhoBC,
        parameterDevice->distributions.f[0],
        boundaryCondition->k,
        boundaryCondition->kN,
        boundaryCondition->numberOfBCnodes,
        parameterDevice->neighborX,
        parameterDevice->neighborY,
        parameterDevice->neighborZ,
        parameterDevice->numberOfNodes,
        parameterDevice->isEvenTimestep,
        static_cast<int>(boundaryCondition->direction),
        macroscopicQuantityFunctor);
    getLastCudaError("PressureNonEquilibriumIncompressible_Device execution failed");
}
//////////////////////////////////////////////////////////////////////////
void PressureNonEquilibriumCompressible(LBMSimulationParameter* parameterDevice, QforDirectionalBoundaryCondition* boundaryCondition)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDevice->numberofthreads, boundaryCondition->numberOfBCnodes);
    dim3 threads(parameterDevice->numberofthreads, 1, 1);

    auto macroscopicQuantityFunctor = [] __device__(const real *const &f /*[27]*/, real &rho, real &vx1, real &vx2, real &vx3) {
        vf::lbm::getCompressibleMacroscopicValues(f, rho, vx1, vx2, vx3);
    };

    PressureNonEquilibrium_Device<decltype(macroscopicQuantityFunctor)><<< grid, threads >>> (
        boundaryCondition->RhoBC,
        parameterDevice->distributions.f[0],
        boundaryCondition->k,
        boundaryCondition->kN,
        boundaryCondition->numberOfBCnodes,
        parameterDevice->neighborX,
        parameterDevice->neighborY,
        parameterDevice->neighborZ,
        parameterDevice->numberOfNodes,
        parameterDevice->isEvenTimestep,
        static_cast<int>(boundaryCondition->direction),
        macroscopicQuantityFunctor);
    getLastCudaError("PressureNonEquilibriumCompressible_Device execution failed");
}

//! \}
