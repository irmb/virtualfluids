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

#include "BoundaryConditions/Stress/Stress_Device.cuh"
#include "Parameter/Parameter.h"

void StressCompressible(Parameter *para,  QforBoundaryConditions* boundaryCondition, const int level)
{
    dim3 grid = vf::cuda::getCudaGrid(  para->getParD(level)->numberofthreads, boundaryCondition->numberOfBCnodes);
    dim3 threads(para->getParD(level)->numberofthreads, 1, 1 );

    StressCompressible_Device<<< grid, threads >>> (
        para->getParD(level)->distributions.f[0],
        boundaryCondition->k,
        boundaryCondition->kN,
        boundaryCondition->q27[0],
        boundaryCondition->numberOfBCnodes,
        para->getParD(level)->omega,
        para->getParD(level)->turbViscosity,
        para->getParD(level)->velocityX,
        para->getParD(level)->velocityY,
        para->getParD(level)->velocityY,
        boundaryCondition->normalX,
        boundaryCondition->normalY,
        boundaryCondition->normalZ,
        boundaryCondition->Vx,
        boundaryCondition->Vy,
        boundaryCondition->Vz,
        boundaryCondition->Vx1,
        boundaryCondition->Vy1,
        boundaryCondition->Vz1,
        para->getParD(level)->wallModel.samplingOffset,
        para->getParD(level)->wallModel.z0,
        para->getHasWallModelMonitor(),
        para->getParD(level)->wallModel.u_star,
        para->getParD(level)->wallModel.Fx,
        para->getParD(level)->wallModel.Fy,
        para->getParD(level)->wallModel.Fz,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("StressCompressible_Device execution failed");
}

//////////////////////////////////////////////////////////////////////////
void StressBounceBackCompressible(Parameter *para,  QforBoundaryConditions* boundaryCondition, const int level)
{
    dim3 grid = vf::cuda::getCudaGrid( para->getParD(level)->numberofthreads, boundaryCondition->numberOfBCnodes);
    dim3 threads(para->getParD(level)->numberofthreads, 1, 1 );

    StressBounceBackCompressible_Device<<< grid, threads >>> (
        para->getParD(level)->distributions.f[0],
        boundaryCondition->k,
        boundaryCondition->kN,
        boundaryCondition->q27[0],
        boundaryCondition->numberOfBCnodes,
        para->getParD(level)->velocityX,
        para->getParD(level)->velocityY,
        para->getParD(level)->velocityY,
        boundaryCondition->normalX,
        boundaryCondition->normalY,
        boundaryCondition->normalZ,
        boundaryCondition->Vx,
        boundaryCondition->Vy,
        boundaryCondition->Vz,
        boundaryCondition->Vx1,
        boundaryCondition->Vy1,
        boundaryCondition->Vz1,
        para->getParD(level)->wallModel.samplingOffset,
        para->getParD(level)->wallModel.z0,
        para->getHasWallModelMonitor(),
        para->getParD(level)->wallModel.u_star,
        para->getParD(level)->wallModel.Fx,
        para->getParD(level)->wallModel.Fy,
        para->getParD(level)->wallModel.Fz,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("StressBounceBackCompressible_Device execution failed");
}

//////////////////////////////////////////////////////////////////////////
void StressBounceBackPressureCompressible(Parameter *para,  QforBoundaryConditions* boundaryCondition, const int level)
{
    dim3 grid = vf::cuda::getCudaGrid( para->getParD(level)->numberofthreads, boundaryCondition->numberOfBCnodes);
    dim3 threads(para->getParD(level)->numberofthreads, 1, 1 );

    StressBounceBackPressureCompressible_Device<<< grid, threads >>> (
        para->getParD(level)->distributions.f[0],
        boundaryCondition->k,
        boundaryCondition->kN,
        boundaryCondition->q27[0],
        boundaryCondition->numberOfBCnodes,
        para->getParD(level)->velocityX,
        para->getParD(level)->velocityY,
        para->getParD(level)->velocityY,
        boundaryCondition->normalX,
        boundaryCondition->normalY,
        boundaryCondition->normalZ,
        boundaryCondition->Vx,
        boundaryCondition->Vy,
        boundaryCondition->Vz,
        boundaryCondition->Vx1,
        boundaryCondition->Vy1,
        boundaryCondition->Vz1,
        para->getParD(level)->wallModel.samplingOffset,
        para->getParD(level)->wallModel.z0,
        para->getHasWallModelMonitor(),
        para->getParD(level)->wallModel.u_star,
        para->getParD(level)->wallModel.Fx,
        para->getParD(level)->wallModel.Fy,
        para->getParD(level)->wallModel.Fz,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("BBStressPressureDevice27 execution failed");
}

//! \}
