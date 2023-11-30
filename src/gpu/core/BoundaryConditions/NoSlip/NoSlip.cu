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
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

#include "LBM/LB.h"
#include <cuda_helper/CudaGrid.h>

#include "BoundaryConditions/NoSlip/NoSlip_Device.cuh"
#include "Parameter/Parameter.h"

void NoSlipBounceBack(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
    dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
    dim3 threads(parameterDevice->numberofthreads, 1, 1 );

    NoSlipBounceBack_Device<<< grid, threads >>> (
        parameterDevice->distributions.f[0],
        boundaryCondition->k,
        boundaryCondition->q27[0],
        boundaryCondition->numberOfBCnodes,
        parameterDevice->neighborX,
        parameterDevice->neighborY,
        parameterDevice->neighborZ,
        parameterDevice->numberOfNodes,
        parameterDevice->isEvenTimestep);
    getLastCudaError("NoSlipBounceBack_Device execution failed");
}

void NoSlipInterpolatedIncompressible(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
    dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
    dim3 threads(parameterDevice->numberofthreads, 1, 1 );

    NoSlipInterpolatedIncompressible_Device<<< grid, threads >>> (
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
      getLastCudaError("NoSlipInterpolatedIncompressible_Device execution failed");
}

void NoSlipInterpolatedCompressible(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
    dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
    dim3 threads(parameterDevice->numberofthreads, 1, 1 );

    NoSlipInterpolatedCompressible_Device<<< grid, threads >>> (
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
    getLastCudaError("NoSlipInterpolatedCompressible_Device execution failed");
}

