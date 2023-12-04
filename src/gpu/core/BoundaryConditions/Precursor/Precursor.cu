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

#include "BoundaryConditions/Precursor/Precursor_Device.cuh"
#include "Parameter/Parameter.h"

void PrecursorNonReflectiveCompressible(
    LBMSimulationParameter* parameterDevice,
    QforPrecursorBoundaryConditions* boundaryCondition,
    real timeRatio,
    real velocityRatio)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(parameterDevice->numberofthreads, boundaryCondition->numberOfBCnodes);

    PrecursorNonReflectiveCompressible_Device<<< grid.grid, grid.threads >>>(
        boundaryCondition->k,
        boundaryCondition->numberOfBCnodes,
        boundaryCondition->numberOfPrecursorNodes,
        boundaryCondition->sizeQ,
        parameterDevice->omega,
        parameterDevice->distributions.f[0],
        boundaryCondition->q27[0],
        parameterDevice->neighborX,
        parameterDevice->neighborY,
        parameterDevice->neighborZ,
        boundaryCondition->planeNeighbor0PP,
        boundaryCondition->planeNeighbor0PM,
        boundaryCondition->planeNeighbor0MP,
        boundaryCondition->planeNeighbor0MM,
        boundaryCondition->weights0PP,
        boundaryCondition->weights0PM,
        boundaryCondition->weights0MP,
        boundaryCondition->weights0MM,
        boundaryCondition->last,
        boundaryCondition->current,
        boundaryCondition->velocityX,
        boundaryCondition->velocityY,
        boundaryCondition->velocityZ,
        timeRatio,
        velocityRatio,
        parameterDevice->numberOfNodes,
        parameterDevice->isEvenTimestep);
    getLastCudaError("PrecursorNonReflectiveCompressible_Device execution failed");
}

void PrecursorDistributions(
    LBMSimulationParameter* parameterDevice,
    QforPrecursorBoundaryConditions* boundaryCondition,
    real timeRatio,
    real velocityRatio)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(parameterDevice->numberofthreads, boundaryCondition->numberOfBCnodes);

    PrecursorDistributions_Device<<< grid.grid, grid.threads >>>(
        boundaryCondition->k,
        boundaryCondition->numberOfBCnodes,
        boundaryCondition->numberOfPrecursorNodes,
        parameterDevice->distributions.f[0],
        parameterDevice->neighborX,
        parameterDevice->neighborY,
        parameterDevice->neighborZ,
        boundaryCondition->planeNeighbor0PP,
        boundaryCondition->planeNeighbor0PM,
        boundaryCondition->planeNeighbor0MP,
        boundaryCondition->planeNeighbor0MM,
        boundaryCondition->weights0PP,
        boundaryCondition->weights0PM,
        boundaryCondition->weights0MP,
        boundaryCondition->weights0MM,
        boundaryCondition->last,
        boundaryCondition->current,
        timeRatio,
        parameterDevice->numberOfNodes,
        parameterDevice->isEvenTimestep);
    getLastCudaError("PrecursorDistributions_Device execution failed");

}

