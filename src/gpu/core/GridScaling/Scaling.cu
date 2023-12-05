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
#include "Scaling.cuh"

#include <helper_cuda.h>
#include <cuda_helper/CudaGrid.h>

#include "lbm/MacroscopicQuantities.h"
#include "lbm/constants/D3Q27.h"

#include "basics/constants/NumericConstants.h"

#include "Utilities/KernelUtilities.h"
#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;


// coarse to fine
template <bool hasTurbulentViscosity>
__global__ void scaleCoarseToFineCompressible_Device(
    real* distributionsCoarse,
    real* distributionsFine,
    uint* neighborXcoarse,
    uint* neighborYcoarse,
    uint* neighborZcoarse,
    uint* neighborXfine,
    uint* neighborYfine,
    uint* neighborZfine,
    unsigned long long numberOfLBnodesCoarse,
    unsigned long long numberOfLBnodesFine,
    bool isEvenTimestep,
    uint* indicesCoarseMMM,
    uint* indicesFineMMM,
    uint numberOfInterfaceNodes,
    real omegaCoarse,
    real omegaFine,
    real* turbulentViscosityCoarse,
    real* turbulentViscosityFine,
    ICellNeigh offsetCF);


__global__ void scaleCoarseToFineAdvectionDiffusion_Device(
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
    ICellNeigh neighborCoarseToFine);


// fine to coarse
template <bool hasTurbulentViscosity>
__global__ void scaleFineToCoarseCompressible_Device(
    real* distributionsCoarse,
    real* distributionsFine,
    uint* neighborXcoarse,
    uint* neighborYcoarse,
    uint* neighborZcoarse,
    uint* neighborXfine,
    uint* neighborYfine,
    uint* neighborZfine,
    unsigned long long numberOfLBnodesCoarse,
    unsigned long long numberOfLBnodesFine,
    bool isEvenTimestep,
    uint* indicesCoarse000,
    uint* indicesFineMMM,
    uint numberOfInterfaceNodes,
    real omegaCoarse,
    real omegaFine,
    real* turbulentViscosityCoarse,
    real* turbulentViscosityFine,
    ICellNeigh offsetFC);


__global__ void scaleFineToCoarseAdvectionDiffusion_Device(
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
    ICellNeigh neighborFineToCoarse);

//////////////////////////////////////////////////////////////////////////

template <bool hasTurbulentViscosity>
void scaleCoarseToFineCompressible(
    LBMSimulationParameter* parameterDeviceC,
    LBMSimulationParameter* parameterDeviceF,
    ICells* coarseToFine,
    ICellNeigh& neighborCoarseToFine,
    CUstream_st* stream)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDeviceC->numberofthreads, coarseToFine->numberOfCells);
    dim3 threads(parameterDeviceC->numberofthreads, 1, 1);

    scaleCoarseToFineCompressible_Device<hasTurbulentViscosity><<<grid, threads, 0, stream>>>(
        parameterDeviceC->distributions.f[0],
        parameterDeviceF->distributions.f[0],
        parameterDeviceC->neighborX,
        parameterDeviceC->neighborY,
        parameterDeviceC->neighborZ,
        parameterDeviceF->neighborX,
        parameterDeviceF->neighborY,
        parameterDeviceF->neighborZ,
        parameterDeviceC->numberOfNodes,
        parameterDeviceF->numberOfNodes,
        parameterDeviceC->isEvenTimestep,
        coarseToFine->coarseCellIndices,
        coarseToFine->fineCellIndices,
        coarseToFine->numberOfCells,
        parameterDeviceC->omega,
        parameterDeviceF->omega,
        parameterDeviceC->turbViscosity,
        parameterDeviceF->turbViscosity,
        neighborCoarseToFine);

    getLastCudaError("scaleCoarseToFineCompressible_Device execution failed");
}
template void scaleCoarseToFineCompressible<true>(
    LBMSimulationParameter* parameterDeviceC,
    LBMSimulationParameter* parameterDeviceF,
    ICells* coarseToFine,
    ICellNeigh& neighborCoarseToFine,
    CUstream_st* stream);
template void scaleCoarseToFineCompressible<false>(
    LBMSimulationParameter* parameterDeviceC,
    LBMSimulationParameter* parameterDeviceF,
    ICells* coarseToFine,
    ICellNeigh& neighborCoarseToFine,
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
    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCoarseToFineAdvectionDiffusion_Device<<<grid.grid, grid.threads>>>(
        DC,
        DF,
        DD27C,
        DD27F,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
        isEvenTimestep,
        posCSWB,
        posFSWB,
        kCF,
        nu,
        diffusivity_fine,
        neighborCoarseToFine);
    getLastCudaError("scaleCoarseToFineAdvectionDiffusion_Device execution failed");
}


template <bool hasTurbulentViscosity>
void scaleFineToCoarseCompressible(
    LBMSimulationParameter* parameterDeviceC,
    LBMSimulationParameter* parameterDeviceF,
    ICells* fineToCoarse,
    ICellNeigh& neighborFineToCoarse,
    CUstream_st* stream)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDeviceC->numberofthreads, fineToCoarse->numberOfCells);
    dim3 threads(parameterDeviceC->numberofthreads, 1, 1);

    scaleFineToCoarseCompressible_Device<hasTurbulentViscosity><<<grid, threads, 0, stream>>>(
        parameterDeviceC->distributions.f[0],
        parameterDeviceF->distributions.f[0],
        parameterDeviceC->neighborX,
        parameterDeviceC->neighborY,
        parameterDeviceC->neighborZ,
        parameterDeviceF->neighborX,
        parameterDeviceF->neighborY,
        parameterDeviceF->neighborZ,
        parameterDeviceC->numberOfNodes,
        parameterDeviceF->numberOfNodes,
        parameterDeviceC->isEvenTimestep,
        fineToCoarse->coarseCellIndices,
        fineToCoarse->fineCellIndices,
        fineToCoarse->numberOfCells,
        parameterDeviceC->omega,
        parameterDeviceF->omega,
        parameterDeviceC->turbViscosity,
        parameterDeviceF->turbViscosity,
        neighborFineToCoarse);

    getLastCudaError("scaleFineToCoarseCompressible_Device execution failed");
}
template void scaleFineToCoarseCompressible<true>(
    LBMSimulationParameter* parameterDeviceC,
    LBMSimulationParameter* parameterDeviceF,
    ICells* fineToCoarse,
    ICellNeigh& neighborFineToCoarse,
    CUstream_st* stream);
template void scaleFineToCoarseCompressible<false>(
    LBMSimulationParameter* parameterDeviceC,
    LBMSimulationParameter* parameterDeviceF,
    ICells* fineToCoarse,
    ICellNeigh& neighborFineToCoarse,
    CUstream_st* stream);


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
    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFineToCoarseAdvectionDiffusion_Device<<<grid.grid, grid.threads>>>(
        DC,
        DF,
        DD27C,
        DD27F,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
        isEvenTimestep,
        posC,
        posFSWB,
        kFC,
        nu,
        diffusivity_coarse,
        neighborFineToCoarse);
    getLastCudaError("scaleFineToCoarseAdvectionDiffusion_Device execution failed");
}
