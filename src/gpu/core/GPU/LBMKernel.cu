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
// includes, cuda
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

#include "LBM/LB.h"
#include <cuda_helper/CudaGrid.h>

// includes, kernels
#include "GPU/GPU_Kernels.cuh"

#include "Parameter/Parameter.h"

//////////////////////////////////////////////////////////////////////////
void QADPressDev7(
    unsigned int numberOfThreads,
    real* DD,
    real* DD7,
    real* temp,
    real* velo,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QADPress7<<< grid.grid, grid.threads >>>(
        DD,
        DD7,
        temp,
        velo,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QADPress7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADPressDev27(
    unsigned int numberOfThreads,
    real* DD,
    real* DD27,
    real* temp,
    real* velo,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QADPress27<<< grid.grid, grid.threads >>>(
        DD,
        DD27,
        temp,
        velo,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QADPress27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADPressNEQNeighborDev27(
    unsigned int numberOfThreads,
    real* DD,
    real* DD27,
    int* k_Q,
    int* k_N,
    int numberOfBCnodes,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QADPressNEQNeighbor27<<< grid.grid, grid.threads >>>(
        DD,
        DD27,
        k_Q,
        k_N,
        numberOfBCnodes,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
       getLastCudaError("QADPressNEQNeighbor27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADVelDev7(
    unsigned int numberOfThreads,
    real* DD,
    real* DD7,
    real* temp,
    real* velo,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QADVel7<<< grid.grid, grid.threads >>> (
        DD,
        DD7,
        temp,
        velo,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QADVel7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADVelDev27(
    unsigned int numberOfThreads,
    real* DD,
    real* DD27,
    real* temp,
    real* velo,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QADVel27<<< grid.grid, grid.threads >>> (
        DD,
        DD27,
        temp,
        velo,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QADVel27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADDev7(
    unsigned int numberOfThreads,
    real* DD,
    real* DD7,
    real* temp,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QAD7<<< grid.grid, grid.threads >>> (
        DD,
        DD7,
        temp,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QAD7 execution failed");
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void ADSlipVelDevComp(
    uint numberOfThreads,
    real * normalX,
    real * normalY,
    real * normalZ,
    real * distributions,
    real * distributionsAD,
    int* QindexArray,
    real * Qarrays,
    uint numberOfBCnodes,
    real omegaDiffusivity,
    uint * neighborX,
    uint * neighborY,
    uint * neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    AD_SlipVelDeviceComp <<< grid.grid, grid.threads >>> (
        normalX,
        normalY,
        normalZ,
        distributions,
        distributionsAD,
        QindexArray,
        Qarrays,
        numberOfBCnodes,
        omegaDiffusivity,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("AD_SlipVelDeviceComp execution failed");
}
//////////////////////////////////////////////////////////////////////////

void QADDirichletDev27(
    unsigned int numberOfThreads,
    real* DD,
    real* DD27,
    real* temp,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QADDirichlet27<<< grid.grid, grid.threads >>> (
        DD,
        DD27,
        temp,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QADDirichletDev27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADBBDev27(
    unsigned int numberOfThreads,
    real* DD,
    real* DD27,
    real* temp,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QADBB27<<< grid.grid, grid.threads >>> (
        DD,
        DD27,
        temp,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QADBB27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QNoSlipADincompDev7(
    unsigned int numberOfThreads,
    real* DD,
    real* DD7,
    real* temp,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QNoSlipADincomp7<<< grid.grid, grid.threads >>> (
        DD,
        DD7,
        temp,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QNoSlipADincomp7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QNoSlipADincompDev27(
    unsigned int numberOfThreads,
    real* DD,
    real* DD27,
    real* temp,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QNoSlipADincomp27<<< grid.grid, grid.threads >>> (
        DD,
        DD27,
        temp,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QNoSlipADincomp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADVeloIncompDev7(
    unsigned int numberOfThreads,
    real* DD,
    real* DD7,
    real* temp,
    real* velo,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QADVeloIncomp7<<< grid.grid, grid.threads >>> (
        DD,
        DD7,
        temp,
        velo,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QADVeloIncomp7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADVeloIncompDev27(
    unsigned int numberOfThreads,
    real* DD,
    real* DD27,
    real* temp,
    real* velo,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QADVeloIncomp27<<< grid.grid, grid.threads >>> (
        DD,
        DD27,
        temp,
        velo,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QADVeloIncomp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADPressIncompDev7(
    unsigned int numberOfThreads,
    real* DD,
    real* DD7,
    real* temp,
    real* velo,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QADPressIncomp7<<< grid.grid, grid.threads >>>(
        DD,
        DD7,
        temp,
        velo,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QADPressIncomp7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADPressIncompDev27(
    unsigned int numberOfThreads,
    real* DD,
    real* DD27,
    real* temp,
    real* velo,
    real diffusivity,
    int* k_Q,
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    QADPressIncomp27<<< grid.grid, grid.threads >>>(
        DD,
        DD27,
        temp,
        velo,
        diffusivity,
        k_Q,
        QQ,
        numberOfBCnodes,
        om1,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("QADPressIncomp27 execution failed");
}
