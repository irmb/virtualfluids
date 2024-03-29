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
//! \addtogroup gpu_Kernel Kernel
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "Kernel/ADKernelManager.h"
#include "Cuda/CudaMemoryManager.h"

#include "BoundaryConditions/AdvectionDiffusion/AdvectionDiffusion.h"
#include "Parameter/Parameter.h"
#include "Kernel/AdvectionDiffusionKernel.h"
#include "PostProcessor/Concentration.cuh"

ADKernelManager::ADKernelManager(SPtr<Parameter> parameter, std::vector<SPtr<AdvectionDiffusionKernel>>& adkernels): para(parameter), adkernels(adkernels){}

////////////////////////////////////////////////////////////////////////////////
void ADKernelManager::setInitialNodeValuesAD(const int level, SPtr<CudaMemoryManager> cudaMemoryManager) const
{
    for (size_t index = 1; index <= para->getParH(level)->numberOfNodes; index++) {
        const real coordX = para->getParH(level)->coordinateX[index];
        const real coordY = para->getParH(level)->coordinateY[index];
        const real coordZ = para->getParH(level)->coordinateZ[index];

        real concentration;

        // call functor object with initial condition
        if (para->getInitialConditionAD()) {
            para->getInitialConditionAD()(coordX, coordY, coordZ, concentration);
        } else {
            concentration = real(0.0);
        }

        para->getParH(level)->concentration[index] = concentration;
    }

    cudaMemoryManager->cudaCopyConcentrationHostToDevice(level);
}

////////////////////////////////////////////////////////////////////////////////
void ADKernelManager::runADcollisionKernel(const int level)const
{
    adkernels[level]->run();
}

void ADKernelManager::runADslipBCKernel(const int level) const{
    if (para->getParD(level)->slipBC.numberOfBCnodes > 1) {
        AdvectionDiffusionSlipVelocityCompressible(
            para->getParD(level)->numberofthreads,
            para->getParD(level)->slipBC.normalX,
            para->getParD(level)->slipBC.normalY,
            para->getParD(level)->slipBC.normalZ,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->distributionsAD.f[0],
            para->getParD(level)->slipBC.k,
            para->getParD(level)->slipBC.q27[0],
            para->getParD(level)->slipBC.numberOfBCnodes,
            para->getParD(level)->omegaDiffusivity,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep);
    }
}

void ADKernelManager::runADgeometryBCKernel(const int level) const
{
    if (para->getParD(level)->geometryBC.numberOfBCnodes > 0) {
            AdvectionDiffusionBounceBack(
                para->getParD(level)->numberofthreads,
                para->getParD(level)->distributions.f[0],
                para->getParD(level)->distributionsAD.f[0],
                para->getParD(level)->AdvectionDiffusionNoSlipBC.concentration,
                para->getParD(level)->diffusivity,
                para->getParD(level)->AdvectionDiffusionNoSlipBC.k,
                para->getParD(level)->geometryBC.q27[0],
                para->getParD(level)->AdvectionDiffusionNoSlipBC.numberOfBcNodes,
                para->getParD(level)->omega,
                para->getParD(level)->neighborX,
                para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ,
                para->getParD(level)->numberOfNodes,
                para->getParD(level)->isEvenTimestep);
    }
}

void ADKernelManager::runADDirichletBCKernel(const int level) const{
    if (para->getParD(level)->AdvectionDiffusionDirichletBC.numberOfBcNodes > 0){
            AdvectionDiffusionDirichlet(
                para->getParD(level)->numberofthreads,
                para->getParD(level)->distributions.f[0],
                para->getParD(level)->distributionsAD.f[0],
                para->getParD(level)->AdvectionDiffusionDirichletBC.concentrationBC,
                para->getParD(level)->diffusivity,
                para->getParD(level)->velocityBC.k,
                para->getParD(level)->velocityBC.q27[0],
                para->getParD(level)->velocityBC.numberOfBCnodes,
                para->getParD(level)->omega,
                para->getParD(level)->neighborX,
                para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ,
                para->getParD(level)->numberOfNodes,
                para->getParD(level)->isEvenTimestep);

    }
}

////////////////////////////////////////////////////////////////////////////////
void ADKernelManager::printAD(const int level, SPtr<CudaMemoryManager> cudaMemoryManager) const
{
    CalcConcentration27(
        para->getParD(level)->numberofthreads,
        para->getParD(level)->concentration,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->distributionsAD.f[0],
        para->getParD(level)->isEvenTimestep);

    cudaMemoryManager->cudaCopyConcentrationDeviceToHost(level);
}

//! \}
