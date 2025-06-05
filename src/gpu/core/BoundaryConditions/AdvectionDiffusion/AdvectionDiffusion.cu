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
// includes, cuda
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>

#include "Calculation/Calculation.h"
#include <cuda_helper/CudaGrid.h>

#include "BoundaryConditions/AdvectionDiffusion/AdvectionDiffusionBounceBack.cuh"
#include "BoundaryConditions/AdvectionDiffusion/AdvectionDiffusionDirichlet.cuh"
#include "BoundaryConditions/AdvectionDiffusion/AdvectionDiffusionNeumann.cuh"
#include "BoundaryConditions/AdvectionDiffusion/AdvectionDiffusionFlux.cuh"
#include "Parameter/Parameter.h"

void AdvectionDiffusionBounceBack(LBMSimulationParameter* parameterDevice,
                                  AdvectionDiffusionNoSlipBoundaryConditions bcParameters)
{
    const vf::cuda::CudaGrid grid(parameterDevice->numberofthreads, bcParameters.numberOfBCnodes);

    AdvectionDiffusionBounceBack_Device<<<grid.grid, grid.threads>>>(
        parameterDevice->distributionsAD.f[0], bcParameters, parameterDevice->neighborX, parameterDevice->neighborY,
        parameterDevice->neighborZ, parameterDevice->numberOfNodes, parameterDevice->isEvenTimestep);
    getLastCudaError("AdvectionDiffusionBounceBack_Device execution failed");
}

void AdvectionDiffusionFluxTurbulentViscosityCompressible(
    LBMSimulationParameter* parameterDevice, AdvectionDiffusionFluxBoundaryConditions bcParameters)
{
    const vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(parameterDevice->numberofthreads, bcParameters.numberOfBCnodes);

    AdvectionDiffusionFlux_Device<
        BoundaryConditionFactory::AdvectionDiffusionFluxBC::FluxTurbulentViscosityCompressible>
        <<<grid.grid, grid.threads>>>(parameterDevice->distributionsAD.f[0], bcParameters,
                                      parameterDevice->velocityX, parameterDevice->velocityY, parameterDevice->velocityZ,
                                      parameterDevice->turbulentDiffusivity, parameterDevice->diffusivity,
                                      parameterDevice->omegaDiffusivity, parameterDevice->neighborX,
                                      parameterDevice->neighborY, parameterDevice->neighborZ, parameterDevice->numberOfNodes,
                                      parameterDevice->isEvenTimestep);
    getLastCudaError("AdvectionDiffusionFlux_Device execution failed");
}

void AdvectionDiffusionFluxCompressible(LBMSimulationParameter* parameterDevice,
                                                AdvectionDiffusionFluxBoundaryConditions bcParameters)
{
    const vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(parameterDevice->numberofthreads, bcParameters.numberOfBCnodes);

    AdvectionDiffusionFlux_Device<
        BoundaryConditionFactory::AdvectionDiffusionFluxBC::FluxCompressible><<<grid.grid, grid.threads>>>(
        parameterDevice->distributionsAD.f[0], bcParameters, parameterDevice->velocityX,
        parameterDevice->velocityY, parameterDevice->velocityZ, parameterDevice->turbulentDiffusivity,
        parameterDevice->diffusivity, parameterDevice->omegaDiffusivity, parameterDevice->neighborX,
        parameterDevice->neighborY, parameterDevice->neighborZ, parameterDevice->numberOfNodes,
        parameterDevice->isEvenTimestep);
    getLastCudaError("AdvectionDiffusionFlux_Device execution failed");
}

void AdvectionDiffusionFluxBounceBack(LBMSimulationParameter* parameterDevice,
                                              AdvectionDiffusionFluxBoundaryConditions bcParameters)
{
    const vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(parameterDevice->numberofthreads, bcParameters.numberOfBCnodes);

    AdvectionDiffusionFlux_Device<BoundaryConditionFactory::AdvectionDiffusionFluxBC::FluxBounceBack>
        <<<grid.grid, grid.threads>>>(parameterDevice->distributionsAD.f[0], bcParameters,
                                      parameterDevice->velocityX, parameterDevice->velocityY, parameterDevice->velocityZ,
                                      parameterDevice->turbulentDiffusivity, parameterDevice->diffusivity,
                                      parameterDevice->omegaDiffusivity, parameterDevice->neighborX,
                                      parameterDevice->neighborY, parameterDevice->neighborZ, parameterDevice->numberOfNodes,
                                      parameterDevice->isEvenTimestep);
    getLastCudaError("AdvectionDiffusionFlux_Device execution failed");
}

void AdvectionDiffusionDirichletAntiBounceBackSlip(LBMSimulationParameter* parameterDevice,
                                                   AdvectionDiffusionDirichletBoundaryConditions bcParameters)
{
    const vf::cuda::CudaGrid grid(parameterDevice->numberofthreads, bcParameters.numberOfBCnodes);

    AdvectionDiffusionDirichlet_Device<BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletAntiBounceBackSlip>
        <<<grid.grid, grid.threads>>>(parameterDevice->distributionsAD.f[0], bcParameters, parameterDevice->neighborX,
                                      parameterDevice->neighborY, parameterDevice->neighborZ, parameterDevice->velocityX,
                                      parameterDevice->velocityY, parameterDevice->velocityZ, parameterDevice->numberOfNodes,
                                      parameterDevice->omegaDiffusivity, parameterDevice->isEvenTimestep);
    getLastCudaError("AdvectionDiffusionDirichlet_Device execution failed");
}

void AdvectionDiffusionDirichletInterpolatedSlip(LBMSimulationParameter* parameterDevice,
                                                 AdvectionDiffusionDirichletBoundaryConditions bcParameters)
{
    const vf::cuda::CudaGrid grid(parameterDevice->numberofthreads, bcParameters.numberOfBCnodes);

    AdvectionDiffusionDirichlet_Device<BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletAntiBounceBackSlip>
        <<<grid.grid, grid.threads>>>(parameterDevice->distributionsAD.f[0], bcParameters, parameterDevice->neighborX,
                                      parameterDevice->neighborY, parameterDevice->neighborZ, parameterDevice->velocityX,
                                      parameterDevice->velocityY, parameterDevice->velocityZ, parameterDevice->numberOfNodes,
                                      parameterDevice->omegaDiffusivity, parameterDevice->isEvenTimestep);
    getLastCudaError("AdvectionDiffusionDirichlet_Device execution failed");
}

void AdvectionDiffusionDirichletAntiBounceBackNoSlip(LBMSimulationParameter* parameterDevice,
                                                     AdvectionDiffusionDirichletBoundaryConditions bcParameters)
{
    const vf::cuda::CudaGrid grid(parameterDevice->numberofthreads, bcParameters.numberOfBCnodes);

    AdvectionDiffusionDirichlet_Device<
        BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletAntiBounceBackNoSlip><<<grid.grid, grid.threads>>>(
        parameterDevice->distributionsAD.f[0], bcParameters, parameterDevice->neighborX, parameterDevice->neighborY,
        parameterDevice->neighborZ, parameterDevice->velocityX, parameterDevice->velocityY, parameterDevice->velocityZ,
        parameterDevice->numberOfNodes, parameterDevice->omegaDiffusivity, parameterDevice->isEvenTimestep);
    getLastCudaError("AdvectionDiffusionDirichlet_Device execution failed");
}

void AdvectionDiffusionDirichletInterpolatedNoSlip(LBMSimulationParameter* parameterDevice,
                                                   AdvectionDiffusionDirichletBoundaryConditions bcParameters)
{
    const vf::cuda::CudaGrid grid(parameterDevice->numberofthreads, bcParameters.numberOfBCnodes);

    AdvectionDiffusionDirichlet_Device<
        BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletAntiBounceBackNoSlip><<<grid.grid, grid.threads>>>(
        parameterDevice->distributionsAD.f[0], bcParameters, parameterDevice->neighborX, parameterDevice->neighborY,
        parameterDevice->neighborZ, parameterDevice->velocityX, parameterDevice->velocityY, parameterDevice->velocityZ,
        parameterDevice->numberOfNodes, parameterDevice->omegaDiffusivity, parameterDevice->isEvenTimestep);
    getLastCudaError("AdvectionDiffusionDirichlet_Device execution failed");
}

void AdvectionDiffusionNeumannAntiBounceBackSlip(LBMSimulationParameter* parameterDevice,
                                                 AdvectionDiffusionNeumannBoundaryConditions bcParameters)
{
    const vf::cuda::CudaGrid grid(parameterDevice->numberofthreads, bcParameters.numberOfBCnodes);

    AdvectionDiffusionNeumann_Device<BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannAntiBounceBackSlip>
        <<<grid.grid, grid.threads>>>(parameterDevice->distributionsAD.f[0], bcParameters, parameterDevice->neighborX,
                                      parameterDevice->neighborY, parameterDevice->neighborZ, parameterDevice->velocityX,
                                      parameterDevice->velocityY, parameterDevice->velocityZ, parameterDevice->numberOfNodes,
                                      parameterDevice->omegaDiffusivity, parameterDevice->isEvenTimestep);
    getLastCudaError("AdvectionDiffusionNeumann_Device execution failed");
}

void AdvectionDiffusionNeumannInterpolatedSlip(LBMSimulationParameter* parameterDevice,
                                               AdvectionDiffusionNeumannBoundaryConditions bcParameters)
{
    const vf::cuda::CudaGrid grid(parameterDevice->numberofthreads, bcParameters.numberOfBCnodes);

    AdvectionDiffusionNeumann_Device<BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannAntiBounceBackSlip>
        <<<grid.grid, grid.threads>>>(parameterDevice->distributionsAD.f[0], bcParameters, parameterDevice->neighborX,
                                      parameterDevice->neighborY, parameterDevice->neighborZ, parameterDevice->velocityX,
                                      parameterDevice->velocityY, parameterDevice->velocityZ, parameterDevice->numberOfNodes,
                                      parameterDevice->omegaDiffusivity, parameterDevice->isEvenTimestep);
    getLastCudaError("AdvectionDiffusionNeumann_Device execution failed");
}

void AdvectionDiffusionNeumannAntiBounceBackNoSlip(LBMSimulationParameter* parameterDevice,
                                                   AdvectionDiffusionNeumannBoundaryConditions bcParameters)
{
    const vf::cuda::CudaGrid grid(parameterDevice->numberofthreads, bcParameters.numberOfBCnodes);

    AdvectionDiffusionNeumann_Device<BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannAntiBounceBackNoSlip>
        <<<grid.grid, grid.threads>>>(parameterDevice->distributionsAD.f[0], bcParameters, parameterDevice->neighborX,
                                      parameterDevice->neighborY, parameterDevice->neighborZ, parameterDevice->velocityX,
                                      parameterDevice->velocityY, parameterDevice->velocityZ, parameterDevice->numberOfNodes,
                                      parameterDevice->omegaDiffusivity, parameterDevice->isEvenTimestep);
    getLastCudaError("AdvectionDiffusionNeumann_Device execution failed");
}

void AdvectionDiffusionNeumannInterpolatedNoSlip(LBMSimulationParameter* parameterDevice,
                                                 AdvectionDiffusionNeumannBoundaryConditions bcParameters)
{
    const vf::cuda::CudaGrid grid(parameterDevice->numberofthreads, bcParameters.numberOfBCnodes);

    AdvectionDiffusionNeumann_Device<BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannAntiBounceBackNoSlip>
        <<<grid.grid, grid.threads>>>(parameterDevice->distributionsAD.f[0], bcParameters, parameterDevice->neighborX,
                                      parameterDevice->neighborY, parameterDevice->neighborZ, parameterDevice->velocityX,
                                      parameterDevice->velocityY, parameterDevice->velocityZ, parameterDevice->numberOfNodes,
                                      parameterDevice->omegaDiffusivity, parameterDevice->isEvenTimestep);
    getLastCudaError("AdvectionDiffusionNeumann_Device execution failed");
}

//! \}
