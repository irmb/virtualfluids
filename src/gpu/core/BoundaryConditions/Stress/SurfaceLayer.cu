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
//! \author Henry Korb
//=======================================================================================
#include <cuda_helper/CudaGrid.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>

#include <basics/DataTypes.h>

#include "BoundaryConditions/BoundaryConditionFactory.h"
#include "BoundaryConditions/Stress/SurfaceLayer.h"
#include "BoundaryConditions/Stress/SurfaceLayer_Device.cuh"
#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"
#include "Stress.h"

namespace vf::gpu {

using StressBC = BoundaryConditionFactory::StressBC;
using HeatFluxBC = BoundaryConditionFactory::SurfaceLayerBC;

TemperatureParameters getTemperatureParameters(LBMSimulationParameter* parameterDevice)
{
    return {
        parameterDevice->concentration,        parameterDevice->diffusivity,          parameterDevice->gravity,
        parameterDevice->turbulentDiffusivity, parameterDevice->referenceTemperature, parameterDevice->distributionsAD.f[0]
    };
}

template <StressBC stressBCType, HeatFluxBC heatFluxBCType, bool delayed>
void SurfaceLayer(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* surfaceLayerBoundaryCondition)
{
    const GridParameter gridParams = getStressBCGridParameter(parameterDevice);
    const TemperatureParameters tempParams = getTemperatureParameters(parameterDevice);

    const vf::cuda::CudaGrid grid(parameterDevice->numberofthreads, surfaceLayerBoundaryCondition->numberOfBCnodes);

    SurfaceLayerDevice27<stressBCType, heatFluxBCType, delayed><<<grid.grid, grid.threads>>>(
        gridParams, *surfaceLayerBoundaryCondition, parameterDevice->surfaceLayerWallModel.momentumParameters,
        parameterDevice->surfaceLayerWallModel.temperatureParameters, tempParams);

    getLastCudaError("SurfaceLayerDevice27 execution failed");
}

void SurfaceLayerBounceBackCompressibleHeatFlux(LBMSimulationParameter* parameterDevice,
                                                QforBoundaryConditions* surfaceLayerBoundaryCondition)
{
    SurfaceLayer<StressBC::StressBounceBackCompressible, HeatFluxBC::SurfaceHeatFlux, false>(parameterDevice,
                                                                                      surfaceLayerBoundaryCondition);
}

void SurfaceLayerBounceBackWithPressureCompressibleHeatFlux(LBMSimulationParameter* parameterDevice,
                                                            QforBoundaryConditions* surfaceLayerBoundaryCondition)
{
    SurfaceLayer<StressBC::StressBounceBackWithPressureCompressible, HeatFluxBC::SurfaceHeatFlux, false>(
        parameterDevice, surfaceLayerBoundaryCondition);
}


void SurfaceLayerInterpolatedCompressibleHeatFlux(LBMSimulationParameter* parameterDevice,
                                                  QforBoundaryConditions* surfaceLayerBoundaryCondition)
{
    SurfaceLayer<StressBC::StressInterpolatedCompressible, HeatFluxBC::SurfaceHeatFlux, false>(parameterDevice,
                                                                                        surfaceLayerBoundaryCondition);
}

void SurfaceLayerBounceBackCompressibleSurfaceTemperature(LBMSimulationParameter* parameterDevice,
                                                          QforBoundaryConditions* surfaceLayerBoundaryCondition)
{
    SurfaceLayer<StressBC::StressBounceBackCompressible, HeatFluxBC::SurfaceTemperature, false>(
        parameterDevice, surfaceLayerBoundaryCondition);
}

void SurfaceLayerBounceBackWithPressureCompressibleSurfaceTemperature(LBMSimulationParameter* parameterDevice,
                                                                      QforBoundaryConditions* surfaceLayerBoundaryCondition)
{
    SurfaceLayer<StressBC::StressBounceBackWithPressureCompressible, HeatFluxBC::SurfaceTemperature, false>(
        parameterDevice, surfaceLayerBoundaryCondition);
}

void SurfaceLayerInterpolatedCompressibleSurfaceTemperature(LBMSimulationParameter* parameterDevice,
                                                            QforBoundaryConditions* surfaceLayerBoundaryCondition)
{
    SurfaceLayer<StressBC::StressInterpolatedCompressible, HeatFluxBC::SurfaceTemperature, false>(
        parameterDevice, surfaceLayerBoundaryCondition);
}

}