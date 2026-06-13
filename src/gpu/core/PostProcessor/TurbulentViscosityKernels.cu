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
//! \addtogroup gpu_PostProcessor PostProcessor
//! \ingroup gpu_core core
//! \{
//! \author Henry Korb, Henrik Asmuth
//======================================================================================

#include "TurbulentViscosityKernels.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cuda_helper/CudaGrid.h>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

#include "cuda_helper/CudaIndexCalculation.h"
#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"

using namespace vf::basics::constant;

namespace vf::gpu {

constexpr real derivative(const real* quantity, uint node, uint neighborP, uint neighborM, bool neighborPisFluid, bool neighborMisFluid)
{
    if (neighborPisFluid && neighborMisFluid)
        return c1o2 * (quantity[neighborP] - quantity[neighborM]);
    if (neighborPisFluid)
        return quantity[neighborP] - quantity[node];
    if (neighborMisFluid)
        return quantity[node] - quantity[neighborM];
    return NAN;
}


__global__ void calculateTurbulentViscosityAMDKernel(const real* vx, const real* vy, const real* vz, real* turbulentViscosity,
                                                const uint* neighborX, const uint* neighborY, const uint* neighborZ,
                                                const uint* neighborMMM, const uint* typeOfGridNode,
                                                unsigned long long numberOfLBnodes, real SGSConstant)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= numberOfLBnodes)
        return;
    if (typeOfGridNode[nodeIndex] != GEO_FLUID)
        return;

    const uint kP00 = neighborX[nodeIndex];
    const uint k0P0 = neighborY[nodeIndex];
    const uint k00P = neighborZ[nodeIndex];
    const uint kMMM = neighborMMM[nodeIndex];
    const uint kM00 = neighborZ[neighborY[kMMM]];
    const uint k0M0 = neighborZ[neighborX[kMMM]];
    const uint k00M = neighborY[neighborX[kMMM]];

    const bool fluidP00 = typeOfGridNode[kP00] == GEO_FLUID;
    const bool fluidM00 = typeOfGridNode[kM00] == GEO_FLUID;
    const bool fluid0P0 = typeOfGridNode[k0P0] == GEO_FLUID;
    const bool fluid0M0 = typeOfGridNode[k0M0] == GEO_FLUID;
    const bool fluid00P = typeOfGridNode[k00P] == GEO_FLUID;
    const bool fluid00M = typeOfGridNode[k00M] == GEO_FLUID;

    const real dvxdx = derivative(vx, nodeIndex, kP00, kM00, fluidP00, fluidM00);
    const real dvydx = derivative(vy, nodeIndex, kP00, kM00, fluidP00, fluidM00);
    const real dvzdx = derivative(vz, nodeIndex, kP00, kM00, fluidP00, fluidM00);

    const real dvxdy = derivative(vx, nodeIndex, k0P0, k0M0, fluid0P0, fluid0M0);
    const real dvydy = derivative(vy, nodeIndex, k0P0, k0M0, fluid0P0, fluid0M0);
    const real dvzdy = derivative(vz, nodeIndex, k0P0, k0M0, fluid0P0, fluid0M0);

    const real dvxdz = derivative(vx, nodeIndex, k00P, k00M, fluid00P, fluid00M);
    const real dvydz = derivative(vy, nodeIndex, k00P, k00M, fluid00P, fluid00M);
    const real dvzdz = derivative(vz, nodeIndex, k00P, k00M, fluid00P, fluid00M);

    turbulentViscosity[nodeIndex] =
        vf::lbm::calcTurbulentViscosityAMD(SGSConstant, dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz);
}

void calculateTurbulentViscosityAMD(Parameter* para, int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, para->getParH(level)->numberOfNodes);
    calculateTurbulentViscosityAMDKernel<<<grid.grid, grid.threads>>>(
        para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
        para->getParD(level)->turbulentViscosity, para->getParD(level)->neighborX, para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ, para->getParD(level)->neighborInverse, para->getParD(level)->typeOfGridNode,
        para->getParD(level)->numberOfNodes, para->getSGSConstant());
    getLastCudaError("calcAMD execution failed");
}

__global__ void calculateTurbulentViscosityAndDiffusivityAMDStratifiedKernel(
    const real* vx, const real* vy, const real* vz, const real* temperature, const real* referenceTemperature, unsigned long long numberOfLBnodes,
    const uint* typeOfGridNode, const uint* neighborX, const uint* neighborY, const uint* neighborZ, const uint* neighborMMM,
    real buoyancyParameter, real SGSConstant, real* turbulentViscosity, real* turbulentDiffusivity)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= numberOfLBnodes)
        return;
    if (typeOfGridNode[nodeIndex] != GEO_FLUID)
        return;

    const uint kP00 = neighborX[nodeIndex];
    const uint k0P0 = neighborY[nodeIndex];
    const uint k00P = neighborZ[nodeIndex];
    const uint kMMM = neighborMMM[nodeIndex];
    const uint kM00 = neighborZ[neighborY[kMMM]];
    const uint k0M0 = neighborZ[neighborX[kMMM]];
    const uint k00M = neighborY[neighborX[kMMM]];

    const bool fluidP00 = typeOfGridNode[kP00] == GEO_FLUID;
    const bool fluidM00 = typeOfGridNode[kM00] == GEO_FLUID;
    const bool fluid0P0 = typeOfGridNode[k0P0] == GEO_FLUID;
    const bool fluid0M0 = typeOfGridNode[k0M0] == GEO_FLUID;
    const bool fluid00P = typeOfGridNode[k00P] == GEO_FLUID;
    const bool fluid00M = typeOfGridNode[k00M] == GEO_FLUID;

    const real dvxdx = derivative(vx, nodeIndex, kP00, kM00, fluidP00, fluidM00);
    const real dvydx = derivative(vy, nodeIndex, kP00, kM00, fluidP00, fluidM00);
    const real dvzdx = derivative(vz, nodeIndex, kP00, kM00, fluidP00, fluidM00);
    const real dthetadx = derivative(temperature, nodeIndex, kP00, kM00, fluidP00, fluidM00);

    const real dvxdy = derivative(vx, nodeIndex, k0P0, k0M0, fluid0P0, fluid0M0);
    const real dvydy = derivative(vy, nodeIndex, k0P0, k0M0, fluid0P0, fluid0M0);
    const real dvzdy = derivative(vz, nodeIndex, k0P0, k0M0, fluid0P0, fluid0M0);
    const real dthetady = derivative(temperature, nodeIndex, k0P0, k0M0, fluid0P0, fluid0M0);

    const real dvxdz = derivative(vx, nodeIndex, k00P, k00M, fluid00P, fluid00M);
    const real dvydz = derivative(vy, nodeIndex, k00P, k00M, fluid00P, fluid00M);
    const real dvzdz = derivative(vz, nodeIndex, k00P, k00M, fluid00P, fluid00M);
    const real dthetadz = derivative(temperature, nodeIndex, k00P, k00M, fluid00P, fluid00M);
    const real dthetaRefdz = derivative(referenceTemperature, nodeIndex, k00P, k00M, fluid00P, fluid00M);

    turbulentViscosity[nodeIndex] =
        vf::lbm::calcTurbulentViscosityAMDStratified(SGSConstant, dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy,
                                                     dvzdz, buoyancyParameter, dthetadx, dthetady, dthetadz - dthetaRefdz);
    turbulentDiffusivity[nodeIndex] = vf::lbm::advection_diffusion::calcTurbulentDiffusivityAMD(
        SGSConstant, dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz, dthetadx, dthetady, dthetadz);
}

void calculateTurbulentViscosityAndDiffusivityAMDStratified(Parameter* para, int level)
{
    auto& parD = para->getParDeviceAsReference(level);
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(parD.numberofthreads, parD.numberOfNodes);
    calculateTurbulentViscosityAndDiffusivityAMDStratifiedKernel<<<grid.grid, grid.threads>>>(
        parD.velocityX, parD.velocityY, parD.velocityZ, parD.concentration, parD.localReferenceTemperature,
        parD.numberOfNodes, parD.typeOfGridNode, parD.neighborX, parD.neighborY, parD.neighborZ, parD.neighborInverse,
        para->getScaledBuoyancyFactor(level), para->getSGSConstant(), parD.turbulentViscosity, parD.turbulentDiffusivity);
    getLastCudaError("calcAMDStratified failed");
}

__global__ void calculateTurbulentDiffusivityMoengKernel(const real* temperature, real* turbulentDiffusivity, const real* turbulentViscosity,
                                                    const uint* neighborX, const uint* neighborY, const uint* neighborZ, const uint* neighborMMM,
                                                    const uint* typeOfGridNode, unsigned long long numberOfLBnodes,
                                                    real buoyancyParameter)
{
    using namespace vf::lbm::advection_diffusion;
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= numberOfLBnodes)
        return;
    if (typeOfGridNode[nodeIndex] != GEO_FLUID)
        return;
    const uint neighborP = neighborZ[nodeIndex];
    const uint neighborM = neighborX[neighborY[neighborMMM[nodeIndex]]];

    const bool fluidP = typeOfGridNode[neighborP] == GEO_FLUID;
    const bool fluidM = typeOfGridNode[neighborM] == GEO_FLUID;

    const real temperatureGradient = derivative(temperature, nodeIndex, neighborP, neighborM, fluidP, fluidM);

    const real turbulentViscosityNode = turbulentViscosity[nodeIndex];

    turbulentDiffusivity[nodeIndex] =
        calcTurbulentDiffusivityMoeng(temperatureGradient, turbulentViscosityNode, buoyancyParameter);
}

void calculateTurbulentDiffusivityMoeng(Parameter* para, int level)
{
    vf::cuda::CudaGrid grid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);
    calculateTurbulentDiffusivityMoengKernel<<<grid.grid, grid.threads>>>(
        para->getParD(level)->concentration, para->getParD(level)->turbulentDiffusivity, para->getParD(level)->turbulentViscosity,
        para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
        para->getParD(level)->neighborInverse, para->getParD(level)->typeOfGridNode, para->getParD(level)->numberOfNodes,
        para->getScaledBuoyancyFactor(level));
    getLastCudaError("calcTurbulentDiffusivityMoeng execution failed");
}

}
//! \}
