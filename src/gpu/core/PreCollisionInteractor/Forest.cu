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
//! \addtogroup gpu_PreCollisionInteractor PreCollisionInteractor
//! \ingroup gpu_core core
//! \{
//! \author Mohammad Mehdi Mohammadi
//======================================================================================
#include "Forest.h"

#include <cstdio>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <basics/constants/NumericConstants.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <cuda_helper/CudaGrid.h>
#include <logger/Logger.h>

#include "Cuda/CudaMemoryManager.h"
#include "Cuda/CudaStreamManager.h"
#include "DataStructureInitializer/GridProvider.h"
#include "Parameter/Parameter.h"
#include "Utilities/GeometryUtils.h"
#include "Utilities/KernelUtilities.h"
#include "cuda_helper/CudaIndexCalculation.h"

using namespace vf::basics::constant;
namespace vf::gpu {

struct GridData
{
    const uint* indices;
    const uint nIndices;
    const real *coordsX, *coordsY, *coordsZ;
    const real *vx, *vy, *vz;
    real *fx, *fy, *fz;
    const real velocityRatio;
    const real forceRatio;
};

struct ForestData
{
    const real dragCoefficient;
    const real* leafAreaDensity;
    real *previousVelocityX, *previousVelocityY, *previousVelocityZ;
};

__global__ void readVelocityApplyBodyForces(GridData gridData, const ForestData forestData)
{
    const uint index = vf::cuda::get1DIndexFrom2DBlock();

    if (index >= gridData.nIndices)
        return;

    const uint gridIndex = gridData.indices[index];

    const real ux = gridData.vx[gridIndex];
    const real uy = gridData.vy[gridIndex];
    const real uz = gridData.vz[gridIndex];

    const real uxAveraged = c1o2 * (ux + forestData.previousVelocityX[index]);
    const real uyAveraged = c1o2 * (uy + forestData.previousVelocityY[index]);
    const real uzAveraged = c1o2 * (uz + forestData.previousVelocityZ[index]);

    forestData.previousVelocityX[index] = ux;
    forestData.previousVelocityY[index] = uy;
    forestData.previousVelocityZ[index] = uz;

    const real velocityMagnitude =
        sqrt(uxAveraged * uxAveraged + uyAveraged * uyAveraged + uzAveraged * uzAveraged) * gridData.velocityRatio;

    const real gridForceX = -velocityMagnitude * forestData.dragCoefficient * forestData.leafAreaDensity[index] * gridData.velocityRatio * uxAveraged;
    const real gridForceY = -velocityMagnitude * forestData.dragCoefficient * forestData.leafAreaDensity[index] * gridData.velocityRatio * uyAveraged;
    const real gridForceZ = -velocityMagnitude * forestData.dragCoefficient * forestData.leafAreaDensity[index] * gridData.velocityRatio * uzAveraged;

    gridData.fx[gridIndex] += gridForceX / gridData.forceRatio;
    gridData.fy[gridIndex] += gridForceY / gridData.forceRatio;
    gridData.fz[gridIndex] += gridForceZ / gridData.forceRatio;
}

void Forest::init()
{
    if (!para->getIsBodyForce())
        throw std::runtime_error("try to allocate Forest but BodyForce is not set in Parameter.");
    this->initForestIndices();
    this->initForestVelocities();
    this->initLeafAreaDensity();
}

void Forest::interact(int level, unsigned int t)
{
    if (level != this->level)
        return;

    const GridData gridData { this->forestIndicesD,
                              this->numberOfIndices,
                              para->getParD(this->level)->coordinateX,
                              para->getParD(this->level)->coordinateY,
                              para->getParD(this->level)->coordinateZ,
                              para->getParD(this->level)->velocityX,
                              para->getParD(this->level)->velocityY,
                              para->getParD(this->level)->velocityZ,
                              para->getParD(this->level)->forceX_SP,
                              para->getParD(this->level)->forceY_SP,
                              para->getParD(this->level)->forceZ_SP,
                              para->getScaledVelocityRatio(this->level),
                              para->getScaledForceRatio(this->level) };

    const ForestData forestData {
                                  this->dragCoefficient,
                                  this->leafAreaDensityD, 
                                  this->forestVelocitiesXPreviousTimestepD,
                                  this->forestVelocitiesYPreviousTimestepD,
                                  this->forestVelocitiesZPreviousTimestepD,
                                  };

    vf::cuda::CudaGrid boxGrid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, this->numberOfIndices);

    readVelocityApplyBodyForces<<<boxGrid.grid, boxGrid.threads>>>(gridData, forestData);
}

Forest::~Forest()
{
    cudaMemoryManager->cudaFreeForestIndices(this);
    cudaMemoryManager->cudaFreeForestVelocities(this);
    cudaMemoryManager->cudaFreeLeafAreaDensity(this);
}

void Forest::initForestIndices()
{
    std::vector<uint> forestIndices;

    std::vector<real> maxCoordX = para->getMaxCoordX();
    std::vector<real> maxCoordY = para->getMaxCoordY();
    std::vector<real> minCoordX = para->getMinCoordX();
    std::vector<real> minCoordY = para->getMinCoordY();
    const real spatialResolution = para->getLengthRatio();

    for (uint index = 1u; index < para->getParH(level)->numberOfNodes; ++index) {
        const bool withinX = (para->getParH(level)->coordinateX[index] > (minCoordX[0] + c2o1) *spatialResolution) &&
                             (para->getParH(level)->coordinateX[index] <= (maxCoordX[0] - c2o1) * spatialResolution);
        const bool withinY = (para->getParH(level)->coordinateY[index] > (minCoordY[0] + c2o1) *spatialResolution) &&
                             (para->getParH(level)->coordinateY[index] <= (maxCoordY[0] - c2o1) * spatialResolution);
        const bool withinZ =
            (para->getParH(level)->coordinateZ[index] > c0o1) && (para->getParH(level)->coordinateZ[index] <= height);

        if (withinX && withinY && withinZ) {
            forestIndices.push_back(index);
        }
    }

    this->numberOfIndices = static_cast<uint>(forestIndices.size());

    cudaMemoryManager->cudaAllocForestIndices(this);
    std::copy(forestIndices.begin(), forestIndices.end(), this->forestIndicesH);
    cudaMemoryManager->cudaCopyForestIndicesHtoD(this);
}

void Forest::initLeafAreaDensity()
{
    std::vector<real> leafAreaDensity;

    real maxLeafAreaDensity = c0o1;
    const real maxHeightLeafAreaDensity = 0.7 * this->height;

    if (this->plantAreaDensityFunction == PlantAreaDensityFunction::TopHeavy)
    {
        const int nSteps = c50o1;
        const real dz = this->height / nSteps;
        real integralShape = c0o1;

        for(int i=0; i<nSteps; ++i)
        {
            const real z = i * dz;
            const real nVal = (z < maxHeightLeafAreaDensity) ? c6o1 : c1o2;
            const real shape = powf((this->height - maxHeightLeafAreaDensity)/(this->height - z), nVal) *
                               expf(nVal * (c1o1 - (this->height - maxHeightLeafAreaDensity)/(this->height - z)));
            integralShape += shape * dz;
        }

        maxLeafAreaDensity = this->plantAreaIndex / integralShape;
    }

    for (uint i = 0; i < this->numberOfIndices; i++) 
    {
        const uint gridIndex = this->forestIndicesH[i];
        const real z = para->getParH(this->level)->coordinateZ[gridIndex];
        
        real variableLeafAreaDensity = c0o1;

        switch (this->plantAreaDensityFunction) 
        {
            case PlantAreaDensityFunction::TopHeavy:
            {
                const real variableLeafAreaDensityExponent = (z < maxHeightLeafAreaDensity) ? c6o1 : c1o2;

                variableLeafAreaDensity = maxLeafAreaDensity *
                                         powf((this->height - maxHeightLeafAreaDensity)/(this->height - z), variableLeafAreaDensityExponent) *
                                         expf(variableLeafAreaDensityExponent * (c1o1 - (this->height - maxHeightLeafAreaDensity)/(this->height - z)));
            } break;
            case PlantAreaDensityFunction::BottomHeavy:
                variableLeafAreaDensity = c2o1 * (this->plantAreaIndex / (this->height * this->height)) * (this->height - z);
                break;
            case PlantAreaDensityFunction::Flat:
                variableLeafAreaDensity = this->plantAreaIndex / this->height;
                break;
        }

        leafAreaDensity.push_back(variableLeafAreaDensity);
    }

        cudaMemoryManager->cudaAllocLeafAreaDensity(this);

        std::copy(leafAreaDensity.begin(), leafAreaDensity.end(), this->leafAreaDensityH);
        cudaMemoryManager->cudaCopyLeafAreaDensityHtoD(this);
}

void Forest::initForestVelocities()
{
    cudaMemoryManager->cudaAllocForestVelocities(this);

    std::vector<real> velocitiesX;
    std::vector<real> velocitiesY;
    std::vector<real> velocitiesZ;

    for (uint i = 0; i < this->numberOfIndices; i++) {
        uint idx = this->forestIndicesH[i];
        velocitiesX.push_back(para->getParH(this->level)->velocityX[idx]);
        velocitiesY.push_back(para->getParH(this->level)->velocityY[idx]);
        velocitiesZ.push_back(para->getParH(this->level)->velocityZ[idx]);
    }

    cudaMemoryManager->cudaCopyForestVelocitiesHtoD(this, velocitiesX, velocitiesY, velocitiesZ);
}

void Forest::getTaggedFluidNodes(GridProvider* gridProvider)
{
    std::vector<uint> indicesInBox(this->forestIndicesH, this->forestIndicesH + this->numberOfIndices);
    gridProvider->tagFluidNodeIndices(indicesInBox, CollisionTemplate::AllFeatures, this->level);
}
}

//! \}
