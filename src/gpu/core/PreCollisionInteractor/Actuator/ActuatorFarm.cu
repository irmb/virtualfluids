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
//! \author Henrik Asmuth, Henry Korb
//======================================================================================
#include "ActuatorFarm.h"
#include "ActuatorFarmInlines.h"

#include <algorithm>
#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <limits>

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
namespace vf::gpu
{

struct GridData
{
    const uint* indices;
    const uint nIndices;
    const real *coordsX, *coordsY, *coordsZ;
    const uint *neighborsX, *neighborsY, *neighborsZ, *neighborsWSB;
    const real *vx, *vy, *vz;
    real *fx, *fy, *fz;
    const real deltaX, velocityRatio, forceRatio;
};

struct TurbineData
{
    const real *posX, *posY, *posZ;
    const uint numberOfTurbines;
    const real smearingWidth;
};

struct ComponentData
{
    const real referenceLength;
    const real rotorBoundingMargin;
    const uint numberOfBladePointsPerTurbine;
    const uint totalNumberOfPoints;
    const real *coordsX, *coordsY, *coordsZ;
    real *velocitiesX, *velocitiesY, *velocitiesZ;
    const real *forcesX, *forcesY, *forcesZ;
    uint* gridIndices;
    const uint numberOfHubPointsPerTurbine;
    const real hubRadius, hubLength, hubPositionOffset;
    const uint numberOfTowerPointsPerTurbine;
    const real towerRadius, towerOffset, maxTowerHeight;
    const bool useVAWTVolume;
    const real vawtRotorHeight;
    const bool flagLocalSmearingWidth;
    const real* localSmearingWidth;
};

__global__ void interpolateVelocities(const GridData gridData, ComponentData componentData)
{
    const uint pointIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (pointIndex >= componentData.totalNumberOfPoints)
        return;

    const real coordX = componentData.coordsX[pointIndex];
    const real coordY = componentData.coordsY[pointIndex];
    const real coordZ = componentData.coordsZ[pointIndex];

    const uint kMMM = findNearestCellBSW(componentData.gridIndices[pointIndex], gridData.coordsX, gridData.coordsY,
                                         gridData.coordsZ, coordX, coordY, coordZ, gridData.neighborsX, gridData.neighborsY,
                                         gridData.neighborsZ, gridData.neighborsWSB);

    componentData.gridIndices[pointIndex] = kMMM;

    uint kPMM, kMPM, kMMP, kPPM, kPMP, kMPP, kPPP;
    getNeighborIndicesOfBSW(kMMM, kPMM, kMPM, kMMP, kPPM, kPMP, kMPP, kPPP, gridData.neighborsX, gridData.neighborsY,
                            gridData.neighborsZ);

    const real distX = (coordX - gridData.coordsX[kMMM]) / gridData.deltaX;
    const real distY = (coordY - gridData.coordsY[kMMM]) / gridData.deltaX;
    const real distZ = (coordZ - gridData.coordsZ[kMMM]) / gridData.deltaX;

    componentData.velocitiesX[pointIndex] =
        trilinearInterpolation(distX, distY, distZ, kMMM, kPMM, kMPM, kMMP, kPPM, kPMP, kMPP, kPPP, gridData.vx) *
        gridData.velocityRatio;
    componentData.velocitiesY[pointIndex] =
        trilinearInterpolation(distX, distY, distZ, kMMM, kPMM, kMPM, kMMP, kPPM, kPMP, kMPP, kPPP, gridData.vy) *
        gridData.velocityRatio;
    componentData.velocitiesZ[pointIndex] =
        trilinearInterpolation(distX, distY, distZ, kMMM, kPMM, kMPM, kMMP, kPPM, kPMP, kMPP, kPPP, gridData.vz) *
        gridData.velocityRatio;
}

template <bool UseLocalSmearingWidth>
__device__ void accumulateSmearingForces(const ComponentData& componentData, real gridCoordX, real gridCoordY,
                                         real gridCoordZ, uint startIndex, uint count, real smearingWidth, real& gridForceX,
                                         real& gridForceY, real& gridForceZ)
{
    for (uint i = 0; i < count; i++) {
        const uint node = startIndex + i;
        const real distX = componentData.coordsX[node] - gridCoordX;
        const real distY = componentData.coordsY[node] - gridCoordY;
        const real distZ = componentData.coordsZ[node] - gridCoordZ;

        if constexpr (UseLocalSmearingWidth)
            smearingWidth = componentData.localSmearingWidth[node];

        const real eta = gaussianSmearing(distX, distY, distZ, smearingWidth);
        gridForceX += componentData.forcesX[node] * eta;
        gridForceY += componentData.forcesY[node] * eta;
        gridForceZ += componentData.forcesZ[node] * eta;
    }
}

template <bool UseVAWT, bool HasHub, bool HasTower, bool UseLocalSmearingWidth>
__global__ void applyBodyForces(GridData gridData, const TurbineData turbineData, const ComponentData componentData)
{
    const uint index = vf::cuda::get1DIndexFrom2DBlock();

    if (index >= gridData.nIndices)
        return;

    const uint gridIndex = gridData.indices[index];
    const real gridCoordX = gridData.coordsX[gridIndex];
    const real gridCoordY = gridData.coordsY[gridIndex];
    const real gridCoordZ = gridData.coordsZ[gridIndex];

    real gridForceX = c0o1;
    real gridForceY = c0o1;
    real gridForceZ = c0o1;

    const uint bladePointsTotal = componentData.numberOfBladePointsPerTurbine * turbineData.numberOfTurbines;

    for (uint turbine = 0; turbine < turbineData.numberOfTurbines; turbine++) {

        const real turbinePosX = turbineData.posX[turbine];
        const real turbinePosY = turbineData.posY[turbine];
        const real turbinePosZ = turbineData.posZ[turbine];

        // --- Rotor bounding-volume check ---------------------
        bool isInRotorVolume;
        if constexpr (UseVAWT) {
            isInRotorVolume = inCylinderVolume<z, true>(
                gridCoordX, gridCoordY, gridCoordZ, turbinePosX, turbinePosY, turbinePosZ, componentData.vawtRotorHeight,
                c1o2 * componentData.referenceLength, componentData.rotorBoundingMargin);
        } else {
            const real distToHubX = gridCoordX - turbinePosX;
            const real distToHubY = gridCoordY - turbinePosY;
            const real distToHubZ = gridCoordZ - turbinePosZ;
            isInRotorVolume =
                inSphereVolume(distToHubX, distToHubY, distToHubZ, componentData.referenceLength, turbineData.smearingWidth);
        }

        // --- Blade forces ---------------------------------
        if (isInRotorVolume)
            accumulateSmearingForces<UseLocalSmearingWidth>(
                componentData, gridCoordX, gridCoordY, gridCoordZ, turbine * componentData.numberOfBladePointsPerTurbine,
                componentData.numberOfBladePointsPerTurbine, turbineData.smearingWidth, gridForceX, gridForceY, gridForceZ);

        if constexpr (HasHub) {
            const real centerX = turbinePosX - componentData.hubPositionOffset;
            const real centerY = turbinePosY;
            const real centerZ = turbinePosZ;

            if (inCylinderVolume<x, false>(gridCoordX, gridCoordY, gridCoordZ, centerX + c1o2 * componentData.hubLength,
                                           centerY, centerZ, componentData.hubLength, componentData.hubRadius,
                                           c3o2 * turbineData.smearingWidth))
                accumulateSmearingForces<UseLocalSmearingWidth>(
                    componentData, gridCoordX, gridCoordY, gridCoordZ,
                    bladePointsTotal + turbine * componentData.numberOfHubPointsPerTurbine,
                    componentData.numberOfHubPointsPerTurbine, turbineData.smearingWidth, gridForceX, gridForceY,
                    gridForceZ);
        }

        if constexpr (HasTower) {
            const real towerCenterX = turbinePosX + componentData.towerOffset;
            const real towerCenterY = turbinePosY;
            const real towerTopZ = turbinePosZ - componentData.hubRadius;
            const uint hubPointsTotal = componentData.numberOfHubPointsPerTurbine * turbineData.numberOfTurbines;

            if (inCylinderVolume<z, false>(gridCoordX, gridCoordY, gridCoordZ, towerCenterX, towerCenterY,
                                           towerTopZ - c1o2 * componentData.maxTowerHeight, componentData.maxTowerHeight,
                                           componentData.towerRadius, c3o2 * turbineData.smearingWidth))
                accumulateSmearingForces<UseLocalSmearingWidth>(
                    componentData, gridCoordX, gridCoordY, gridCoordZ,
                    bladePointsTotal + hubPointsTotal + turbine * componentData.numberOfTowerPointsPerTurbine,
                    componentData.numberOfTowerPointsPerTurbine, turbineData.smearingWidth, gridForceX, gridForceY,
                    gridForceZ);
        }
    }

    gridData.fx[gridIndex] += gridForceX / gridData.forceRatio;
    gridData.fy[gridIndex] += gridForceY / gridData.forceRatio;
    gridData.fz[gridIndex] += gridForceZ / gridData.forceRatio;
}

static void launchApplyBodyForces(bool useVAWT, bool hasHub, bool hasTower, bool UseLocalSmearingWidth,
                                  const vf::cuda::CudaGrid& sphereGrid, cudaStream_t stream, const GridData& gridData,
                                  const TurbineData& turbineData, const ComponentData& componentData)
{

    if (useVAWT) {
        if (UseLocalSmearingWidth)
            applyBodyForces<true, false, false, true><<<sphereGrid.grid, sphereGrid.threads, 0, stream>>>(gridData, turbineData, componentData);
        else
            applyBodyForces<true, false, false, false><<<sphereGrid.grid, sphereGrid.threads, 0, stream>>>(gridData, turbineData, componentData);
    } else {
        if (hasHub && hasTower)
            applyBodyForces<false, true, true, false><<<sphereGrid.grid, sphereGrid.threads, 0, stream>>>(gridData, turbineData, componentData);
        else if (hasHub)
            applyBodyForces<false, true, false, false><<<sphereGrid.grid, sphereGrid.threads, 0, stream>>>(gridData, turbineData, componentData);
        else if (hasTower)
            applyBodyForces<false, false, true, false><<<sphereGrid.grid, sphereGrid.threads, 0, stream>>>(gridData, turbineData, componentData);
        else
            applyBodyForces<false, false, false, false><<<sphereGrid.grid, sphereGrid.threads, 0, stream>>>(gridData, turbineData, componentData);
    }
}

void ActuatorFarm::init()
{
    if (!para->getIsBodyForce())
        throw std::runtime_error("try to allocate ActuatorFarm but BodyForce is not set in Parameter.");
    if (para->getDensityRatio() == c0o1)
        throw std::runtime_error("Parameter::densityRatio is zero. Non-zero density ratio needed for actuator farm.");
    if (!this->hasVAWTRotorVolume() && this->flagLocalSmearingWidth)
        throw std::runtime_error("ActuatorFarm::flagLocalSmearingWidth is only allowed for VAWTs.");
    if (this->hasVAWTRotorVolume() && this->numberOfHubPointsPerTurbine)
        throw std::runtime_error("ActuatorFarm::numberOfHubPointsPerTurbine is not defined for VAWTs.");
    if (this->hasVAWTRotorVolume() && this->numberOfTowerPointsPerTurbine)
        throw std::runtime_error("ActuatorFarm::numberOfTowerPointsPerTurbine is not defined for VAWTs.");
    this->initTurbineGeometries();
    this->initCoords();
    this->initIndices();
    this->initVelocities();
    this->initForces();
    this->initBoundingVolumes();

    this->streamIndex = para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::ActuatorFarm);
}

void ActuatorFarm::interact(int level, uint t)
{
    if (level != this->level)
        return;
    cudaStream_t stream = para->getStreamManager()->getStream(CudaStreamIndex::ActuatorFarm, this->streamIndex);

    if (useHostArrays)
        cudaMemoryManager->cudaCopyCoordsHtoD(this);

    if (this->writeOutput && ((t - this->tStartOut) % this->tOut == 0) && t >= this->tStartOut) {
        if (!useHostArrays) {
            cudaMemoryManager->cudaCopyCoordsDtoH(this);
            cudaMemoryManager->cudaCopyVelocitiesDtoH(this);
            cudaMemoryManager->cudaCopyForcesDtoH(this);
        }
        this->write(this->getFilename(t));
    }

    const GridData gridData { this->boundingVolumeIndicesD,
                              this->numberOfIndices,
                              para->getParD(this->level)->coordinateX,
                              para->getParD(this->level)->coordinateY,
                              para->getParD(this->level)->coordinateZ,
                              para->getParD(this->level)->neighborX,
                              para->getParD(this->level)->neighborY,
                              para->getParD(this->level)->neighborZ,
                              para->getParD(this->level)->neighborInverse,
                              para->getParD(this->level)->velocityX,
                              para->getParD(this->level)->velocityY,
                              para->getParD(this->level)->velocityZ,
                              para->getParD(this->level)->forceX_SP,
                              para->getParD(this->level)->forceY_SP,
                              para->getParD(this->level)->forceZ_SP,
                              para->getScaledLengthRatio(level),
                              para->getScaledVelocityRatio(level),
                              para->getScaledForceRatio(level) };

    const TurbineData turbineData { this->turbinePosXD, this->turbinePosYD, this->turbinePosZD, this->numberOfTurbines,
                                    this->smearingWidth };

    const ComponentData componentData {
        this->diameter,
        this->getVAWTRotorBoundingMargin(),
        this->numberOfBladePointsPerTurbine,
        this->getTotalNumberOfPoints(),
        this->coordsXDCurrentTimestep,
        this->coordsYDCurrentTimestep,
        this->coordsZDCurrentTimestep,
        this->velocitiesXDCurrentTimestep,
        this->velocitiesYDCurrentTimestep,
        this->velocitiesZDCurrentTimestep,
        this->forcesXDCurrentTimestep,
        this->forcesYDCurrentTimestep,
        this->forcesZDCurrentTimestep,
        this->indicesD,
        this->numberOfHubPointsPerTurbine,
        this->hubRadius,
        this->hubLength,
        this->hubPositionOffset,
        this->numberOfTowerPointsPerTurbine,
        this->towerRadius,
        this->towerOffset,
        this->maxTowerHeight,
        this->hasVAWTRotorVolume(),
        this->vawtRotorHeight,
        this->flagLocalSmearingWidth,
        this->getAllBladeLocalSmearingWidthDevice(),
    };
    vf::cuda::CudaGrid grid(para->getParH(level)->numberofthreads, this->getTotalNumberOfPoints());
    interpolateVelocities<<<grid.grid, grid.threads, 0, stream>>>(gridData, componentData);
    cudaStreamSynchronize(stream);

    if (useHostArrays)
        cudaMemoryManager->cudaCopyVelocitiesDtoH(this);

    const uint subIterationTimestep = para->getTimeStep(level, t, false);
    const real deltaT = para->getScaledTimeRatio(level);
    const real time = subIterationTimestep * deltaT;
    this->updateForcesAndCoordinates(time, deltaT);
    this->swapDeviceArrays();

    if (useHostArrays) {
        cudaMemoryManager->cudaCopyForcesHtoD(this);
    }
    vf::cuda::CudaGrid sphereGrid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, this->numberOfIndices);
    launchApplyBodyForces(this->hasVAWTRotorVolume(), this->numberOfHubPointsPerTurbine > 0,
                          this->numberOfTowerPointsPerTurbine > 0, this->flagLocalSmearingWidth, sphereGrid, stream,
                          gridData, turbineData, componentData);
    cudaStreamSynchronize(stream);
}

ActuatorFarm::~ActuatorFarm()
{
    cudaMemoryManager->cudaFreeBladeGeometries(this);
    cudaMemoryManager->cudaFreeCoords(this);
    cudaMemoryManager->cudaFreeVelocities(this);
    cudaMemoryManager->cudaFreeForces(this);
    cudaMemoryManager->cudaFreeIndices(this);
    cudaMemoryManager->cudaFreeBoundingVolumeIndices(this);
}

void ActuatorFarm::getTaggedFluidNodes(GridProvider* gridProvider)
{
    std::vector<uint> indicesInBoundingVolumes(this->boundingVolumeIndicesH,
                                               this->boundingVolumeIndicesH + this->numberOfIndices);
    gridProvider->tagFluidNodeIndices(indicesInBoundingVolumes, CollisionTemplate::AllFeatures, this->level);
}

void ActuatorFarm::initTurbineGeometries()
{
    cudaMemoryManager->cudaAllocBladeGeometries(this);

    std::copy(initialTurbinePositionsX.begin(), initialTurbinePositionsX.end(), turbinePosXH);
    std::copy(initialTurbinePositionsY.begin(), initialTurbinePositionsY.end(), turbinePosYH);
    std::copy(initialTurbinePositionsZ.begin(), initialTurbinePositionsZ.end(), turbinePosZH);

    cudaMemoryManager->cudaCopyBladeGeometriesHtoD(this);
}

void ActuatorFarm::initCoords()
{
    cudaMemoryManager->cudaAllocCoords(this);

    for (uint turbine = 0; turbine < this->numberOfTurbines; turbine++) {
        for (uint blade = 0; blade < this->numberOfBlades; blade++) {
            const real localAzimuth = this->azimuths[turbine] + blade * c2Pi / static_cast<real>(this->numberOfBlades);

            for (uint bladePoint = 0; bladePoint < this->numberOfPointsPerBlade; bladePoint++) {
                const uint node = calcPointIndexInBladeArrays({ turbine, blade, bladePoint }, this->numberOfPointsPerBlade,
                                                              this->numberOfBlades);

                real x, y, z;
                rotateFromBladeToGlobal(c0o1, c0o1, this->bladeRadii[bladePoint], x, y, z, localAzimuth);
                getAllBladeCoordsX()[node] = x + this->turbinePosXH[turbine];
                getAllBladeCoordsY()[node] = y + this->turbinePosYH[turbine];
                getAllBladeCoordsZ()[node] = z + this->turbinePosZH[turbine];
            }
        }
    }

    if (numberOfHubPoints > 0) {
        uint pointIndex = 0;
        for (uint turbine = 0; turbine < numberOfTurbines; turbine++)
            generateHubAxisPoints(turbine, pointIndex);
    }

    if (numberOfTowerPoints > 0) {
        uint pointIndex = 0;
        for (uint turbine = 0; turbine < numberOfTurbines; turbine++)
            generateTowerAxisPoints(turbine, pointIndex);
    }

    cudaMemoryManager->cudaCopyCoordsHtoD(this);
    std::swap(this->coordsXDCurrentTimestep, this->coordsXDPreviousTimestep);
    std::swap(this->coordsYDCurrentTimestep, this->coordsYDPreviousTimestep);
    std::swap(this->coordsZDCurrentTimestep, this->coordsZDPreviousTimestep);
    cudaMemoryManager->cudaCopyCoordsHtoD(this);
}

void ActuatorFarm::initVelocities()
{
    cudaMemoryManager->cudaAllocVelocities(this);

    const uint totalPoints = getTotalNumberOfPoints();
    std::fill_n(velocitiesXH, totalPoints, c0o1);
    std::fill_n(velocitiesYH, totalPoints, c0o1);
    std::fill_n(velocitiesZH, totalPoints, c0o1);

    cudaMemoryManager->cudaCopyVelocitiesHtoD(this);
    std::swap(this->velocitiesXDCurrentTimestep, this->velocitiesXDPreviousTimestep);
    std::swap(this->velocitiesYDCurrentTimestep, this->velocitiesYDPreviousTimestep);
    std::swap(this->velocitiesZDCurrentTimestep, this->velocitiesZDPreviousTimestep);
    cudaMemoryManager->cudaCopyVelocitiesHtoD(this);
}

void ActuatorFarm::initForces()
{
    cudaMemoryManager->cudaAllocForces(this);

    const uint totalPoints = getTotalNumberOfPoints();
    const bool requiresLocalSmearingWidth = this->requiresLocalSmearingWidth();
    std::fill_n(forcesXH, totalPoints, c0o1);
    std::fill_n(forcesYH, totalPoints, c0o1);
    std::fill_n(forcesZH, totalPoints, c0o1);

    if (requiresLocalSmearingWidth)
        std::fill_n(localSmearingWidthH, totalPoints, this->smearingWidth);

    cudaMemoryManager->cudaCopyForcesHtoD(this);
    std::swap(this->forcesXDCurrentTimestep, this->forcesXDPreviousTimestep);
    std::swap(this->forcesYDCurrentTimestep, this->forcesYDPreviousTimestep);
    std::swap(this->forcesZDCurrentTimestep, this->forcesZDPreviousTimestep);
    if (requiresLocalSmearingWidth)
        std::swap(this->localSmearingWidthDCurrentTimestep, this->localSmearingWidthDPreviousTimestep);
    cudaMemoryManager->cudaCopyForcesHtoD(this);
}

void ActuatorFarm::initIndices()
{
    cudaMemoryManager->cudaAllocIndices(this);

    std::fill_n(indicesH, getTotalNumberOfPoints(), 1);

    cudaMemoryManager->cudaCopyIndicesHtoD(this);
}

void ActuatorFarm::initBoundingVolumes()
{
    std::vector<uint> nodesInBoundingVolumes;

    // --- Rotor Bounding Volume ---
    const bool useVAWTVolume = this->hasVAWTRotorVolume();
    const real rotorBoundingSmearingWidth = this->getRotorBoundingSmearingWidth();
    const real rotorBoundingMargin = this->getVAWTRotorBoundingMargin();
    const real rotorBoundingVolumeRadius =
        useVAWTVolume ? c1o2 * this->diameter + rotorBoundingMargin
                      : getRotorBoundingVolumeRadius(this->diameter, rotorBoundingSmearingWidth, false);
    const real rotorBoundingVolumeHeight = this->vawtRotorHeight + c2o1 * rotorBoundingMargin;
    const real deltaX = para->getScaledLengthRatio(level);
    const real effectiveRotorRadius = std::max(c0o1, rotorBoundingVolumeRadius - deltaX);
    const real rotorBoundingVolumeInnerRadius = std::max(c0o1, c1o2 * this->diameter - rotorBoundingMargin);
    const real effectiveRotorInnerRadius =
        rotorBoundingVolumeInnerRadius > c0o1 ? std::max(c0o1, rotorBoundingVolumeInnerRadius + deltaX) : c0o1;
    const real effectiveRotorHeight = std::max(c0o1, rotorBoundingVolumeHeight - c2o1 * deltaX);
    const uint minimumNumberOfNodesPerRotorBoundingVolume =
        useVAWTVolume ? std::max(uint(1), uint(cPi *
                                               std::max(c0o1, std::pow(effectiveRotorRadius, c2o1) -
                                                                  std::pow(effectiveRotorInnerRadius, c2o1)) *
                                               effectiveRotorHeight / std::pow(deltaX, c3o1)))
                      : uint(c4o3 * cPi * std::pow(effectiveRotorRadius, c3o1) / std::pow(deltaX, c3o1));

    for (uint turbine = 0; turbine < this->numberOfTurbines; turbine++) {

        const real turbinePosX = this->turbinePosXH[turbine];
        const real turbinePosY = this->turbinePosYH[turbine];
        const real turbinePosZ = this->turbinePosZH[turbine];

        uint nodesInThisTurbineBoundingVolume = 0;

        for (size_t pos = 1; pos <= para->getParH(this->level)->numberOfNodes; pos++) {
            const real nodeX = para->getParH(this->level)->coordinateX[pos];
            const real nodeY = para->getParH(this->level)->coordinateY[pos];
            const real nodeZ = para->getParH(this->level)->coordinateZ[pos];

            bool inRotorVolume = false;
            if (useVAWTVolume) {
                inRotorVolume = inCylinderVolume<z, true>(nodeX, nodeY, nodeZ, turbinePosX, turbinePosY, turbinePosZ,
                                                          this->vawtRotorHeight, c1o2 * this->diameter, rotorBoundingMargin);
            } else {
                const real distX = nodeX - turbinePosX;
                const real distY = nodeY - turbinePosY;
                const real distZ = nodeZ - turbinePosZ;
                inRotorVolume = inSphereVolume(distX, distY, distZ, this->diameter, rotorBoundingSmearingWidth);
            }

            if (inRotorVolume) {
                nodesInBoundingVolumes.push_back(uint(pos));
                nodesInThisTurbineBoundingVolume++;
            }
        }

        if (nodesInThisTurbineBoundingVolume < minimumNumberOfNodesPerRotorBoundingVolume) {
            VF_LOG_CRITICAL("Found only {} nodes in rotor bounding volume of turbine no. {}, expected at least {}!",
                            nodesInThisTurbineBoundingVolume, turbine, minimumNumberOfNodesPerRotorBoundingVolume);
            throw std::runtime_error(
                "ActuatorFarm::initBoundingVolumes: Turbine rotor bounding volume partially out of domain.");
        }
    }

    // --- Hub Bounding Volume ---
    if (numberOfHubPoints > 0) {
        for (uint turbine = 0; turbine < numberOfTurbines; turbine++) {
            const real centerX = turbinePosXH[turbine] - hubPositionOffset;
            const real centerY = turbinePosYH[turbine];
            const real centerZ = turbinePosZH[turbine];

            for (size_t pos = 1; pos <= para->getParH(level)->numberOfNodes; pos++) {
                const real nodeX = para->getParH(level)->coordinateX[pos];
                const real nodeY = para->getParH(level)->coordinateY[pos];
                const real nodeZ = para->getParH(level)->coordinateZ[pos];

                if (inCylinderVolume<x, false>(nodeX, nodeY, nodeZ, centerX + c1o2 * hubLength, centerY, centerZ, hubLength,
                                               hubRadius, c3o2 * this->smearingWidth)) {
                    nodesInBoundingVolumes.push_back((uint)pos);
                }
            }
        }
    }

    // --- Tower Bounding Volume ---
    if (numberOfTowerPoints > 0) {
        maxTowerHeight = *std::max_element(towerHeights.begin(), towerHeights.end());
        for (uint turbine = 0; turbine < numberOfTurbines; turbine++) {
            const real towerCenterX = turbinePosXH[turbine] + towerOffset;
            const real towerCenterY = turbinePosYH[turbine];
            const real towerTopZ = turbinePosZH[turbine] - hubRadius;
            const real towerBottomZ = towerTopZ - towerHeights[turbine];

            real localMinZ = std::numeric_limits<real>::max();
            for (size_t pos = 1; pos <= para->getParH(level)->numberOfNodes; pos++) {
                const real nodeX = para->getParH(level)->coordinateX[pos];
                const real nodeY = para->getParH(level)->coordinateY[pos];
                const real nodeZ = para->getParH(level)->coordinateZ[pos];

                if (std::abs(nodeX - towerCenterX) < deltaX && std::abs(nodeY - towerCenterY) < deltaX)
                    localMinZ = std::min(localMinZ, nodeZ);

                if (inCylinderVolume<z, false>(nodeX, nodeY, nodeZ, towerCenterX, towerCenterY,
                                               towerTopZ - c1o2 * towerHeights[turbine], towerHeights[turbine], towerRadius,
                                               c3o2 * smearingWidth)) {
                    nodesInBoundingVolumes.push_back((uint)pos);
                }
            }

            if (towerBottomZ < localMinZ)
                throw std::runtime_error("ActuatorFarm::initBoundingVolumes: Tower of turbine " + std::to_string(turbine) +
                                         " extends below the domain (tower bottom z = " + std::to_string(towerBottomZ) +
                                         ", domain min z = " + std::to_string(localMinZ) +
                                         "). Reduce the tower height so that it does not extend below the domain.");
        }
    }

    std::sort(nodesInBoundingVolumes.begin(), nodesInBoundingVolumes.end());
    nodesInBoundingVolumes.erase(std::unique(nodesInBoundingVolumes.begin(), nodesInBoundingVolumes.end()),
                                 nodesInBoundingVolumes.end());

    this->numberOfIndices = uint(nodesInBoundingVolumes.size());

    cudaMemoryManager->cudaAllocBoundingVolumeIndices(this);
    std::copy(nodesInBoundingVolumes.begin(), nodesInBoundingVolumes.end(), this->boundingVolumeIndicesH);
    cudaMemoryManager->cudaCopyBoundingVolumeIndicesHtoD(this);
}

void ActuatorFarm::setAllBladeCoords(const real* bladeCoordsX, const real* bladeCoordsY, const real* bladeCoordsZ) const
{
    std::copy_n(bladeCoordsX, this->numberOfBladePoints, this->getAllBladeCoordsX());
    std::copy_n(bladeCoordsY, this->numberOfBladePoints, this->getAllBladeCoordsY());
    std::copy_n(bladeCoordsZ, this->numberOfBladePoints, this->getAllBladeCoordsZ());
}

void ActuatorFarm::setAllBladeVelocities(const real* bladeVelocitiesX, const real* bladeVelocitiesY,
                                         const real* bladeVelocitiesZ) const
{
    std::copy_n(bladeVelocitiesX, this->numberOfBladePoints, this->getAllBladeVelocitiesX());
    std::copy_n(bladeVelocitiesY, this->numberOfBladePoints, this->getAllBladeVelocitiesY());
    std::copy_n(bladeVelocitiesZ, this->numberOfBladePoints, this->getAllBladeVelocitiesZ());
}

void ActuatorFarm::setAllBladeForces(const real* bladeForcesX, const real* bladeForcesY, const real* bladeForcesZ) const
{
    std::copy_n(bladeForcesX, this->numberOfBladePoints, this->getAllBladeForcesX());
    std::copy_n(bladeForcesY, this->numberOfBladePoints, this->getAllBladeForcesY());
    std::copy_n(bladeForcesZ, this->numberOfBladePoints, this->getAllBladeForcesZ());
}

void ActuatorFarm::setTurbineBladeCoords(size_t turbine, const real* bladeCoordsX, const real* bladeCoordsY,
                                         const real* bladeCoordsZ) const
{
    std::copy_n(bladeCoordsX, this->numberOfBladePointsPerTurbine,
                &this->getAllBladeCoordsX()[turbine * this->numberOfBladePointsPerTurbine]);
    std::copy_n(bladeCoordsY, this->numberOfBladePointsPerTurbine,
                &this->getAllBladeCoordsY()[turbine * this->numberOfBladePointsPerTurbine]);
    std::copy_n(bladeCoordsZ, this->numberOfBladePointsPerTurbine,
                &this->getAllBladeCoordsZ()[turbine * this->numberOfBladePointsPerTurbine]);
}

void ActuatorFarm::setTurbineBladeVelocities(size_t turbine, const real* bladeVelocitiesX, const real* bladeVelocitiesY,
                                             const real* bladeVelocitiesZ) const
{
    std::copy_n(bladeVelocitiesX, this->numberOfBladePointsPerTurbine,
                &this->getAllBladeVelocitiesX()[turbine * this->numberOfBladePointsPerTurbine]);
    std::copy_n(bladeVelocitiesY, this->numberOfBladePointsPerTurbine,
                &this->getAllBladeVelocitiesY()[turbine * this->numberOfBladePointsPerTurbine]);
    std::copy_n(bladeVelocitiesZ, this->numberOfBladePointsPerTurbine,
                &this->getAllBladeVelocitiesZ()[turbine * this->numberOfBladePointsPerTurbine]);
}

void ActuatorFarm::setTurbineBladeForces(size_t turbine, const real* bladeForcesX, const real* bladeForcesY,
                                         const real* bladeForcesZ) const
{
    std::copy_n(bladeForcesX, this->numberOfBladePointsPerTurbine,
                &this->getAllBladeForcesX()[turbine * this->numberOfBladePointsPerTurbine]);
    std::copy_n(bladeForcesY, this->numberOfBladePointsPerTurbine,
                &this->getAllBladeForcesY()[turbine * this->numberOfBladePointsPerTurbine]);
    std::copy_n(bladeForcesZ, this->numberOfBladePointsPerTurbine,
                &this->getAllBladeForcesZ()[turbine * this->numberOfBladePointsPerTurbine]);
}

void ActuatorFarm::swapDeviceArrays()
{
    std::swap(this->coordsXDPreviousTimestep, this->coordsXDCurrentTimestep);
    std::swap(this->coordsYDPreviousTimestep, this->coordsYDCurrentTimestep);
    std::swap(this->coordsZDPreviousTimestep, this->coordsZDCurrentTimestep);

    std::swap(this->velocitiesXDPreviousTimestep, this->velocitiesXDCurrentTimestep);
    std::swap(this->velocitiesYDPreviousTimestep, this->velocitiesYDCurrentTimestep);
    std::swap(this->velocitiesZDPreviousTimestep, this->velocitiesZDCurrentTimestep);

    std::swap(this->forcesXDPreviousTimestep, this->forcesXDCurrentTimestep);
    std::swap(this->forcesYDPreviousTimestep, this->forcesYDCurrentTimestep);
    std::swap(this->forcesZDPreviousTimestep, this->forcesZDCurrentTimestep);
    if (this->requiresLocalSmearingWidth())
        std::swap(this->localSmearingWidthDPreviousTimestep, this->localSmearingWidthDCurrentTimestep);
}

void ActuatorFarm::generateHubAxisPoints(uint turbineIndex, uint& pointIndex)
{
    const real centerX = turbinePosXH[turbineIndex] - hubPositionOffset;
    const real centerY = turbinePosYH[turbineIndex];
    const real centerZ = turbinePosZH[turbineIndex];

    const real segmentLength = hubLength / static_cast<real>(numberOfHubPointsPerTurbine);

    for (uint i = 0; i < numberOfHubPointsPerTurbine; i++) {
        getAllHubCoordsX()[pointIndex] = centerX + static_cast<real>(i) * segmentLength;
        getAllHubCoordsY()[pointIndex] = centerY;
        getAllHubCoordsZ()[pointIndex] = centerZ;
        pointIndex++;
    }
}

void ActuatorFarm::generateTowerAxisPoints(uint turbineIndex, uint& pointIndex)
{
    const real centerX = turbinePosXH[turbineIndex] + towerOffset;
    const real centerY = turbinePosYH[turbineIndex];
    const real towerTop = turbinePosZH[turbineIndex] - hubRadius;
    const real towerHeight = this->towerHeights[turbineIndex];

    const real segmentHeight = towerHeight / static_cast<real>(numberOfTowerPointsPerTurbine);

    for (uint i = 0; i < numberOfTowerPointsPerTurbine; i++) {
        getAllTowerCoordsX()[pointIndex] = centerX;
        getAllTowerCoordsY()[pointIndex] = centerY;
        getAllTowerCoordsZ()[pointIndex] = towerTop - (static_cast<real>(i) + c1o2) * segmentHeight;
        pointIndex++;
    }
}

std::string ActuatorFarm::getFilename(uint t) const
{
    return para->getOutputPath() + this->outputName + "_ID_" + std::to_string(para->getMyProcessID()) + "_t_" +
           std::to_string(t);
}

void ActuatorFarm::write(const std::string& filename) const
{
    const uint totalPoints = getTotalNumberOfPoints();

    std::vector<std::string> dataNames = { "VelocitiesX", "VelocitiesY", "VelocitiesZ", "ForcesX", "ForcesY",
                                           "ForcesZ",     "Blade",       "Hub",         "Tower" };
    std::vector<UbTupleFloat3> nodes(totalPoints);
    std::vector<std::vector<double>> nodeData(dataNames.size());
    for (auto& data : nodeData)
        data.resize(totalPoints, 0.0);

    for (uint i = 0; i < totalPoints; i++) {
        nodes[i] = UbTupleFloat3(this->coordsXH[i], this->coordsYH[i], this->coordsZH[i]);
        nodeData[0][i] = this->velocitiesXH[i];
        nodeData[1][i] = this->velocitiesYH[i];
        nodeData[2][i] = this->velocitiesZH[i];
        nodeData[3][i] = this->forcesXH[i];
        nodeData[4][i] = this->forcesYH[i];
        nodeData[5][i] = this->forcesZH[i];
    }

    // Boolean flags for Paraview Threshold filtering:
    // Select scalar "Blade"/"Hub"/"Tower" and threshold [1, 1]
    for (uint i = 0; i < numberOfBladePoints; i++)
        nodeData[6][i] = 1.0;
    for (uint i = 0; i < numberOfHubPoints; i++)
        nodeData[7][numberOfBladePoints + i] = 1.0;
    for (uint i = 0; i < numberOfTowerPoints; i++)
        nodeData[8][numberOfBladePoints + numberOfHubPoints + i] = 1.0;

    this->appendOutputData(dataNames, nodeData);

    WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(filename, nodes, dataNames, nodeData);
}

} // namespace vf::gpu

//! \}
