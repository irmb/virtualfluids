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
//! \addtogroup gpu_Samplers Samplers
//! \ingroup gpu_core core
//! \{
#include "WallModelProbe.h"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/reduce.h>

#include <basics/DataTypes.h>

#include "Cuda/CudaMemoryManager.h"
#include "DataStructureInitializer/GridProvider.h"
#include "Parameter/Parameter.h"
#include "Utilities.h"
#include "basics/constants/NumericConstants.h"

using namespace vf::basics::constant;
using valueIterator = thrust::device_vector<real>::iterator;
using indexIterator = thrust::device_vector<uint>::iterator;

///////////////////////////////////////////////////////////////////////////////////

std::vector<std::string> WallModelProbe::getVariableNames()
{
    std::vector<std::string> variableNames;
    variableNames.emplace_back("velocityNode_spatioTemporalMean");
    variableNames.emplace_back("velocitySample_spatialMean");
    variableNames.emplace_back("frictionVelocity_spatialMean");
    variableNames.emplace_back("momentumExchangeX_spatialMean");
    variableNames.emplace_back("momentumExchangeY_spatialMean");
    variableNames.emplace_back("momentumExchangeZ_spatialMean");
    if (evaluatePressureGradient) {
        variableNames.emplace_back("pressureGradientX_spatialMean");
        variableNames.emplace_back("pressureGradientY_spatialMean");
        variableNames.emplace_back("pressureGradientZ_spatialMean");
    }
    if (sampleSurfaceLayer) {
        variableNames.emplace_back("temperatureNode_spatialMean");
        variableNames.emplace_back("temperatureSample_spatialMean");
        variableNames.emplace_back("frictionTemperature_spatialMean");
        variableNames.emplace_back("surfaceHeatFlux_spatialMean");
        variableNames.emplace_back("surfaceTemperature_spatialMean");
    }
    if (computeTemporalAverages) {
    variableNames.emplace_back("velocityNode_spatioTemporalMean");
    variableNames.emplace_back("velocitySample_spatioTemporalMean");
    variableNames.emplace_back("frictionVelocity_spatioTemporalMean");
    variableNames.emplace_back("momentumExchangeX_spatioTemporalMean");
    variableNames.emplace_back("momentumExchangeY_spatioTemporalMean");
    variableNames.emplace_back("momentumExchangeZ_spatioTemporalMean");
        if (evaluatePressureGradient) {
            variableNames.emplace_back("pressureGradientX_spatioTemporalMean");
            variableNames.emplace_back("pressureGradientY_spatioTemporalMean");
            variableNames.emplace_back("pressureGradientZ_spatioTemporalMean");
        }
        if (sampleSurfaceLayer) {
            variableNames.emplace_back("temperatureNode_spatioTemporalMean");
            variableNames.emplace_back("temperatureSample_spatioTemporalMean");
            variableNames.emplace_back("frictionTemperature_spatioTemporalMean");
            variableNames.emplace_back("surfaceHeatFlux_spatioTemporalMean");
            variableNames.emplace_back("surfaceTemperature_spatioTemporalMean");
        }
    }

    return variableNames;
}

///////////////////////////////////////////////////////////////////////////////////
void WallModelProbe::init()
{

    std::vector<std::string> variableNames = getVariableNames();

    const uint numberOfQuantities = getNumberOfInstantaneousQuantities();

    const real x[1] { 0 };
    const real y[1] { 0 };
    const real z[1] { 0 };

    for (int level = 0; level <= para->getMaxLevel(); level++) {
        const std::string fileName = outputPath + makeTimeseriesFileName(probeName, level, para->getMyProcessID());
        writeTimeseriesFileHeader(fileName, 1, variableNames, x, y, z);
        const uint numberOfFluidNodes = evaluatePressureGradient ? countFluidNodes(level) : 0;
        levelData.emplace_back(fileName, numberOfFluidNodes);
        levelData.back().averagedData.emplace_back(numberOfQuantities, 0);
    }
}

void WallModelProbe::sample(int level, uint t)
{
    const uint tLevel = para->getTimeStep(level, t, false);
    const bool isCoarseTimestep = tLevel % t == 0;

    auto* data = &levelData[level];

    const uint tAfterStartAvg = t - tStartSampling;
    const uint tAfterStartOut = t - tStartWritingOutput;

    if (t >= tStartSampling && ((tAfterStartAvg % tBetweenSamples == 0 && isCoarseTimestep) || sampleEveryTimestep)) {
        calculateQuantities(data, tLevel, level);
    }

    if (t >= tStartWritingOutput && isCoarseTimestep && tAfterStartOut % tBetweenWriting == 0) {
        write(level);
    }
}

///////////////////////////////////////////////////////////////////////////////////

template <typename T>
T computeMean(T* devicePointer, uint numberOfPoints, T conversionFactor)
{
    thrust::device_ptr<T> thrustPointer = thrust::device_pointer_cast(devicePointer);
    return thrust::reduce(thrustPointer, thrustPointer + numberOfPoints) / real(numberOfPoints) * conversionFactor;
}

struct isValidNode
{
    __host__ __device__ real operator()(thrust::tuple<real, uint> x)
    {
        return thrust::get<1>(x) == GEO_FLUID ? thrust::get<0>(x) : c0o1;
    }
};

template <typename T>
T computeIndexBasedMean(T* devicePointer, uint* typeOfGridNode, uint numberOfNodes, uint numberOfFluidNodes,
                        T conversionFactor)
{
    thrust::device_ptr<T> thrustPointer = thrust::device_pointer_cast(devicePointer);
    thrust::device_ptr<uint> typePointer = thrust::device_pointer_cast(typeOfGridNode);
    auto begin = thrust::make_zip_iterator(thrust::make_tuple(thrustPointer, typePointer));
    auto end = thrust::make_zip_iterator(thrust::make_tuple(thrustPointer + numberOfNodes, typePointer + numberOfNodes));
    auto iter_begin = thrust::make_transform_iterator(begin, isValidNode());
    auto iter_end = thrust::make_transform_iterator(end, isValidNode());

    return thrust::reduce(iter_begin, iter_end) / real(numberOfFluidNodes) * conversionFactor;
}

template <typename T>
void computeAndSaveMean(T* devicePointer, uint numberOfPoints, std::vector<T>& quantityArray, T conversionFactor)
{
    quantityArray.push_back(computeMean(devicePointer, numberOfPoints, conversionFactor));
}

template <typename T>
void computeAndSaveIndexBasedMean(T* devicePointer, uint* typeOfGridNode, uint numberOfNodes, uint numberOfFluidNodes,
                                  std::vector<real>& quantitiesArray, T conversionFactor)
{
    quantitiesArray.push_back(
        computeIndexBasedMean(devicePointer, typeOfGridNode, numberOfNodes, numberOfFluidNodes, conversionFactor));
}

template <typename T>
void computeTemporalAverage(std::vector<T>& quantityArray, T oldMean, T currentValue, real invNumberOfAverages)
{
    quantityArray.push_back(computeNewTimeAverage(oldMean, currentValue, invNumberOfAverages));
}

uint WallModelProbe::countFluidNodes(int level)
{
    uint* typePointer = para->getParH(level)->typeOfGridNode;
    return std::count(typePointer, typePointer + para->getParH(level)->numberOfNodes, GEO_FLUID);
}

void WallModelProbe::calculateQuantities(WallModelProbe::LevelData* data, uint t, int level)
{

    auto parD = para->getParD(level);
    const uint nPoints = sampleSurfaceLayer ? parD->surfaceLayerBC.numberOfBCnodes : parD->stressBC.numberOfBCnodes;
    if (nPoints < 1)
        return;
    const uint numberOfQuantities = getNumberOfInstantaneousQuantities();

    data->timestepTime.push_back(t * para->getScaledTimeRatio(level));

    std::vector<real>& newInstantaneous =  data->instantaneousData.emplace_back();
    newInstantaneous.reserve(numberOfQuantities);

    const real velocityFactor = para->getScaledVelocityRatio(level);
    const real stressFactor = para->getScaledStressRatio(level);
    const real forceFactor = para->getScaledForceRatio(level);

    const auto& momentumWallModel =
    sampleSurfaceLayer ? parD->surfaceLayerWallModel.momentumParameters : parD->momentumWallModel;


    computeAndSaveMean(momentumWallModel.velocityMagnitudeNode, nPoints, newInstantaneous, velocityFactor);
    computeAndSaveMean(momentumWallModel.velocityMagnitudeSample, nPoints, newInstantaneous, velocityFactor);

    computeAndSaveMean(momentumWallModel.frictionVelocity, nPoints, newInstantaneous, velocityFactor);

    computeAndSaveMean(momentumWallModel.forceX, nPoints, newInstantaneous, outputStress ? stressFactor : forceFactor);
    computeAndSaveMean(momentumWallModel.forceY, nPoints, newInstantaneous, outputStress ? stressFactor : forceFactor);
    computeAndSaveMean(momentumWallModel.forceZ, nPoints, newInstantaneous, outputStress ? stressFactor : forceFactor);

    if (this->evaluatePressureGradient) {
        computeAndSaveIndexBasedMean(parD->forceX_SP, parD->typeOfGridNode, parD->numberOfNodes,
                                     data->numberOfFluidNodes, newInstantaneous, forceFactor);
        computeAndSaveIndexBasedMean(parD->forceY_SP, parD->typeOfGridNode, parD->numberOfNodes,
                                     data->numberOfFluidNodes, newInstantaneous, forceFactor);
        computeAndSaveIndexBasedMean(parD->forceZ_SP, parD->typeOfGridNode, parD->numberOfNodes,
                                     data->numberOfFluidNodes, newInstantaneous, forceFactor);
    }

    if (sampleSurfaceLayer) {
        const auto& temperatureWallModel = parD->surfaceLayerWallModel.temperatureParameters;
        computeAndSaveMean(temperatureWallModel.temperatureSample, nPoints, newInstantaneous, c1o1);
        computeAndSaveMean(temperatureWallModel.temperatureNode, nPoints, newInstantaneous, c1o1);
        computeAndSaveMean(temperatureWallModel.temperatureScale, nPoints, newInstantaneous, c1o1);
        computeAndSaveMean(temperatureWallModel.surfaceHeatFlux, nPoints, newInstantaneous, velocityFactor);
        computeAndSaveMean(temperatureWallModel.surfaceTemperature, nPoints, newInstantaneous, c1o1);
    }

    if (computeTemporalAverages) {
        if (t > tStartTemporalAveraging) {
            const real inverseNumberOfAveragedValues = c1o1 / real(data->numberOfAveragedValues + 1);
            std::vector<real>& oldAverages = data->averagedData.back();
            std::vector<real> newAverages;
            newAverages.reserve(numberOfQuantities);
            for (uint i = 0; i < numberOfQuantities; i++)
                computeTemporalAverage(newAverages, oldAverages[i], newInstantaneous[i], inverseNumberOfAveragedValues);
            data->averagedData.push_back(newAverages);
            data->numberOfAveragedValues++;
        } else {
            std::vector newAverages(numberOfQuantities, c0o1);
            data->averagedData.push_back(newAverages);
        }
    }
}

void WallModelProbe::write(int level)
{
    auto* data = &levelData[level];

    if (data->timestepTime.empty())
        return;

    std::vector<std::vector<real>> dataToWrite;
    for (size_t i = 0; i < data->timestepTime.size(); i++) {
        std::vector<real> row;
        row.reserve(data->instantaneousData[i].size() + (computeTemporalAverages ? data->averagedData[i + 1].size() : 0) +
                    1);
        row.push_back(data->timestepTime[i]);
        std::copy(data->instantaneousData[i].begin(), data->instantaneousData[i].end(), std::back_inserter(row));
        if (computeTemporalAverages)
            std::copy(data->averagedData[i + 1].begin(), data->averagedData[i + 1].end(), std::back_inserter(row));
        dataToWrite.push_back(row);
    }

    appendDataToTimeseriesFile(data->timeseriesFileName, dataToWrite);

    data->timestepTime.clear();
    data->instantaneousData.clear();
    if (computeTemporalAverages) {
        auto lastTimestep = data->averagedData.back();
        data->averagedData.clear();
        data->averagedData.push_back(lastTimestep);
    }
}

//! \}
