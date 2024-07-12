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
//! \addtogroup gpu_Samplers
//! \ingroup gpu_core core
//! \{
//! \author Henry Korb, Henrik Asmuth
//=======================================================================================

#include "Probe.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

#include "gpu/core/Calculation/Calculation.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/Output/FilePartCalculator.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Utilities/GeometryUtils.h"
#include "gpu/cuda_helper/CudaGrid.h"
#include "gpu/cuda_helper/CudaIndexCalculation.h"

#include "Utilities.h"

using namespace vf::basics::constant;

__host__ __device__ int calcArrayIndex(int node, int nNodes, int timestep, int nTimesteps, int quantity)
{
    return node + nNodes * (timestep + nTimesteps * quantity);
}

uint calcOldTimestep(uint currentTimestep, uint lastTimestepInOldSeries)
{
    return currentTimestep > 0 ? currentTimestep - 1 : lastTimestepInOldSeries;
}

real* getStatisticArray(Probe::ProbeData probeData, Probe::Statistic statistic)
{
    switch (statistic) {
        case Probe::Statistic::Instantaneous:
            return probeData.instantaneous;
        case Probe::Statistic::Means:
            return probeData.means;
        case Probe::Statistic::Variances:
            return probeData.variances;
        default:
            throw std::runtime_error("getStatisticArray: Statistic unavailable!");
    }
}

__host__ __device__ real computeMean(real oldMean, real newValue, real inverseCount)
{
    return oldMean + (newValue - oldMean) * inverseCount;
}

__host__ __device__ real computeAndSaveMean(real* quantityArray, real oldValue, uint index, real currentValue, real invCount)
{
    const real newValue = computeMean(oldValue, currentValue, invCount);
    quantityArray[index] = newValue;
    return newValue;
}

__host__ __device__ real computeVariance(real oldVariance, real oldMean, real newMean, real currentValue,
                                         uint numberOfAveragedValues, real inverseCount)
{
    return (numberOfAveragedValues * oldVariance + (currentValue - oldMean) * (currentValue - newMean)) * inverseCount;
}

__host__ __device__ real computeAndSaveVariance(real* quantityArray, real oldVariance, uint indexNew, real currentValue,
                                                real oldMean, real newMean, uint numberOfAveragedValues, real inverseCount)
{
    const real newVariance =
        computeVariance(oldVariance, oldMean, newMean, currentValue, numberOfAveragedValues, inverseCount);
    quantityArray[indexNew] = newVariance;
    return newVariance;
}

__forceinline__ __device__ void computeStatistics(uint nAveragedValues, uint currentTimestep, uint lastTimestep,
                                                  uint nPoints, uint nodeIndex, bool computeInstant, bool computeMean,
                                                  bool computeVariance, uint iQuantity, real currentValue,
                                                  const Probe::ProbeData& probeData, real invCount)
{
    const uint indexCurrent = calcArrayIndex(nodeIndex, nPoints, currentTimestep, probeData.numberOfTimesteps, iQuantity);
    const uint indexLast = calcArrayIndex(nodeIndex, nPoints, lastTimestep, probeData.numberOfTimesteps, iQuantity);
    if (computeInstant)
        probeData.instantaneous[indexCurrent] = currentValue;
    if (computeMean) {
        const real meanOld = probeData.means[indexLast];
        const real meanNew = computeAndSaveMean(probeData.means, meanOld, indexCurrent, currentValue, invCount);
        if (nAveragedValues == 0)
            return;
        if (computeVariance) {
            const real varianceOld = probeData.variances[indexLast];
            computeAndSaveVariance(probeData.variances, varianceOld, indexCurrent, currentValue, meanOld, meanNew,
                                   nAveragedValues, invCount);
        }
    }
}

__global__ void calculateQuantitiesKernel(uint numberOfAveragedValues, Probe::GridParams gridParams,
                                          Probe::ProbeData probeData, uint currentTimestep, uint lastTimestep)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= probeData.numberOfPoints)
        return;

    const uint gridNodeIndex = probeData.indices[nodeIndex];

    const real invCount = c1o1 / real(numberOfAveragedValues + 1);
    const real nPoints = probeData.numberOfPoints;

    computeStatistics(numberOfAveragedValues, currentTimestep, lastTimestep, nPoints, nodeIndex,
                      probeData.computeInstantaneous, probeData.computeMeans, probeData.computeVariances, 0,
                      gridParams.velocityX[gridNodeIndex], probeData, invCount);
    computeStatistics(numberOfAveragedValues, currentTimestep, lastTimestep, nPoints, nodeIndex,
                      probeData.computeInstantaneous, probeData.computeMeans, probeData.computeVariances, 1,
                      gridParams.velocityY[gridNodeIndex], probeData, invCount);
    computeStatistics(numberOfAveragedValues, currentTimestep, lastTimestep, nPoints, nodeIndex,
                      probeData.computeInstantaneous, probeData.computeMeans, probeData.computeVariances, 2,
                      gridParams.velocityZ[gridNodeIndex], probeData, invCount);
    computeStatistics(numberOfAveragedValues, currentTimestep, lastTimestep, nPoints, nodeIndex,
                      probeData.computeInstantaneous, probeData.computeMeans, probeData.computeVariances, 3,
                      gridParams.density[gridNodeIndex], probeData, invCount);
}

std::vector<Probe::PostProcessingVariable> Probe::getPostProcessingVariables(Statistic statistic, int level) const
{
    const real velocityRatio = para->getScaledVelocityRatio(level);
    const real stressRatio = para->getScaledStressRatio(level);
    const real densityRatio = para->getScaledDensityRatio(level);
    std::vector<PostProcessingVariable> postProcessingVariables;
    switch (statistic) {
        case Statistic::Instantaneous:
            postProcessingVariables.emplace_back("vx", velocityRatio);
            postProcessingVariables.emplace_back("vy", velocityRatio);
            postProcessingVariables.emplace_back("vz", velocityRatio);
            postProcessingVariables.emplace_back("rho", densityRatio);
            break;
        case Statistic::Means:
            postProcessingVariables.emplace_back("vx_mean", velocityRatio);
            postProcessingVariables.emplace_back("vy_mean", velocityRatio);
            postProcessingVariables.emplace_back("vz_mean", velocityRatio);
            postProcessingVariables.emplace_back("rho_mean", densityRatio);
            break;
        case Statistic::Variances:
            postProcessingVariables.emplace_back("vx_var", stressRatio);
            postProcessingVariables.emplace_back("vy_var", stressRatio);
            postProcessingVariables.emplace_back("vz_var", stressRatio);
            postProcessingVariables.emplace_back("rho_var", densityRatio);
            break;

        default:
            throw std::runtime_error("Probe::getPostProcessingVariables: Statistic unavailable!");
            break;
    }
    return postProcessingVariables;
}

std::vector<Probe::PostProcessingVariable> Probe::getAllPostProcessingVariables(int level) const
{
    std::vector<PostProcessingVariable> postProcessingVariables;
    if (enableComputationInstantaneous) {
        for (auto p : getPostProcessingVariables(Statistic::Instantaneous, level))
            postProcessingVariables.push_back(p);
    }
    if (enableComputationMeans) {
        for (auto p : getPostProcessingVariables(Statistic::Means, level))
            postProcessingVariables.push_back(p);
    }
    if (enableComputationVariances) {
        for (auto p : getPostProcessingVariables(Statistic::Variances, level))
            postProcessingVariables.push_back(p);
    }
    return postProcessingVariables;
}

void Probe::init()
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        this->addLevelData(level);

        if (this->outputTimeSeries) {
            auto levelData = levelDatas[level];
            const std::string fileName = makeTimeseriesFileName(probeName, level, para->getMyProcessID());
            auto variableNames = getVarNames();
            writeTimeseriesFileHeader(fileName, levelData.probeDataH.numberOfPoints, variableNames,
                                      levelData.coordinatesX.data(), levelData.coordinatesY.data(),
                                      levelData.coordinatesZ.data());
            timeseriesFileNames.push_back(fileName);
        }
    }
}

void Probe::addLevelData(int level)
{
    std::vector<uint> indices;
    std::vector<real> coordinatesX, coordinatesY, coordinatesZ;

    const real* nodeCoordinatesX = para->getParH(level)->coordinateX;
    const real* nodeCoordinatesY = para->getParH(level)->coordinateY;
    const real* nodeCoordinatesZ = para->getParH(level)->coordinateZ;
    const real deltaX = nodeCoordinatesX[para->getParH(level)->neighborX[1]] - nodeCoordinatesX[1];
    for (unsigned long long pos = 1; pos < para->getParH(level)->numberOfNodes; pos++) {
        const real pointCoordX = nodeCoordinatesX[pos];
        const real pointCoordY = nodeCoordinatesY[pos];
        const real pointCoordZ = nodeCoordinatesZ[pos];
        for (auto point : points) {
            const real distX = point.x - pointCoordX;
            const real distY = point.y - pointCoordY;
            const real distZ = point.z - pointCoordZ;
            if (distX <= deltaX && distY <= deltaX && distZ <= deltaX && distX > c0o1 && distY > c0o1 && distZ > c0o1 &&
                isValidProbePoint(pos, para.get(), level)) {
                indices.push_back(static_cast<uint>(pos));
                coordinatesX.push_back(pointCoordX);
                coordinatesY.push_back(pointCoordY);
                coordinatesZ.push_back(pointCoordZ);
                continue;
            }
        }
        for (auto plane : planes) {
            const real distanceX = pointCoordX - plane.startX;
            const real distanceY = pointCoordY - plane.startY;
            const real distanceZ = pointCoordZ - plane.startZ;

            if (distanceX <= plane.length && distanceY <= plane.width && distanceZ <= plane.height && distanceX >= c0o1 &&
                distanceY >= c0o1 && distanceZ >= c0o1 && isValidProbePoint(pos, para.get(), level)) {
                indices.push_back(static_cast<uint>(pos));
                coordinatesX.push_back(pointCoordX);
                coordinatesY.push_back(pointCoordY);
                coordinatesZ.push_back(pointCoordZ);
                continue;
            }
        }
    }

    const uint numberOfQuantities = static_cast<uint>(getPostProcessingVariables(Statistic::Instantaneous, 0).size());

    ProbeData probeData(enableComputationInstantaneous, enableComputationMeans, enableComputationVariances,
                        static_cast<uint>(indices.size()), numberOfQuantities, getNumberOfTimestepsInTimeseries(level));

    const uint sizeData = probeData.numberOfPoints * probeData.numberOfTimesteps * numberOfQuantities;

    levelDatas.emplace_back(probeData, probeData, coordinatesX, coordinatesY, coordinatesZ);

    cudaMemoryManager->cudaAllocProbeData(this, level);

    std::copy(indices.begin(), indices.end(), levelDatas[level].probeDataH.indices);

    if (enableComputationInstantaneous)
        std::fill_n(levelDatas[level].probeDataH.instantaneous, sizeData, c0o1);
    if (enableComputationMeans)
        std::fill_n(levelDatas[level].probeDataH.means, sizeData, c0o1);
    if (enableComputationVariances)
        std::fill_n(levelDatas[level].probeDataH.variances, sizeData, c0o1);

    cudaMemoryManager->cudaCopyProbeDataHtoD(this, level);
}

void Probe::sample(int level, uint t)
{
    auto levelData = &levelDatas[level];
    if (levelData->probeDataH.numberOfPoints == 0)
        return;
    const uint tLevel = para->getTimeStep(level, t, false);

    //! if averageEveryTimestep the probe will be evaluated in every sub-timestep of each respective level
    //! else, the probe will only be evaluated in each synchronous time step tBetweenAverages

    const uint levelFactor = exp2(level);

    const uint tStartAvgLevel = this->tStartAveraging * levelFactor;
    const uint tAfterStartAvg = tLevel - tStartAvgLevel;
    const uint tAvgLevel = this->tBetweenAverages * levelFactor;
    const bool averageThisTimestep = this->averageEveryTimestep || (tAfterStartAvg % tAvgLevel == 0);

    const uint tStartOutLevel = this->tStartWritingOutput * levelFactor;
    const uint tAfterStartOut = tLevel - tStartOutLevel;
    const uint tOutLevel = this->tBetweenWriting * levelFactor;
    const bool outputThisTimestep = tAfterStartOut % tOutLevel == 0;

    auto gridParams = getGridParams(para->getParD(level).get());

    const vf::cuda::CudaGrid grid(para->getParD(level)->numberofthreads, levelData->probeDataD.numberOfPoints);

    if ((t > this->tStartAveraging) && averageThisTimestep) {
        if (outputTimeSeries) {
            const uint lastTimestep = calcOldTimestep(levelData->timeseriesParams.currentTimestep,
                                                      levelData->timeseriesParams.lastTimestepInOldTimeseries);
            const uint currentTimestep = levelData->timeseriesParams.currentTimestep;
            calculateQuantitiesKernel<<<grid.grid, grid.threads>>>(levelData->numberOfAveragedValues, gridParams,
                                                                   levelData->probeDataD, currentTimestep, lastTimestep);
            if (tLevel >= tStartOutLevel)
                levelData->timeseriesParams.currentTimestep++;
        } else {
            calculateQuantitiesKernel<<<grid.grid, grid.threads>>>(levelData->numberOfAveragedValues, gridParams,
                                                                   levelData->probeDataD, 0, 0);
        }
        levelData->numberOfAveragedValues++;
    }
    //! output only in synchronous timesteps
    if ((t > this->tStartWritingOutput) && outputThisTimestep) {
        cudaMemoryManager->cudaCopyProbeDataDtoH(this, level);
        if (outputTimeSeries) {
            this->appendTimeseriesFile(level, t);
            levelData->timeseriesParams.lastTimestepInOldTimeseries =
                levelData->timeseriesParams.currentTimestep > 0 ? levelData->timeseriesParams.currentTimestep - 1 : 0;
            levelData->timeseriesParams.currentTimestep = 0;
        } else {
            this->writeGridFiles(level, t);
            if (level == 0)
                this->writeParallelFile(t);
        }
    }
}

Probe::~Probe()
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        cudaMemoryManager->cudaFreeProbeData(this, level);
    }
}

void Probe::addStatistic(Statistic variable)
{
    switch (variable) {
        case Statistic::Instantaneous:
            enableComputationInstantaneous = true;
            break;
        case Statistic::Means:
            enableComputationMeans = true;
            break;
        case Statistic::Variances:
            enableComputationVariances = true;
            enableComputationMeans = true;
            break;
        default:
            throw std::runtime_error("Probe::addStatistic: Statistic unavailable!");
            break;
    }
}

void Probe::writeGridFiles(int level, int t)
{
    const int tWrite = this->fileNameLU ? t : (t - tStartWritingOutput) / this->tBetweenWriting;

    const uint numberOfParts = levelDatas[level].probeDataH.numberOfPoints / FilePartCalculator::limitOfNodesForVTK + 1;

    std::vector<std::string> fnames;
    for (uint i = 1; i <= numberOfParts; i++) {
        this->writeGridFile(level, tWrite, i);
    }
}

void Probe::writeParallelFile(int t)
{
    const int t_write = this->fileNameLU ? t : t / this->tBetweenWriting;
    const std::string filename = this->outputPath + makeParallelFileName(probeName, para->getMyProcessID(), t_write);

    std::vector<std::string> nodedatanames = this->getVarNames();
    std::vector<std::string> cellNames;

    getWriter()->writeParallelFile(filename, fileNamesForCollectionFile, nodedatanames, cellNames);

    this->fileNamesForCollectionFile.clear();
}

void Probe::appendStatisticToNodeData(Statistic statistic, uint startPos, uint endPos, uint timestep, int level,
                                      std::vector<std::vector<double>>& nodedata)
{
    auto levelData = &levelDatas[level];
    const uint numberOfNodes = levelData->probeDataH.numberOfPoints;
    const real* data = getStatisticArray(levelData->probeDataH, statistic);
    std::vector<PostProcessingVariable> postProcessingVariables = this->getPostProcessingVariables(statistic, level);
    for (uint arr = 0; arr < uint(postProcessingVariables.size()); arr++) {
        std::vector<double> quantityData(numberOfNodes);
        const real coeff = postProcessingVariables[arr].conversionFactor;
        const int startIndex =
            calcArrayIndex(startPos, numberOfNodes, timestep, levelData->probeDataH.numberOfTimesteps, arr);
        for (uint idx = 0; idx < endPos - startPos; idx++) {
            quantityData[idx] = double(data[startIndex + idx] * coeff);
        }
        nodedata.push_back(quantityData);
    }
}

void Probe::writeGridFile(int level, int t, uint part)
{
    const std::string fname = this->outputPath + makeGridFileName(probeName, level, para->getMyProcessID(), t, part);

    std::vector<UbTupleFloat3> nodes;
    std::vector<std::string> nodedatanames = this->getVarNames();

    std::vector<std::vector<double>> nodedata;

    auto levelData = &levelDatas[level];

    const uint startpos = (part - 1) * FilePartCalculator::limitOfNodesForVTK;
    const uint sizeOfNodes =
        std::min(FilePartCalculator::limitOfNodesForVTK, levelData->probeDataH.numberOfPoints - startpos);
    const uint endpos = startpos + sizeOfNodes;

    //////////////////////////////////////////////////////////////////////////
    nodes.resize(sizeOfNodes);

    for (uint pos = startpos; pos < endpos; pos++) {
        nodes[pos - startpos] = makeUbTuple(float(levelData->coordinatesX[pos]), float(levelData->coordinatesY[pos]),
                                            float(levelData->coordinatesZ[pos]));
    }

    if (enableComputationInstantaneous)
        appendStatisticToNodeData(Statistic::Instantaneous, startpos, endpos, 0, level, nodedata);
    if (enableComputationMeans)
        appendStatisticToNodeData(Statistic::Means, startpos, endpos, level, 0, nodedata);
    if (enableComputationVariances)
        appendStatisticToNodeData(Statistic::Variances, startpos, endpos, level, 0, nodedata);
    std::string fullName = getWriter()->writeNodesWithNodeData(fname, nodes, nodedatanames, nodedata);
    this->fileNamesForCollectionFile.push_back(fullName.substr(fullName.find_last_of('/') + 1));
}

void Probe::appendStatisticToTimestepData(int timestep, std::vector<real>& timestepData, Statistic statistic, int level)
{
    std::vector<PostProcessingVariable> variables = this->getPostProcessingVariables(statistic, level);
    auto probeData = levelDatas[level].probeDataH;
    const real* data = getStatisticArray(probeData, statistic);

    for (uint variable = 0; variable < uint(variables.size()); variable++) {
        const real conversionFactor = variables[variable].conversionFactor;
        const uint startIndex = calcArrayIndex(0, probeData.numberOfPoints, timestep, probeData.numberOfTimesteps, variable);

        for (uint point = 0; point < probeData.numberOfPoints; point++) {
            timestepData.push_back(data[startIndex + point] * conversionFactor);
        }
    }
}

std::vector<real> Probe::getTimestepData(real time, int timestep, int level)
{
    std::vector<real> timestepData;
    timestepData.push_back(time);

    if (enableComputationInstantaneous)
        appendStatisticToTimestepData(timestep, timestepData, Statistic::Instantaneous, level);
    if (enableComputationMeans)
        appendStatisticToTimestepData(timestep, timestepData, Statistic::Means, level);
    if (enableComputationVariances)
        appendStatisticToTimestepData(timestep, timestepData, Statistic::Variances, level);
    return timestepData;
}

void Probe::appendTimeseriesFile(int level, int t)
{
    const uint tAvg_level = this->tBetweenAverages == 1 ? this->tBetweenAverages : this->tBetweenAverages * exp2(-level);
    const real deltaT = para->getTimeRatio() * tAvg_level;
    auto levelData = levelDatas[level];

    const real tStart = (t - this->tBetweenWriting) * para->getTimeRatio();

    std::vector<std::vector<real>> timestepData;

    for (uint timestep = 0; timestep < levelData.timeseriesParams.currentTimestep; timestep++) {
        const real time = tStart + timestep * deltaT;
        timestepData.push_back(this->getTimestepData(time, timestep, level));
    }
    appendDataToTimeseriesFile(this->timeseriesFileNames[level], timestepData);
    levelData.timeseriesParams.lastTimestepInOldTimeseries = levelData.timeseriesParams.currentTimestep;
    levelData.timeseriesParams.currentTimestep = 0;
}

std::vector<std::string> Probe::getVarNames()
{
    std::vector<std::string> varNames;
    for (auto variable : getAllPostProcessingVariables(0))
        varNames.push_back(variable.name);
    return varNames;
}

void Probe::getTaggedFluidNodes(GridProvider* gridProvider)
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        auto probeData = levelDatas[level].probeDataH;
        std::vector<uint> probeIndices(probeData.indices, probeData.indices + probeData.numberOfPoints);
        gridProvider->tagFluidNodeIndices(probeIndices, CollisionTemplate::WriteMacroVars, level);
    }
}

Probe::GridParams Probe::getGridParams(LBMSimulationParameter* para)
{
    return {
        para->velocityX,
        para->velocityY,
        para->velocityZ,
        para->rho,
    };
}

bool isCoarseInterpolationCell(unsigned long long pointIndex, Parameter* para, int level)
{
    if (level == para->getMaxLevel())
        return false;
    auto interpolationCells = para->getParH(level)->fineToCoarse;
    for (uint i = 0; i < interpolationCells.numberOfCells; i++) {
        if (interpolationCells.coarseCellIndices[i] == pointIndex) {
            return true;
        }
    }
    return false;
}

bool isFineInterpolationCell(unsigned long long pointIndex, Parameter* para, int level)
{
    if (level == 0)
        return false;
    auto interpolationCells = para->getParH(level - 1)->coarseToFine;
    const uint* neighborX = para->getParH(level)->neighborX;
    const uint* neighborY = para->getParH(level)->neighborY;
    const uint* neighborZ = para->getParH(level)->neighborZ;
    for (uint i = 0; i < interpolationCells.numberOfCells; i++) {
        const uint kMMM = interpolationCells.fineCellIndices[i];
        uint kPMM, kMPM, kMMP, kPPM, kPMP, kMPP, kPPP;
        getNeighborIndicesOfBSW(kMMM, kPMM, kMPM, kMMP, kPPM, kPMP, kMPP, kPPP, neighborX, neighborY, neighborZ);
        if (kMMM == pointIndex || kPMM == pointIndex || kMPM == pointIndex || kMMP == pointIndex || kPPM == pointIndex ||
            kPMP == pointIndex || kMPP == pointIndex || kPPP == pointIndex) {
            return true;
        }
    }
    return false;
}

bool isValidProbePoint(unsigned long long pointIndex, Parameter* para, int level)
{
    return GEO_FLUID == para->getParH(level)->typeOfGridNode[pointIndex] &&
           !isCoarseInterpolationCell(pointIndex, para, level) && !isFineInterpolationCell(pointIndex, para, level);
}
//! \}
