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
#include "PlanarAverageProbe.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>

#include <helper_cuda.h>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/reduce.h>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <basics/utilities/UbTuple.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/Output/FilePartCalculator.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/PostProcessor/MacroscopicQuantities.cuh"
#include "gpu/core/Samplers/Utilities.h"
#include "gpu/core/Utilities/KernelUtilities.h"
#include "gpu/cuda_helper/CudaGrid.h"
#include "gpu/cuda_helper/CudaIndexCalculation.h"

using namespace vf::basics::constant;

using valIterator = thrust::device_vector<real>::iterator;
using indIterator = thrust::device_vector<unsigned long long>::iterator;
using permIterator = thrust::permutation_iterator<valIterator, indIterator>;
using iterPair = std::pair<permIterator, permIterator>;

iterPair getPermutationIterators(real* values, unsigned long long* indices, uint numberOfIndices)
{
    auto val_pointer = thrust::device_pointer_cast(values);
    auto indices_pointer = thrust::device_pointer_cast(indices);
    permIterator iter_begin(val_pointer, indices_pointer);
    permIterator iter_end(val_pointer, indices_pointer + numberOfIndices);
    return std::make_pair(iter_begin, iter_end);
}

struct shiftIndex
{
    const uint* neighborNormal;
    __host__ __device__ unsigned long long operator()(unsigned long long& pointIndices)
    {
        return neighborNormal[pointIndices];
    }
};

struct covariance
{
    const real mean_x, mean_y;
    covariance(real mean_x, real mean_y) : mean_x(mean_x), mean_y(mean_y)
    {
    }

    template <typename Tuple>
    __host__ __device__ real operator()(const Tuple& t) const
    {
        return (thrust::get<0>(t) - mean_x) * (thrust::get<1>(t) - mean_y);
    }
};

struct skewness
{
    const real mean;
    skewness(real mean) : mean(mean)
    {
    }
    __host__ __device__ real operator()(const real& x) const
    {
        return (x - mean) * (x - mean) * (x - mean);
    }
};

struct flatness
{
    const real mean;
    flatness(real mean) : mean(mean)
    {
    }
    __host__ __device__ real operator()(const real& x) const
    {
        return (x - mean) * (x - mean) * (x - mean) * (x - mean);
    }
};

struct Means
{
    real vx, vy, vz;
};

struct Covariances
{
    real vxvx, vyvy, vzvz, vxvy, vxvz, vyvz;
};

struct Skewnesses
{
    real Sx, Sy, Sz;
};

struct Flatnesses
{
    real Fx, Fy, Fz;
};

///////////////////////////////////////////////////////////////////////////////////

__global__ void moveIndicesInPosNormalDir(unsigned long long* pointIndices, uint nPoints, const uint* neighborNormal)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= nPoints)
        return;

    pointIndices[nodeIndex] = neighborNormal[pointIndices[nodeIndex]];
}

std::string getName(std::string quantity, bool timeAverage)
{
    if (timeAverage) {
        return quantity + "_spatioTemporalMean";
    }
    return quantity + "_spatialMean";
}

///////////////////////////////////////////////////////////////////////////////////
std::vector<std::string> PlanarAverageProbe::getVariableNames(PlanarAverageProbe::Statistic statistic,
                                                              bool namesForTimeAverages)
{

    std::vector<std::string> variableNames;
    switch (statistic) {
        case Statistic::Means:
            variableNames.emplace_back(getName("vx", namesForTimeAverages));
            variableNames.emplace_back(getName("vy", namesForTimeAverages));
            variableNames.emplace_back(getName("vz", namesForTimeAverages));
            if (para->getUseTurbulentViscosity()) {
                variableNames.emplace_back(getName("EddyViscosity", namesForTimeAverages));
            }
            break;
        case Statistic::Covariances:
            variableNames.emplace_back(getName("vxvx", namesForTimeAverages));
            variableNames.emplace_back(getName("vyvy", namesForTimeAverages));
            variableNames.emplace_back(getName("vzvz", namesForTimeAverages));
            variableNames.emplace_back(getName("vxvy", namesForTimeAverages));
            variableNames.emplace_back(getName("vxvz", namesForTimeAverages));
            variableNames.emplace_back(getName("vyvz", namesForTimeAverages));
            break;
        case Statistic::Skewness:
            variableNames.emplace_back(getName("SkewnessX", namesForTimeAverages));
            variableNames.emplace_back(getName("SkewnessY", namesForTimeAverages));
            variableNames.emplace_back(getName("SkewnessZ", namesForTimeAverages));
            break;
        case Statistic::Flatness:
            variableNames.emplace_back(getName("FlatnessX", namesForTimeAverages));
            variableNames.emplace_back(getName("FlatnessY", namesForTimeAverages));
            variableNames.emplace_back(getName("FlatnessZ", namesForTimeAverages));
            break;

        default:
            throw std::runtime_error("PlanarAverageProbe::getVariableNames: Statistic unavailable!");
            break;
    }

    return variableNames;
}

void PlanarAverageProbe::addStatistic(PlanarAverageProbe::Statistic statistic)
{
    if (!isStatisticIn(statistic, statistics))
        statistics.push_back(statistic);
}

///////////////////////////////////////////////////////////////////////////////////
void PlanarAverageProbe::init()
{
    size_t numberOfVariables = 0;
    for (auto statistic : statistics) {
        numberOfVariables += getVariableNames(statistic, false).size();
    }

    for (int level = 0; level <= para->getMaxLevel(); level++) {
        auto data = &levelData.emplace_back();
        std::vector<unsigned long long> indices = findIndicesInPlane(level);
        findCoordinatesForPlanes(level, data->coordinateX, data->coordinateY, data->coordinateZ);
        data->numberOfPlanes = static_cast<uint>(data->coordinateX.size());
        data->numberOfPointsPerPlane = static_cast<uint>(indices.size());
        cudaMemoryManager->cudaAllocPlanarAverageProbeIndices(this, level);
        std::copy(indices.begin(), indices.end(), data->indicesOfFirstPlaneH);
        cudaMemoryManager->cudaCopyPlanarAverageProbeIndicesHtoD(this, level);
        data->instantaneous.resize(data->numberOfPlanes, std::vector<real>(numberOfVariables, c0o1));
        if (computeTimeAverages) {
            data->timeAverages.resize(data->numberOfPlanes, std::vector<real>(numberOfVariables, c0o1));
        }
    }
}

void PlanarAverageProbe::sample(int level, uint t)
{
    if (t < tStartAveraging)
        return;

    const uint tLevel = para->getTimeStep(level, t, true);
    const uint levelFactor = exp2(level);
    const uint tAfterStartAveraging = tLevel - tStartAveraging * levelFactor;
    const uint tAfterStartWriting = tLevel - tStartWritingOutput * levelFactor;
    const bool doTimeAverages = computeTimeAverages && tLevel > (tStartTemporalAveraging * levelFactor);

    if (tAfterStartAveraging % (tBetweenAverages * levelFactor) == 0) {
        calculateQuantities(level, doTimeAverages);
    }

    if (tAfterStartWriting > 0 && tAfterStartWriting % (tBetweenWriting * levelFactor) == 0) {
        const int tWrite = this->nameFilesWithFileCount ? (t - tStartWritingOutput) / this->tBetweenWriting : t;

        writeGridFile(level, tWrite);
        if (level == 0)
            writeParallelFile(tWrite);
    }
}

PlanarAverageProbe::~PlanarAverageProbe()
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        cudaMemoryManager->cudaFreePlanarAverageProbeIndices(this, level);
    }
}

std::vector<unsigned long long> PlanarAverageProbe::findIndicesInPlane(int level)
{
    std::vector<unsigned long long> indices;
    auto param = para->getParH(level);

    const real* coordinatesPlaneNormal = [&] {
        switch (planeNormal) {
            case PlaneNormal::x:
                return param->coordinateX;
            case PlaneNormal::y:
                return param->coordinateY;
            case PlaneNormal::z:
                return param->coordinateZ;
            default:
                throw std::runtime_error("PlaneNormal not defined!");
        }
    }();

    const auto firstIndex = param->neighborZ[param->neighborY[param->neighborX[1]]];

    const real coordFirstPlane = coordinatesPlaneNormal[firstIndex];

    for (unsigned long long node = 1; node < param->numberOfNodes; node++) {
        if (coordinatesPlaneNormal[node] == coordFirstPlane && param->typeOfGridNode[node] == GEO_FLUID)
            indices.push_back(node);
    }

    return indices;
}

void PlanarAverageProbe::findCoordinatesForPlanes(int level, std::vector<real>& coordinateX, std::vector<real>& coordinateY,
                                                  std::vector<real>& coordinateZ)
{
    const unsigned long long startIndex =
        para->getParH(level)->neighborZ[para->getParH(level)->neighborY[para->getParH(level)->neighborX[1]]];
    unsigned long long nextIndex = startIndex;
    do {
        switch (planeNormal) {
            case PlaneNormal::x:
                coordinateX.push_back(para->getParH(level)->coordinateX[nextIndex]);
                coordinateY.push_back(999999);
                coordinateZ.push_back(999999);
                nextIndex = para->getParH(level)->neighborX[nextIndex];
                break;
            case PlaneNormal::y:
                coordinateX.push_back(999999);
                coordinateY.push_back(para->getParH(level)->coordinateY[nextIndex]);
                coordinateZ.push_back(999999);
                nextIndex = para->getParH(level)->neighborY[nextIndex];
                break;
            case PlaneNormal::z:
                coordinateX.push_back(999999);
                coordinateY.push_back(999999);
                coordinateZ.push_back(para->getParH(level)->coordinateZ[nextIndex]);
                nextIndex = para->getParH(level)->neighborZ[nextIndex];
                break;
            default:
                throw std::runtime_error("PlaneNormal not defined!");
        }
    } while (GEO_FLUID == para->getParH(level)->typeOfGridNode[nextIndex] && startIndex != nextIndex);
}

real computeMean(iterPair x, real invNPointsPerPlane)
{
    const real sum = thrust::reduce(std::get<0>(x), std::get<1>(x));
    return sum * invNPointsPerPlane;
}

Means computeMeans(iterPair vx, iterPair vy, iterPair vz, real invNPointsPerPlane)
{
    Means means;
    means.vx = computeMean(vx, invNPointsPerPlane);
    means.vy = computeMean(vy, invNPointsPerPlane);
    means.vz = computeMean(vz, invNPointsPerPlane);
    return means;
}

real computeCovariance(iterPair x, iterPair y, real mean_x, real mean_y, real invNPointsPerPlane)
{
    auto begin = thrust::make_zip_iterator(thrust::make_tuple(x.first, y.first));
    auto end = thrust::make_zip_iterator(thrust::make_tuple(x.second, y.second));
    return thrust::transform_reduce(begin, end, covariance(mean_x, mean_y), c0o1, thrust::plus<real>()) * invNPointsPerPlane;
}

Covariances computeCovariances(iterPair vx, iterPair vy, iterPair vz, Means means, real invNPointsPerPlane)
{
    Covariances covariances;

    covariances.vxvx = computeCovariance(vx, vx, means.vx, means.vx, invNPointsPerPlane);
    covariances.vyvy = computeCovariance(vy, vy, means.vy, means.vy, invNPointsPerPlane);
    covariances.vzvz = computeCovariance(vz, vz, means.vz, means.vz, invNPointsPerPlane);
    covariances.vxvy = computeCovariance(vx, vy, means.vx, means.vy, invNPointsPerPlane);
    covariances.vxvz = computeCovariance(vx, vz, means.vx, means.vz, invNPointsPerPlane);
    covariances.vyvz = computeCovariance(vy, vz, means.vy, means.vz, invNPointsPerPlane);

    return covariances;
}

real computeSkewness(iterPair x, real mean, real covariance, real invNPointsPerPlane)
{
    return thrust::transform_reduce(x.first, x.second, skewness(mean), c0o1, thrust::plus<real>()) * invNPointsPerPlane *
           pow(covariance, -1.5f);
}

Skewnesses computeSkewnesses(Means means, Covariances covariances, iterPair vx, iterPair vy, iterPair vz,
                             real invNPointsPerPlane)
{
    Skewnesses skewnesses;

    skewnesses.Sx = computeSkewness(vx, means.vx, covariances.vxvx, invNPointsPerPlane);
    skewnesses.Sy = computeSkewness(vy, means.vy, covariances.vyvy, invNPointsPerPlane);
    skewnesses.Sz = computeSkewness(vz, means.vz, covariances.vzvz, invNPointsPerPlane);

    return skewnesses;
}

real computeFlatness(iterPair x, real mean, real covariance, real invNPointsPerPlane)
{
    return thrust::transform_reduce(x.first, x.second, flatness(mean), c0o1, thrust::plus<real>()) * invNPointsPerPlane *
           pow(covariance, -c2o1);
}

Flatnesses computeFlatnesses(iterPair vx, iterPair vy, iterPair vz, Means means, Covariances covariances,
                             real invNPointsPerPlane)
{
    Flatnesses flatnesses;

    flatnesses.Fx = computeFlatness(vx, means.vx, covariances.vxvx, invNPointsPerPlane);
    flatnesses.Fy = computeFlatness(vy, means.vy, covariances.vyvy, invNPointsPerPlane);
    flatnesses.Fz = computeFlatness(vz, means.vz, covariances.vzvz, invNPointsPerPlane);

    return flatnesses;
}

std::vector<real> PlanarAverageProbe::computePlaneStatistics(int level)
{
    auto data = &levelData[level];
    auto parameter = para->getParD(level);
    const real invNPointsPerPlane = c1o1 / static_cast<real>(data->numberOfPointsPerPlane);

    const auto velocityX =
        getPermutationIterators(parameter->velocityX, data->indicesOfFirstPlaneD, data->numberOfPointsPerPlane);
    const auto velocityY =
        getPermutationIterators(parameter->velocityY, data->indicesOfFirstPlaneD, data->numberOfPointsPerPlane);
    const auto velocityZ =
        getPermutationIterators(parameter->velocityZ, data->indicesOfFirstPlaneD, data->numberOfPointsPerPlane);
    const auto turbulentViscosity =
        getPermutationIterators(parameter->turbViscosity, data->indicesOfFirstPlaneD, data->numberOfPointsPerPlane);

    std::vector<real> averages;

    if (!isStatisticIn(PlanarAverageProbe::Statistic::Means, statistics))
        return averages;

    const real velocityRatio = para->getScaledVelocityRatio(level);
    const real viscosityRatio = para->getScaledViscosityRatio(level);
    const real stressRatio = para->getScaledStressRatio(level);

    const auto means = computeMeans(velocityX, velocityY, velocityZ, invNPointsPerPlane);
    averages.push_back(means.vx * velocityRatio);
    averages.push_back(means.vy * velocityRatio);
    averages.push_back(means.vz * velocityRatio);
    if (para->getUseTurbulentViscosity())
        averages.push_back(computeMean(turbulentViscosity, invNPointsPerPlane) * viscosityRatio);

    if (!isStatisticIn(PlanarAverageProbe::Statistic::Covariances, statistics))
        return averages;

    const auto covariances = computeCovariances(velocityX, velocityY, velocityZ, means, invNPointsPerPlane);
    averages.push_back(covariances.vxvx * stressRatio);
    averages.push_back(covariances.vyvy * stressRatio);
    averages.push_back(covariances.vzvz * stressRatio);
    averages.push_back(covariances.vxvy * stressRatio);
    averages.push_back(covariances.vxvz * stressRatio);
    averages.push_back(covariances.vyvz * stressRatio);

    if (!isStatisticIn(PlanarAverageProbe::Statistic::Skewness, statistics))
        return averages;

    const auto skewnesses = computeSkewnesses(means, covariances, velocityX, velocityY, velocityZ, invNPointsPerPlane);
    averages.push_back(skewnesses.Sx);
    averages.push_back(skewnesses.Sy);
    averages.push_back(skewnesses.Sz);

    if (!isStatisticIn(PlanarAverageProbe::Statistic::Flatness, statistics))
        return averages;

    const auto flatnesses = computeFlatnesses(velocityX, velocityY, velocityZ, means, covariances, invNPointsPerPlane);
    averages.push_back(flatnesses.Fx);
    averages.push_back(flatnesses.Fy);
    averages.push_back(flatnesses.Fz);

    return averages;
}

std::vector<real> computeNewTimeAverages(std::vector<real>& oldAverages, std::vector<real>& instantaneous,
                                         real invNumberOfTimesteps)
{
    std::vector<real> newAverages;
    newAverages.reserve(oldAverages.size());
    for (uint i = 0; i < oldAverages.size(); i++) {
        newAverages.push_back(computeNewTimeAverage(oldAverages[i], instantaneous[i], invNumberOfTimesteps));
    }
    return newAverages;
}

void PlanarAverageProbe::calculateQuantities(int level, bool doTimeAverages)
{
    auto data = &levelData[level];
    auto parameter = para->getParD(level);
    calculateMacroscopicQuantitiesCompressible(
        parameter->velocityX, parameter->velocityY, parameter->velocityZ, parameter->rho, parameter->pressure,
        parameter->typeOfGridNode, parameter->neighborX, parameter->neighborY, parameter->neighborZ,
        parameter->numberOfNodes, parameter->numberofthreads, parameter->distributions.f[0], parameter->isEvenTimestep);
    cudaDeviceSynchronize();

    const real invNumberOfTimesteps = c1o1 / static_cast<real>(data->numberOfTimestepsInTimeAverage + 1);

    const auto indices = thrust::device_pointer_cast(data->indicesOfFirstPlaneD);
    const auto shiftOp = shiftIndex { getNeighborIndicesInPlaneNormal(level) };

    for (uint plane = 0; plane < data->numberOfPlanes; plane++) {
        data->instantaneous[plane] = computePlaneStatistics(level);
        if (doTimeAverages)
            data->timeAverages[plane] =
                computeNewTimeAverages(data->timeAverages[plane], data->instantaneous[plane], invNumberOfTimesteps);

        thrust::transform(indices, indices + data->numberOfPointsPerPlane, indices, shiftOp);
    }
    cudaMemoryManager->cudaCopyPlanarAverageProbeIndicesHtoD(this, level);
    getLastCudaError("PlanarAverageProbe::calculateQuantities execution failed");
}

const uint* PlanarAverageProbe::getNeighborIndicesInPlaneNormal(int level)
{
    switch (planeNormal) {
        case PlaneNormal::x:
            return para->getParD(level)->neighborX;
        case PlaneNormal::y:
            return para->getParD(level)->neighborY;
        case PlaneNormal::z:
            return para->getParD(level)->neighborZ;
    };

    throw std::runtime_error("PlaneNormal not defined!");
}

void PlanarAverageProbe::writeParallelFile(uint tWrite)
{
    const std::string filename = this->outputPath + makeParallelFileName(probeName, para->getMyProcessID(), tWrite);

    std::vector<std::string> nodedatanames = this->getAllVariableNames();
    std::vector<std::string> cellNames;

    WbWriterVtkXmlBinary::getInstance()->writeParallelFile(filename, fileNamesForCollectionFile, nodedatanames, cellNames);

    this->fileNamesForCollectionFile.clear();
}

void PlanarAverageProbe::copyDataToNodedata(std::vector<std::vector<real>>& data, std::vector<std::vector<double>>& nodedata)
{
    int arrayIndex = 0;

    for (auto statistic : statistics) {
        for (auto variable : this->getVariableNames(statistic, false)) {
            std::vector<double> array;
            array.reserve(data.size());
            for (auto& plane : data) {
                array.emplace_back(plane[arrayIndex]);
            }
            nodedata.push_back(array);
            arrayIndex++;
        }
    }
}

void PlanarAverageProbe::writeGridFile(int level, uint tWrite)
{
    auto data = &levelData[level];
    const std::string fname = outputPath + makeGridFileName(probeName, level, para->getMyProcessID(), tWrite, 1);

    std::vector<UbTupleFloat3> nodes;
    nodes.reserve(data->numberOfPlanes);
    for (uint pos = 0; pos < data->numberOfPlanes; pos++) {
        nodes.push_back(
            makeUbTuple(float(data->coordinateX[pos]), float(data->coordinateY[pos]), float(data->coordinateZ[pos])));
    }

    std::vector<std::string> nodedatanames = this->getAllVariableNames();

    std::vector<std::vector<double>> nodedata;
    copyDataToNodedata(data->instantaneous, nodedata);
    if (computeTimeAverages)
        copyDataToNodedata(data->timeAverages, nodedata);

    std::string fullName =
        WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(fname, nodes, nodedatanames, nodedata);
    this->fileNamesForCollectionFile.push_back(fullName.substr(fullName.find_last_of('/') + 1));
}

std::vector<std::string> PlanarAverageProbe::getAllVariableNames()
{
    std::vector<std::string> varNames;
    for (auto statistic : statistics) {
        for (auto variable : getVariableNames(statistic, false)) {
            varNames.push_back(variable);
        }
    }
    if (this->computeTimeAverages) {
        for (auto statistic : statistics) {
            for (auto variable : getVariableNames(statistic, true))
                varNames.push_back(variable);
        }
    }
    return varNames;
}

bool isStatisticIn(PlanarAverageProbe::Statistic statistic, std::vector<PlanarAverageProbe::Statistic> statistics)
{
    return std::find(statistics.begin(), statistics.end(), statistic) != statistics.end();
}

//! \}