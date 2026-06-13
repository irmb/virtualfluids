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
#include "Calculation/Calculation.h"
#include "PlanarAverageProbe.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/inner_product.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/reduce.h>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <basics/geometry3d/Axis.h>
#include <basics/utilities/UbTuple.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/PostProcessor/MacroscopicQuantities.cuh"
#include "gpu/core/Samplers/Utilities.h"

using namespace vf::basics::constant;
namespace vf::gpu {

using valIterator = thrust::device_vector<real>::iterator;
using indexPointer = thrust::device_ptr<uint>;
using indIterator = thrust::device_vector<uint>::iterator;
using permIterator = thrust::permutation_iterator<valIterator, indIterator>;
using iterPair = std::pair<permIterator, permIterator>;

iterPair getIteratorPair(real* values, indexPointer indices, uint nNodes)
{
    auto val_pointer = thrust::device_pointer_cast(values);
    permIterator iter_begin(val_pointer, indices);
    permIterator iter_end(val_pointer, indices + nNodes);
    return std::make_pair(iter_begin, iter_end);
}

struct covariance
{
    const real mean_x, mean_y;
    covariance(real mean_x, real mean_y) : mean_x(mean_x), mean_y(mean_y)
    {
    }

    template <typename Tuple>
    constexpr real operator()(const Tuple& t) const
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
    constexpr real operator()(const real& x) const
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
    constexpr real operator()(const real& x) const
    {
        return (x - mean) * (x - mean) * (x - mean) * (x - mean);
    }
};

struct Means
{
    // phi is the scalar
    real vx, vy, vz, phi;
};

struct Covariances
{
    real vxvx, vyvy, vzvz, vxvy, vxvz, vyvz, phiphi, vxphi, vyphi, vzphi;
};

struct Skewnesses
{
    real Sx, Sy, Sz, Sphi;
};

struct Flatnesses
{
    real Fx, Fy, Fz, Fphi;
};

std::string getName(const std::string& quantity, bool timeAverage)
{
    if (timeAverage) {
        return quantity + "_spatioTemporalMean";
    }
    return quantity + "_spatialMean";
}

///////////////////////////////////////////////////////////////////////////////////
std::vector<std::string> PlanarAverageProbe::getVariableNames(Statistic statistic, bool namesForTimeAverages) const
{

    std::vector<std::string> variableNames;
    switch (statistic) {
        case Statistic::Means:
            variableNames.emplace_back(getName("vx", namesForTimeAverages));
            variableNames.emplace_back(getName("vy", namesForTimeAverages));
            variableNames.emplace_back(getName("vz", namesForTimeAverages));
            if (para->getUseTurbulentViscosity())
                variableNames.emplace_back(getName("EddyViscosity", namesForTimeAverages));
            if (sampleSubgridScaleFluxes) {
                variableNames.emplace_back(getName("SGSFluxXX", namesForTimeAverages));
                variableNames.emplace_back(getName("SGSFluxXY", namesForTimeAverages));
                variableNames.emplace_back(getName("SGSFluxXZ", namesForTimeAverages));
                variableNames.emplace_back(getName("SGSFluxYY", namesForTimeAverages));
                variableNames.emplace_back(getName("SGSFluxYZ", namesForTimeAverages));
                variableNames.emplace_back(getName("SGSFluxZZ", namesForTimeAverages));
            }
            if (sampleScalar)
                variableNames.emplace_back(getName("phi", namesForTimeAverages));
            if (sampleScalar && para->getUseTurbulentViscosity())
                variableNames.emplace_back(getName("EddyDiffusivity", namesForTimeAverages));
            if (sampleScalar && sampleSubgridScaleFluxes){
                variableNames.emplace_back(getName("SGSFluxPhiX", namesForTimeAverages));
                variableNames.emplace_back(getName("SGSFluxPhiY", namesForTimeAverages));
                variableNames.emplace_back(getName("SGSFluxPhiZ", namesForTimeAverages));
            }
            break;
        case Statistic::Covariances:
            variableNames.emplace_back(getName("vxvx", namesForTimeAverages));
            variableNames.emplace_back(getName("vyvy", namesForTimeAverages));
            variableNames.emplace_back(getName("vzvz", namesForTimeAverages));
            variableNames.emplace_back(getName("vxvy", namesForTimeAverages));
            variableNames.emplace_back(getName("vxvz", namesForTimeAverages));
            variableNames.emplace_back(getName("vyvz", namesForTimeAverages));
            if (sampleScalar) {
                variableNames.emplace_back(getName("phiphi", namesForTimeAverages));
                variableNames.emplace_back(getName("vxphi", namesForTimeAverages));
                variableNames.emplace_back(getName("vyphi", namesForTimeAverages));
                variableNames.emplace_back(getName("vzphi", namesForTimeAverages));
            }
            break;
        case Statistic::Skewness:
            variableNames.emplace_back(getName("SkewnessX", namesForTimeAverages));
            variableNames.emplace_back(getName("SkewnessY", namesForTimeAverages));
            variableNames.emplace_back(getName("SkewnessZ", namesForTimeAverages));
            if (sampleScalar)
                variableNames.emplace_back(getName("SkewnessPhi", namesForTimeAverages));
            break;
        case Statistic::Flatness:
            variableNames.emplace_back(getName("FlatnessX", namesForTimeAverages));
            variableNames.emplace_back(getName("FlatnessY", namesForTimeAverages));
            variableNames.emplace_back(getName("FlatnessZ", namesForTimeAverages));
            if (sampleScalar)
                variableNames.emplace_back(getName("FlatnessPhi", namesForTimeAverages));
            break;

        default:
            throw std::runtime_error("PlanarAverageProbe::getVariableNames: Statistic unavailable!");
            break;
    }

    return variableNames;
}

void PlanarAverageProbe::addStatistic(Statistic statistic)
{
    if (!isStatisticIn(statistic, statistics))
        statistics.push_back(statistic);
}

///////////////////////////////////////////////////////////////////////////////////
void PlanarAverageProbe::init()
{
    if (sampleScalar && !this->para->getDiffOn())
        throw std::runtime_error("PlaneAverageProbe: Scalar can only be sampled if diff is on!");
    if (sampleSubgridScaleFluxes && !this->para->getUseTurbulentViscosity())
        throw std::runtime_error("PlanarAverageProbe: sampleSubgridScaleFluxes is true but useTurbulentViscosity of "
                                 "parameter is false.");

    size_t numberOfVariables = 0;
    for (auto statistic : statistics) {
        numberOfVariables += getVariableNames(statistic, false).size();
    }

    for (int level = 0; level <= para->getMaxLevel(); level++) {
        auto* data = &levelData.emplace_back();

        findCoordinatesForPlanes(level, data->coordinateX, data->coordinateY, data->coordinateZ,
                                 data->numberOfNodesPerPlane);
        data->numberOfPlanes = static_cast<uint>(data->coordinateX.size());
        data->maxNumberOfPointsPerPlane =
            std::accumulate(data->numberOfNodesPerPlane.begin(), data->numberOfNodesPerPlane.end(), uint {},
                            [](uint acc, uint nNodesPerPlane) { return std::max(acc, nNodesPerPlane); });
        cudaMemoryManager->cudaAllocPlanarAverageProbeIndices(this, level);
        if (sampleSubgridScaleFluxes)
            cudaMemoryManager->cudaAllocPlanarAverageProbeSubgridScaleFluxes(this, level);
        data->instantaneous.resize(data->numberOfPlanes, std::vector<real>(numberOfVariables, c0o1));
        if (computeTimeAverages) {
            data->timeAverages.resize(data->numberOfPlanes, std::vector<real>(numberOfVariables, c0o1));
        }
    }
}

void PlanarAverageProbe::sample(int level, uint t)
{
    if (t < tStartSampling)
        return;

    const uint tLevel = para->getTimeStep(level, t, true);
    const uint levelFactor = std::exp2(level);
    const uint tAfterStartAveraging = tLevel - tStartSampling * levelFactor;
    const uint tAfterStartWriting = tLevel - tStartWritingOutput * levelFactor;
    const bool doTimeAverages = computeTimeAverages && tLevel > (tStartTemporalAveraging * levelFactor);

    if (tAfterStartAveraging % (tBetweenSamples * levelFactor) == 0) {
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

struct isCoordAndFluid
{
    const real* coords;
    const uint* typeOfGridNode;
    const real coord;
    __host__ __device__ bool operator()(const unsigned long long index) const
    {
        return coords[index] == coord && typeOfGridNode[index] == GEO_FLUID;
    }
};

void PlanarAverageProbe::findCoordinatesForPlanes(int level, std::vector<real>& coordinateX, std::vector<real>& coordinateY,
                                                  std::vector<real>& coordinateZ, std::vector<uint>& numberOfNodesPerPlane)
{
    const real* coords = getPlaneNormalCoordinatesH(level);
    const uint* typeOfGridNode = para->getParH(level)->typeOfGridNode;

    std::vector<real> planeCoords;
    real coordOfCurrentPlane = coords[1];
    uint count = 1;
    for (unsigned long long index = 1; index < para->getParH(level)->numberOfNodes; index++) {
        if (isCoordAndFluid { coords, typeOfGridNode, coordOfCurrentPlane }(index)) {
            count++;
            continue;
        }
        if (typeOfGridNode[index] != GEO_FLUID)
            continue;

        const real coordOfLastPlane = coordOfCurrentPlane;
        coordOfCurrentPlane = coords[index];

        if (count == 1) // skip empty planes
            continue;

        numberOfNodesPerPlane.emplace_back(count);
        planeCoords.emplace_back(coordOfLastPlane);
        count = 1;
    }

    switch (planeNormal) {
        case x:
            std::copy(planeCoords.begin(), planeCoords.end(), std::back_inserter(coordinateX));
            std::fill_n(std::back_inserter(coordinateY), planeCoords.size(), 999999);
            std::fill_n(std::back_inserter(coordinateZ), planeCoords.size(), 999999);
            break;
        case y:
            std::fill_n(std::back_inserter(coordinateX), planeCoords.size(), 999999);
            std::copy(planeCoords.begin(), planeCoords.end(), std::back_inserter(coordinateY));
            std::fill_n(std::back_inserter(coordinateZ), planeCoords.size(), 999999);
            break;
        case z:
            std::fill_n(std::back_inserter(coordinateX), planeCoords.size(), 999999);
            std::fill_n(std::back_inserter(coordinateY), planeCoords.size(), 999999);
            std::copy(planeCoords.begin(), planeCoords.end(), std::back_inserter(coordinateZ));
            break;
    }
}

real computeMean(iterPair x, uint nNodes)
{
    return thrust::reduce(std::get<0>(x), std::get<1>(x)) / nNodes;
}

real computeMean(real* x, indexPointer indices, uint nNodes)
{
    return computeMean(getIteratorPair(x, indices, nNodes), nNodes);
}

real computeMean(real* x, uint nNodes)
{
    return thrust::reduce(thrust::device, x, x + nNodes) / nNodes;
}

Means computeMeans(iterPair vx, iterPair vy, iterPair vz, iterPair phi, uint nNodes, bool sampleScalar)
{
    Means means;
    means.vx = computeMean(vx, nNodes);
    means.vy = computeMean(vy, nNodes);
    means.vz = computeMean(vz, nNodes);
    if (sampleScalar)
        means.phi = computeMean(phi, nNodes);
    return means;
}

real computeCovariance(iterPair x, iterPair y, real mean_x, real mean_y, uint nNodes)
{
    auto begin = thrust::make_zip_iterator(thrust::make_tuple(x.first, y.first));
    auto end = thrust::make_zip_iterator(thrust::make_tuple(x.second, y.second));
    return thrust::transform_reduce(begin, end, covariance(mean_x, mean_y), c0o1, thrust::plus<real>()) / nNodes;
}

Covariances computeCovariances(iterPair vx, iterPair vy, iterPair vz, iterPair phi, Means means, uint nNodes,
                               bool sampleScalar)
{
    Covariances covariances;

    covariances.vxvx = computeCovariance(vx, vx, means.vx, means.vx, nNodes);
    covariances.vyvy = computeCovariance(vy, vy, means.vy, means.vy, nNodes);
    covariances.vzvz = computeCovariance(vz, vz, means.vz, means.vz, nNodes);
    covariances.vxvy = computeCovariance(vx, vy, means.vx, means.vy, nNodes);
    covariances.vxvz = computeCovariance(vx, vz, means.vx, means.vz, nNodes);
    covariances.vyvz = computeCovariance(vy, vz, means.vy, means.vz, nNodes);
    if (sampleScalar) {
        covariances.phiphi = computeCovariance(phi, phi, means.phi, means.phi, nNodes);
        covariances.vxphi = computeCovariance(vx, phi, means.vx, means.phi, nNodes);
        covariances.vyphi = computeCovariance(vy, phi, means.vy, means.phi, nNodes);
        covariances.vzphi = computeCovariance(vz, phi, means.vz, means.phi, nNodes);
    }

    return covariances;
}

real computeSkewness(iterPair x, real mean, real covariance, uint nNodes)
{
    return thrust::transform_reduce(x.first, x.second, skewness(mean), c0o1, thrust::plus<real>()) /
           (nNodes * std::pow(covariance, c3o2));
}

Skewnesses computeSkewnesses(Means means, Covariances covariances, iterPair vx, iterPair vy, iterPair vz, iterPair phi,
                             uint nNodes, bool sampleScalar)
{
    Skewnesses skewnesses;

    skewnesses.Sx = computeSkewness(vx, means.vx, covariances.vxvx, nNodes);
    skewnesses.Sy = computeSkewness(vy, means.vy, covariances.vyvy, nNodes);
    skewnesses.Sz = computeSkewness(vz, means.vz, covariances.vzvz, nNodes);
    if (sampleScalar)
        skewnesses.Sphi = computeSkewness(phi, means.phi, covariances.phiphi, nNodes);

    return skewnesses;
}

real computeFlatness(iterPair x, real mean, real covariance, uint nNodes)
{
    return thrust::transform_reduce(x.first, x.second, flatness(mean), c0o1, thrust::plus<real>()) /
           (nNodes * std::pow(covariance, 2));
}

Flatnesses computeFlatnesses(iterPair vx, iterPair vy, iterPair vz, iterPair phi, Means means, Covariances covariances,
                             uint nNodes, bool sampleScalar)
{
    Flatnesses flatnesses;

    flatnesses.Fx = computeFlatness(vx, means.vx, covariances.vxvx, nNodes);
    flatnesses.Fy = computeFlatness(vy, means.vy, covariances.vyvy, nNodes);
    flatnesses.Fz = computeFlatness(vz, means.vz, covariances.vzvz, nNodes);
    if (sampleScalar)
        flatnesses.Fphi = computeFlatness(phi, means.phi, covariances.phiphi, nNodes);
    return flatnesses;
}

std::vector<real> PlanarAverageProbe::computePlaneStatistics(int level, uint nNodes)
{
    auto* data = &levelData[level];
    auto parameter = para->getParD(level);
    indexPointer indices = thrust::device_pointer_cast<uint>(data->indicesD);

    const auto velocityX = getIteratorPair(parameter->velocityX, indices, nNodes);
    const auto velocityY = getIteratorPair(parameter->velocityY, indices, nNodes);
    const auto velocityZ = getIteratorPair(parameter->velocityZ, indices, nNodes);

    const auto phi = getIteratorPair(parameter->concentration, indices, nNodes);

    std::vector<real> averages;

    if (!isStatisticIn(Statistic::Means, statistics))
        return averages;

    const real velocityRatio = para->getScaledVelocityRatio(level);
    const real viscosityRatio = para->getScaledViscosityRatio(level);
    const real stressRatio = para->getScaledStressRatio(level);

    const auto means = computeMeans(velocityX, velocityY, velocityZ, phi, nNodes, sampleScalar);
    averages.push_back(means.vx * velocityRatio);
    averages.push_back(means.vy * velocityRatio);
    averages.push_back(means.vz * velocityRatio);
    if (para->getUseTurbulentViscosity())
        averages.push_back(computeMean(parameter->turbulentViscosity, indices, nNodes) * viscosityRatio);
    if (sampleSubgridScaleFluxes) {
        averages.push_back(computeMean(data->subgridScaleFluxXX, nNodes) * stressRatio);
        averages.push_back(computeMean(data->subgridScaleFluxXY, nNodes) * stressRatio);
        averages.push_back(computeMean(data->subgridScaleFluxXZ, nNodes) * stressRatio);
        averages.push_back(computeMean(data->subgridScaleFluxYY, nNodes) * stressRatio);
        averages.push_back(computeMean(data->subgridScaleFluxYZ, nNodes) * stressRatio);
        averages.push_back(computeMean(data->subgridScaleFluxZZ, nNodes) * stressRatio);
    }
    if (sampleScalar)
        averages.push_back(means.phi);
    if (sampleScalar && para->getUseTurbulentViscosity())
        averages.push_back(computeMean(parameter->turbulentDiffusivity, indices, nNodes) * viscosityRatio);
    if (sampleScalar && sampleSubgridScaleFluxes) {
        averages.push_back(computeMean(data->subgridScaleFluxPhiX, nNodes) * velocityRatio);
        averages.push_back(computeMean(data->subgridScaleFluxPhiY, nNodes) * velocityRatio);
        averages.push_back(computeMean(data->subgridScaleFluxPhiZ, nNodes) * velocityRatio);
    }

    if (!isStatisticIn(Statistic::Covariances, statistics))
        return averages;

    const auto covariances = computeCovariances(velocityX, velocityY, velocityZ, phi, means, nNodes, sampleScalar);
    averages.push_back(covariances.vxvx * stressRatio);
    averages.push_back(covariances.vyvy * stressRatio);
    averages.push_back(covariances.vzvz * stressRatio);
    averages.push_back(covariances.vxvy * stressRatio);
    averages.push_back(covariances.vxvz * stressRatio);
    averages.push_back(covariances.vyvz * stressRatio);

    if (sampleScalar) {
        averages.push_back(covariances.phiphi);
        averages.push_back(covariances.vxphi * velocityRatio);
        averages.push_back(covariances.vyphi * velocityRatio);
        averages.push_back(covariances.vzphi * velocityRatio);
    }

    if (!isStatisticIn(Statistic::Skewness, statistics))
        return averages;

    const auto skewnesses =
        computeSkewnesses(means, covariances, velocityX, velocityY, velocityZ, phi, nNodes, sampleScalar);
    averages.push_back(skewnesses.Sx);
    averages.push_back(skewnesses.Sy);
    averages.push_back(skewnesses.Sz);
    if (sampleScalar)
        averages.push_back(skewnesses.Sphi);

    if (!isStatisticIn(Statistic::Flatness, statistics))
        return averages;

    const auto flatnesses =
        computeFlatnesses(velocityX, velocityY, velocityZ, phi, means, covariances, nNodes, sampleScalar);
    averages.push_back(flatnesses.Fx);
    averages.push_back(flatnesses.Fy);
    averages.push_back(flatnesses.Fz);
    if (sampleScalar)
        averages.push_back(flatnesses.Fphi);

    return averages;
}

std::vector<real> computeNewTimeAverages(std::vector<real>& oldAverages, std::vector<real>& instantaneous,
                                         uint numberOfTimesteps)
{
    std::vector<real> newAverages;
    newAverages.reserve(oldAverages.size());
    for (uint i = 0; i < oldAverages.size(); i++) {
        newAverages.push_back(computeNewTimeAverage(oldAverages[i], instantaneous[i], numberOfTimesteps));
    }
    return newAverages;
}

void PlanarAverageProbe::calculateQuantities(int level, bool doTimeAverages)
{
    auto* data = &levelData[level];
    auto parameter = para->getParD(level);
    calculateMacroscopicQuantitiesCompressible(
        parameter->velocityX, parameter->velocityY, parameter->velocityZ, parameter->rho, parameter->pressure,
        parameter->typeOfGridNode, parameter->neighborX, parameter->neighborY, parameter->neighborZ,
        parameter->numberOfNodes, parameter->numberofthreads, parameter->distributions.f[0], parameter->isEvenTimestep);

    auto counting = thrust::make_counting_iterator<unsigned long long>(1);

    for (uint plane = 0; plane < data->numberOfPlanes; plane++) {
        thrust::copy_if(
            thrust::device, counting, counting + parameter->numberOfNodes - counting[0], data->indicesD,
            isCoordAndFluid { getPlaneNormalCoordinatesD(level), parameter->typeOfGridNode, data->coordinateZ[plane] });
        const uint numberOfNodesInPlane = data->numberOfNodesPerPlane[plane];
        if (sampleSubgridScaleFluxes)
            calculateSubGridScaleFluxesCompressible(
                data->indicesD, numberOfNodesInPlane, data->subgridScaleFluxXX, data->subgridScaleFluxXY,
                data->subgridScaleFluxXZ, data->subgridScaleFluxYY, data->subgridScaleFluxYZ, data->subgridScaleFluxZZ,
                data->subgridScaleFluxPhiX, data->subgridScaleFluxPhiY, data->subgridScaleFluxPhiZ,
                parameter->typeOfGridNode, parameter->velocityX, parameter->velocityY, parameter->velocityZ,
                parameter->concentration, parameter->turbulentViscosity, parameter->turbulentDiffusivity,
                parameter->neighborX, parameter->neighborY, parameter->neighborZ, parameter->numberOfNodes,
                parameter->distributions.f[0], parameter->distributionsAD.f[0], parameter->omega, parameter->omegaDiffusivity, parameter->numberofthreads,
                parameter->isEvenTimestep, sampleScalar);
        data->instantaneous[plane] = computePlaneStatistics(level, numberOfNodesInPlane);
        if (doTimeAverages)
            data->timeAverages[plane] = computeNewTimeAverages(data->timeAverages[plane], data->instantaneous[plane],
                                                               data->numberOfTimestepsInTimeAverage + 1U);
    }
    if (doTimeAverages)
        data->numberOfTimestepsInTimeAverage++;
}

const real* PlanarAverageProbe::getPlaneNormalCoordinatesH(int level)
{
    switch (planeNormal) {
        case Axis::x:
            return para->getParH(level)->coordinateX;
        case Axis::y:
            return para->getParH(level)->coordinateY;
        case Axis::z:
            return para->getParH(level)->coordinateZ;
    };
    throw std::runtime_error("PlaneNormal not defined");
}
const real* PlanarAverageProbe::getPlaneNormalCoordinatesD(int level)
{
    switch (planeNormal) {
        case Axis::x:
            return para->getParD(level)->coordinateX;
        case Axis::y:
            return para->getParD(level)->coordinateY;
        case Axis::z:
            return para->getParD(level)->coordinateZ;
    };
    throw std::runtime_error("PlaneNormal not defined");
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
    auto* data = &levelData[level];
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
        for (const auto& variable : getVariableNames(statistic, false)) {
            varNames.push_back(variable);
        }
    }
    if (this->computeTimeAverages) {
        for (auto statistic : statistics) {
            for (const auto& variable : getVariableNames(statistic, true))
                varNames.push_back(variable);
        }
    }
    return varNames;
}

bool isStatisticIn(PlanarAverageProbe::Statistic statistic, std::vector<PlanarAverageProbe::Statistic> statistics)
{
    return std::find(statistics.begin(), statistics.end(), statistic) != statistics.end();
}

}
//! \}