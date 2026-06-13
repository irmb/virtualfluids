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
//! \author Henrik Asmuth, Henry Korb
//! \date 13/05/2022
//=======================================================================================

#ifndef PlanarAverageProbe_H
#define PlanarAverageProbe_H

#include "Sampler.h"

#include <stdexcept>
#include <string>
#include <vector>

#include <basics/DataTypes.h>
#include <basics/geometry3d/Axis.h>
#include <logger/Logger.h>

namespace vf::gpu {

class Parameter;
class CudaMemoryManager;

//! \brief Computes spatial statistics across x, y or z-normal planes defined by planeNormal.
//! The planes include all points of the domain at each respective position along that normal direction.
//! The spatial statistics can additionally be averaged in time.
//! The name phi is used to denote the scalar field.
class PlanarAverageProbe : public Sampler
{
public:
    enum class Statistic {
        Means,
        Covariances,
        Skewness,
        Flatness,
    };
    struct LevelData;

public:
    //! \param tStartSampling The first timestep at which the probe samples planar averages.
    //! \param tStartTemporalAveraging The first timestep at which temporal averaging starts.
    //! \param tBetweenSamples The number of timesteps between samples.
    //! \param tStartWritingOutput The first timestep at which the probe writes output.
    //! \param tBetweenWriting The number of timesteps between writing output.
    //! \param planeNormal The normal direction of the planes along which the probe samples.
    //! \param computeTimeAverages If true, the probe computes time averages.
    //! \param sampleScalar If true, the probe samples statistics related to the scalar.
    //! \param sampleSubgridScaleFluxes If true, the probe samples sub-gridscale fluxes.
    PlanarAverageProbe(SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaMemoryManager, const std::string& outputPath,
                       const std::string& probeName, uint tStartSampling, uint tStartTemporalAveraging, uint tBetweenSamples,
                       uint tStartWritingOutput, uint tBetweenWriting, Axis planeNormal, bool computeTimeAverages,
                       bool sampleScalar = false, bool sampleSubgridScaleFluxes = false)
        : para(std::move(para)), cudaMemoryManager(std::move(cudaMemoryManager)), tStartSampling(tStartSampling),
          tStartTemporalAveraging(tStartTemporalAveraging), tBetweenSamples(tBetweenSamples),
          tStartWritingOutput(tStartWritingOutput), tBetweenWriting(tBetweenWriting),
          computeTimeAverages(computeTimeAverages), planeNormal(planeNormal), sampleScalar(sampleScalar),
          sampleSubgridScaleFluxes(sampleSubgridScaleFluxes), Sampler(outputPath, probeName)
    {
        if (tBetweenSamples == 0)
            throw std::runtime_error("PlanarAverageProbe: tBetweenSamples is 0! tBetweenSamples must be larger than 0");
        if (tBetweenWriting == 0)
            throw std::runtime_error("PlanarAverageProbe: tBetweenWriting is 0! tBetweenWriting must be larger than 0");
        if (tStartTemporalAveraging < tStartSampling && computeTimeAverages)
            throw std::runtime_error("PlaneAverageProbe: tStartTemporalAveraging must be larger than tStartSampling!");

        VF_LOG_INFO(
            "Created planar averaging probe, output path: " + outputPath + ", probe name: " + probeName +
                " start sampling: {}, time steps between sampling: {}, start writing: {}, time steps between writing: {}",
            tStartSampling, tBetweenSamples, tStartWritingOutput, tBetweenWriting);
        if(sampleScalar)
            VF_LOG_INFO("sampling scalar.");
        if(sampleSubgridScaleFluxes)
            VF_LOG_INFO("sampling sub-gridscale fluxes.");
    }
    ~PlanarAverageProbe() override;

    void init() override;
    void sample(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* /**/) override {};
    LevelData* getLevelData(int level)
    {
        return &levelData[level];
    }
    void addAllAvailableStatistics()
    {
        addStatistic(Statistic::Means);
        addStatistic(Statistic::Covariances);
        addStatistic(Statistic::Skewness);
        addStatistic(Statistic::Flatness);
    }
    void addStatistic(Statistic statistic);
    void setFileNameToNOut()
    {
        nameFilesWithFileCount = true;
    }
    bool getSampleScalar() const
    {
        return sampleScalar;
    }

private:
    std::vector<std::string> getVariableNames(Statistic statistic, bool namesForTimeAverages) const;
    void copyDataToNodedata(std::vector<std::vector<real>>& data, std::vector<std::vector<double>>& nodeData);
    void calculateQuantities(int level, bool doTimeAverages);
    void findCoordinatesForPlanes(int level, std::vector<real>& coordinateX, std::vector<real>& coordinateY,
                                  std::vector<real>& coordinateZ, std::vector<uint>& numberOfNodesPerPlane);
    std::vector<real> computePlaneStatistics(int level, uint nNodes);

    std::vector<std::string> getAllVariableNames();
    const uint* getNeighborIndicesInPlaneNormal(int level);
    const real* getPlaneNormalCoordinatesH(int level);
    const real* getPlaneNormalCoordinatesD(int level);
    void writeGridFile(int level, uint tWrite);
    void writeParallelFile(uint tWrite);

private:
    SPtr<Parameter> para;
    SPtr<CudaMemoryManager> cudaMemoryManager;
    const uint tStartSampling, tStartTemporalAveraging, tBetweenSamples, tStartWritingOutput, tBetweenWriting;
    const bool computeTimeAverages, sampleScalar, sampleSubgridScaleFluxes;
    bool nameFilesWithFileCount = false;
    const Axis planeNormal;
    std::vector<Statistic> statistics;
    std::vector<LevelData> levelData;
    std::vector<std::string> fileNamesForCollectionFile;
};

bool isStatisticIn(PlanarAverageProbe::Statistic statistic, std::vector<PlanarAverageProbe::Statistic> statistics);

struct PlanarAverageProbe::LevelData
{
    uint* indicesD;
    real *subgridScaleFluxXX, *subgridScaleFluxXY, *subgridScaleFluxXZ, *subgridScaleFluxYY, *subgridScaleFluxYZ,
        *subgridScaleFluxZZ, *subgridScaleFluxPhiX, *subgridScaleFluxPhiY, *subgridScaleFluxPhiZ;
    std::vector<uint> numberOfNodesPerPlane;
    uint numberOfPlanes {}, maxNumberOfPointsPerPlane {}, numberOfTimestepsInTimeAverage {};
    std::vector<real> coordinateX, coordinateY, coordinateZ;
    std::vector<std::vector<real>> instantaneous;
    std::vector<std::vector<real>> timeAverages;
};

}

#endif
//! \}
