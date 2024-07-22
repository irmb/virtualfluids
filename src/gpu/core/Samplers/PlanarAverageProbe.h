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
//! \author Henrik Asmuth
//! \date 13/05/2022
//! \brief Probe computing statistics across planes spanning the entire domain
//!

//!
//=======================================================================================

#ifndef PlanarAverageProbe_H
#define PlanarAverageProbe_H

#include "Sampler.h"

#include <string>
#include <vector>

#include <basics/DataTypes.h>
#include <basics/geometry3d/Axis.h>

class Parameter;
class CudaMemoryManager;

//! \brief Computes spatial statistics across x, y or z-normal planes defined by planeNormal.
//! The planes include all points of the domain at each respective position along that normal direction.
//! The spatial statistics can additionally be averaged in time.
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
    PlanarAverageProbe(SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaMemoryManager, const std::string outputPath,
                       const std::string probeName, uint tStartAveraging, uint tStartTemporalAveraging,
                       uint tBetweenAverages, uint tStartWritingOutput, uint tBetweenWriting, Axis planeNormal,
                       bool computeTimeAverages, bool computeStatisticsOfConcentration);
    ~PlanarAverageProbe();

    void init() override;
    void sample(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* gridProvider) override {};
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

private:
    std::vector<std::string> getVariableNames(Statistic statistic, bool namesForTimeAverages) const;
    void copyDataToNodedata(std::vector<std::vector<real>>& data, std::vector<std::vector<double>>& nodeData);
    void calculateQuantities(int level, bool doTimeAverages);
    std::vector<unsigned long long> findIndicesInPlane(int level);
    void findCoordinatesForPlanes(int level, std::vector<real>& coordinateX, std::vector<real>& coordinateY,
                                  std::vector<real>& coordinateZ);
    std::vector<real> computePlaneStatistics(int level);

    std::vector<std::string> getAllVariableNames();
    const uint* getNeighborIndicesInPlaneNormal(int level);
    void writeGridFile(int level, uint tWrite);
    void writeParallelFile(uint tWrite);

private:
    SPtr<Parameter> para;
    SPtr<CudaMemoryManager> cudaMemoryManager;
    uint tStartAveraging, tStartTemporalAveraging, tBetweenAverages, tStartWritingOutput, tBetweenWriting;
    const bool computeTimeAverages, computeStatisticsOfConcentration;
    bool nameFilesWithFileCount = false;
    Axis planeNormal;
    std::vector<Statistic> statistics;
    std::vector<LevelData> levelData;
    std::vector<std::string> fileNamesForCollectionFile;
};

bool isStatisticIn(PlanarAverageProbe::Statistic statistic, std::vector<PlanarAverageProbe::Statistic> statistics);

struct PlanarAverageProbe::LevelData
{
    unsigned long long *indicesOfFirstPlaneH, *indicesOfFirstPlaneD;
    uint numberOfPlanes {}, numberOfPointsPerPlane {}, numberOfTimestepsInTimeAverage {};
    std::vector<real> coordinateX, coordinateY, coordinateZ;
    std::vector<std::vector<real>> instantaneous;
    std::vector<std::vector<real>> timeAverages;
};

#endif
//! \}
