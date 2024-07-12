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
//!
//!
//=======================================================================================

#ifndef WallModelProbe_H
#define WallModelProbe_H

#include <string>
#include <vector>

#include <basics/PointerDefinitions.h>
#include <logger/Logger.h>

#include "Sampler.h"

///////////////////////////////////////////////////////////////////////////////////

class Parameter;
class CudaMemoryManager;

struct WallModelProbeLevelData
{
    uint numberOfAveragedValues, numberOfFluidNodes;
    std::string timeseriesFileName;
    std::vector<std::vector<real>> instantaneousData, averagedData;
    std::vector<real> timestepTime;
    bool firstWrite = true;
    WallModelProbeLevelData(std::string fileName, uint numberOfFluidNodes)
        : timeseriesFileName(fileName), numberOfFluidNodes(numberOfFluidNodes)
    {
    }
};

//! \brief Probe computing statistics of all relevant wall model quantities used in the StressBC kernels
//! Computes spatial statistics for all grid points of the StressBC
//! The spatial statistics can additionally be averaged in time.
class WallModelProbe : public Sampler
{
public:
    WallModelProbe(SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaMemoryManager, const std::string outputPath,
                   const std::string probeName, uint tStartAveraging, uint tStartTemporalAveraging, uint tBetweenAverages,
                   uint tStartWritingOutput, uint tBetweenWriting, bool averageEveryTimestep, bool computeTemporalAverages,
                   bool outputStress, bool evaluatePressureGradient)
        : para(para), cudaMemoryManager(cudaMemoryManager), tStartAveraging(tStartAveraging),
          tStartTemporalAveraging(tStartTemporalAveraging), tStartWritingOutput(tStartWritingOutput),
          tBetweenAverages(tBetweenAverages), tBetweenWriting(tBetweenWriting), averageEveryTimestep(averageEveryTimestep),
          computeTemporalAverages(computeTemporalAverages), outputStress(outputStress),
          evaluatePressureGradient(evaluatePressureGradient), Sampler(outputPath, probeName)
    {
        if (tStartTemporalAveraging < tStartAveraging)
            throw std::runtime_error("WallModelProbe: tStartTemporalAveraging must be larger than tStartAveraging!");
        if (averageEveryTimestep)
            VF_LOG_INFO("WallModelProbe: averageEveryTimestep is true, ignoring tBetweenAverages");
        if (tBetweenAverages == 0 && !averageEveryTimestep)
            throw std::runtime_error("WallModelProbe: tBetweenAverages must be larger than 0!");
    }

    ~WallModelProbe() = default;

    void init() override;
    void sample(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* gridProvider) override {};

private:
    std::vector<std::string> getVariableNames();
    int getNumberOfInstantaneousQuantities()
    {
        return evaluatePressureGradient ? 13 : 10;
    }
    void calculateQuantities(WallModelProbeLevelData* levelData, uint t, int level);
    void write(int level);
    uint countFluidNodes(int level);

private:
    SPtr<Parameter> para;
    SPtr<CudaMemoryManager> cudaMemoryManager;
    uint tStartAveraging, tStartTemporalAveraging, tBetweenAverages, tStartWritingOutput, tBetweenWriting;
    bool outputStress = false;
    bool evaluatePressureGradient = false;
    bool computeTemporalAverages = false;
    bool averageEveryTimestep = false;
    std::vector<WallModelProbeLevelData> levelData;
};

#endif
//! \}
