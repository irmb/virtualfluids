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

#ifndef WallModelProbe_H
#define WallModelProbe_H

#include <string>
#include <vector>

#include <basics/PointerDefinitions.h>
#include <logger/Logger.h>

#include "Sampler.h"

namespace vf::gpu {

class Parameter;
class CudaMemoryManager;

//! \brief Probe computing statistics of all relevant wall model quantities used in the StressBC kernels
//! Computes spatial statistics for all grid points of the StressBC
//! The spatial statistics can additionally be averaged in time.
class WallModelProbe : public Sampler
{
public:
    WallModelProbe(SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaMemoryManager, const std::string& outputPath,
                   const std::string& probeName, uint tStartSampling, uint tStartTemporalAveraging, uint tBetweenSamples,
                   uint tStartWritingOutput, uint tBetweenWriting, bool sampleEveryTimestep, bool computeTemporalAverages,
                   bool outputStress, bool evaluatePressureGradient, bool sampleSurfaceLayer)
        : para(std::move(para)), cudaMemoryManager(std::move(cudaMemoryManager)), tStartSampling(tStartSampling),
          tStartTemporalAveraging(tStartTemporalAveraging), tStartWritingOutput(tStartWritingOutput),
          tBetweenSamples(tBetweenSamples), tBetweenWriting(tBetweenWriting), sampleEveryTimestep(sampleEveryTimestep),
          computeTemporalAverages(computeTemporalAverages), outputStress(outputStress),
          evaluatePressureGradient(evaluatePressureGradient), sampleSurfaceLayer(sampleSurfaceLayer), Sampler(outputPath, probeName)
    {
        if(tBetweenSamples == 0)
            throw std::runtime_error("WallModelProbe: tBetweenSamples is 0! If you want to sample every time step set sampleEveryTimestep to true");
        if(tBetweenWriting == 0)
            throw std::runtime_error("WallModelProbe: tBetweenWriting is 0! tBetweenWriting must be larger than 0");
        if (tStartTemporalAveraging < tStartSampling)
            throw std::runtime_error("WallModelProbe: tStartTemporalAveraging must be larger than tStartSampling!");
        if (sampleEveryTimestep)
            VF_LOG_INFO("WallModelProbe: sampleEveryTimestep is true, ignoring tBetweenSamples");
        if (tBetweenSamples == 0 && !sampleEveryTimestep)
            throw std::runtime_error("WallModelProbe: tBetweenSamples must be larger than 0!");
        VF_LOG_INFO(
        "Created wall model probe, output path: " + outputPath + ", probe name: " + probeName +
            " start sampling: {}, time steps between sampling: {}, start writing: {}, time steps between writing: {}",
        tStartSampling, tBetweenSamples, tStartWritingOutput, tBetweenWriting);
    }

    ~WallModelProbe() override = default;

    void init() override;
    void sample(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* /*gridProvider*/) override {};

    struct LevelData;

private:
    std::vector<std::string> getVariableNames();
    uint getNumberOfInstantaneousQuantities() const
    {
        uint numberOfQuantities = 6;
        if(evaluatePressureGradient) numberOfQuantities += 3;
        if(sampleSurfaceLayer) numberOfQuantities += 4;
        return numberOfQuantities;
    }
    void calculateQuantities(LevelData* levelData, uint t, int level);
    void write(int level);
    uint countFluidNodes(int level);

private:
    SPtr<Parameter> para;
    SPtr<CudaMemoryManager> cudaMemoryManager;
    const uint tStartSampling, tStartTemporalAveraging, tBetweenSamples, tStartWritingOutput, tBetweenWriting;
    const bool outputStress = false;
    const bool evaluatePressureGradient = false;
    const bool computeTemporalAverages = false;
    const bool sampleEveryTimestep = false;
    const bool sampleSurfaceLayer = false;
    std::vector<LevelData> levelData;
};

struct WallModelProbe::LevelData
{
    uint numberOfAveragedValues {}, numberOfFluidNodes {};
    std::string timeseriesFileName;
    std::vector<std::vector<real>> instantaneousData, averagedData;
    std::vector<real> timestepTime;
    bool firstWrite = true;
    LevelData(std::string fileName, uint numberOfFluidNodes)
        : timeseriesFileName(std::move(fileName)), numberOfFluidNodes(numberOfFluidNodes)
    {
    }
};

}

#endif
//! \}
