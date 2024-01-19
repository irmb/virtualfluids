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
//! \author Henrik Asmuth
//! \date 13/05/2022
//=======================================================================================

#ifndef WallModelProbe_H
#define WallModelProbe_H

#include <basics/PointerDefinitions.h>

#include "Probe.h"

///////////////////////////////////////////////////////////////////////////////////

//! \brief Probe computing statistics of all relevant wall model quantities used in the StressBC kernels
//!
//! Computes spatial statistics for all grid points of the StressBC 
//! The spatial statistics can additionally be averaged in time.
//!
class WallModelProbe : public Probe
{
public: 
    WallModelProbe(
        const std::string probeName,
        const std::string outputPath,
        uint tStartAvg,
        uint tStartTmpAvg,
        uint tAvg,
        uint tStartOut,
        uint tOut
    ):  Probe(probeName, 
            outputPath,
            tStartAvg,
            tStartTmpAvg,
            tAvg,
            tStartOut, 
            tOut,
            false,
            true)
    {
        if (tStartTmpAvg<tStartAvg)   throw std::runtime_error("Probe: tStartTmpAvg must be larger than tStartAvg!");
    }

    ~WallModelProbe() = default;


    void setForceOutputToStress(bool _outputStress){ this->outputStress = _outputStress; }
    void setEvaluatePressureGradient(bool _evalPressGrad){ this->evaluatePressureGradient = _evalPressGrad; }

private:
    bool isAvailableStatistic(Statistic _variable) override;

    std::vector<PostProcessingVariable> getPostProcessingVariables(Statistic variable) override;

    void findPoints(std::vector<int>& probeIndices_level,
                    std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                    std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                    int level) override;
    void calculateQuantities(SPtr<ProbeStruct> probeStruct, uint t, int level) override;
    void getTaggedFluidNodes(GridProvider* gridProvider) override {};
    uint getNumberOfTimestepsInTimeseries(int level) override;

private:
    bool outputStress = false; //!> if true, output wall force is converted to a stress 
    bool evaluatePressureGradient = false; //!> if true, mean global pressure gradient will also be evaluated
};

#endif
//! \}
