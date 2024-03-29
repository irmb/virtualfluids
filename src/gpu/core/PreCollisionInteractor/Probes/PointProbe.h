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
//! \author Henry Korb, Henrik Asmuth
//! \date 13/05/2022
//=======================================================================================

#ifndef PointProbe_H
#define PointProbe_H

#include "Probe.h"

//! \brief Probe computing statistics for a set of points in space
//!
//! The set of points can be defined by providing a list or on an x-normal plane (the latter being somewhat redundant with PlaneProbe)
//! All statistics are temporal.
//!
class PointProbe: public Probe
{
public:
    PointProbe(
        const std::string probeName,
        const std::string outputPath,
        uint tStartAvg,
        uint tAvg,
        uint tStartOut,
        uint tOut,
        bool outputTimeseries = false
    ): Probe(probeName, 
             outputPath,
             tStartAvg, 
             0,
             tAvg,
             tStartOut, 
             tOut,
             true,
             outputTimeseries)
    {}

    ~PointProbe() = default;

    void addProbePoint(real pointCoordX, real pointCoordY, real pointCoordZ);
    void addProbePointsFromList(std::vector<real>& _pointCoordsX, std::vector<real>& _pointCoordsY, std::vector<real>& _pointCoordsZ);
    void getTaggedFluidNodes(GridProvider* gridProvider) override;
    
private:
    bool isAvailableStatistic(Statistic _variable) override;

    std::vector<PostProcessingVariable> getPostProcessingVariables(Statistic variable) override;

    void findPoints(std::vector<int>& probeIndices_level,
                    std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                    std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                    int level) override;

    void calculateQuantities(SPtr<ProbeStruct> probeStruct, uint t, int level) override;

private:
    std::vector<real> pointCoordsX, pointCoordsY, pointCoordsZ;
    uint getNumberOfTimestepsInTimeseries(int level) override
    {
        (void)para;
        return outputTimeSeries ? tOut * exp2(level) : 1;
    }
};

#endif
//! \}
