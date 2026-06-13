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
//! \addtogroup gpu_Output Output
//! \ingroup gpu_core core
//! \{
//! \author Soeren Peters
//=======================================================================================
#include "PerformanceMeasurement.h"

#include <logger/Logger.h>

#include "Parameter/Parameter.h"

namespace vf::gpu {

PerformanceMeasurement::PerformanceMeasurement(const Parameter& para) : timestepStart(para.getTimestepStart())
{
    for (int level = para.getCoarse(); level <= para.getFine(); level++) {
        totalNumberOfNodesCorrected += para.getParH(level)->numberOfNodes * pow(2., level);
        totalNumberOfNodes += para.getParH(level)->numberOfNodes;
    }
}

double PerformanceMeasurement::getNups() const
{
    return nups;
}

double PerformanceMeasurement::totalRuntimeInSeconds() const
{
    return totalTime;
}

void PerformanceMeasurement::log(vf::basics::Timer& timer, uint timestep, vf::parallel::Communicator& communicator)
{
    totalTime += timer.getTimeInSeconds();
    const uint numberOfTimeSteps = timestep - timestepStart;

    nups = numberOfTimeSteps * totalNumberOfNodesCorrected / totalTime;
    const real bandwidth = (27.0 + 1.0) * 4.0 * numberOfTimeSteps * totalNumberOfNodes / (totalTime * 1.0E9);

    const double meganups = nups / 1.0E6;

    if (this->firstOutput && communicator.isRoot()) // only display the legend once
    {
        VF_LOG_INFO("PID \t --- Average performance ---  Processing time (ms) \t Nups in Mio \t Bandwidth in GB/sec");
        this->firstOutput = false;
    }

    VF_LOG_INFO(" {} \t --- Average performance --- {:>8.1f}/ {:<8.1f} \t   {:5.1f} \t       {:4.1f}",
                communicator.getProcessID(), timer.getTimeInSeconds() * 1000, totalTime * 1000, meganups, bandwidth);

    // When using multiple GPUs, sum the nups of all processes
    if (communicator.getNumberOfProcesses() > 1) {
        const double nupsSum = communicator.reduceSum(nups);
        if (communicator.isRoot())
            VF_LOG_INFO("Sum of all {} processes: Nups in Mio: {:.1f}", communicator.getNumberOfProcesses(), nupsSum*1e-6);
    }
}

}

//! \}
