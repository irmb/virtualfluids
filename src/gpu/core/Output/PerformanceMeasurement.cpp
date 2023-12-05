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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \author Soeren Peters
//=======================================================================================
#include "PerformanceMeasurement.h"

#include <logger/Logger.h>

#include "Parameter/Parameter.h"

void PerformanceMeasurement::print(vf::basics::Timer& timer, uint timestep, Parameter* para,
                                   vf::parallel::Communicator& communicator)
{
    totalTime += timer.getTimeInSeconds();
    real fnups = 0.0;
    real bandwidth = 0.0;

    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        fnups += 1000.0 * (timestep - para->getTimestepStart()) * para->getParH(lev)->numberOfNodes * pow(2., lev) /
                 totalTime * 1000;
        bandwidth += (27.0 + 1.0) * 4.0 * 1000.0 * (timestep - para->getTimestepStart()) *
                     para->getParH(lev)->numberOfNodes / totalTime;
    }

    if (this->firstOutput && communicator.getProcessID() == 0) // only display the legend once
    {
        VF_LOG_INFO("PID \t --- Average performance ---  Processing time (ms) \t Nups in Mio \t Bandwidth in GB/sec");
        this->firstOutput = false;
    }

    VF_LOG_INFO(" {} \t --- Average performance --- {:>8.1f}/ {:<8.1f} \t   {:5.1f} \t       {:4.1f}",
                communicator.getProcessID(), timer.getTimeInSeconds() * 1000, totalTime * 1000, fnups, bandwidth);

    // When using multiple GPUs, sum the nups of all processes
    if (communicator.getNumberOfProcesses() > 1) {
        double nupsSum = communicator.reduceSum(fnups);
        if (communicator.getProcessID() == 0)
            VF_LOG_INFO("Sum of all {} processes: Nups in Mio: {:.1f}", communicator.getNumberOfProcesses(), nupsSum);
    }
}