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
//! \file TimeseriesSimulationObserver.cpp
//! \ingroup SimulationObservers
//! \author Sonja Uphoff
//=======================================================================================

#include "TimeseriesSimulationObserver.h"

#include <fstream>

#include <parallel/Communicator.h>
#include "Grid3D.h"
#include "IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "UbScheduler.h"

TimeseriesSimulationObserver::TimeseriesSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<IntegrateValuesHelper> h1,
                                             const std::string &path, std::shared_ptr<vf::parallel::Communicator> comm)
    : SimulationObserver(grid, s), h1(h1), path(path), comm(comm)
{
    if (comm->getProcessID() == comm->getRoot()) {
        std::ofstream ostr;
        // fname = path+"/timeseries/timeseries"+UbSystem::toString(grid->getTimeStep())+".csv";
        fname = path + ".csv";
        UBLOG(logINFO, "TimeseriesWriterSimulationObserver::fname:" << fname);
        ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
        if (!ostr) {
            ostr.clear();
            std::string file_path = UbSystem::getPathFromString(fname);
            if (file_path.size() > 0) {
                UbSystem::makeDirectory(file_path);
                ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }
        ostr << "step;rho;vx;vy;vz;volume\n";
        ostr.close();
        UBLOG(logINFO, "TimeseriesWriterSimulationObserver::Constructor:end");
    }
}
//////////////////////////////////////////////////////////////////////////
TimeseriesSimulationObserver::~TimeseriesSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
void TimeseriesSimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        collectData(step);
}
//////////////////////////////////////////////////////////////////////////
void TimeseriesSimulationObserver::collectData(real step)
{
    h1->calculateMQ();

    UBLOG(logDEBUG3, "TimeseriesWriterSimulationObserver::update:" << step);

    if (comm->getProcessID() == comm->getRoot()) {
        int istep = static_cast<int>(step);
        std::ofstream ostr;
        real cellsVolume = h1->getCellsVolume();

        real rho    = (h1->getRho()) / cellsVolume;
        real vx     = (h1->getVx1()) / cellsVolume;
        real vy     = (h1->getVx2()) / cellsVolume;
        real vz     = (h1->getVx3()) / cellsVolume;
        real volume = cellsVolume;

        ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
        if (!ostr) {
            ostr.clear();
            std::string path = UbSystem::getPathFromString(fname);
            if (path.size() > 0) {
                UbSystem::makeDirectory(path);
                ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }

        ostr << istep << ";" << rho << ";" << vx << ";" << vy << ";" << vz << ";" << volume << "\n";
        ostr.close();
    }
}
