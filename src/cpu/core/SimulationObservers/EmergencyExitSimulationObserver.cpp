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
//! \file EmergencyExitSimulationObserver.h
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================
#include "EmergencyExitSimulationObserver.h"
#include <parallel/Communicator.h>
#include "Grid3D.h"
#include "MPIIORestartSimulationObserver.h"
#include "UbLogger.h"
#include "UbScheduler.h"
#include <basics/utilities/UbFileInputASCII.h>
#include <basics/utilities/UbFileOutputASCII.h>

EmergencyExitSimulationObserver::EmergencyExitSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                                   SPtr<MPIIORestartSimulationObserver> rp, std::shared_ptr<vf::parallel::Communicator> comm)
    : SimulationObserver(grid, s), path(path), rp(rp), comm(comm)
{
    this->path = path + "/exit";
    metafile   = this->path + "/stop.txt";
    if (comm->getProcessID() == comm->getRoot()) {
        // checkMetafile();
        writeMetafile(false);
    }
    comm->barrier();
}
//////////////////////////////////////////////////////////////////////////
EmergencyExitSimulationObserver::~EmergencyExitSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
void EmergencyExitSimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        collectData(step);

    UBLOG(logDEBUG3, "EmergencyExitSimulationObserver::update:" << step);
}

void EmergencyExitSimulationObserver::collectData(real step)
{
    if (readMetafile()) {
        rp->update((int)step);
        if (comm->getProcessID() == comm->getRoot())
            UBLOG(logINFO, "EmergencyExitSimulationObserver save step: " << step);
        comm->barrier();
        exit(EXIT_SUCCESS);
    }
}
//////////////////////////////////////////////////////////////////////////
void EmergencyExitSimulationObserver::writeMetafile(int /*status*/)
{
    UbFileOutputASCII out(metafile);
    out.writeBool(false);
}
//////////////////////////////////////////////////////////////////////////
bool EmergencyExitSimulationObserver::readMetafile()
{
    UbFileInputASCII in(metafile);
    return in.readBool();
}
//////////////////////////////////////////////////////////////////////////
void EmergencyExitSimulationObserver::checkMetafile()
{
    std::ifstream file(metafile.c_str());
    if (!file.is_open()) {
        writeMetafile(false);
        return;
    }
    file.close();
}
