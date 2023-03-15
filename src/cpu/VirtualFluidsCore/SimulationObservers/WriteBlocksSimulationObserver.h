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
//! \file WriteBlocksSimulationObserver.h
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef WriteBlocksSimulationObserver_H_
#define WriteBlocksSimulationObserver_H_

#include <PointerDefinitions.h>
#include <string>

#include "SimulationObserver.h"

namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;
class WbWriter;

//! \class WriteBlocksSimulationObserver
//! \brief A class writes a block grid to a VTK-file
class WriteBlocksSimulationObserver : public SimulationObserver
{
public:
    //! \brief Construct WriteBlocksSimulationObserver object.
    //! \pre The Grid3D and UbScheduler objects must exist.
    //! \param grid is observable Grid3D object
    //! \param s is UbScheduler object for scheduling of observer
    //! \param path is path of folder for output
    //! \param writer is WbWriter object
    //! \param comm is Communicator object
    WriteBlocksSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, WbWriter *const writer,
                           std::shared_ptr<vf::mpi::Communicator> comm);
    ~WriteBlocksSimulationObserver() override;

    void process(real step) override;

protected:
    //! Collect data for VTK-file
    //! \param step is a time step
    void collectData(real step);

    std::string path;
    WbWriter *writer;
    std::shared_ptr<vf::mpi::Communicator> comm;
};

#endif
