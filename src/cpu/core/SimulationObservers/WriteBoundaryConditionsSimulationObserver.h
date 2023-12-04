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
//! \file WriteBoundaryConditionsSimulationObserver.h
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef WriteBoundaryConditionsSimulationObserver_H
#define WriteBoundaryConditionsSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "SimulationObserver.h"
#include "UbTuple.h"

namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;
class WbWriter;
class Block3D;
class LBMUnitConverter;

//! \brief A class writes boundary conditions information to a VTK-file
class WriteBoundaryConditionsSimulationObserver : public SimulationObserver
{
public:
    WriteBoundaryConditionsSimulationObserver();
    //! \brief Construct WriteBoundaryConditionsSimulationObserver object
    //! \pre The Grid3D and UbScheduler objects must exist
    //! \param grid is observable Grid3D object
    //! \param s is UbScheduler object for scheduling of observer
    //! \param path is path of folder for output
    //! \param writer is WbWriter object
    //! \param comm is Communicator object
    WriteBoundaryConditionsSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                       WbWriter *const writer, std::shared_ptr<vf::parallel::Communicator> comm);
    ~WriteBoundaryConditionsSimulationObserver() override = default;

    void update(real step) override;

protected:
    //! Collect data for VTK-file
    //! \param step is a time step
    void collectData(real step);
    void addDataGeo(SPtr<Block3D> block);
    void clearData();

private:
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleUInt8> cells;
    std::vector<std::string> datanames;
    std::vector<std::vector<double>> data;
    std::string path;
    WbWriter *writer;
    std::vector<std::vector<SPtr<Block3D>>> blockVector;
    int minInitLevel;
    int maxInitLevel;
    int gridRank;
    std::shared_ptr<vf::parallel::Communicator> comm;
};
#endif
