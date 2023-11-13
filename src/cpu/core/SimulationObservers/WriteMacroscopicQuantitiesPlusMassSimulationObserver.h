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
//! \file WriteMacroscopicQuantitiesPlusMassSimulationObserver.h
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef WriteMacroscopicQuantitiesPlusMassSimulationObserver_H
#define WriteMacroscopicQuantitiesPlusMassSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "SimulationObserver.h"
#include "LBMSystem.h"
#include "UbTuple.h"

namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;
class LBMUnitConverter;
class WbWriter;
class Block3D;

//! \brief A class writes macroscopic quantities information to a VTK-file
class WriteMacroscopicQuantitiesPlusMassSimulationObserver : public SimulationObserver
{
public:
    WriteMacroscopicQuantitiesPlusMassSimulationObserver();
    //! \brief Construct WriteMacroscopicQuantitiesPlusMassSimulationObserver object
    //! \pre The Grid3D and UbScheduler objects must exist
    //! \param grid is observable Grid3D object
    //! \param s is UbScheduler object for scheduling of observer
    //! \param path is path of folder for output
    //! \param writer is WbWriter object
    //! \param conv is LBMUnitConverter object
    //! \param comm is Communicator object
    WriteMacroscopicQuantitiesPlusMassSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                          WbWriter *const writer, SPtr<LBMUnitConverter> conv, std::shared_ptr<vf::parallel::Communicator> comm);
    ~WriteMacroscopicQuantitiesPlusMassSimulationObserver() override = default;

    void update(real step) override;

protected:
    //! Collect data for VTK-file
    //! \param step is a time step
    void collectData(real step);
    //! Collect data for VTK-file
    //! \param block is a time step
    void addDataMQ(SPtr<Block3D> block);
    void clearData();

private:
    void init();
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleUInt8> cells;
    std::vector<std::string> datanames;
    std::vector<std::vector<real>> data;
    std::string path;
    WbWriter *writer;
    SPtr<LBMUnitConverter> conv;
    std::vector<std::vector<SPtr<Block3D>>> blockVector;
    int minInitLevel;
    int maxInitLevel;
    int gridRank;
    std::shared_ptr<vf::parallel::Communicator> comm;

    using CalcMacrosFct = void (*)(const real *const &, real &, real &, real &, real &);
    CalcMacrosFct calcMacros;
};

#endif
