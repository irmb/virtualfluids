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
//! \file WriteGbObjectsSimulationObserver.h
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef WriteGbObjectsSimulationObserver_h__
#define WriteGbObjectsSimulationObserver_h__

#include "SimulationObserver.h"
#include "UbTuple.h"

#include <vector>

class GbObject3D;
namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;
class WbWriter;

//! \brief     Writes geometry objects as VTK unstructured grid.
//! \details   Writes geometry objects as VTK unstructured grid. Use addGbObject() for add a GbObjects.
//! \author    Konstantin Kutscher
//! \date      December 2018

class WriteGbObjectsSimulationObserver : public SimulationObserver
{
public:
    WriteGbObjectsSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, WbWriter *const writer,
                              std::shared_ptr<vf::parallel::Communicator> comm);
    ~WriteGbObjectsSimulationObserver() override;
    //! calls collectData.
    void update(real step) override;
    //! adds geometry object
    void addGbObject(SPtr<GbObject3D> object);

protected:
    void collectData(real step);

private:
    std::vector<SPtr<GbObject3D>> objects;
    std::string path;
    WbWriter *writer;
    std::shared_ptr<vf::parallel::Communicator> comm;
};

#endif // WriteGbObjectsSimulationObserver_h__