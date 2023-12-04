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
//! \file QCriterionSimulationObserver.h
//! \ingroup SimulationObservers
//! \author Sonja Uphoff
//=======================================================================================

#ifndef QCriterionSimulationObserver_H
#define QCriterionSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "SimulationObserver.h"
#include "LBMSystem.h"
#include "UbTuple.h"

namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;
class WbWriter;
class Block3D;

//! \brief  Computes the value Q with which vortices can be visualized as isocontours to Q=0, writes to .vtk, For
//! uniform, serial setups only! \details writes at given time intervals specified in scheduler (s)
//!          Processing: paraview, take isolines of entry for Q-criterion vortex detection
//!             Q-Criterion: Visualize Vorteces as regions where Vorticity is larger than strain rate (Hunt, 1988)
//! \author  Sonja Uphoff

class QCriterionSimulationObserver : public SimulationObserver
{
public:
    QCriterionSimulationObserver(SPtr<Grid3D> grid, const std::string &path, WbWriter *const writer, SPtr<UbScheduler> s,
                          std::shared_ptr<vf::parallel::Communicator> comm);
    //! Make update if timestep is write-timestep specified in SPtr<UbScheduler> s
    void update(real step) override;

protected:
    //! Prepare data and write in .vtk file
    void collectData(real step);
    //! Q is computed for all points in a block. Data for writing is added to data and cell vectors.
    void addData(const SPtr<Block3D> block);
    //! After writing to .vtk-file, all vectors are reset
    void clearData();
    //! Computes macroscopic velocities
    void computeVelocity(real *f, real *v, bool compressible);
    //! Computes average and RMS values of macroscopic quantities
    void getNeighborVelocities(int offx, int offy, int offz, int ix1, int ix2, int ix3, const SPtr<Block3D> block,
                               real *vE, real *vW);

private:
    void init();
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleUInt8> cells;
    std::vector<std::string> datanames; // only one entry for QKrit-SimulationObserver: Q
    std::vector<std::vector<double>> data;
    std::vector<std::vector<SPtr<Block3D>>> blockVector;
    int minInitLevel; // go through all levels for block vector of current process from minInitLevel to maxInitLevel
    int maxInitLevel;
    int gridRank; // comm-Rank des aktuellen prozesses
    std::string path;
    WbWriter *writer;
    std::shared_ptr<vf::parallel::Communicator> comm;
    enum Values { xdir = 0, ydir = 1, zdir = 2 }; // labels for the different components
};

#endif
