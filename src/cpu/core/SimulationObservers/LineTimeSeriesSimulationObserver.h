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
//! \file LineTimeSeriesSimulationObserver.h
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef LineTimeSeriesSimulationObserver_h__
#define LineTimeSeriesSimulationObserver_h__

#include <PointerDefinitions.h>
#include <string>

#include <mpi.h>

#include "SimulationObserver.h"
#include "LBMSystem.h"

namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;
class GbLine3D;

//! \brief  Writes to .csv file time series for a line in x1 direction.
//! \details It can be used to compute for given time range  the time averaged two-point correlations for a line. <br>
//!  \f$ R_{ij}(x_{a},x{b},t) = <u_{i}(x_{a},t)u_{j}(x_{a}+r,t)> \f$   <br>
//
//! \author  Konstantin Kutscher

class LineTimeSeriesSimulationObserver : public SimulationObserver
{
public:
    enum Direction { X1, X2, X3 };

public:
    LineTimeSeriesSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, SPtr<GbLine3D> line,
                              int level, std::shared_ptr<vf::parallel::Communicator> comm);
    ~LineTimeSeriesSimulationObserver() override = default;

    void update(real step) override;
    void writeLine(const std::string &path);

protected:
    void collectData();

private:
    std::string path;
    std::string fname;
    bool root;
    SPtr<GbLine3D> line;
    // function pointer
    using CalcMacrosFct = void (*)(const real *const &, real &, real &, real &, real &);
    CalcMacrosFct calcMacros;
    int blocknx;
    int blockix1;
    int blockix2;
    int blockix3;
    int level;
    int ix1;
    int ix2;
    int ix3;
    int length;
    MPI_Comm mpi_comm;
    int numOfProc;
    int gridRank;
    Direction dir;
};
#endif // LineTimeSeriesSimulationObserver_h__
