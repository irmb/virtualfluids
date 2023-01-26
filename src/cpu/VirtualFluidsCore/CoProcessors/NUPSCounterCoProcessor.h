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
//! \file NUPSCounterCoProcessor.h
//! \ingroup CoProcessors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef NUPSCOUNTERCoProcessor_H_
#define NUPSCOUNTERCoProcessor_H_

#include <PointerDefinitions.h>

#include "CoProcessor.h"
#include "basics/utilities/UbTiming.h"

namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;

//! \class NUPSCounterCoProcessor
//! \brief A class calculates Nodal Updates Per Second (NUPS)
class NUPSCounterCoProcessor : public CoProcessor
{
public:
    //! \brief Construct NUPSCounterCoProcessor object for grid object and scheduler object.
    //! \pre The Grid3D and UbScheduler objects must exist.
    //! \param grid is observable Grid3D object
    //! \param s is UbScheduler object for scheduling of observer
    //! \param numOfThreads is number of threads
    //! \param comm is Communicator object
    NUPSCounterCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, int numOfThreads, std::shared_ptr<vf::mpi::Communicator> comm);
    ~NUPSCounterCoProcessor() override;

    void process(real step) override;

protected:
    //! Collect data for calculation of NUPS
    //! \param step is a time step
    void collectData(real step);
    UbTimer timer;
    int numOfThreads;
    real numberOfNodes;
    real numberOfBlocks;
    real nup;
    real nup_t;
    real nupsStep;
    std::shared_ptr<vf::mpi::Communicator> comm;
};

#endif
