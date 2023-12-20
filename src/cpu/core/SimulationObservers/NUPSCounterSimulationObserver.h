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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_SimulationObservers SimulationObservers
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef NUPSCOUNTERSimulationObserver_H_
#define NUPSCOUNTERSimulationObserver_H_

#include <PointerDefinitions.h>

#include "SimulationObserver.h"

#include <basics/Timer/Timer.h>

namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;

//! \class NUPSCounterSimulationObserver
//! \brief A class calculates Nodal Updates Per Second (NUPS)
class NUPSCounterSimulationObserver : public SimulationObserver
{
public:
    //! \brief Construct NUPSCounterSimulationObserver object for grid object and scheduler object.
    //! \pre The Grid3D and UbScheduler objects must exist.
    //! \param grid is observable Grid3D object
    //! \param s is UbScheduler object for scheduling of observer
    //! \param numOfThreads is number of threads
    //! \param comm is Communicator object
    NUPSCounterSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, int numOfThreads, std::shared_ptr<vf::parallel::Communicator> comm);
    ~NUPSCounterSimulationObserver() override;

    void update(real step) override;

protected:
    //! Collect data for calculation of NUPS
    //! \param step is a time step
    void collectData(real step);
    vf::basics::Timer timer;
    int numOfThreads;
    real numberOfNodes;
    real numberOfBlocks;
    real nup;
    real nup_t;
    real nupsStep;
    std::shared_ptr<vf::parallel::Communicator> comm;
};

#endif

//! \}
