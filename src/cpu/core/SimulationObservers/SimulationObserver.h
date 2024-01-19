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

#ifndef SimulationObserver_H
#define SimulationObserver_H

#include <PointerDefinitions.h>
#include "lbm/constants/D3Q27.h"

class Grid3D;
class UbScheduler;

//! \class SimulationObserver
//! \brief An abstract class implements observer design pettern
class SimulationObserver
{
public:
    //! Class default constructor
    SimulationObserver();
    //! \brief Construct SimulationObserver object for grid object and scheduler object.
    //! \pre The Grid3D and UbScheduler objects must exist.
    //! \param grid is observable Grid3D object
    //! \param s is UbScheduler object for scheduling of observer
    //! \details
    //! Class SimulationObserver implements the observer design pettern. SimulationObserver object is observer.  Grid3D object is
    //! observable.
    SimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s);
    //! Class destructor
    virtual ~SimulationObserver();
    //! \brief Updates observer
    //! \param step is the actual time step
    virtual void update(real step) = 0;

protected:
    SPtr<Grid3D> grid;
    SPtr<UbScheduler> scheduler;
};
#endif

//! \}
