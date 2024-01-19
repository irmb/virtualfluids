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
#ifndef TimeDependentBCSimulationObserver_H
#define TimeDependentBCSimulationObserver_H

#include <PointerDefinitions.h>
#include <vector>

#include "SimulationObserver.h"

class Interactor3D;
class Grid3D;

//! \brief The class update interactors depend of time step.
//! \details TimeDependentBCSimulationObserver update every time step information in BCs throw Interactors
//! \author Sonja Uphoff, Kostyantyn Kucher
class TimeDependentBCSimulationObserver : public SimulationObserver
{
public:
    TimeDependentBCSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s);
    ~TimeDependentBCSimulationObserver() override;

    void update(real step) override;

    //! add interactors to SimulationObserver
    void addInteractor(SPtr<Interactor3D> interactor);

private:
    std::vector<SPtr<Interactor3D>> interactors;
};

#endif /* TimeDependentBCSimulationObserver_H */

//! \}
