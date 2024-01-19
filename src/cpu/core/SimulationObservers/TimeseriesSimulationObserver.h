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
//! \author Sonja Uphoff
//=======================================================================================


#ifndef TimeseriesSimulationObserver_H
#define TimeseriesSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>

#include "SimulationObserver.h"

namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;
class IntegrateValuesHelper;

//! \brief     Writes timeseries of density and velocity to a file.
//! \details   Uses Integrate values helper, scheduler must be set in testcase.
//! \author    Sonja Uphoff
//! \date      May 2013

class TimeseriesSimulationObserver : public SimulationObserver
{
public:
    TimeseriesSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<IntegrateValuesHelper> h1,
                          const std::string &path, std::shared_ptr<vf::parallel::Communicator> comm);
    ~TimeseriesSimulationObserver() override;

    //! calls collectData.
    void update(real step) override;

protected:
    void collectData(real step);

    //! object that can compute spacial average values in 3D-subdomain.
    SPtr<IntegrateValuesHelper> h1;
    std::shared_ptr<vf::parallel::Communicator> comm;

private:
    std::string path; //! output filename, e.g.  pathname + "/steps/timeseries"
    std::string fname;
};

#endif /* TimeseriesSimulationObserver_H */

//! \}
