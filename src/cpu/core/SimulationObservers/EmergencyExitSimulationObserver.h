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

#ifndef EmergencyExitSimulationObserver_H
#define EmergencyExitSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>

#include "SimulationObserver.h"

class MPIIORestartSimulationObserver;
namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;

class EmergencyExitSimulationObserver : public SimulationObserver
{
public:
    EmergencyExitSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                             SPtr<MPIIORestartSimulationObserver> rp, std::shared_ptr<vf::parallel::Communicator> comm);
    ~EmergencyExitSimulationObserver() override;

    void update(real step) override;

protected:
    void collectData(real step);
    void writeMetafile(int status);
    bool readMetafile();
    void checkMetafile();

private:
    std::string path;
    std::shared_ptr<vf::parallel::Communicator> comm;
    SPtr<MPIIORestartSimulationObserver> rp;
    std::string metafile;
};

#endif /* EmergencyExitSimulationObserver_H */

//! \}
