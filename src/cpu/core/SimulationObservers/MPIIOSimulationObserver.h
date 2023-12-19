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
//! \author Alena Karanchuk
//=======================================================================================
#ifndef _MPIIOSimulationObserver_H_
#define _MPIIOSimulationObserver_H_

#include "SimulationObserver.h"
#include <PointerDefinitions.h>
#include <mpi.h>
#include <string>

class Grid3D;
class UbScheduler;
namespace vf::parallel {class Communicator;}

//! \class MPIWriteBlocksBESimulationObserver
//! \brief Writes the grid each timestep into the files and reads the grip from the files before regenerating
class MPIIOSimulationObserver : public SimulationObserver
{
public:
    MPIIOSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, std::shared_ptr<vf::parallel::Communicator> comm);
    ~MPIIOSimulationObserver() override;

    //! Each timestep writes the grid into the files
    void update(real step) override = 0;

    //! Writes the blocks of the grid into the file cpBlocks.bin
    void writeBlocks(int step);

    //! Reads the blocks of the grid from the file cpBlocks.bin
    void readBlocks(int step);

    //! The function truncates the data files
    void clearAllFiles(int step);

    //! The function write a time step of last check point
    void writeCpTimeStep(int step);
    //! The function read a time step of last check point
    int readCpTimeStep();

protected:
    std::string path;
    std::shared_ptr<vf::parallel::Communicator> comm;
    MPI_Datatype gridParamType, block3dType, dataSetParamType, boundCondType, arrayPresenceType;
};
#endif // ! _MPIIOSimulationObserver_H_
#define _MPIIOSimulationObserver_H_

//! \}
