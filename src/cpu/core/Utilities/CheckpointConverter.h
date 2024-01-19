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
//! \addtogroup cpu_Utilities Utilities
//! \ingroup cpu_core core
//! \{
//! \author Alena Karanchuk
//=======================================================================================

#ifndef _UTILITACONVERTOR_H_
#define _UTILITACONVERTOR_H_

#include "MPIIODataStructures.h"
#include <PointerDefinitions.h>
#include <mpi.h>
#include <string>
#include <vector>

class Grid3D;
namespace vf::parallel {class Communicator;}

//! \class UtilConvertor
//! \brief Converts timestep data from MPIIORestartSimulationObserver format into MPIIOMigrationSimulationObserver format
class CheckpointConverter
{
public:
    CheckpointConverter(SPtr<Grid3D> grid, const std::string &path, std::shared_ptr<vf::parallel::Communicator> comm);
    virtual ~CheckpointConverter();

    void convert(int step, int procCount);
    void convertBlocks(int step, int procCount);
    void convertDataSet(int step, int procCount);
    void convertBC(int step, int procCount);
    void convert___Array(int step, int procCount, std::string filenameR, std::string filenameW);

protected:
    std::string path;
    std::shared_ptr<vf::parallel::Communicator> comm;
    SPtr<Grid3D> grid;

private:
    MPI_Datatype gridParamType, block3dType;
    MPI_Datatype dataSetParamType, dataSetTypeRead, dataSetTypeWrite;
    MPI_Datatype boundCondType, boundCondType1000;

    MPIIODataStructures::boundCondParam boundCondParamStr;
};

#endif

//! \}
