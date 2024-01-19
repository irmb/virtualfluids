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
#ifndef _MPIIOMigrationSimulationObserver_H_
#define _MPIIOMigrationSimulationObserver_H_

#include <mpi.h>
#include <string>

#include "MPIIOSimulationObserver.h"
#include "MPIIODataStructures.h"

class Grid3D;
class UbScheduler;
namespace vf::parallel {class Communicator;}
class BCSet;
class LBMKernel;
class Grid3DVisitor;

//! \class MPIWriteBlocksSimulationObserver
//! \brief Writes the grid each timestep into the files and reads the grip from the files before regenerating
class MPIIOMigrationSimulationObserver : public MPIIOSimulationObserver
{
public:
    enum Arrays {
        AverageDensity      = 1,
        AverageVelocity     = 2,
        AverageFluktuations = 3,
        AverageTriple       = 4,
        ShearStressVal      = 5,
        RelaxationFactor    = 6,
        PhaseField1         = 7,
        PhaseField2         = 8,
        PressureField       = 9
    };

    MPIIOMigrationSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<Grid3DVisitor> mV, const std::string &path, std::shared_ptr<vf::parallel::Communicator> comm);
    ~MPIIOMigrationSimulationObserver() override;
    //! Each timestep writes the grid into the files
    void update(real step) override;
    //! Reads the grid from the files before grid reconstruction
    void restart(int step);
    //! Writes the blocks of the grid into the file cpBlocks.bin
    void writeBlocks(int step);
    //! Writes the datasets of the blocks into the file cpDataSet.bin
    void writeDataSet(int step);
    void write4DArray(int step, Arrays arrType, std::string fname);
    void write3DArray(int step, Arrays arrType, std::string fname);
    //   void writeAverageDensityArray(int step);
    //   void writeAverageVelocityArray(int step);
    //   void writeAverageFluktuationsArray(int step);
    //   void writeAverageTripleArray(int step);
    //   void writeShearStressValArray(int step);
    //   void writeRelaxationFactor(int step);
    //! Writes the boundary conditions of the blocks into the file cpBC.bin
    void writeBoundaryConds(int step);

    //! Reads the blocks of the grid from the file cpBlocks.bin
    void readBlocks(int step);
    //! Reads the datasets of the blocks from the file cpDataSet.bin
    void readDataSet(int step);
    void readArray(int step, Arrays arrType, std::string fname);
    //   void readAverageDensityArray(int step);
    //   void readAverageVelocityArray(int step);
    //   void readAverageFluktuationsArray(int step);
    //   void readAverageTripleArray(int step);
    //   void readShearStressValArray(int step);
    //   void readRelaxationFactor(int step);
    //! Reads the boundary conditions of the blocks from the file cpBC.bin
    void readBoundaryConds(int step);
    //! The function sets LBMKernel
    void setLBMKernel(SPtr<LBMKernel> kernel);
    //! The function sets BCSet
    void setBCSet(SPtr<BCSet> BCSet);
    //! The function truncates the data files
    void clearAllFiles(int step);
    // void setNu(real nu);

private:
    // MPI_Datatype gridParamType, block3dType;
    MPI_Datatype dataSetType, dataSetSmallType, dataSetDoubleType;
    MPI_Datatype boundCondParamType, boundCondTypeAdd, bcindexmatrixType;

    MPIIODataStructures::boundCondParam boundCondParamStr;
    SPtr<LBMKernel> lbmKernel;
    SPtr<BCSet> bcSet;
    SPtr<Grid3DVisitor> metisVisitor;

};

#endif

//! \}
