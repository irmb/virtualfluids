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
//! \file MPIWriteBlocksBESimulationObserver.h
//! \ingroup SimulationObservers
//! \author Alena Karanchuk
//=======================================================================================
#ifndef _MPIIOMigrationBESimulationObserver_H_
#define _MPIIOMigrationBESimulationObserver_H_

#include <mpi.h>
#include <string>
#include <vector>

#include "MPIIOSimulationObserver.h"
#include "MPIIODataStructures.h"

class Grid3D;
class UbScheduler;
namespace vf::parallel {class Communicator;}
class BCSet;
class LBMKernel;
class Grid3DVisitor;

//! \class MPIWriteBlocksBESimulationObserver
//! \brief Writes the grid each timestep into the files and reads the grip from the files before regenerating
class MPIIOMigrationBESimulationObserver : public MPIIOSimulationObserver
{
    enum Arrays {
        AverageDensity      = 1,
        AverageVelocity     = 2,
        AverageFluktuations = 3,
        AverageTriple       = 4,
        ShearStressVal      = 5,
        RelaxationFactor    = 6,
        PhaseField1         = 7,
        PhaseField2         = 8,
        PressureField = 9
    };

public:
    MPIIOMigrationBESimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<Grid3DVisitor> mV, const std::string &path,
                                std::shared_ptr<vf::parallel::Communicator> comm);
    ~MPIIOMigrationBESimulationObserver() override;
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
    // void writeAverageDensityArray(int step);
    // void writeAverageVelocityArray(int step);
    // void writeAverageFluktuationsArray(int step);
    // void writeAverageTripleArray(int step);
    // void writeShearStressValArray(int step);
    // void writeRelaxationFactor(int step);
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
    void setNu(real nu);
    void setNuLG(real cfL, real cfG);
    void setDensityRatio(real dr);

    void blocksExchange(int tagN, int ind1, int ind2, int doubleCountInBlock, std::vector<real> &pV,
                        std::vector<real> *rawDataReceive);

private:
    // MPI_Datatype gridParamType, block3dType;
    //   MPI_Datatype dataSetType, dataSetSmallType;
    MPI_Datatype dataSetDoubleType;
    //   MPI_Datatype boundCondParamType, boundCondTypeAdd;
    MPI_Datatype bcindexmatrixType;
    MPI_Datatype sendBlockDoubleType, sendBlockIntType;

    MPIIODataStructures::boundCondParam boundCondParamStr;
    SPtr<LBMKernel> lbmKernel;
    SPtr<BCSet> bcSet;
    SPtr<Grid3DVisitor> metisVisitor;
    real nue;
    real nuL;
    real nuG;
    real densityRatio;

};

#endif
