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
//! \file MPIIODataStructures.h
//! \ingroup Parallel
//! \author Alena Karanchuk
//=======================================================================================

#ifndef _MPI_STRUCTURES_H_
#define _MPI_STRUCTURES_H_

#include "lbm/constants/D3Q27.h"

namespace MPIIODataStructures
{
//! \struct GridParam
//! \brief Structure describes parameters of the grid
//! \details The structure is nessasary to restore the grid correctly
struct GridParam {
    double trafoParams[33]; // not float!!!
    double deltaX;          // not float!!!
    int blockNx1;
    int blockNx2;
    int blockNx3;
    int nx1;
    int nx2;
    int nx3;
    bool periodicX1;
    bool periodicX2;
    bool periodicX3;
    bool active;
    bool transformation;
};

//! \struct Block3d
//! \brief Structure contains information of the block
//! \details The structure is used to write the data describing the block in the grid when saving the grid
//! and to read it when restoring the grid
struct Block3d {
    int x1;
    int x2;
    int x3;
    int bundle;
    int rank;
    int lrank;
    int part;
    int globalID;
    int localID;
    int level;
    int interpolationFlagCF;
    int interpolationFlagFC;
    int counter;
    bool active;
};

//! \struct dataSetParam
//! \brief Structure describes parameters of the dataSet that are equal in all blocks
//! \details The structure used to store some parameters needed to restore dataSet arrays
struct dataSetParam {
    int nx1;
    int nx2;
    int nx3;
    int nx[4]; // nx1, nx2, nx3, nx4
};

//! \struct DataSetRestart
//! \brief Structure describes parameters of the dataSet in MPIIORestartSimulationObserver format
//! \details The structure is used when reading from the file
struct DataSetRestart {
    double collFactor; // not float!!!
    double deltaT;     // not float!!!
    double collFactorL; // for Multiphase model  // not float!!!
    double collFactorG; // for Multiphase model // not float!!!
    double densityRatio;// for Multiphase model // not float!!!
    int x1;
    int x2;
    int x3;
    int level;
    int ghostLayerWidth;
    bool compressible;
    bool withForcing;
};

//! \struct DataSetMigration
//! \brief Structure describes parameters of the dataSet in MPIIOMigrationSimulationObserver format
//! \details The structure is used to find the needed block in the grid when restoring a dataSet
struct DataSetMigration {
    double collFactor;  // not float!!!
    double deltaT;      // not float!!!
    double collFactorL; // for Multiphase model
    double collFactorG; // for Multiphase model
    double densityRatio;// for Multiphase model
    int globalID;
    int ghostLayerWidth;
    bool compressible;
    bool withForcing;
};

//! \struct DataSetSmallRead
//! \brief Structure describes parameters of the DataSetSmall in MPIIORestartSimulationObserver format
//! \details The structure is used when reading from the file
struct DataSetSmallRestart {
    int x1;
    int x2;
    int x3;
    int level;
};
//! \struct dataSetSmall
//! \brief Structure containes information identifying the block in MPIIOMigrationSimulationObserver format
//! \details The structure is used to find the needed block in the grid when restoring a dataSet arrays
struct DataSetSmallMigration {
    int globalID;
};

//! \struct BoundaryCondition
//! \brief Structure containes information about boundary conditions of the block
//! \details The structure is used to write data describing boundary conditions of the blocks when saving the grid
//! and to read it when restoring the grid
struct BoundaryCondition {
    long long noslipBoundaryFlags; //    MPI_LONG_LONG
    long long slipBoundaryFlags;
    long long velocityBoundaryFlags;
    long long densityBoundaryFlags;
    long long wallModelBoundaryFlags;

    float bcVelocityX1;    //  not double!!!
    float bcVelocityX2;    //  not double!!!
    float bcVelocityX3;    //  not double!!!
    float bcDensity;       //  not double!!!
    float bcPhaseField;    //  not double!!!

    float nx1, nx2, nx3;    //  not double!!!
    float  q[26];           //  not double!!!

    char bcStrategyKey;
};

//! \struct boundCondParam
//! \brief Structure describes parameters of the boundaryConditions that are equal in all blocks
//! \details The structure used to store some parameters needed to restore boundaryConditions arrays
struct boundCondParam {
    int nx1;
    int nx2;
    int nx3;
    int bcindexmatrixCount; // how many bcindexmatrix-values in one (any) block
};

//! \struct BCAddRead
//! \brief Structure describes parameters of the BCAdd in MPIIORestartSimulationObserver format
//! \details The structure is used when reading from the file
struct BCAddRestart {
    int x1; //    to find the right block
    int x2;
    int x3;
    int level;
    int boundCond_count;      //    how many BoundaryCondition-structures are in this block
    int indexContainer_count; // how many indexContainer-values are in this block
};

//! \struct BCAdd
//! \brief Structure containes information identifying the block and some parameters of the arrays
//! \of boundary conditions that are equal in all blocks in MPIIOMigrationSimulationObserver format
//! \details The structure is used to find the needed block in the grid when restoring a dataSet
//! and to set common parameters
struct BCAddMigration {
    int globalID;
    int boundCond_count;      //    how many BoundaryCondition-structures are in this block
    int indexContainer_count; // how many indexContainer-values are in this block
};

struct DSArraysPresence {
    bool isAverageDensityArrayPresent;
    bool isAverageVelocityArrayPresent;
    bool isAverageFluktuationsArrayPresent;
    bool isAverageTripleArrayPresent;
    bool isShearStressValArrayPresent;
    bool isRelaxationFactorPresent;
    bool isPhaseField1Present;
    bool isPhaseField2Present;
    bool isPressureFieldPresent;
};
} // namespace MPIIODataStructures
#endif