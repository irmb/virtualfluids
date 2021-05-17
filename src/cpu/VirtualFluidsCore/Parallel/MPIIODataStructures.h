#ifndef _MPI_STRUCTURES_H_
#define _MPI_STRUCTURES_H_

namespace MPIIODataStructures
{
//! \struct GridParam
//! \brief Structure describes parameters of the grid
//! \details The structure is nessasary to restore the grid correctly
struct GridParam {
    double trafoParams[33];
    double deltaX;
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
//! \brief Structure describes parameters of the dataSet in MPIIORestartCoProcessor format
//! \details The structure is used when reading from the file
struct DataSetRestart {
    double collFactor;
    double deltaT;
    double collFactorL; // for Multiphase model
    double collFactorG; // for Multiphase model
    double densityRatio;// for Multiphase model
    int x1;
    int x2;
    int x3;
    int level;
    int ghostLayerWidth;
    bool compressible;
    bool withForcing;
};

//! \struct DataSetMigration
//! \brief Structure describes parameters of the dataSet in MPIIOMigrationCoProcessor format
//! \details The structure is used to find the needed block in the grid when restoring a dataSet
struct DataSetMigration {
    double collFactor;
    double deltaT;
    double collFactorL; // for Multiphase model
    double collFactorG; // for Multiphase model
    double densityRatio;// for Multiphase model
    int globalID;
    int ghostLayerWidth;
    bool compressible;
    bool withForcing;
};

//! \struct DataSetSmallRead
//! \brief Structure describes parameters of the DataSetSmall in MPIIORestartCoProcessor format
//! \details The structure is used when reading from the file
struct DataSetSmallRestart {
    int x1;
    int x2;
    int x3;
    int level;
};
//! \struct dataSetSmall
//! \brief Structure containes information identifying the block in MPIIOMigrationCoProcessor format
//! \details The structure is used to find the needed block in the grid when restoring a dataSet arrays
struct DataSetSmallMigration {
    int globalID;
};

//! \struct BoundaryCondition
//! \brief Structure containes information about boundary conditions of the block
//! \details The structure is used to write data describing boundary conditions of the blocks when saving the grid
//! and to read it when restoring the grid
struct BoundaryCondition {
    long long noslipBoundaryFlags; //	MPI_LONG_LONG
    long long slipBoundaryFlags;
    long long velocityBoundaryFlags;
    long long densityBoundaryFlags;
    long long wallModelBoundaryFlags;

    float bcVelocityX1;
    float bcVelocityX2;
    float bcVelocityX3;
    float bcDensity;

    float bcLodiDensity;
    float bcLodiVelocityX1;
    float bcLodiVelocityX2;
    float bcLodiVelocityX3;
    float bcLodiLentgh;

    float nx1, nx2, nx3;
    float q[26];

    char algorithmType;
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
//! \brief Structure describes parameters of the BCAdd in MPIIORestartCoProcessor format
//! \details The structure is used when reading from the file
struct BCAddRestart {
    int x1; //	to find the right block
    int x2;
    int x3;
    int level;
    int boundCond_count;      //	how many BoundaryCondition-structures are in this block
    int indexContainer_count; // how many indexContainer-values are in this block
};

//! \struct BCAdd
//! \brief Structure containes information identifying the block and some parameters of the arrays
//! \of boundary conditions that are equal in all blocks in MPIIOMigrationCoProcessor format
//! \details The structure is used to find the needed block in the grid when restoring a dataSet
//! and to set common parameters
struct BCAddMigration {
    int globalID;
    int boundCond_count;      //	how many BoundaryCondition-structures are in this block
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
};
} // namespace MPIIODataStructures
#endif