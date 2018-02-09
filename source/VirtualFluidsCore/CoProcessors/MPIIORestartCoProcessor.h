#ifndef _MPIIORestartCoProcessor_H_
#define _MPIIORestartCoProcessor_H_

#include <mpi.h>
#include <PointerDefinitions.h>
#include <string>

#include "CoProcessor.h"

class Grid3D;
class UbScheduler;
class Communicator;
class BCProcessor;
class LBMKernel;

//! \class MPIWriteBlocksCoProcessor 
//! \brief Writes the grid each timestep into the files and reads the grip from the files before regenerating  
class MPIIORestartCoProcessor : public CoProcessor
{
   //! \struct GridParam
   //! \brief Structure describes parameters of the grid
   //! \details The structure is nessasary to restore the grid correctly
   struct GridParam
   {
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
   struct Block3d
   {
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
   //! \brief Structure describes parameters of one array of the dataSet that are equal in all blocks
   //! \details The structure used to store some parameters needed to restore dataSet arrays
   struct dataSetParam
   {
      int nx1;
      int nx2;
      int nx3;
      int nx[4]; // nx1, nx2, nx3, nx4
   };

   //! \struct dataSet
   //! \brief Structure containes information identifying the block 
   //! \details The structure is used to find the needed block in the grid when restoring a dataSet
   struct DataSet
   {
      double collFactor;
      double deltaT;
      int x1;
      int x2;
      int x3;
      int level;
      int ghostLayerWidth;
      bool compressible;
      bool withForcing;
   };

   //! \struct dataSetSmall
   //! \brief Structure containes information identifying the block 
   //! \details The structure is used to find the needed block in the grid when restoring a dataSet arrays
   struct DataSetSmall
   {
      int x1;
      int x2;
      int x3;
      int level;
   };

   //! \struct BoundaryCondition
   //! \brief Structure containes information about boundary conditions of the block 
   //! \details The structure is used to write data describing boundary conditions of the blocks when saving the grid 
   //! and to read it when restoring the grid
   struct BoundaryCondition
   {
      long long noslipBoundaryFlags;	//	MPI_LONG_LONG
      long long slipBoundaryFlags;
      long long velocityBoundaryFlags;
      long long densityBoundaryFlags;
      long long wallModelBoundaryFlags;

      float  bcVelocityX1;
      float  bcVelocityX2;
      float  bcVelocityX3;
      float  bcDensity;

      float  bcLodiDensity;
      float  bcLodiVelocityX1;
      float  bcLodiVelocityX2;
      float  bcLodiVelocityX3;
      float  bcLodiLentgh;

      float  nx1, nx2, nx3;
      float q[26];					//	MPI_FLOAT

      char algorithmType;
   };

   //! \struct boundCondParam
   //! \brief Structure describes parameters of the boundaryConditions that are equal in all blocks
   //! \details The structure used to store some parameters needed to restore boundaryConditions arrays
   struct boundCondParam
   {
      int nx1;
      int nx2;
      int nx3;
      int bcindexmatrixCount;	// how many bcindexmatrix-values in one (any) block 
   };

   //! \struct BCAdd
   //! \brief Structure containes information identifying the block 
   //! and some parameters of the arrays of boundary conditions that are equal in all blocks
   //! \details The structure is used to find the needed block in the grid when restoring a dataSet
   //! and to set common parameters
   struct BCAdd
   {
      int x1;		//	to find the right block
      int x2;
      int x3;
      int level;
      int boundCond_count;		//	how many BoundaryCondition-structures are in this block
      int indexContainer_count;	// how many indexContainer-values are in this block
   };

   struct DSArraysPresence
   {
      bool isAverageDensityArrayPresent;
      bool isAverageVelocityArrayPresent;
      bool isAverageFluktuationsArrayPresent;
      bool isAverageTripleArrayPresent;
      bool isShearStressValArrayPresent;
      bool isRelaxationFactorPresent;
   };

public:
   MPIIORestartCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, SPtr<Communicator> comm);
   virtual ~MPIIORestartCoProcessor();
   //! Each timestep writes the grid into the files
   void process(double step);
   //! Reads the grid from the files before grid reconstruction
   void restart(int step);
   //! Writes the blocks of the grid into the file cpBlocks.bin
   void writeBlocks(int step);
   //! Writes the datasets of the blocks into the file cpDataSet.bin
   void writeDataSet(int step);
   void writeAverageDensityArray(int step);
   void writeAverageVelocityArray(int step);
   void writeAverageFluktuationsArray(int step);
   void writeAverageTripleArray(int step);
   void writeShearStressValArray(int step);
   void writeRelaxationFactor(int step);
   //! Writes the boundary conditions of the blocks into the file cpBC.bin
   void writeBoundaryConds(int step);

   //! Reads the blocks of the grid from the file cpBlocks.bin
   void readBlocks(int step);
   //! Reads the datasets of the blocks from the file cpDataSet.bin
   void readDataSet(int step);
   void readAverageDensityArray(int step);
   void readAverageVelocityArray(int step);
   void readAverageFluktuationsArray(int step);
   void readAverageTripleArray(int step);
   void readShearStressValArray(int step);
   void readRelaxationFactor(int step);
   //! Reads the boundary conditions of the blocks from the file cpBC.bin
   void readBoundaryConds(int step);
   //! The function sets LBMKernel
   void setLBMKernel(SPtr<LBMKernel> kernel);
   //!The function sets BCProcessor
   void setBCProcessor(SPtr<BCProcessor> bcProcessor);
   //!The function truncates the data files
   void clearAllFiles(int step);

protected:
   std::string path;
   SPtr<Communicator> comm;

private:
   MPI_Datatype gridParamType, block3dType, arrayPresenceType;
   MPI_Datatype dataSetParamType, dataSetType, dataSetSmallType, dataSetDoubleType;
   MPI_Datatype boundCondParamType, boundCondType, boundCondType1000, boundCondTypeAdd, bcindexmatrixType;

   boundCondParam boundCondParamStr;
   SPtr<LBMKernel> lbmKernel;
   SPtr<BCProcessor> bcProcessor;
};

#endif 
