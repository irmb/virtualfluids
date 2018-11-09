#ifndef _UTILITACONVERTOR_H_
#define _UTILITACONVERTOR_H_

#include <mpi.h>
#include <PointerDefinitions.h>
#include <string>
#include <vector>

class Grid3D;
class Communicator;

//! \class UtilConvertor 
//! \brief Converts timestep data from MPIIORestartCoProcessor format into MPIIOMigrationCoProcessor format  
class CheckpointConverter
{
   //! \struct GridParam
   //! \brief Structure describes parameters of the grid
   //! \details The structure is necessary to restore the grid correctly
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
   //! \brief Structure describes parameters of the dataSet that are equal in all blocks
   //! \details The structure used to store some parameters needed to restore dataSet arrays
   struct dataSetParam
   {
      int nx1;
      int nx2;
      int nx3;
      int nx[4]; //nx1, nx2, nx3, nx4
   };

   //! \struct DataSetRead
   //! \brief Structure describes parameters of the dataSet in MPIIORestartCoProcessor format
   //! \details The structure is used when reading from the file
   struct DataSetRead
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

   //! \struct DataSetWrite
   //! \brief Structure describes parameters of the dataSet in MPIIOMigrationCoProcessor format
   //! \details The structure is used when writing to the file
   struct DataSetWrite
   {
      double collFactor;
      double deltaT;
      int globalID;
      int ghostLayerWidth;
      bool compressible;
      bool withForcing;
   };

   //! \struct DataSetSmallRead
   //! \brief Structure describes parameters of the DataSetSmall in MPIIORestartCoProcessor format
   //! \details The structure is used when reading from the file
   struct DataSetSmallRead
   {
      int x1;
      int x2;
      int x3;
      int level;
   };

   //! \struct DataSetSmallWrite
   //! \brief Structure describes parameters of the DataSetSmall in MPIIOMigrationCoProcessor format
   //! \details The structure is used when writing to the file
   struct DataSetSmallWrite
   {
      int globalID;
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

   //! \struct BCAddRead
   //! \brief Structure describes parameters of the BCAdd in MPIIORestartCoProcessor format
   //! \details The structure is used when reading from the file
   struct BCAddRead
   {
      int x1;		//	to find the right block
      int x2;
      int x3;
      int level;
      int boundCond_count;		//	how many BoundaryCondition-structures are in this block
      int indexContainer_count;	// how many indexContainer-values are in this block
   };

   //! \struct BCAddWrite
   //! \brief Structure describes parameters of the BCAdd in MPIIOMigrationCoProcessor format
   //! \details The structure is used when writing to the file
   struct BCAddWrite
   {
      int globalID;
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
   CheckpointConverter(SPtr<Grid3D> grid, const std::string& path, SPtr<Communicator> comm);
   virtual ~CheckpointConverter();

   void convert(int step, int procCount);
   void convertBlocks(int step, int procCount);
   void convertDataSet(int step, int procCount);
   void convertBC(int step, int procCount);
   void convert___Array(int step, int procCount, std::string filenameR, std::string filenameW);

protected:
   std::string path;
   SPtr<Communicator> comm;
   SPtr<Grid3D> grid;

private:
   MPI_Datatype gridParamType, block3dType;
   MPI_Datatype dataSetParamType, dataSetDoubleType, dataSetTypeRead, dataSetTypeWrite, dataSetSmallTypeRead, dataSetSmallTypeWrite;
   MPI_Datatype boundCondType, boundCondType1000, boundCondParamType, boundCondTypeAddRead, boundCondTypeAddWrite, bcindexmatrixType;

   boundCondParam boundCondParamStr;
};

#endif 
