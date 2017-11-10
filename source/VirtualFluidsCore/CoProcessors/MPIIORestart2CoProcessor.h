#ifndef _MPIIORestart2CoProcessor_H_
#define _MPIIORestart2CoProcessor_H_

#include "mpi.h"

#include "CoProcessor.h"
#include "Communicator.h"
#include "WbWriter.h"

#include <boost/shared_ptr.hpp>

class MPIIORestart2CoProcessor;
typedef boost::shared_ptr<MPIIORestart2CoProcessor> MPIIORestart2CoProcessorPtr;

//! \class MPIWriteBlocksCoProcessor 
//! \brief Writes the grid each timestep into the files and reads the grip from the files before regenerating  
class MPIIORestart2CoProcessor: public CoProcessor
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

   //! \struct blockParam
   //! \brief Structure describes parameters of the block that are equal in all blocks
   //! \details The structure used to store some parameters needed to restore dataSet arrays
   struct BlockParam
    {
       int nx1;   //	to find the right block
       int nx2;
       int nx3;
       int nx[9][4]; // 9 arrays x (nx1, nx2, nx3, nx4)
       int doubleCountInBlock;   // how many double-values are in all arrays dataSet in one (any) block
       int bcindexmatrix_count;	// how many bcindexmatrix-values are in one (any) block 
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

   //! \struct dataSet
   //! \brief Structure containes information identifying the block 
   //! \details The structure is used to find the needed block in the grid when restoring a dataSet
   struct DataSet
	{
      double collFactor;
      double deltaT;
      int globalID;
      int ghostLayerWidth;
      bool compressible;
      bool withForcing;
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
		
		float  nx1,nx2,nx3;
		float q[26];					//	MPI_FLOAT

      char algorithmType;
   };

   //! \struct BCAdd
   //! \brief Structure containes information identifying the block 
   //! and some parameters of the arrays of boundary conditions that are equal in all blocks
   //! \details The structure is used to find the needed block in the grid when restoring a dataSet
   //! and to set common parameters
   struct BCAdd
	{
      int globalID;
		//int x1;		//	to find the right block
		//int x2;		
		//int x3;		
		//int level;	
      int boundCond_count;		//	how many BoundaryCondition-structures are in this block
      int indexContainer_count;	// how many indexContainer-values are in this block
   };

public:
   MPIIORestart2CoProcessor(Grid3DPtr grid, UbSchedulerPtr s, const std::string& path, CommunicatorPtr comm);
   virtual ~MPIIORestart2CoProcessor();
   //! Each timestep writes the grid into the files
   void process(double step);
   //! Reads the grid from the files before grid reconstruction
   void restart(int step);
   //! Writes the blocks of the grid into the file outputBlocks.bin
   void writeBlocks(int step);
   //! Writes the datasets of the blocks into the file outputDataSet.bin
   void writeDataSet(int step);
   //! Writes the boundary conditions of the blocks into the file outputBoundCond.bin
   void writeBoundaryConds(int step);
   //! Reads the blocks of the grid from the file outputBlocks.bin
   void readBlocks(int step);
   //! Reads the datasets of the blocks from the file outputDataSet.bin
   void readDataSet(int step);
   //! Reads the boundary conditions of the blocks from the file outputBoundCond.bin
   void readBoundaryConds(int step);
   //! The function sets number of ranks that read simultaneously 
   void setChunk(int val);
   //! The function sets LBMKernel
   void setLBMKernel(LBMKernelPtr kernel);
   //!The function sets BCProcessor
   void setBCProcessor(BCProcessorPtr bcProcessor);
   //!The function truncates the data files
   void clearAllFiles(int step);

protected:
   std::string path;
   CommunicatorPtr comm;
   bool mpiTypeFreeFlag;

private:
	MPI_Datatype gridParamType, blockParamType, block3dType, dataSetType, dataSetDoubleType, boundCondType, boundCondType1000, boundCondTypeAdd, bcindexmatrixType;
   BlockParam blockParamStr;
   int chunk;
   LBMKernelPtr lbmKernel;
   BCProcessorPtr bcProcessor;

};

#endif 
