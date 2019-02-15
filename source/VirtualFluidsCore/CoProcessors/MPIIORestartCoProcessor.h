#ifndef _MPIIORestartCoProcessor_H_
#define _MPIIORestartCoProcessor_H_

#include <mpi.h>
#include <PointerDefinitions.h>
#include <string>

#include "CoProcessor.h"
#include "MPIIODataStructures.h"

class Grid3D;
class UbScheduler;
class Communicator;
class BCProcessor;
class LBMKernel;

//! \class MPIWriteBlocksCoProcessor 
//! \brief Writes the grid each timestep into the files and reads the grip from the files before regenerating  
class MPIIORestartCoProcessor : public CoProcessor
{
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
   //!The function write a time step of last check point
   void writeCpTimeStep(int step);
   //!The function read a time step of last check point
   int readCpTimeStep();

protected:
   std::string path;
   SPtr<Communicator> comm;

private:
   MPI_Datatype gridParamType, block3dType, arrayPresenceType;
   MPI_Datatype dataSetParamType, dataSetType, dataSetSmallType, dataSetDoubleType;
   MPI_Datatype boundCondParamType, boundCondType, boundCondType1000, boundCondTypeAdd, bcindexmatrixType;

   MPIIODataStructures::boundCondParam boundCondParamStr;
   SPtr<LBMKernel> lbmKernel;
   SPtr<BCProcessor> bcProcessor;
};

#endif 
