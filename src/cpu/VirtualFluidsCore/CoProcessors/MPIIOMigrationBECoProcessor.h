#ifndef _MPIIOMigrationBECoProcessor_H_
#define _MPIIOMigrationBECoProcessor_H_

#include <mpi.h>
#include <string>
#include <vector>

#include "MPIIOCoProcessor.h"
#include "MPIIODataStructures.h"

class Grid3D;
class UbScheduler;
class Communicator;
class BCProcessor;
class LBMKernel;

//! \class MPIWriteBlocksBECoProcessor 
//! \brief Writes the grid each timestep into the files and reads the grip from the files before regenerating  
class MPIIOMigrationBECoProcessor : public MPIIOCoProcessor
{
  enum Arrays { AverageDensity = 1, AverageVelocity = 2, AverageFluktuations = 3, AverageTriple = 4, ShearStressVal = 5, RelaxationFactor = 6};

public:

   MPIIOMigrationBECoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, SPtr<Communicator> comm);
   virtual ~MPIIOMigrationBECoProcessor();
   //! Each timestep writes the grid into the files
   void process(double step);
   //! Reads the grid from the files before grid reconstruction
   void restart(int step);
   //! Writes the blocks of the grid into the file cpBlocks.bin
   void writeBlocks(int step);
   //! Writes the datasets of the blocks into the file cpDataSet.bin
   void writeDataSet(int step);
   void write4DArray(int step, Arrays arrType, std::string fname);
   void write3DArray(int step, Arrays arrType, std::string fname);
  // void writeAverageDensityArray(int step);
   //void writeAverageVelocityArray(int step);
   //void writeAverageFluktuationsArray(int step);
   //void writeAverageTripleArray(int step);
   //void writeShearStressValArray(int step);
   //void writeRelaxationFactor(int step);
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
   //!The function sets BCProcessor
   void setBCProcessor(SPtr<BCProcessor> bcProcessor);
   //!The function truncates the data files
   void clearAllFiles(int step);
   void setNu(double nu);

   void blocksExchange(int tagN, int ind1, int ind2, int doubleCountInBlock, std::vector<double>& pV, std::vector<double>* rawDataReceive);

protected:
   //std::string path;
   //SPtr<Communicator> comm;

private:
   //MPI_Datatype gridParamType, block3dType;
//   MPI_Datatype dataSetType, dataSetSmallType;
   MPI_Datatype dataSetDoubleType;
//   MPI_Datatype boundCondParamType, boundCondTypeAdd;
   MPI_Datatype bcindexmatrixType;
   MPI_Datatype sendBlockDoubleType, sendBlockIntType;

   MPIIODataStructures::boundCondParam boundCondParamStr;
   SPtr<LBMKernel> lbmKernel;
   SPtr<BCProcessor> bcProcessor;
   double nue;

};

#endif 
