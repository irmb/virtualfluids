#ifndef  _MPIIOCoProcessor_H_
#define _MPIIOCoProcessor_H_

#include "CoProcessor.h"
#include <PointerDefinitions.h>
#include <string>
#include <mpi.h>

class Grid3D;
class UbScheduler;
class Communicator;

//! \class MPIWriteBlocksBECoProcessor 
//! \brief Writes the grid each timestep into the files and reads the grip from the files before regenerating  
class MPIIOCoProcessor : public CoProcessor
{
public:
   MPIIOCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, SPtr<Communicator> comm); 
   virtual ~MPIIOCoProcessor();

   //! Each timestep writes the grid into the files
   virtual void process(double step) = 0;

   //! Writes the blocks of the grid into the file cpBlocks.bin
   void writeBlocks(int step);

   //! Reads the blocks of the grid from the file cpBlocks.bin
   void readBlocks(int step);

   //!The function truncates the data files
   void clearAllFiles(int step);
protected:
   std::string path;
   SPtr<Communicator> comm;
   MPI_Datatype gridParamType, block3dType;
};
#endif // ! _MPIIOCoProcessor_H_
#define _MPIIOCoProcessor_H_
