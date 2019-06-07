#include "RenumberGridVisitor.h"
#include "Grid3DSystem.h"
#include "Grid3D.h"
#include "Block3D.h"
//#include <mpi.h>

RenumberGridVisitor::RenumberGridVisitor(SPtr<Communicator> com)
 : comm(com)
{
}

//////////////////////////////////////////////////////////////////////////
void RenumberGridVisitor::visit(SPtr<Grid3D> grid)
{
   int counter = 0;

   //UBLOG(logDEBUG5, "RenumberGridVisitor::visit() - start");
   std::vector<SPtr<Block3D>> blocks;
   int gridRank = grid->getRank();
   int size;
   //MPI_Comm_size(MPI_COMM_WORLD, &size);
   size = comm->getNumberOfProcesses();

   int minInitLevel = grid->getCoarsestInitializedLevel();
   int maxInitLevel = grid->getFinestInitializedLevel();

   Grid3D::BlockIDMap& blockIdMap = grid->getBlockIDs();
   blockIdMap.clear();

   for(int rank = 0; rank < size; rank++)
   {
      for (int level = minInitLevel; level <= maxInitLevel; level++)
      {
         std::vector<SPtr<Block3D>> blockVector;
         grid->getBlocks(level, blockVector);
         for (SPtr<Block3D> block : blockVector)
         {
            if(block->getRank() == rank)
            { 
               block->setGlobalID(counter);
               blockIdMap.insert(std::make_pair(counter, block));
               counter++;
            }
         }
      }
   }

   //UBLOG(logDEBUG5, "RenumberGridVisitor::visit() - end");
}
