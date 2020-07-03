#include "PQueuePartitioningGridVisitor.h"

#include <vector>

#include "PriorityQueueDecompositor.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "UbLogger.h"


PQueuePartitioningGridVisitor::PQueuePartitioningGridVisitor(int numOfParts) : numOfParts(numOfParts)
{

}
//////////////////////////////////////////////////////////////////////////
void PQueuePartitioningGridVisitor::visit(SPtr<Grid3D> grid)
{
   UBLOG(logDEBUG5, "PQueuePartitioningGridVisitor::visit() - start");
   std::vector<SPtr<Block3D>> blocks;
   std::vector<int> weights;
   std::vector< std::vector <SPtr<Block3D>> > parts;
   int gridRank = grid->getRank();

   int minInitLevel = grid->getCoarsestInitializedLevel();
   int maxInitLevel = grid->getFinestInitializedLevel();

   for(int level = minInitLevel; level<=maxInitLevel; level++)
   {
       std::vector<SPtr<Block3D>> blockVector;
      grid->getBlocks(level, gridRank, true, blockVector);
      for(SPtr<Block3D> block : blockVector)
      {
         if (block)
         {
            blocks.push_back(block);
            weights.push_back(block->getNumberOfLocalConnectors()*(1<<block->getLevel()));
         }
      }
   }
   PriorityQueueDecompositor <SPtr<Block3D>> dec = PriorityQueueDecompositor <SPtr<Block3D>> (blocks, weights, numOfParts);
   dec.getDecomposition(parts);

   int i = 0;
   for(std::vector<SPtr<Block3D>> p : parts)
   {
      for(SPtr<Block3D> block : p)
      {
         block->setPart(i);
      }
      i++;
      
   }
   UBLOG(logDEBUG5, "PQueuePartitioningGridVisitor::visit() - end");
}
