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
void PQueuePartitioningGridVisitor::visit(Grid3DPtr grid)
{
   UBLOG(logDEBUG5, "PQueuePartitioningGridVisitor::visit() - start");
   std::vector<Block3DPtr> blocks;
   std::vector<int> weights;
   std::vector< std::vector <Block3DPtr> > parts;
   int gridRank = grid->getRank();

   int minInitLevel = grid->getCoarsestInitializedLevel();
   int maxInitLevel = grid->getFinestInitializedLevel();

   for(int level = minInitLevel; level<=maxInitLevel; level++)
   {
       std::vector<Block3DPtr> blockVector;
      grid->getBlocks(level, gridRank, true, blockVector);
      for(Block3DPtr block : blockVector)
      {
         if (block)
         {
            blocks.push_back(block);
            weights.push_back(block->getNumberOfLocalConnectors()*(1<<block->getLevel()));
         }
      }
   }
   PriorityQueueDecompositor <Block3DPtr> dec = PriorityQueueDecompositor <Block3DPtr> (blocks, weights, numOfParts);
   dec.getDecomposition(parts);

   int i = 0;
   for(std::vector<Block3DPtr> p : parts)
   {
      for(Block3DPtr block : p)
      {
         block->setPart(i);
      }
      i++;
      
   }
   UBLOG(logDEBUG5, "PQueuePartitioningGridVisitor::visit() - end");
}
