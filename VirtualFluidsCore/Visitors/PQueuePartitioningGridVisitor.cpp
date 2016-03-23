#include "PQueuePartitioningGridVisitor.h"
#include "PriorityQueueDecompositor.h"
#include <vector>
#include <boost/foreach.hpp>

using namespace std;

PQueuePartitioningGridVisitor::PQueuePartitioningGridVisitor(int numOfParts) : numOfParts(numOfParts)
{

}
//////////////////////////////////////////////////////////////////////////
void PQueuePartitioningGridVisitor::visit(Grid3DPtr grid)
{
   UBLOG(logDEBUG5, "PQueuePartitioningGridVisitor::visit() - start");
   vector<Block3DPtr> blocks;
   vector<int> weights;
   vector< vector <Block3DPtr> > parts;
   int gridRank = grid->getRank();

   int minInitLevel = grid->getCoarsestInitializedLevel();
   int maxInitLevel = grid->getFinestInitializedLevel();

   for(int level = minInitLevel; level<=maxInitLevel; level++)
   {
      vector<Block3DPtr> blockVector;
      grid->getBlocks(level, gridRank, true, blockVector);
      BOOST_FOREACH(Block3DPtr block, blockVector)
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
   BOOST_FOREACH(vector<Block3DPtr> p, parts)
   {
      BOOST_FOREACH(Block3DPtr block, p)
      {
         block->setPart(i);
      }
      i++;
      
   }
   UBLOG(logDEBUG5, "PQueuePartitioningGridVisitor::visit() - end");
}
