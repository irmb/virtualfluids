#ifndef InitDistributionsWithCoarseGridBlockVisitor_h__
#define InitDistributionsWithCoarseGridBlockVisitor_h__

#include "Block3DVisitor.h"
#include <mpi.h>

class InitDistributionsWithInterpolationGridVisitor : public Grid3DVisitor
{
public:
   InitDistributionsWithInterpolationGridVisitor(Grid3DPtr oldGrid, InterpolationProcessorPtr iProcessor, LBMReal nu);
   ~InitDistributionsWithInterpolationGridVisitor();
   void visit(Grid3DPtr grid);
private:
   void copyLocalBlock(Block3DPtr oldBlock, Block3DPtr newBlock);
   void copyRemoteBlock(Block3DPtr oldBlock, Block3DPtr newBlock);
   void interpolateLocalBlockCoarseToFine(Block3DPtr oldBlock, Block3DPtr newBlock);
   void interpolateRemoteBlockCoarseToFine(Block3DPtr oldBlock, Block3DPtr newBlock);
   void interpolateLocalBlockFineToCoarse(Block3DPtr oldBlock, Block3DPtr newBlock);
   void interpolateRemoteBlockFineToCoarse(Block3DPtr oldBlock, Block3DPtr newBlock);

   Grid3DPtr newGrid;
   Grid3DPtr oldGrid;
   LBMReal nu;

   InterpolationProcessorPtr iProcessor;
};

#endif // InitDistributionsWithVelocityProfileBlockVisitor_h__
