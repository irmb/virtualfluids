#ifndef InitDistributionsWithCoarseGridBlockVisitor_h__
#define InitDistributionsWithCoarseGridBlockVisitor_h__

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"
#include "LBMSystem.h"

class Grid3D;
class Block3D;
class InterpolationProcessor;

class InitDistributionsWithInterpolationGridVisitor : public Block3DVisitor
{
public:
   InitDistributionsWithInterpolationGridVisitor(SPtr<Grid3D> oldGrid, SPtr<InterpolationProcessor> iProcessor, LBMReal nu);
   ~InitDistributionsWithInterpolationGridVisitor();
   void visit(SPtr<Grid3D> grid);
private:
   void copyLocalBlock(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock);
   void copyRemoteBlock(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock);
   void interpolateLocalBlockCoarseToFine(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock);
   void interpolateRemoteBlockCoarseToFine(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock);
   void interpolateLocalBlockFineToCoarse(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock);
   void interpolateRemoteBlockFineToCoarse(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock);

   SPtr<Grid3D> newGrid;
   SPtr<Grid3D> oldGrid;
   LBMReal nu;

   SPtr<InterpolationProcessor> iProcessor;
};

#endif // InitDistributionsWithVelocityProfileBlockVisitor_h__
