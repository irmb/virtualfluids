#ifndef InitDistributionsWithCoarseGridBlockVisitor_h__
#define InitDistributionsWithCoarseGridBlockVisitor_h__

#include <memory>

#include "Block3DVisitor.h"
#include "LBMSystem.h"

class Grid3D;
class Block3D;
class InterpolationProcessor;

class InitDistributionsWithInterpolationGridVisitor : public Block3DVisitor
{
public:
   InitDistributionsWithInterpolationGridVisitor(std::shared_ptr<Grid3D> oldGrid, std::shared_ptr<InterpolationProcessor> iProcessor, LBMReal nu);
   ~InitDistributionsWithInterpolationGridVisitor();
   void visit(std::shared_ptr<Grid3D> grid);
private:
   void copyLocalBlock(std::shared_ptr<Block3D> oldBlock, std::shared_ptr<Block3D> newBlock);
   void copyRemoteBlock(std::shared_ptr<Block3D> oldBlock, std::shared_ptr<Block3D> newBlock);
   void interpolateLocalBlockCoarseToFine(std::shared_ptr<Block3D> oldBlock, std::shared_ptr<Block3D> newBlock);
   void interpolateRemoteBlockCoarseToFine(std::shared_ptr<Block3D> oldBlock, std::shared_ptr<Block3D> newBlock);
   void interpolateLocalBlockFineToCoarse(std::shared_ptr<Block3D> oldBlock, std::shared_ptr<Block3D> newBlock);
   void interpolateRemoteBlockFineToCoarse(std::shared_ptr<Block3D> oldBlock, std::shared_ptr<Block3D> newBlock);

   std::shared_ptr<Grid3D> newGrid;
   std::shared_ptr<Grid3D> oldGrid;
   LBMReal nu;

   std::shared_ptr<InterpolationProcessor> iProcessor;
};

#endif // InitDistributionsWithVelocityProfileBlockVisitor_h__
