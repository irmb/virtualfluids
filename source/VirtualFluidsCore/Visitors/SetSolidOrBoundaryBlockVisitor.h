#ifndef SET_SOLID_OR_TRANS_BLOCK_VISITOR_H
#define SET_SOLID_OR_TRANS_BLOCK_VISITOR_H

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;
class Interactor3D;

class SetSolidOrBoundaryBlockVisitor : public Block3DVisitor
{
public:
   enum BlockType { SOLID, BC };
public:
   SetSolidOrBoundaryBlockVisitor(SPtr<Interactor3D> interactor, BlockType type);
   virtual ~SetSolidOrBoundaryBlockVisitor() {}

   virtual void visit(SPtr<Grid3D> grid, SPtr<Block3D> block);

private:
   SPtr<Interactor3D> interactor;
   BlockType type;
};

#endif

