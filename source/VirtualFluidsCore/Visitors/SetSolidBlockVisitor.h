#ifndef SET_SOLID_OR_TRANS_BLOCK_VISITOR_H
#define SET_SOLID_OR_TRANS_BLOCK_VISITOR_H

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;
class Interactor3D;

enum class BlockType { SOLID, BC };

class SetSolidBlockVisitor : public Block3DVisitor
{
public:
   SetSolidBlockVisitor(SPtr<Interactor3D> interactor, BlockType type);
   virtual ~SetSolidBlockVisitor() {}

   virtual void visit(SPtr<Grid3D> grid, SPtr<Block3D> block);

private:
    SPtr<Interactor3D> interactor;
   BlockType type;
};

#endif

