#ifndef SET_SOLID_OR_TRANS_BLOCK_VISITOR_H
#define SET_SOLID_OR_TRANS_BLOCK_VISITOR_H

#include <memory>

#include "Block3DVisitor.h"

#include <VirtualFluidsDefinitions.h>

class Grid3D;
class Block3D;
class Interactor3D;

enum class BlockType { SOLID, BC };

class VF_PUBLIC SetSolidBlockVisitor : public Block3DVisitor
{
public:
   SetSolidBlockVisitor(std::shared_ptr<Interactor3D> interactor, BlockType type);
   virtual ~SetSolidBlockVisitor() {}

   virtual void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block);

private:
    std::shared_ptr<Interactor3D> interactor;
   BlockType type;
};

#endif

