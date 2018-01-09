#ifndef SetUndefinedNodesBlockVisitor_h
#define SetUndefinedNodesBlockVisitor_h

#include <memory>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;
class BCArray3D;

class SetUndefinedNodesBlockVisitor : public Block3DVisitor
{
public:
   SetUndefinedNodesBlockVisitor();

   virtual ~SetUndefinedNodesBlockVisitor() {}

   void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;

private:
   void setNodesUndefined( int startix1, int endix1, int startix2, int endix2, int startix3, int endix3, std::shared_ptr<BCArray3D> bcMatix );
};
#endif
