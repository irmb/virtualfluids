#ifndef SetUndefinedNodesBlockVisitor_h
#define SetUndefinedNodesBlockVisitor_h

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;
class BCArray3D;

class SetUndefinedNodesBlockVisitor : public Block3DVisitor
{
public:
   SetUndefinedNodesBlockVisitor(bool twoTypeOfConnectorsCheck = true);

   virtual ~SetUndefinedNodesBlockVisitor() {}

   void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;
protected:
   void setNodesUndefined( int startix1, int endix1, int startix2, int endix2, int startix3, int endix3, SPtr<BCArray3D> bcMatix );
private:
   bool twoTypeOfConnectorsCheck;

};
#endif
