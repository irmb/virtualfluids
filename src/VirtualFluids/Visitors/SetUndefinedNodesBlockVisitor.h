#ifndef SetUndefinedNodesBlockVisitor_h
#define SetUndefinedNodesBlockVisitor_h

#include "Block3DVisitor.h"
#include "LBMKernel.h"
#include "BCArray3D.h"
#include "BoundaryConditions.h"
#include "BCProcessor.h"

class SetUndefinedNodesBlockVisitor : public Block3DVisitor
{
public:
   SetUndefinedNodesBlockVisitor();

   virtual ~SetUndefinedNodesBlockVisitor() {}

   void visit(Grid3DPtr grid, Block3DPtr block);

private:
   void setNodesUndefined( int startix1, int endix1, int startix2, int endix2, int startix3, int endix3, BCArray3D& bcMatix );
};
#endif
