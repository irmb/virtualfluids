#ifndef D3Q27SetUndefinedNodesBlockVisitor_h
#define D3Q27SetUndefinedNodesBlockVisitor_h

#include "Block3DVisitor.h"
#include "LBMKernel3D.h"
#include "BCArray3D.h"
#include "D3Q27BoundaryCondition.h"
#include "D3Q27ETBCProcessor.h"

class D3Q27SetUndefinedNodesBlockVisitor : public Block3DVisitor
{
public:
   D3Q27SetUndefinedNodesBlockVisitor();

   virtual ~D3Q27SetUndefinedNodesBlockVisitor() {}

   void visit(Grid3DPtr grid, Block3DPtr block);

private:
   void setNodesUndefined( int startix1, int endix1, int startix2, int endix2, int startix3, int endix3, BCArray3D<D3Q27BoundaryCondition>& bcMatix );
};
#endif
