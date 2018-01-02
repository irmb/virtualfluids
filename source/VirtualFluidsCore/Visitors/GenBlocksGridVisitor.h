#ifndef GenBlocksGridVisitor_h
#define GenBlocksGridVisitor_h

#include "Grid3DVisitor.h" 
#include <numerics/geometry3d/GbObject3D.h>

class GenBlocksGridVisitor : public Grid3DVisitor
{
public:
   GenBlocksGridVisitor(GbObject3DPtr boundingBox);
   virtual ~GenBlocksGridVisitor(){}

   void visit(Grid3DPtr grid);

private:
   UbTupleInt3 minInd, maxInd;
   GbObject3DPtr boundingBox;
   void fillExtentWithBlocks(Grid3DPtr grid);
   void genBlocks(Grid3DPtr grid);
};

#endif 
