#ifndef GenBlocksGridVisitor_h
#define GenBlocksGridVisitor_h

#include "Grid3DVisitor.h" 
#include <numerics/geometry3d/GbObject3D.h>

class GenBlocksGridVisitor : public Grid3DVisitor
{
public:
   GenBlocksGridVisitor(GbObject3DPtr boundingBox);
   GenBlocksGridVisitor(int nx1, int nx2, int nx3);
   virtual ~GenBlocksGridVisitor(){}

   void visit(Grid3DPtr grid);

private:
   double orgX1, orgX2, orgX3;
   UbTupleInt3 minInd, maxInd;
   int nx1, nx2, nx3;
   bool withDeltaX;
   GbObject3DPtr boundingBox;
   void fillExtentWithBlocks(Grid3DPtr grid);
   void findOrigin(Grid3DPtr grid);
   void genBlocks(Grid3DPtr grid);
};

#endif 
