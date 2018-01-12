#ifndef GenBlocksGridVisitor_h
#define GenBlocksGridVisitor_h

#include <PointerDefinitions.h>

#include <basics/utilities/UbTuple.h>

#include "Grid3DVisitor.h" 

class GbObject3D;
class Grid3D;

class GenBlocksGridVisitor : public Grid3DVisitor
{
public:
   GenBlocksGridVisitor(SPtr<GbObject3D> boundingBox);
   virtual ~GenBlocksGridVisitor(){}

   void visit(SPtr<Grid3D> grid);

private:
   UbTupleInt3 minInd, maxInd;
   SPtr<GbObject3D> boundingBox;
   void fillExtentWithBlocks(SPtr<Grid3D> grid);
   void genBlocks(SPtr<Grid3D> grid);
};

#endif 
