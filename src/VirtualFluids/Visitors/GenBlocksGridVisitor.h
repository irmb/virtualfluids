#ifndef GenBlocksGridVisitor_h
#define GenBlocksGridVisitor_h

#include <memory>

#include <VirtualFluidsBasics/basics/utilities/UbTuple.h>

#include "Grid3DVisitor.h" 

class GbObject3D;
class Grid3D;

class GenBlocksGridVisitor : public Grid3DVisitor
{
public:
   GenBlocksGridVisitor(std::shared_ptr<GbObject3D> boundingBox);
   virtual ~GenBlocksGridVisitor(){}

   void visit(std::shared_ptr<Grid3D> grid);

private:
   UbTupleInt3 minInd, maxInd;
   std::shared_ptr<GbObject3D> boundingBox;
   void fillExtentWithBlocks(std::shared_ptr<Grid3D> grid);
   void genBlocks(std::shared_ptr<Grid3D> grid);
};

#endif 
