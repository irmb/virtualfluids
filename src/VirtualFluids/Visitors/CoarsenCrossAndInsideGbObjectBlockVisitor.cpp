#include "CoarsenCrossAndInsideGbObjectBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include <numerics/geometry3d/GbObject3D.h>

CoarsenCrossAndInsideGbObjectBlockVisitor::CoarsenCrossAndInsideGbObjectBlockVisitor()
   : Block3DVisitor(), notActive(true)
{
}
//////////////////////////////////////////////////////////////////////////
CoarsenCrossAndInsideGbObjectBlockVisitor::CoarsenCrossAndInsideGbObjectBlockVisitor(GbObject3DPtr geoObject, int fineLevel, int coarseLevel)
   : Block3DVisitor(fineLevel, fineLevel), geoObject(geoObject), notActive(true), coarseLevel(coarseLevel)
{

}
//////////////////////////////////////////////////////////////////////////
CoarsenCrossAndInsideGbObjectBlockVisitor::~CoarsenCrossAndInsideGbObjectBlockVisitor()
{
}
//////////////////////////////////////////////////////////////////////////
void CoarsenCrossAndInsideGbObjectBlockVisitor::visit(const Grid3DPtr grid, Block3DPtr block)
{
   int fineLevel = block->getLevel();
   if (notActive && block->isNotActive()) return;
   if (fineLevel>this->getStopLevel()) return;

   UbTupleDouble3 coords = grid->getBlockWorldCoordinates(block);
   UbTupleDouble3 deltas = grid->getBlockLengths(block);
   if (geoObject->isCellInsideOrCuttingGbObject3D(val<1>(coords)
      , val<2>(coords)
      , val<3>(coords)
      , val<1>(coords)+val<1>(deltas)
      , val<2>(coords)+val<2>(deltas)
      , val<3>(coords)+val<3>(deltas)))
   {
      grid->collapseBlock(block->getX1(), block->getX2(), block->getX3(), fineLevel, coarseLevel);
   }

   return;
}
//////////////////////////////////////////////////////////////////////////
