#include "RefineCrossAndInsideGbObjectBlockVisitor.h"

RefineCrossAndInsideGbObjectBlockVisitor::RefineCrossAndInsideGbObjectBlockVisitor()
:  Block3DVisitor() , notActive(true)
{
}
//////////////////////////////////////////////////////////////////////////
RefineCrossAndInsideGbObjectBlockVisitor::RefineCrossAndInsideGbObjectBlockVisitor(GbObject3DPtr geoObject, int refineLevel)
   : Block3DVisitor(0,refineLevel-1), geoObject(geoObject), notActive(true)
{

}
//////////////////////////////////////////////////////////////////////////
RefineCrossAndInsideGbObjectBlockVisitor::~RefineCrossAndInsideGbObjectBlockVisitor()
{
}
//////////////////////////////////////////////////////////////////////////
void RefineCrossAndInsideGbObjectBlockVisitor::visit(const Grid3DPtr grid, Block3DPtr block)
{
   int level = block->getLevel();
   if( notActive && block->isNotActive() ) return;
   if( level > this->getStopLevel() ) return;

   UbTupleDouble3 coords = grid->getBlockWorldCoordinates(block);
   UbTupleDouble3 deltas = grid->getBlockLengths(block);
   if(geoObject->isCellInsideOrCuttingGbObject3D(  val<1>(coords) 
      , val<2>(coords)
      , val<3>(coords)
      , val<1>(coords)+val<1>(deltas)
      , val<2>(coords)+val<2>(deltas)
      , val<3>(coords)+val<3>(deltas)) ) 
   {
      grid->expandBlock(block->getX1(),block->getX2(),block->getX3(),level); 
   } 

   return;
}
//////////////////////////////////////////////////////////////////////////
