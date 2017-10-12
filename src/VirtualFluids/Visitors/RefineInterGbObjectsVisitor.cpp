#include "RefineInterGbObjectsVisitor.h"

RefineInterGbObjectsBlockVisitor::RefineInterGbObjectsBlockVisitor() 
   : Block3DVisitor(-1, -1)
{
}
//////////////////////////////////////////////////////////////////////////
RefineInterGbObjectsBlockVisitor::RefineInterGbObjectsBlockVisitor(GbObject3DPtr includeGbObject3D, GbObject3DPtr excludeGbObject3D, int startlevel, int stoplevel)
   : Block3DVisitor(startlevel, stoplevel)
{
   this->includeGbObjects3D.push_back(includeGbObject3D);
   this->excludeGbObjects3D.push_back(excludeGbObject3D);
}
//////////////////////////////////////////////////////////////////////////
RefineInterGbObjectsBlockVisitor::RefineInterGbObjectsBlockVisitor(std::vector<GbObject3DPtr> includeGbObjects3D, std::vector<GbObject3DPtr> excludeGbObjects3D, int startlevel, int stoplevel)
   : Block3DVisitor(startlevel, stoplevel)
{
   this->includeGbObjects3D = includeGbObjects3D;
   this->excludeGbObjects3D = excludeGbObjects3D;
}
//////////////////////////////////////////////////////////////////////////
void RefineInterGbObjectsBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   UbTupleDouble3 coords = grid->getBlockWorldCoordinates(block);
   UbTupleDouble3 delta  = grid->getBlockLengths(block);

   double cellMinX1 = val<1>(coords);
   double cellMinX2 = val<2>(coords);
   double cellMinX3 = val<3>(coords);
   double cellMaxX1 = val<1>(coords)+val<1>(delta);
   double cellMaxX2 = val<2>(coords)+val<2>(delta);
   double cellMaxX3 = val<3>(coords)+val<3>(delta);

   bool insideInclude = false;
   for(size_t i=0; i<includeGbObjects3D.size(); i++)
   {
      if(   includeGbObjects3D[i]->isCellInsideOrCuttingGbObject3D(cellMinX1,cellMinX2,cellMinX3,cellMaxX1,cellMaxX2,cellMaxX3) )
      {
         insideInclude = true;
         break;
      }
   }

   bool insideExclude = false;
   for(size_t e=0; e<excludeGbObjects3D.size(); e++)
   {
      if(excludeGbObjects3D[e]->isCellInsideGbObject3D(cellMinX1, cellMinX2, cellMinX3, cellMaxX1, cellMaxX2, cellMaxX3)) 
      {
         insideExclude = true;
         break;
      }
   }

   if(insideInclude && !insideExclude)         
   {
      int ix1, ix2, ix3, level;
      ix1 = block->getX1();
      ix2 = block->getX2();
      ix3 = block->getX3();
      level = block->getLevel();
      grid->expandBlock(ix1,ix2,ix3,level); 
   }
}
