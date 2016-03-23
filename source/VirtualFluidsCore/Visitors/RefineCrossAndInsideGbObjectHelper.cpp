#include "RefineCrossAndInsideGbObjectHelper.h"
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"
#include "RatioBlockVisitor.h"
#include "RatioSmoothBlockVisitor.h"
#include "OverlapBlockVisitor.h"
#include "SetInterpolationDirsBlockVisitor.h"
#include <D3Q27System.h>


RefineCrossAndInsideGbObjectHelper::RefineCrossAndInsideGbObjectHelper(Grid3DPtr grid, int maxRefineLevel) :
                                    grid(grid),
                                    maxRefineLevel(maxRefineLevel)
{
}
//////////////////////////////////////////////////////////////////////////
RefineCrossAndInsideGbObjectHelper::~RefineCrossAndInsideGbObjectHelper(void)
{
}
//////////////////////////////////////////////////////////////////////////
void RefineCrossAndInsideGbObjectHelper::refine()
{
   UBLOG(logDEBUG5,"RefineCrossAndInsideGbObjectHelper: refine - start");	
   
   int size = (int)objects.size();

   for (int i = 0; i < size; i++)
   {
      RefineCrossAndInsideGbObjectBlockVisitor refVisitor(objects[i], levels[i]);
      grid->accept(refVisitor);
   }

   RatioBlockVisitor ratioVisitor(maxRefineLevel);
   grid->accept(ratioVisitor);

   RatioSmoothBlockVisitor ratioSmoothVisitor(maxRefineLevel);
   grid->accept(ratioSmoothVisitor);

   OverlapBlockVisitor overlapVisitor(maxRefineLevel);
   grid->accept(overlapVisitor);

   std::vector<int> dirs;

   for (int i=D3Q27System::E; i<D3Q27System::ENDDIR; i++)
   {
      dirs.push_back(i);
   }
   SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
   grid->accept(interDirsVisitor);
   UBLOG(logDEBUG5,"RefineCrossAndInsideGbObjectHelper: refine - end");	
}
//////////////////////////////////////////////////////////////////////////
void RefineCrossAndInsideGbObjectHelper::addGbObject( GbObject3DPtr object, int refineLevel )
{
   objects.push_back(object);
   levels.push_back(refineLevel);
}
