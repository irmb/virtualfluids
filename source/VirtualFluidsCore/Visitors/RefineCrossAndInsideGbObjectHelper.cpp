#include "RefineCrossAndInsideGbObjectHelper.h"
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"
#include "RatioBlockVisitor.h"
#include "RatioSmoothBlockVisitor.h"
#include "CheckRatioBlockVisitor.h"
#include "OverlapBlockVisitor.h"
#include "SetInterpolationDirsBlockVisitor.h"
#include <D3Q27System.h>


RefineCrossAndInsideGbObjectHelper::RefineCrossAndInsideGbObjectHelper(Grid3DPtr grid, int maxRefineLevel, CommunicatorPtr comm) :
                                    grid(grid),
                                    maxRefineLevel(maxRefineLevel),
                                    comm(comm)
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
   
   if (comm->isRoot())
   {
      int size = (int)objects.size();

      for (int i = 0; i<size; i++)
      {
         RefineCrossAndInsideGbObjectBlockVisitor refVisitor(objects[i], levels[i]);
         grid->accept(refVisitor);
      }

      //RatioBlockVisitor ratioVisitor(maxRefineLevel);
      //grid->accept(ratioVisitor);

      //RatioSmoothBlockVisitor ratioSmoothVisitor(maxRefineLevel);
      //grid->accept(ratioSmoothVisitor);

      RatioBlockVisitor ratioVisitor(maxRefineLevel);
      CheckRatioBlockVisitor checkRatio(maxRefineLevel);
      int count = 0;

      do {
         grid->accept(ratioVisitor);
         checkRatio.resetState();
         grid->accept(checkRatio);
         UBLOG(logINFO, "count = "<<count++<<" state = "<<checkRatio.getState());
      } while (!checkRatio.getState());


      OverlapBlockVisitor overlapVisitor(maxRefineLevel, false);
      grid->accept(overlapVisitor);
   }

   grid->updateDistributedBlocks(comm);

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
