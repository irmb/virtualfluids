#include "RefineAroundGbObjectHelper.h"
#include "RatioBlockVisitor.h"
#include "RatioSmoothBlockVisitor.h"
#include "OverlapBlockVisitor.h"
#include "SetInterpolationDirsBlockVisitor.h"
#include <D3Q27System.h>
#include <Grid3D.h>
#include <D3Q27TriFaceMeshInteractor.h>
#include "Communicator.h"

RefineAroundGbObjectHelper::RefineAroundGbObjectHelper(Grid3DPtr grid, int refineLevel, D3Q27TriFaceMeshInteractorPtr objectIter, double startDistance, double stopDistance, CommunicatorPtr comm) :
   grid(grid),
   refineLevel(refineLevel),
   objectIter(objectIter),
   startDistance(startDistance), 
   stopDistance(stopDistance),
   comm(comm)
{
}
//////////////////////////////////////////////////////////////////////////
RefineAroundGbObjectHelper::~RefineAroundGbObjectHelper(void)
{
}
//////////////////////////////////////////////////////////////////////////
void RefineAroundGbObjectHelper::refine()
{
   UBLOG(logDEBUG5,"RefineCrossAndInsideGbObjectHelper: refine - start");	

   int rank = grid->getRank();
   grid->setRank(0);

   objectIter->refineBlockGridToLevel(refineLevel, startDistance, stopDistance);

   RatioBlockVisitor ratioVisitor(refineLevel);
   grid->accept(ratioVisitor);

   RatioSmoothBlockVisitor ratioSmoothVisitor(refineLevel);
   grid->accept(ratioSmoothVisitor);

   OverlapBlockVisitor overlapVisitor(refineLevel, false);
   grid->accept(overlapVisitor);

   std::vector<int> dirs;
   for (int i=D3Q27System::E; i<=D3Q27System::TS; i++)
   {
      dirs.push_back(i);
   }
   SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
   grid->accept(interDirsVisitor);

   grid->setRank(rank);

   UBLOG(logDEBUG5,"RefineCrossAndInsideGbObjectHelper: refine - end");	
}

