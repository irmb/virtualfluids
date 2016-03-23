#include "D3Q27ETForThinWallBCProcessor.h"

#include "D3Q27EsoTwist3DSplittedVector.h"

D3Q27ETForThinWallBCProcessor::D3Q27ETForThinWallBCProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
D3Q27ETForThinWallBCProcessor::D3Q27ETForThinWallBCProcessor(LBMKernel3DPtr kernel) : D3Q27ETBCProcessor(kernel)
{
   //distributionsTemp = EsoTwist3DPtr(new D3Q27EsoTwist3DSplittedVector(maxX1, maxX2, maxX2, -999.0));
}
//////////////////////////////////////////////////////////////////////////
D3Q27ETForThinWallBCProcessor::~D3Q27ETForThinWallBCProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
BCProcessorPtr D3Q27ETForThinWallBCProcessor::clone(LBMKernel3DPtr kernel)
{
   BCProcessorPtr bcProcessor(new D3Q27ETForThinWallBCProcessor(kernel));
   return bcProcessor;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETForThinWallBCProcessor::applyBC()
{
   D3Q27ETBCProcessor::applyBC();

   LBMReal fReturn;

   for (int x3 = minX3; x3 < maxX3; x3++)
   {
      for (int x2 = minX2; x2 < maxX2; x2++)
      {
         for (int x1 = minX1; x1 < maxX1; x1++)
         {
            if (!bcArray.isSolid(x1, x2, x3) && !bcArray.isUndefined(x1, x2, x3))
            {
               if ((bcPtr=bcArray.getBC(x1, x2, x3)) != NULL)
               {
                  if (bcPtr->hasNoSlipBoundary())
                  {
                     for (int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
                     {
                        if (bcPtr->hasNoSlipBoundaryFlag(fdir))
                        {
                           switch (bcPtr->getNoSlipSecondaryOption(fdir))
                           {
                           case 2:
                           {
                              //quadratic bounce back with for thin walls
                              const int invDir = D3Q27System::INVDIR[fdir];
                              fReturn = distributionsTemp->getDistributionInvForDirection(x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
                              distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
                           }
                           break;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETForThinWallBCProcessor::applyPostCollisionBC()
{
   D3Q27ETBCProcessor::applyPostCollisionBC();

   BOOST_FOREACH(BoundaryConditionPtr bc, postBC)
   {
      if (bc->getType() == BoundaryCondition::NoSlip)
      {
         bc->apply();
      }
   }
}


