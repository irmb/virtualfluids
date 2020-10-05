#include "ThinWallBCProcessor.h"

#include "ThinWallNoSlipBCAlgorithm.h"

#include "LBMKernel.h"

//////////////////////////////////////////////////////////////////////////
ThinWallBCProcessor::ThinWallBCProcessor(SPtr<ILBMKernel> kernel) : BCProcessor(kernel)
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCProcessor> ThinWallBCProcessor::clone(SPtr<ILBMKernel> kernel)
{
   SPtr<BCProcessor> bcProcessor(new ThinWallBCProcessor(kernel));
   return bcProcessor;
}
//////////////////////////////////////////////////////////////////////////
void ThinWallBCProcessor::applyPostCollisionBC()
{
   BCProcessor::applyPostCollisionBC();

   for(SPtr<BCAlgorithm> bc : postBC)
   {
      if (bc->getType() == BCAlgorithm::ThinWallNoSlipBCAlgorithm)
      {
         dynamicPointerCast<ThinWallNoSlipBCAlgorithm>(bc)->setPass(2); 
         bc->applyBC();
         dynamicPointerCast<ThinWallNoSlipBCAlgorithm>(bc)->setPass(1);
      }
   }
}


