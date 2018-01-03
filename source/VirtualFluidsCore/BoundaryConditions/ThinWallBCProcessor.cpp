#include "ThinWallBCProcessor.h"

#include "ThinWallNoSlipBCAlgorithm.h"

#include "LBMKernel.h"

ThinWallBCProcessor::ThinWallBCProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
ThinWallBCProcessor::ThinWallBCProcessor(LBMKernelPtr kernel) : BCProcessor(kernel)
{

}
//////////////////////////////////////////////////////////////////////////
ThinWallBCProcessor::~ThinWallBCProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
BCProcessorPtr ThinWallBCProcessor::clone(LBMKernelPtr kernel)
{
   BCProcessorPtr bcProcessor(new ThinWallBCProcessor(kernel));
   return bcProcessor;
}
//////////////////////////////////////////////////////////////////////////
void ThinWallBCProcessor::applyPostCollisionBC()
{
   BCProcessor::applyPostCollisionBC();

   for(BCAlgorithmPtr bc : postBC)
   {
      if (bc->getType() == BCAlgorithm::ThinWallNoSlipBCAlgorithm)
      {
         std::dynamic_pointer_cast<ThinWallNoSlipBCAlgorithm>(bc)->setPass(2); 
         bc->applyBC();
         std::dynamic_pointer_cast<ThinWallNoSlipBCAlgorithm>(bc)->setPass(1);
      }
   }
}


