#include "ThinWallBCProcessor.h"

#include "ThinWallNoSlipBCAlgorithm.h"

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

   BOOST_FOREACH(BCAlgorithmPtr bc, postBC)
   {
      if (bc->getType() == BCAlgorithm::ThinWallNoSlipBCAlgorithm)
      {
         boost::dynamic_pointer_cast<ThinWallNoSlipBCAlgorithm>(bc)->setPass(2); 
         bc->apply();
         boost::dynamic_pointer_cast<ThinWallNoSlipBCAlgorithm>(bc)->setPass(1);
      }
   }
}


