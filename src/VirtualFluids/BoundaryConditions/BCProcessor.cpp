#include "BCProcessor.h"
#include "D3Q27EsoTwist3DSplittedVector.h"

BCProcessor::BCProcessor()
{
   
}
//////////////////////////////////////////////////////////////////////////
BCProcessor::BCProcessor(LBMKernelPtr kernel)
{
   DistributionArray3DPtr distributions = boost::dynamic_pointer_cast<EsoTwist3D>(kernel->getDataSet()->getFdistributions());
   bcArray.resize( distributions->getNX1(), distributions->getNX2(), distributions->getNX3(), BCArray3D::FLUID);
}
//////////////////////////////////////////////////////////////////////////
BCProcessor::~BCProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
BCProcessorPtr BCProcessor::clone(LBMKernelPtr kernel)
{
   BCProcessorPtr bcProcessor(new BCProcessor(kernel));
   return bcProcessor;
}
//////////////////////////////////////////////////////////////////////////
BCArray3D& BCProcessor::getBCArray()
{ 
   return this->bcArray; 
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::addBC(BCAlgorithmPtr bc)
{
   if (bc->isPreCollision())
   {
      preBC.push_back(bc);
   }
   else
   {
      postBC.push_back(bc);
   }
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::applyPreCollisionBC()
{
   BOOST_FOREACH(BCAlgorithmPtr bc, preBC)
   {
      bc->apply();
   }
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::applyPostCollisionBC()
{
   BOOST_FOREACH(BCAlgorithmPtr bc, postBC)
   {
      bc->apply();
   }
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::clearBC()
{
   preBC.clear();
   postBC.clear();
}

