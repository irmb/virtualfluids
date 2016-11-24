#include "D3Q27ETBCProcessor.h"
#include "D3Q27EsoTwist3DSplittedVector.h"

D3Q27ETBCProcessor::D3Q27ETBCProcessor()
{
   
}
//////////////////////////////////////////////////////////////////////////
D3Q27ETBCProcessor::D3Q27ETBCProcessor(LBMKernel3DPtr kernel)
{
   DistributionArray3DPtr distributions = boost::dynamic_pointer_cast<EsoTwist3D>(kernel->getDataSet()->getFdistributions());
   bcArray.resize( distributions->getNX1(), distributions->getNX2(), distributions->getNX3(), BCArray3D<D3Q27BoundaryCondition>::FLUID);
}
//////////////////////////////////////////////////////////////////////////
D3Q27ETBCProcessor::~D3Q27ETBCProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
BCProcessorPtr D3Q27ETBCProcessor::clone(LBMKernel3DPtr kernel)
{
   BCProcessorPtr bcProcessor(new D3Q27ETBCProcessor(kernel));
   BoundaryConditionPtr velocity = this->getBC(BoundaryCondition::Velocity);
   BoundaryConditionPtr density  = this->getBC(BoundaryCondition::Density);
   BoundaryConditionPtr noSlip   = this->getBC(BoundaryCondition::NoSlip);
   BoundaryConditionPtr slip     = this->getBC(BoundaryCondition::Slip);
   if (velocity)bcProcessor->addBC(velocity->clone());
   if (density) bcProcessor->addBC(density->clone());
   if (noSlip)  bcProcessor->addBC(noSlip->clone());
   if (slip)    bcProcessor->addBC(slip->clone());
   return bcProcessor;
}
//////////////////////////////////////////////////////////////////////////
BCArray3D<D3Q27BoundaryCondition>& D3Q27ETBCProcessor::getBCArray()
{ 
   return this->bcArray; 
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETBCProcessor::addBC(BoundaryConditionPtr bc)
{
   BoundaryCondition::Type type = bc->getType();
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
BoundaryConditionPtr D3Q27ETBCProcessor::getBC(BoundaryCondition::Type type)
{
   BOOST_FOREACH(BoundaryConditionPtr bc, preBC)
   {
      if (bc->getType() == type)
      {
         return bc;
      }
   }

   BOOST_FOREACH(BoundaryConditionPtr bc, postBC)
   {
      if (bc->getType() == type)
      {
         return bc;
      }
   }

   return BoundaryConditionPtr();
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETBCProcessor::applyPreCollisionBC()
{
   BOOST_FOREACH(BoundaryConditionPtr bc, preBC)
   {
      bc->apply();
   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETBCProcessor::applyPostCollisionBC()
{
   BOOST_FOREACH(BoundaryConditionPtr bc, postBC)
   {
      bc->apply();
   }
}
