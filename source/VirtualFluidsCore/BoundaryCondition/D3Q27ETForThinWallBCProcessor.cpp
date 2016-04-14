#include "D3Q27ETForThinWallBCProcessor.h"

#include "ThinWallNoSlipBoundaryCondition.h"

D3Q27ETForThinWallBCProcessor::D3Q27ETForThinWallBCProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
D3Q27ETForThinWallBCProcessor::D3Q27ETForThinWallBCProcessor(LBMKernel3DPtr kernel) : D3Q27ETBCProcessor(kernel)
{

}
//////////////////////////////////////////////////////////////////////////
D3Q27ETForThinWallBCProcessor::~D3Q27ETForThinWallBCProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
BCProcessorPtr D3Q27ETForThinWallBCProcessor::clone(LBMKernel3DPtr kernel)
{
   BCProcessorPtr bcProcessor(new D3Q27ETForThinWallBCProcessor(kernel));
   BoundaryConditionPtr velocity = this->getBC(BoundaryCondition::Velocity);
   BoundaryConditionPtr density = this->getBC(BoundaryCondition::Density);
   BoundaryConditionPtr noSlip = this->getBC(BoundaryCondition::NoSlip);
   BoundaryConditionPtr slip = this->getBC(BoundaryCondition::Slip);
   if (velocity)bcProcessor->addBC(velocity->clone());
   if (density) bcProcessor->addBC(density->clone());
   if (noSlip)  bcProcessor->addBC(noSlip->clone());
   if (slip)    bcProcessor->addBC(slip->clone());
   return bcProcessor;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETForThinWallBCProcessor::applyPostCollisionBC()
{
   D3Q27ETBCProcessor::applyPostCollisionBC();

   BOOST_FOREACH(BoundaryConditionPtr bc, postBC)
   {
      if (bc->getType() == BoundaryCondition::NoSlip)
      {
         boost::dynamic_pointer_cast<ThinWallNoSlipBoundaryCondition>(bc)->setPass(2); 
         bc->apply();
         boost::dynamic_pointer_cast<ThinWallNoSlipBoundaryCondition>(bc)->setPass(1);
      }
   }
}


