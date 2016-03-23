#include "ThinWallNoSlipBoundaryCondition.h"

#include "D3Q27EsoTwist3DSplittedVector.h"


ThinWallNoSlipBoundaryCondition::ThinWallNoSlipBoundaryCondition()
{
   BoundaryCondition::type = BoundaryCondition::NoSlip;
   BoundaryCondition::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
ThinWallNoSlipBoundaryCondition::~ThinWallNoSlipBoundaryCondition()
{

}
//////////////////////////////////////////////////////////////////////////
BoundaryConditionPtr ThinWallNoSlipBoundaryCondition::clone()
{
   BoundaryConditionPtr bc(new ThinWallNoSlipBoundaryCondition());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void ThinWallNoSlipBoundaryCondition::applyBC()
{
   LBMReal fReturn;

   for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
   {
      if (bcPtr->hasNoSlipBoundaryFlag(fdir))
      {
         //quadratic bounce back with for thin walls
         const int invDir = D3Q27System::INVDIR[fdir];
         fReturn = distributionsTemp->getDistributionInvForDirection(x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
         distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void ThinWallNoSlipBoundaryCondition::addDistributions(EsoTwist3DPtr distributions)
{
   this->distributions = distributions;
   distributionsTemp = EsoTwist3DPtr(new D3Q27EsoTwist3DSplittedVector(distributions->getNX1(), distributions->getNX2(), distributions->getNX3(), -999.0));
}

