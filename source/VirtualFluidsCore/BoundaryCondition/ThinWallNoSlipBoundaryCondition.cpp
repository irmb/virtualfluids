#include "ThinWallNoSlipBoundaryCondition.h"

#include "D3Q27EsoTwist3DSplittedVector.h"


ThinWallNoSlipBoundaryCondition::ThinWallNoSlipBoundaryCondition()
{
   BoundaryCondition::type = BoundaryCondition::NoSlip;
   BoundaryCondition::preCollision = false;
   pass = 1;
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
   LBMReal f[D3Q27System::ENDF + 1];
   LBMReal feq[D3Q27System::ENDF + 1];
   distributions->getDistributionInv(f, x1, x2, x3);
   LBMReal rho, vx1, vx2, vx3;
   calcMacrosFct(f, rho, vx1, vx2, vx3);
   calcFeqFct(feq, rho, vx1, vx2, vx3);

   LBMReal fReturn;

   for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
   {
      if (bcPtr->hasNoSlipBoundaryFlag(fdir))
      {
         const int invDir = D3Q27System::INVDIR[fdir];
         if (pass == 1)
         {
            LBMReal q = bcPtr->getQ(invDir);
            LBMReal fReturn = ((1.0 - q) / (1.0 + q))*0.5*(f[invDir] - f[fdir] + (f[invDir] + f[fdir] - collFactor*(feq[fdir] + feq[invDir])) / (1.0 - collFactor));
            distributionsTemp->setDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);
         }
         else
         {
            //quadratic bounce back with for thin walls
            fReturn = distributionsTemp->getDistributionInvForDirection(x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);
            distributions->setDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void ThinWallNoSlipBoundaryCondition::addDistributions(DistributionArray3DPtr distributions)
{
   this->distributions = distributions;
   distributionsTemp = EsoTwist3DPtr(new D3Q27EsoTwist3DSplittedVector(distributions->getNX1(), distributions->getNX2(), distributions->getNX3(), -999.0));
}
//////////////////////////////////////////////////////////////////////////
void ThinWallNoSlipBoundaryCondition::setPass(int pass)
{
   this->pass = pass;
}
