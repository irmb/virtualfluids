#include "NoSlipBoundaryCondition.h"

NoSlipBoundaryCondition::NoSlipBoundaryCondition()
{
   BoundaryCondition::type = BoundaryCondition::NoSlip;
   BoundaryCondition::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
NoSlipBoundaryCondition::~NoSlipBoundaryCondition()
{

}
//////////////////////////////////////////////////////////////////////////
BoundaryConditionPtr NoSlipBoundaryCondition::clone()
{
   BoundaryConditionPtr bc(new NoSlipBoundaryCondition());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void NoSlipBoundaryCondition::applyBC()
{
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal feq[D3Q27System::ENDF+1];
   distributions->getDistributionInv(f, x1, x2, x3);
   LBMReal rho, vx1, vx2, vx3;
   calcMacrosFct(f, rho, vx1, vx2, vx3);
   calcFeqFct(feq, rho, vx1, vx2, vx3);

   for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
   {
      if (bcPtr->hasNoSlipBoundaryFlag(fdir))
      {
         //quadratic bounce back
         const int invDir = D3Q27System::INVDIR[fdir];
         LBMReal q = bcPtr->getQ(invDir);
         LBMReal fReturn = ((1.0-q)/(1.0+q))*((f[invDir]-feq[invDir])/(1.0-collFactor)+feq[invDir])+((q/(1.0+q))*(f[invDir]+f[fdir]));
         distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
      }
   }
}
