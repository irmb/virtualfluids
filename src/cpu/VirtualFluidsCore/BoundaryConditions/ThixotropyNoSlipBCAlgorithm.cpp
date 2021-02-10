#include "ThixotropyNoSlipBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

//////////////////////////////////////////////////////////////////////////
void ThixotropyNoSlipBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyNoSlipBCAlgorithm::applyBC()
{
   LBMReal f[D3Q27System::ENDF + 1];
   LBMReal feq[D3Q27System::ENDF + 1];
   distributions->getDistribution(f, x1, x2, x3);
   LBMReal rho, vx1, vx2, vx3;
   calcMacrosFct(f, rho, vx1, vx2, vx3);
   calcFeqFct(feq, rho, vx1, vx2, vx3);

   LBMReal shearRate = D3Q27System::getShearRate(f, collFactor);
   LBMReal collFactorF = getThyxotropyCollFactor(collFactor, shearRate, rho);

   for (int fDir = D3Q27System::FSTARTDIR; fDir <= D3Q27System::FENDDIR; fDir++)
   {
      if (bcPtr->hasNoSlipBoundaryFlag(fDir))
      {
         //quadratic bounce back
         const int invDir = D3Q27System::INVDIR[fDir];
         LBMReal q = bcPtr->getQ(invDir);
         LBMReal fReturn =(f[invDir] + q * f[fDir] + q * collFactorF * (feq[invDir] - f[invDir] + feq[fDir] - f[fDir])) / (1.0 + q);
         distributions->setDistributionInvForDirection(fReturn, x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], invDir);
      }
   }
}