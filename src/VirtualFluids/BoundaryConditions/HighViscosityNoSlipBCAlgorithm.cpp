#include "HighViscosityNoSlipBCAlgorithm.h"

HighViscosityNoSlipBCAlgorithm::HighViscosityNoSlipBCAlgorithm()
{
   BCAlgorithm::type = BCAlgorithm::HighViscosityNoSlipBCAlgorithm;
   BCAlgorithm::preCollision = true;
}
//////////////////////////////////////////////////////////////////////////
HighViscosityNoSlipBCAlgorithm::~HighViscosityNoSlipBCAlgorithm()
{
}
//////////////////////////////////////////////////////////////////////////
BCAlgorithmPtr HighViscosityNoSlipBCAlgorithm::clone()
{
   BCAlgorithmPtr bc(new HighViscosityNoSlipBCAlgorithm());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void HighViscosityNoSlipBCAlgorithm::addDistributions(DistributionArray3DPtr distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void HighViscosityNoSlipBCAlgorithm::applyBC()
{
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal feq[D3Q27System::ENDF+1];
   distributions->getDistribution(f, x1, x2, x3);
   LBMReal rho, vx1, vx2, vx3;
   calcMacrosFct(f, rho, vx1, vx2, vx3);
   calcFeqFct(feq, rho, vx1, vx2, vx3);

   for (int fDir = D3Q27System::FSTARTDIR; fDir<=D3Q27System::FENDDIR; fDir++)
   {
      if (bcPtr->hasNoSlipBoundaryFlag(fDir))
      {
         //quadratic bounce back
         const int invDir = D3Q27System::INVDIR[fDir];
         LBMReal q = bcPtr->getQ(invDir);
         LBMReal fReturn = (f[invDir]+q*f[fDir]+q*collFactor*(feq[invDir]-f[invDir]+feq[fDir]-f[fDir]))/(1.0+q);
         distributions->setDistributionInvForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], invDir);
      }
   }
}

