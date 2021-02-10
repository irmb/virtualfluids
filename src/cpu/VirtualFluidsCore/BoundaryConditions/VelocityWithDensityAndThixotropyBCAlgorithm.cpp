#include "VelocityWithDensityAndThixotropyBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BCArray3D.h"

VelocityWithDensityAndThixotropyBCAlgorithm::VelocityWithDensityAndThixotropyBCAlgorithm()
{
   BCAlgorithm::type = BCAlgorithm::VelocityWithDensityAndThixotropyBCAlgorithm;
   BCAlgorithm::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
VelocityWithDensityAndThixotropyBCAlgorithm::~VelocityWithDensityAndThixotropyBCAlgorithm()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> VelocityWithDensityAndThixotropyBCAlgorithm::clone()
{
   SPtr<BCAlgorithm> bc(new VelocityWithDensityAndThixotropyBCAlgorithm());
   dynamicPointerCast<VelocityWithDensityAndThixotropyBCAlgorithm>(bc)->setLambdaBC(lambdaBC);
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void VelocityWithDensityAndThixotropyBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void VelocityWithDensityAndThixotropyBCAlgorithm::addDistributionsH(SPtr<DistributionArray3D> distributions)
{
   this->distributionsH = distributions;
}
//////////////////////////////////////////////////////////////////////////
void VelocityWithDensityAndThixotropyBCAlgorithm::applyBC()
{
   //velocity bc for non reflecting pressure bc
   LBMReal f[D3Q27System::ENDF+1];
   distributions->getDistributionInv(f, x1, x2, x3);
   
   LBMReal h[D3Q27System::ENDF + 1];
   distributionsH->getDistributionInv(h, x1, x2, x3);

   LBMReal rho, vx1, vx2, vx3, drho;
   calcMacrosFct(f, drho, vx1, vx2, vx3);
   
   rho = 1.0+drho*compressibleFactor;
  
   ///////////////////////////////////////////////////////////////////
   // Thixotropy
   LBMReal lambda = D3Q27System::getDensity(h);

   int nx1 = x1;
   int nx2 = x2;
   int nx3 = x3;
   int direction = -1;

   //flag points in direction of fluid
   if (bcPtr->hasVelocityBoundaryFlag(D3Q27System::E)) { nx1 -= 1; direction = D3Q27System::E; }
   else if (bcPtr->hasVelocityBoundaryFlag(D3Q27System::W)) { nx1 += 1; direction = D3Q27System::W; }
   else if (bcPtr->hasVelocityBoundaryFlag(D3Q27System::N)) { nx2 -= 1; direction = D3Q27System::N; }
   else if (bcPtr->hasVelocityBoundaryFlag(D3Q27System::S)) { nx2 += 1; direction = D3Q27System::S; }
   else if (bcPtr->hasVelocityBoundaryFlag(D3Q27System::T)) { nx3 -= 1; direction = D3Q27System::T; }
   else if (bcPtr->hasVelocityBoundaryFlag(D3Q27System::B)) { nx3 += 1; direction = D3Q27System::B; }
   else	 UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on velocity boundary..."));

   for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++)
   {
      int nX1 = x1 + D3Q27System::DX1[fdir];
      int nX2 = x2 + D3Q27System::DX2[fdir];
      int nX3 = x3 + D3Q27System::DX3[fdir];

      int minX1 = 0;
      int minX2 = 0;
      int minX3 = 0;

      int maxX1 = (int)bcArray->getNX1();
      int maxX2 = (int)bcArray->getNX2();
      int maxX3 = (int)bcArray->getNX3();

      if (minX1 <= nX1 && maxX1 > nX1 && minX2 <= nX2 && maxX2 > nX2 && minX3 <= nX3 && maxX3 > nX3)
      {
         if (bcArray->isSolid(nX1,nX2,nX3))
         {
            const int invDir = D3Q27System::INVDIR[fdir];
            LBMReal velocity = bcPtr->getBoundaryVelocity(fdir);

            LBMReal fReturn = (f[fdir] + f[invDir] - velocity*rho) / 2.0 - drho*D3Q27System::WEIGTH[invDir];
            distributions->setDistributionForDirection(fReturn, nX1, nX2, nX3, invDir);
         }
      }
      
      if (bcPtr->hasVelocityBoundaryFlag(fdir))
      {
         LBMReal htemp = D3Q27System::getCompFeqForDirection(fdir, lambda, vx1, vx2, vx3);
         htemp = D3Q27System::getCompFeqForDirection(fdir, lambdaBC, vx1, vx2, vx3) + h[fdir] - htemp;
         distributionsH->setDistributionForDirection(htemp, nx1, nx2, nx3, fdir);
      }
   }
}
