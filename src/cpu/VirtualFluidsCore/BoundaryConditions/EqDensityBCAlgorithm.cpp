#include "EqDensityBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"


EqDensityBCAlgorithm::EqDensityBCAlgorithm()
{
   BCAlgorithm::type = BCAlgorithm::EqDensityBCAlgorithm;
   BCAlgorithm::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
EqDensityBCAlgorithm::~EqDensityBCAlgorithm()
= default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> EqDensityBCAlgorithm::clone()
{
   SPtr<BCAlgorithm> bc(new EqDensityBCAlgorithm());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void EqDensityBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void EqDensityBCAlgorithm::applyBC()
{
   LBMReal f[D3Q27System::ENDF+1];

   distributions->getDistributionInv(f, x1, x2, x3);
   int nx1 = x1;
   int nx2 = x2;
   int nx3 = x3;

   //flag points in direction of fluid
   if      (bcPtr->hasDensityBoundaryFlag(D3Q27System::E)) { nx1 -= 1;}
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::W)) { nx1 += 1;}
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::N)) { nx2 -= 1;}
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::S)) { nx2 += 1;}
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::T)) { nx3 -= 1;}
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::B)) { nx3 += 1;}
   else UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on density boundary..."));

   LBMReal rho, vx1, vx2, vx3;
   calcMacrosFct(f, rho, vx1, vx2, vx3);
   LBMReal rhoBC = bcPtr->getBoundaryDensity();
   for (int fdir = D3Q27System::STARTF; fdir<=D3Q27System::ENDF; fdir++)
   {
      if (bcPtr->hasDensityBoundaryFlag(fdir))
      {
         //Ehsan: 15.2.2013:
         LBMReal ftemp = calcFeqsForDirFct(fdir, rhoBC, vx1, vx2, vx3);
         distributions->setDistributionForDirection(ftemp, nx1, nx2, nx3, fdir);
      }
   }

}

