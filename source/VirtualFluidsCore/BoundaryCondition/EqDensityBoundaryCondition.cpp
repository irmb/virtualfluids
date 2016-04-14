#include "EqDensityBoundaryCondition.h"
#include <boost/pointer_cast.hpp>

EqDensityBoundaryCondition::EqDensityBoundaryCondition()
{
   BoundaryCondition::type = BoundaryCondition::Density;
   BoundaryCondition::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
EqDensityBoundaryCondition::~EqDensityBoundaryCondition()
{
}
//////////////////////////////////////////////////////////////////////////
BoundaryConditionPtr EqDensityBoundaryCondition::clone()
{
   BoundaryConditionPtr bc(new EqDensityBoundaryCondition());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void EqDensityBoundaryCondition::addDistributions(DistributionArray3DPtr distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void EqDensityBoundaryCondition::applyBC()
{
   LBMReal f[D3Q27System::ENDF+1];

   distributions->getDistributionInv(f, x1, x2, x3);
   int nx1 = x1;
   int nx2 = x2;
   int nx3 = x3;
   int direction = -1;

   //flag points in direction of fluid
   if (bcPtr->hasDensityBoundaryFlag(D3Q27System::E)) { nx1 -= 1; direction = D3Q27System::E; }
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::W)) { nx1 += 1; direction = D3Q27System::W; }
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::N)) { nx2 -= 1; direction = D3Q27System::N; }
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::S)) { nx2 += 1; direction = D3Q27System::S; }
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::T)) { nx3 -= 1; direction = D3Q27System::T; }
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::B)) { nx3 += 1; direction = D3Q27System::B; }
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

