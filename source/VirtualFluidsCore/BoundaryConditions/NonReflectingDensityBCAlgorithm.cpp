#include "NonReflectingDensityBCAlgorithm.h"
#include <boost/pointer_cast.hpp>
#include "D3Q27System.h"

NonReflectingDensityBCAlgorithm::NonReflectingDensityBCAlgorithm()
{
   BCAlgorithm::type = BCAlgorithm::NonReflectingDensityBCAlgorithm;
   BCAlgorithm::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
NonReflectingDensityBCAlgorithm::~NonReflectingDensityBCAlgorithm()
{
}
//////////////////////////////////////////////////////////////////////////
BCAlgorithmPtr NonReflectingDensityBCAlgorithm::clone()
{
   BCAlgorithmPtr bc(new NonReflectingDensityBCAlgorithm());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void NonReflectingDensityBCAlgorithm::addDistributions(DistributionArray3DPtr distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void NonReflectingDensityBCAlgorithm::applyBC()
{
   distributions->swap();
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal ftemp[D3Q27System::ENDF+1];
   distributions->getDistribution(f, x1, x2, x3);
   int nx1 = x1;
   int nx2 = x2;
   int nx3 = x3;
   int direction = -1;

   //flag points in direction of fluid
   if      (bcPtr->hasDensityBoundaryFlag(D3Q27System::E)) { nx1 += 1; direction = D3Q27System::E; }
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::W)) { nx1 -= 1; direction = D3Q27System::W; }
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::N)) { nx2 += 1; direction = D3Q27System::N; }
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::S)) { nx2 -= 1; direction = D3Q27System::S; }
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::T)) { nx3 += 1; direction = D3Q27System::T; }
   else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::B)) { nx3 -= 1; direction = D3Q27System::B; }
   else UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on density boundary..."));

   //#ifdef _DEBUG
   //   if (nx1<0 || nx1>maxX1) UB_THROW(UbException(UB_EXARGS, "nx1<0 || nx1>=lengthX1"));
   //   if (nx2<0 || nx2>maxX2) UB_THROW(UbException(UB_EXARGS, "nx2<0 || nx2>=lengthX2"));
   //   if (nx3<0 || nx3>maxX3) UB_THROW(UbException(UB_EXARGS, "nx3<0 || nx3>=lengthX3"));
   //#endif

   distributions->getDistribution(ftemp, nx1, nx2, nx3);
   LBMReal rho, vx1, vx2, vx3;
   calcMacrosFct(f, rho, vx1, vx2, vx3);

   distributions->swap();

   double dim = 0.01;

   switch (direction)
   {
   case D3Q27System::E:
      f[D3Q27System::E]   = ftemp[D3Q27System::E]   * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[D3Q27System::E]  ;// - rho*dim*D3Q27System::WEIGTH[D3Q27System::E];
      f[D3Q27System::NE]  = ftemp[D3Q27System::NE]  * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[D3Q27System::NE] ;// - rho*dim*D3Q27System::WEIGTH[D3Q27System::NE];
      f[D3Q27System::SE]  = ftemp[D3Q27System::SE]  * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[D3Q27System::SE] ;// - rho*dim*D3Q27System::WEIGTH[D3Q27System::SE];
      f[D3Q27System::TE]  = ftemp[D3Q27System::TE]  * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[D3Q27System::TE] ;// - rho*dim*D3Q27System::WEIGTH[D3Q27System::TE];
      f[D3Q27System::BE]  = ftemp[D3Q27System::BE]  * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[D3Q27System::BE] ;// - rho*dim*D3Q27System::WEIGTH[D3Q27System::BE];
      f[D3Q27System::TNE] = ftemp[D3Q27System::TNE] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[D3Q27System::TNE];// - rho*dim*D3Q27System::WEIGTH[D3Q27System::TNE];
      f[D3Q27System::TSE] = ftemp[D3Q27System::TSE] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[D3Q27System::TSE];// - rho*dim*D3Q27System::WEIGTH[D3Q27System::TSE];
      f[D3Q27System::BNE] = ftemp[D3Q27System::BNE] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[D3Q27System::BNE];// - rho*dim*D3Q27System::WEIGTH[D3Q27System::BNE];
      f[D3Q27System::BSE] = ftemp[D3Q27System::BSE] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[D3Q27System::BSE];// - rho*dim*D3Q27System::WEIGTH[D3Q27System::BSE];

      distributions->setDistributionInvForDirection(f[D3Q27System::E],   x1, x2, x3, D3Q27System::E);
      distributions->setDistributionInvForDirection(f[D3Q27System::NE],  x1, x2, x3, D3Q27System::NE);
      distributions->setDistributionInvForDirection(f[D3Q27System::SE],  x1, x2, x3, D3Q27System::SE);
      distributions->setDistributionInvForDirection(f[D3Q27System::TE],  x1, x2, x3, D3Q27System::TE);
      distributions->setDistributionInvForDirection(f[D3Q27System::BE],  x1, x2, x3, D3Q27System::BE);
      distributions->setDistributionInvForDirection(f[D3Q27System::TNE], x1, x2, x3, D3Q27System::TNE);
      distributions->setDistributionInvForDirection(f[D3Q27System::TSE], x1, x2, x3, D3Q27System::TSE);
      distributions->setDistributionInvForDirection(f[D3Q27System::BNE], x1, x2, x3, D3Q27System::BNE);
      distributions->setDistributionInvForDirection(f[D3Q27System::BSE], x1, x2, x3, D3Q27System::BSE);
      break;
   case D3Q27System::W:
      f[D3Q27System::W]   = ftemp[D3Q27System::W]   * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::W]    ;// - rho*dim*D3Q27System::WEIGTH[D3Q27System::W];
      f[D3Q27System::NW]  = ftemp[D3Q27System::NW]  * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::NW]   ;// - rho*dim*D3Q27System::WEIGTH[D3Q27System::NW];
      f[D3Q27System::SW]  = ftemp[D3Q27System::SW]  * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::SW]   ;// - rho*dim*D3Q27System::WEIGTH[D3Q27System::SW];
      f[D3Q27System::TW]  = ftemp[D3Q27System::TW]  * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::TW]   ;// - rho*dim*D3Q27System::WEIGTH[D3Q27System::TW];
      f[D3Q27System::BW]  = ftemp[D3Q27System::BW]  * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::BW]   ;// - rho*dim*D3Q27System::WEIGTH[D3Q27System::BW];
      f[D3Q27System::TNW] = ftemp[D3Q27System::TNW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::TNW];// - rho*dim*D3Q27System::WEIGTH[D3Q27System::TNW];
      f[D3Q27System::TSW] = ftemp[D3Q27System::TSW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::TSW];// - rho*dim*D3Q27System::WEIGTH[D3Q27System::TSW];
      f[D3Q27System::BNW] = ftemp[D3Q27System::BNW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::BNW];// - rho*dim*D3Q27System::WEIGTH[D3Q27System::BNW];
      f[D3Q27System::BSW] = ftemp[D3Q27System::BSW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::BSW];// - rho*dim*D3Q27System::WEIGTH[D3Q27System::BSW];

      distributions->setDistributionInvForDirection(f[D3Q27System::W],   x1, x2, x3, D3Q27System::W);
      distributions->setDistributionInvForDirection(f[D3Q27System::NW],  x1, x2, x3, D3Q27System::NW);
      distributions->setDistributionInvForDirection(f[D3Q27System::SW],  x1, x2, x3, D3Q27System::SW);
      distributions->setDistributionInvForDirection(f[D3Q27System::TW],  x1, x2, x3, D3Q27System::TW);
      distributions->setDistributionInvForDirection(f[D3Q27System::BW],  x1, x2, x3, D3Q27System::BW);
      distributions->setDistributionInvForDirection(f[D3Q27System::TNW], x1, x2, x3, D3Q27System::TNW);
      distributions->setDistributionInvForDirection(f[D3Q27System::TSW], x1, x2, x3, D3Q27System::TSW);
      distributions->setDistributionInvForDirection(f[D3Q27System::BNW], x1, x2, x3, D3Q27System::BNW);
      distributions->setDistributionInvForDirection(f[D3Q27System::BSW], x1, x2, x3, D3Q27System::BSW);
      break;
   case D3Q27System::N:
      f[D3Q27System::N]   = ftemp[D3Q27System::N]   * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[D3Q27System::N]   ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::N];
      f[D3Q27System::NE]  = ftemp[D3Q27System::NE]  * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[D3Q27System::NE]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::NE];
      f[D3Q27System::NW]  = ftemp[D3Q27System::NW]  * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[D3Q27System::NW]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::NW];
      f[D3Q27System::TN]  = ftemp[D3Q27System::TN]  * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[D3Q27System::TN]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TN];
      f[D3Q27System::BN]  = ftemp[D3Q27System::BN]  * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[D3Q27System::BN]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BN];
      f[D3Q27System::TNE] = ftemp[D3Q27System::TNE] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[D3Q27System::TNE] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TNE];
      f[D3Q27System::TNW] = ftemp[D3Q27System::TNW] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[D3Q27System::TNW] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TNW];
      f[D3Q27System::BNE] = ftemp[D3Q27System::BNE] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[D3Q27System::BNE] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BNE];
      f[D3Q27System::BNW] = ftemp[D3Q27System::BNW] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[D3Q27System::BNW] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BNW];

      distributions->setDistributionInvForDirection(f[D3Q27System::N], x1, x2, x3, D3Q27System::N);
      distributions->setDistributionInvForDirection(f[D3Q27System::NE], x1, x2, x3, D3Q27System::NE);
      distributions->setDistributionInvForDirection(f[D3Q27System::NW], x1, x2, x3, D3Q27System::NW);
      distributions->setDistributionInvForDirection(f[D3Q27System::TN], x1, x2, x3, D3Q27System::TN);
      distributions->setDistributionInvForDirection(f[D3Q27System::BN], x1, x2, x3, D3Q27System::BN);
      distributions->setDistributionInvForDirection(f[D3Q27System::TNE], x1, x2, x3, D3Q27System::TNE);
      distributions->setDistributionInvForDirection(f[D3Q27System::TNW], x1, x2, x3, D3Q27System::TNW);
      distributions->setDistributionInvForDirection(f[D3Q27System::BNE], x1, x2, x3, D3Q27System::BNE);
      distributions->setDistributionInvForDirection(f[D3Q27System::BNW], x1, x2, x3, D3Q27System::BNW);
      break;
   case D3Q27System::S:
      f[D3Q27System::S]   = ftemp[D3Q27System::S]   * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[D3Q27System::S]   ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::S];
      f[D3Q27System::SE]  = ftemp[D3Q27System::SE]  * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[D3Q27System::SE]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::SE];
      f[D3Q27System::SW]  = ftemp[D3Q27System::SW]  * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[D3Q27System::SW]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::SW];
      f[D3Q27System::TS]  = ftemp[D3Q27System::TS]  * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[D3Q27System::TS]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TS];
      f[D3Q27System::BS]  = ftemp[D3Q27System::BS]  * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[D3Q27System::BS]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BS];
      f[D3Q27System::TSE] = ftemp[D3Q27System::TSE] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[D3Q27System::TSE] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TSE];
      f[D3Q27System::TSW] = ftemp[D3Q27System::TSW] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[D3Q27System::TSW] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TSW];
      f[D3Q27System::BSE] = ftemp[D3Q27System::BSE] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[D3Q27System::BSE] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BSE];
      f[D3Q27System::BSW] = ftemp[D3Q27System::BSW] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[D3Q27System::BSW] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BSW];

      distributions->setDistributionInvForDirection(f[D3Q27System::S], x1, x2, x3, D3Q27System::S);
      distributions->setDistributionInvForDirection(f[D3Q27System::SE], x1, x2, x3, D3Q27System::SE);
      distributions->setDistributionInvForDirection(f[D3Q27System::SW], x1, x2, x3, D3Q27System::SW);
      distributions->setDistributionInvForDirection(f[D3Q27System::TS], x1, x2, x3, D3Q27System::TS);
      distributions->setDistributionInvForDirection(f[D3Q27System::BS], x1, x2, x3, D3Q27System::BS);
      distributions->setDistributionInvForDirection(f[D3Q27System::TSE], x1, x2, x3, D3Q27System::TSE);
      distributions->setDistributionInvForDirection(f[D3Q27System::TSW], x1, x2, x3, D3Q27System::TSW);
      distributions->setDistributionInvForDirection(f[D3Q27System::BSE], x1, x2, x3, D3Q27System::BSE);
      distributions->setDistributionInvForDirection(f[D3Q27System::BSW], x1, x2, x3, D3Q27System::BSW);
      break;
   case D3Q27System::T:
      f[D3Q27System::T]   = ftemp[D3Q27System::T]   * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[D3Q27System::T]   ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::T];
      f[D3Q27System::TE]  = ftemp[D3Q27System::TE]  * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[D3Q27System::TE]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TE];
      f[D3Q27System::TW]  = ftemp[D3Q27System::TW]  * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[D3Q27System::TW]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TW];
      f[D3Q27System::TN]  = ftemp[D3Q27System::TN]  * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[D3Q27System::TN]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TN];
      f[D3Q27System::TS]  = ftemp[D3Q27System::TS]  * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[D3Q27System::TS]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TS];
      f[D3Q27System::TNE] = ftemp[D3Q27System::TNE] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[D3Q27System::TNE] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TNE];
      f[D3Q27System::TNW] = ftemp[D3Q27System::TNW] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[D3Q27System::TNW] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TNW];
      f[D3Q27System::TSE] = ftemp[D3Q27System::TSE] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[D3Q27System::TSE] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TSE];
      f[D3Q27System::TSW] = ftemp[D3Q27System::TSW] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[D3Q27System::TSW] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::TSW];

      distributions->setDistributionInvForDirection(f[D3Q27System::T], x1, x2, x3, D3Q27System::T);
      distributions->setDistributionInvForDirection(f[D3Q27System::TE], x1, x2, x3, D3Q27System::TE);
      distributions->setDistributionInvForDirection(f[D3Q27System::TW], x1, x2, x3, D3Q27System::TW);
      distributions->setDistributionInvForDirection(f[D3Q27System::TN], x1, x2, x3, D3Q27System::TN);
      distributions->setDistributionInvForDirection(f[D3Q27System::TS], x1, x2, x3, D3Q27System::TS);
      distributions->setDistributionInvForDirection(f[D3Q27System::TNE], x1, x2, x3, D3Q27System::TNE);
      distributions->setDistributionInvForDirection(f[D3Q27System::TNW], x1, x2, x3, D3Q27System::TNW);
      distributions->setDistributionInvForDirection(f[D3Q27System::TSE], x1, x2, x3, D3Q27System::TSE);
      distributions->setDistributionInvForDirection(f[D3Q27System::TSW], x1, x2, x3, D3Q27System::TSW);
      break;
   case D3Q27System::B:
      f[D3Q27System::B]   = ftemp[D3Q27System::B]   * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[D3Q27System::B]   ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::B];
      f[D3Q27System::BE]  = ftemp[D3Q27System::BE]  * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[D3Q27System::BE]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BE];
      f[D3Q27System::BW]  = ftemp[D3Q27System::BW]  * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[D3Q27System::BW]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BW];
      f[D3Q27System::BN]  = ftemp[D3Q27System::BN]  * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[D3Q27System::BN]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BN];
      f[D3Q27System::BS]  = ftemp[D3Q27System::BS]  * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[D3Q27System::BS]  ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BS];
      f[D3Q27System::BNE] = ftemp[D3Q27System::BNE] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[D3Q27System::BNE] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BNE];
      f[D3Q27System::BNW] = ftemp[D3Q27System::BNW] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[D3Q27System::BNW] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BNW];
      f[D3Q27System::BSE] = ftemp[D3Q27System::BSE] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[D3Q27System::BSE] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BSE];
      f[D3Q27System::BSW] = ftemp[D3Q27System::BSW] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[D3Q27System::BSW] ;//- rho*dim*D3Q27System::WEIGTH[D3Q27System::BSW];

      distributions->setDistributionInvForDirection(f[D3Q27System::B], x1, x2, x3, D3Q27System::B);
      distributions->setDistributionInvForDirection(f[D3Q27System::BE], x1, x2, x3, D3Q27System::BE);
      distributions->setDistributionInvForDirection(f[D3Q27System::BW], x1, x2, x3, D3Q27System::BW);
      distributions->setDistributionInvForDirection(f[D3Q27System::BN], x1, x2, x3, D3Q27System::BN);
      distributions->setDistributionInvForDirection(f[D3Q27System::BS], x1, x2, x3, D3Q27System::BS);
      distributions->setDistributionInvForDirection(f[D3Q27System::BNE], x1, x2, x3, D3Q27System::BNE);
      distributions->setDistributionInvForDirection(f[D3Q27System::BNW], x1, x2, x3, D3Q27System::BNW);
      distributions->setDistributionInvForDirection(f[D3Q27System::BSE], x1, x2, x3, D3Q27System::BSE);
      distributions->setDistributionInvForDirection(f[D3Q27System::BSW], x1, x2, x3, D3Q27System::BSW);
      break;
   default:
      UB_THROW(UbException(UB_EXARGS, "It isn't implemented non reflecting density boundary for this direction!"));
   }
}

