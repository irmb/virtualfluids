#include "SlipBCAlgorithm.h"

SlipBCAlgorithm::SlipBCAlgorithm()
{
   BCAlgorithm::type = BCAlgorithm::SlipBCAlgorithm;
   BCAlgorithm::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
SlipBCAlgorithm::~SlipBCAlgorithm()
{

}
//////////////////////////////////////////////////////////////////////////
BCAlgorithmPtr SlipBCAlgorithm::clone()
{
   BCAlgorithmPtr bc(new SlipBCAlgorithm());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void SlipBCAlgorithm::addDistributions(DistributionArray3DPtr distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void SlipBCAlgorithm::applyBC()
{
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal feq[D3Q27System::ENDF+1];
   distributions->getDistributionInv(f, x1, x2, x3);
   LBMReal rho, vx1, vx2, vx3, drho;
   calcMacrosFct(f, drho, vx1, vx2, vx3);
   calcFeqFct(feq, drho, vx1, vx2, vx3);

   UbTupleFloat3 normale = bcPtr->getNormalVector();
   LBMReal amp = vx1*val<1>(normale)+vx2*val<2>(normale)+vx3*val<3>(normale);

   vx1 = vx1 - amp * val<1>(normale); //normale zeigt von struktur weg!
   vx2 = vx2 - amp * val<2>(normale); //normale zeigt von struktur weg!
   vx3 = vx3 - amp * val<3>(normale); //normale zeigt von struktur weg!

   rho = 1.0+drho*compressibleFactor;

   for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
   {
      if (bcPtr->hasSlipBoundaryFlag(fdir))
      {
         //quadratic bounce back
         const int invDir = D3Q27System::INVDIR[fdir];
         LBMReal q = bcPtr->getQ(invDir);// m+m q=0 stabiler
         //vx3=0;
         LBMReal velocity = 0.0;
         switch (invDir)
         {
         case D3Q27System::E: velocity = (UbMath::c4o9*(+vx1)); break;      //(2/cs^2)(=6)*rho_0(=1 bei imkompr)*wi*u*ei mit cs=1/sqrt(3)
         case D3Q27System::W: velocity = (UbMath::c4o9*(-vx1)); break;      //z.B. aus paper manfred MRT LB models in three dimensions (2002)   
         case D3Q27System::N: velocity = (UbMath::c4o9*(+vx2)); break;
         case D3Q27System::S: velocity = (UbMath::c4o9*(-vx2)); break;
         case D3Q27System::T: velocity = (UbMath::c4o9*(+vx3)); break;
         case D3Q27System::B: velocity = (UbMath::c4o9*(-vx3)); break;
         case D3Q27System::NE: velocity = (UbMath::c1o9*(+vx1+vx2)); break;
         case D3Q27System::SW: velocity = (UbMath::c1o9*(-vx1-vx2)); break;
         case D3Q27System::SE: velocity = (UbMath::c1o9*(+vx1-vx2)); break;
         case D3Q27System::NW: velocity = (UbMath::c1o9*(-vx1+vx2)); break;
         case D3Q27System::TE: velocity = (UbMath::c1o9*(+vx1             +vx3)); break;
         case D3Q27System::BW: velocity = (UbMath::c1o9*(-vx1             -vx3)); break;
         case D3Q27System::BE: velocity = (UbMath::c1o9*(+vx1             -vx3)); break;
         case D3Q27System::TW: velocity = (UbMath::c1o9*(-vx1             +vx3)); break;
         case D3Q27System::TN: velocity = (UbMath::c1o9*(+vx2+vx3)); break;
         case D3Q27System::BS: velocity = (UbMath::c1o9*(-vx2-vx3)); break;
         case D3Q27System::BN: velocity = (UbMath::c1o9*(+vx2-vx3)); break;
         case D3Q27System::TS: velocity = (UbMath::c1o9*(-vx2+vx3)); break;
         case D3Q27System::TNE: velocity = (UbMath::c1o36*(+vx1+vx2+vx3)); break;
         case D3Q27System::BSW: velocity = (UbMath::c1o36*(-vx1-vx2-vx3)); break;
         case D3Q27System::BNE: velocity = (UbMath::c1o36*(+vx1+vx2-vx3)); break;
         case D3Q27System::TSW: velocity = (UbMath::c1o36*(-vx1-vx2+vx3)); break;
         case D3Q27System::TSE: velocity = (UbMath::c1o36*(+vx1-vx2+vx3)); break;
         case D3Q27System::BNW: velocity = (UbMath::c1o36*(-vx1+vx2-vx3)); break;
         case D3Q27System::BSE: velocity = (UbMath::c1o36*(+vx1-vx2-vx3)); break;
         case D3Q27System::TNW: velocity = (UbMath::c1o36*(-vx1+vx2+vx3)); break;
         default: throw UbException(UB_EXARGS, "unknown error");
         }
         LBMReal fReturn = ((1.0-q)/(1.0+q))*((f[invDir]-feq[invDir])/(1.0-collFactor)+feq[invDir])+((q*(f[invDir]+f[fdir])-velocity*rho)/(1.0+q));
         distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
      }
   }
}