#include "SlipBCStrategy.h"
#include "BoundaryConditions.h"
#include "DistributionArray3D.h"

SlipBCStrategy::SlipBCStrategy()
{
    BCStrategy::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
SlipBCStrategy::~SlipBCStrategy() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> SlipBCStrategy::clone()
{
    SPtr<BCStrategy> bc(new SlipBCStrategy());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void SlipBCStrategy::addDistributions(SPtr<DistributionArray3D> distributions) { this->distributions = distributions; }
//////////////////////////////////////////////////////////////////////////
void SlipBCStrategy::applyBC()
{
    using namespace vf::lbm::dir;

    real f[D3Q27System::ENDF + 1];
    real feq[D3Q27System::ENDF + 1];
    distributions->getDistributionInv(f, x1, x2, x3);
    real rho, vx1, vx2, vx3, drho;
    calcMacrosFct(f, drho, vx1, vx2, vx3);
    calcFeqFct(feq, drho, vx1, vx2, vx3);

    UbTupleFloat3 normale = bcPtr->getNormalVector();
    real amp            = vx1 * val<1>(normale) + vx2 * val<2>(normale) + vx3 * val<3>(normale);

    vx1 = vx1 - amp * val<1>(normale); // normale zeigt von struktur weg!
    vx2 = vx2 - amp * val<2>(normale); // normale zeigt von struktur weg!
    vx3 = vx3 - amp * val<3>(normale); // normale zeigt von struktur weg!

    rho = vf::basics::constant::c1o1 + drho * compressibleFactor;

   for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
   {
      if (bcPtr->hasSlipBoundaryFlag(fdir))
      {
         //quadratic bounce back
         const int invDir = D3Q27System::INVDIR[fdir];
         real q = bcPtr->getQ(invDir);// m+m q=0 stabiler
         //vx3=0;
         real velocity = vf::basics::constant::c0o1;
         switch (invDir)
         {
         case DIR_P00: velocity = (vf::basics::constant::c4o9*(+vx1)); break;      //(2/cs^2)(=6)*rho_0(=1 bei imkompr)*wi*u*ei mit cs=1/sqrt(3)
         case DIR_M00: velocity = (vf::basics::constant::c4o9*(-vx1)); break;      //z.B. aus paper manfred MRT LB models in three dimensions (2002)   
         case DIR_0P0: velocity = (vf::basics::constant::c4o9*(+vx2)); break;
         case DIR_0M0: velocity = (vf::basics::constant::c4o9*(-vx2)); break;
         case DIR_00P: velocity = (vf::basics::constant::c4o9*(+vx3)); break;
         case DIR_00M: velocity = (vf::basics::constant::c4o9*(-vx3)); break;
         case DIR_PP0: velocity = (vf::basics::constant::c1o9*(+vx1+vx2)); break;
         case DIR_MM0: velocity = (vf::basics::constant::c1o9*(-vx1-vx2)); break;
         case DIR_PM0: velocity = (vf::basics::constant::c1o9*(+vx1-vx2)); break;
         case DIR_MP0: velocity = (vf::basics::constant::c1o9*(-vx1+vx2)); break;
         case DIR_P0P: velocity = (vf::basics::constant::c1o9*(+vx1+vx3)); break;
         case DIR_M0M: velocity = (vf::basics::constant::c1o9*(-vx1-vx3)); break;
         case DIR_P0M: velocity = (vf::basics::constant::c1o9*(+vx1-vx3)); break;
         case DIR_M0P: velocity = (vf::basics::constant::c1o9*(-vx1+vx3)); break;
         case DIR_0PP: velocity = (vf::basics::constant::c1o9*(+vx2+vx3)); break;
         case DIR_0MM: velocity = (vf::basics::constant::c1o9*(-vx2-vx3)); break;
         case DIR_0PM: velocity = (vf::basics::constant::c1o9*(+vx2-vx3)); break;
         case DIR_0MP: velocity = (vf::basics::constant::c1o9*(-vx2+vx3)); break;
         case DIR_PPP: velocity = (vf::basics::constant::c1o36*(+vx1+vx2+vx3)); break;
         case DIR_MMM: velocity = (vf::basics::constant::c1o36*(-vx1-vx2-vx3)); break;
         case DIR_PPM: velocity = (vf::basics::constant::c1o36*(+vx1+vx2-vx3)); break;
         case DIR_MMP: velocity = (vf::basics::constant::c1o36*(-vx1-vx2+vx3)); break;
         case DIR_PMP: velocity = (vf::basics::constant::c1o36*(+vx1-vx2+vx3)); break;
         case DIR_MPM: velocity = (vf::basics::constant::c1o36*(-vx1+vx2-vx3)); break;
         case DIR_PMM: velocity = (vf::basics::constant::c1o36*(+vx1-vx2-vx3)); break;
         case DIR_MPP: velocity = (vf::basics::constant::c1o36*(-vx1+vx2+vx3)); break;
         default: throw UbException(UB_EXARGS, "unknown error");
         }
         real fReturn = ((vf::basics::constant::c1o1-q)/(vf::basics::constant::c1o1+q))*((f[invDir]-feq[invDir])/(vf::basics::constant::c1o1-collFactor)+feq[invDir])+((q*(f[invDir]+f[fdir])-velocity*rho)/(vf::basics::constant::c1o1+q));
         distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
      }
   }
}