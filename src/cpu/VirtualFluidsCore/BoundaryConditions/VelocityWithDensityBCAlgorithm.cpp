#include "VelocityWithDensityBCAlgorithm.h"
#include "BCArray3D.h"
#include "DistributionArray3D.h"

VelocityWithDensityBCAlgorithm::VelocityWithDensityBCAlgorithm()
{
    BCAlgorithm::type         = BCAlgorithm::VelocityWithDensityBCAlgorithm;
    BCAlgorithm::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
VelocityWithDensityBCAlgorithm::~VelocityWithDensityBCAlgorithm() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> VelocityWithDensityBCAlgorithm::clone()
{
    SPtr<BCAlgorithm> bc(new VelocityWithDensityBCAlgorithm());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void VelocityWithDensityBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void VelocityWithDensityBCAlgorithm::applyBC()
{
    // velocity bc for non reflecting pressure bc
    LBMReal f[D3Q27System::ENDF + 1];
    LBMReal feq[D3Q27System::ENDF + 1];
    distributions->getDistributionInv(f, x1, x2, x3);
    LBMReal rho, vx1, vx2, vx3, drho;
    calcMacrosFct(f, drho, vx1, vx2, vx3);
    calcFeqFct(feq, drho, vx1, vx2, vx3);

    rho = 1.0 + drho * compressibleFactor;

    for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
        // if (bcPtr->hasVelocityBoundaryFlag(fdir))
        //{
        //   const int invDir = D3Q27System::INVDIR[fdir];
        //   LBMReal q = bcPtr->getQ(invDir);// m+m q=0 stabiler
        //   LBMReal velocity = bcPtr->getBoundaryVelocity(invDir);
        //   //normal velocity bc: LBMReal fReturn =
        //   ((1.0-q)/(1.0+q))*((f[invDir]-feq[invDir])/(1.0-collFactor)+feq[invDir])+((q*(f[invDir]+f[fdir])-velocity*rho)/(1.0+q));
        //   //LBMReal fReturn = ((1.0 - q) / (1.0 + q))*((f[invDir] - feq[invDir]) / (1.0 - collFactor) + feq[invDir])
        //   + ((q*(f[invDir] + f[fdir]) - velocity) / (1.0 + q))-drho*D3Q27System::WEIGTH[invDir];
        //   //LBMReal fReturn = ((1.0 - q) / (1.0 + q))*((f[invDir] - feq[invDir]) / (1.0 - collFactor) + feq[invDir])
        //   + ((q*(f[invDir] + f[fdir]) - velocity*rho) / (1.0 + q))-drho*D3Q27System::WEIGTH[invDir]; LBMReal fReturn
        //   = ((1.0 - q) / (1.0 + q))*((f[invDir] - feq[invDir]*collFactor) / (1.0 - collFactor)) + ((q*(f[invDir] +
        //   f[fdir]) - velocity*rho) / (1.0 + q))-drho*D3Q27System::WEIGTH[invDir];
        //
        //   distributions->setDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir], x2 +
        //   D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);
        //}

        int nX1 = x1 + D3Q27System::DX1[fdir];
        int nX2 = x2 + D3Q27System::DX2[fdir];
        int nX3 = x3 + D3Q27System::DX3[fdir];

        int minX1 = 0;
        int minX2 = 0;
        int minX3 = 0;

        int maxX1 = (int)bcArray->getNX1();
        int maxX2 = (int)bcArray->getNX2();
        int maxX3 = (int)bcArray->getNX3();

        if (minX1 <= nX1 && maxX1 > nX1 && minX2 <= nX2 && maxX2 > nX2 && minX3 <= nX3 && maxX3 > nX3) {
            if (bcArray->isSolid(nX1, nX2, nX3)) {
                const int invDir = D3Q27System::INVDIR[fdir];
                //            LBMReal q =1.0;// bcPtr->getQ(invDir);// m+m q=0 stabiler
                LBMReal velocity = bcPtr->getBoundaryVelocity(fdir);
                //            LBMReal fReturn = ((1.0 - q) / (1.0 + q))*((f[fdir] - feq[fdir]*collFactor) / (1.0 -
                //            collFactor)) + ((q*(f[fdir] + f[invDir]) - velocity*rho) / (1.0 +
                //            q))-drho*D3Q27System::WEIGTH[invDir];

                // if q=1
                // LBMReal fReturn = ((q*(f[fdir] + f[invDir]) - velocity*rho) / (1.0 +
                // q))-drho*D3Q27System::WEIGTH[invDir];
                LBMReal fReturn = (f[fdir] + f[invDir] - velocity * rho) / 2.0 - drho * D3Q27System::WEIGTH[invDir];

                distributions->setDistributionForDirection(fReturn, nX1, nX2, nX3, invDir);
            }
        }
    }
}
