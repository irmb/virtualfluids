#include "ThinWallNoSlipBCAlgorithm.h"

#include "BoundaryConditions.h"
#include "D3Q27EsoTwist3DSplittedVector.h"

ThinWallNoSlipBCAlgorithm::ThinWallNoSlipBCAlgorithm()
{
    BCAlgorithm::type         = BCAlgorithm::ThinWallNoSlipBCAlgorithm;
    BCAlgorithm::preCollision = false;
    pass                      = 1;
}
//////////////////////////////////////////////////////////////////////////
ThinWallNoSlipBCAlgorithm::~ThinWallNoSlipBCAlgorithm() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> ThinWallNoSlipBCAlgorithm::clone()
{
    SPtr<BCAlgorithm> bc(new ThinWallNoSlipBCAlgorithm());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void ThinWallNoSlipBCAlgorithm::applyBC()
{
    LBMReal f[D3Q27System::ENDF + 1];
    LBMReal feq[D3Q27System::ENDF + 1];
    distributions->getDistributionInv(f, x1, x2, x3);
    LBMReal rho, vx1, vx2, vx3;
    calcMacrosFct(f, rho, vx1, vx2, vx3);
    calcFeqFct(feq, rho, vx1, vx2, vx3);

    LBMReal fReturn;

    for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
        if (bcPtr->hasNoSlipBoundaryFlag(fdir)) {
            const int invDir = D3Q27System::INVDIR[fdir];
            if (pass == 1) {
                LBMReal q = bcPtr->getQ(invDir);
                fReturn   = ((1.0 - q) / (1.0 + q)) * 0.5 *
                          (f[invDir] - f[fdir] +
                           (f[invDir] + f[fdir] - collFactor * (feq[fdir] + feq[invDir])) / (1.0 - collFactor));
                // distributionsTemp->setDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir], x2 +
                // D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);
                fTemp[fdir] = fReturn;
            } else {
                // quadratic bounce back with for thin walls
                // fReturn = distributionsTemp->getDistributionInvForDirection(x1 + D3Q27System::DX1[invDir], x2 +
                // D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);
                fReturn = fTemp[fdir];
                distributions->setDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir],
                                                           x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir],
                                                           fdir);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void ThinWallNoSlipBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void ThinWallNoSlipBCAlgorithm::setPass(int pass) { this->pass = pass; }
