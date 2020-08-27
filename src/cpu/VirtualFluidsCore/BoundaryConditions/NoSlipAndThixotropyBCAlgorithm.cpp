#include "NoSlipAndThixotropyBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

NoSlipAndThixotropyBCAlgorithm::NoSlipAndThixotropyBCAlgorithm()
{
	BCAlgorithm::type = BCAlgorithm::NoSlipAndThixotropyBCAlgorithm;
	BCAlgorithm::preCollision = false;
	BCAlgorithm::thixotropy = true;
	
}
//////////////////////////////////////////////////////////////////////////
NoSlipAndThixotropyBCAlgorithm::~NoSlipAndThixotropyBCAlgorithm()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> NoSlipAndThixotropyBCAlgorithm::clone()
{
	SPtr<BCAlgorithm> bc(new NoSlipAndThixotropyBCAlgorithm());
	return bc;
}
//////////////////////////////////////////////////////////////////////////
void NoSlipAndThixotropyBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
	this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
//void NoSlipAndThixotropyBCAlgorithm::addDistributionsF(SPtr<DistributionArray3D> distributions)
//{
//	this->distributionsf = distributions;
//}
//////////////////////////////////////////////////////////////////////////
void NoSlipAndThixotropyBCAlgorithm::addDistributionsH(SPtr<DistributionArray3D> distributions)
{
	this->distributionsH = distributions;
}
//////////////////////////////////////////////////////////////////////////
void NoSlipAndThixotropyBCAlgorithm::applyBC()
{
	LBMReal f[D3Q27System::ENDF + 1];
	LBMReal feq[D3Q27System::ENDF + 1];
	LBMReal h[D3Q27System::ENDF + 1];
	LBMReal heq[D3Q27System::ENDF + 1];
	distributions->getDistributionInv(f, x1, x2, x3);
	distributionsH->getDistributionInv(h, x1, x2, x3);
	LBMReal rho, vx1, vx2, vx3;//, concentration, fl1, fl2, fl3, m100;
	calcMacrosFct(f, rho, vx1, vx2, vx3);
	calcFeqFct(feq, rho, vx1, vx2, vx3);

	//calcDiffusionMacrosFctPost(h, concentration, fl1, fl2, fl3, m100, collFactor);
	LBMReal lambda = D3Q27System::getDensity(h);
	D3Q27System::calcCompFeq(heq, lambda, 0., 0., 0.);

	for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++)
	{
		if (bcPtr->hasNoSlipBoundaryFlag(fdir))
		{
			//quadratic bounce back
			const int invDir = D3Q27System::INVDIR[fdir];
			LBMReal q = bcPtr->getQ(invDir);
			LBMReal fReturnf = ((1.0 - q) / (1.0 + q))*((f[invDir] - feq[invDir]) / (1.0 - collFactor) + feq[invDir]) + ((q / (1.0 + q))*(f[invDir] + f[fdir]));
			LBMReal fReturnh = ((1.0 - q) / (1.0 + q))*((h[invDir] - heq[invDir]) / (1.0 - collFactor) + heq[invDir]) + ((q / (1.0 + q))*(h[invDir] + h[fdir]));

			distributions->setDistributionForDirection(fReturnf, x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);
			distributionsH->setDistributionForDirection(fReturnh, x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);

		}
	}
}
