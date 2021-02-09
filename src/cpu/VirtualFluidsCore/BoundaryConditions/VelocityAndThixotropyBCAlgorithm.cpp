#include "VelocityAndThixotropyBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

VelocityAndThixotropyBCAlgorithm::VelocityAndThixotropyBCAlgorithm()
{
	BCAlgorithm::type = BCAlgorithm::VelocityAndThixotropyBCAlgorithm;
	BCAlgorithm::preCollision = false;
	BCAlgorithm::thixotropy = true;
	lambdaBC = 0.0;
}
//////////////////////////////////////////////////////////////////////////
VelocityAndThixotropyBCAlgorithm::~VelocityAndThixotropyBCAlgorithm()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> VelocityAndThixotropyBCAlgorithm::clone()
{
	SPtr<BCAlgorithm> bc(new VelocityAndThixotropyBCAlgorithm());
	dynamicPointerCast<VelocityAndThixotropyBCAlgorithm>(bc)->setLambdaBC(lambdaBC);
	return bc;
}
//////////////////////////////////////////////////////////////////////////
void VelocityAndThixotropyBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
	this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
//void VelocityAndThixotropyBCAlgorithm::addDistributionsF(DistributionArray3DPtr distributions)
//{
//	this->distributionsf = distributions;
//}
//////////////////////////////////////////////////////////////////////////
void VelocityAndThixotropyBCAlgorithm::addDistributionsH(SPtr<DistributionArray3D> distributions)
{
	this->distributionsH = distributions;
}
//////////////////////////////////////////////////////////////////////////
void VelocityAndThixotropyBCAlgorithm::applyBC()
{
	LBMReal f[D3Q27System::ENDF + 1];
	LBMReal feq[D3Q27System::ENDF + 1];
	LBMReal h[D3Q27System::ENDF + 1];

	distributions->getDistributionInv(f, x1, x2, x3);
	distributionsH->getDistributionInv(h, x1, x2, x3);
	
	LBMReal rho, vx1, vx2, vx3, drho;
	calcMacrosFct(f, drho, vx1, vx2, vx3);
	calcFeqFct(feq, drho, vx1, vx2, vx3);

	rho = 1.0 + drho * compressibleFactor;

	//calcDiffusionMacrosFctPost(h, concentration, fl1, fl2, fl3, m100, collFactor);
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

	//lambdaBC = bcPtr->getBoundaryThixotropy();

	//LBMReal rhoBC = bcPtr->getBoundaryDensity();

	//for (int fdir = D3Q27System::STARTF; fdir <= D3Q27System::ENDF; fdir++)
	//{
	//	if (bcPtr->hasDensityBoundaryFlag(fdir))
	//	{
	//		LBMReal ftemp = calcFeqsForDirFct(fdir, rho, vx1, vx2, vx3);
	//		ftemp = calcFeqsForDirFct(fdir, rhoBC, vx1, vx2, vx3) + f[fdir] - ftemp;
	//		distributions->setDistributionForDirection(ftemp, nx1, nx2, nx3, fdir);

	//		LBMReal htemp = D3Q27System::getCompFeqForDirection(fdir, lambda, vx1, vx2, vx3);
	//		htemp = D3Q27System::getCompFeqForDirection(fdir,lambdaBC, vx1, vx2, vx3) + h[fdir] - htemp;
	//		distributionsH->setDistributionForDirection(htemp, nx1, nx2, nx3, fdir);
	//	}
	//}

	for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++)
	{
		if (bcPtr->hasVelocityBoundaryFlag(fdir))
		{
			const int invDir = D3Q27System::INVDIR[fdir];
			LBMReal q = bcPtr->getQ(invDir);// m+m q=0 stabiler
			LBMReal velocity = bcPtr->getBoundaryVelocity(invDir);
			LBMReal fReturn = ((1.0 - q) / (1.0 + q)) * ((f[invDir] - feq[invDir]) / (1.0 - collFactor) + feq[invDir]) + ((q * (f[invDir] + f[fdir]) - velocity * rho) / (1.0 + q));
			distributions->setDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);

			LBMReal htemp = D3Q27System::getCompFeqForDirection(fdir, lambda, vx1, vx2, vx3);
			htemp = D3Q27System::getCompFeqForDirection(fdir, lambdaBC, vx1, vx2, vx3) + h[fdir] - htemp;
			distributionsH->setDistributionForDirection(htemp, nx1, nx2, nx3, fdir);
		}
	}
}
