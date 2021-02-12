#include "DensityAndThixotropyBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

DensityAndThixotropyBCAlgorithm::DensityAndThixotropyBCAlgorithm()
{
	BCAlgorithm::type = BCAlgorithm::DensityAndThixotropyBCAlgorithm;
	BCAlgorithm::preCollision = false;
	BCAlgorithm::thixotropy = true;
	lambdaBC = 0.0;
}
//////////////////////////////////////////////////////////////////////////
DensityAndThixotropyBCAlgorithm::~DensityAndThixotropyBCAlgorithm()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> DensityAndThixotropyBCAlgorithm::clone()
{
	SPtr<BCAlgorithm> bc(new DensityAndThixotropyBCAlgorithm());
	dynamicPointerCast<DensityAndThixotropyBCAlgorithm>(bc)->setLambdaBC(lambdaBC);
	return bc;
}
//////////////////////////////////////////////////////////////////////////
void DensityAndThixotropyBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
	this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
//void DensityAndThixotropyBCAlgorithm::addDistributionsF(DistributionArray3DPtr distributions)
//{
//	this->distributionsf = distributions;
//}
//////////////////////////////////////////////////////////////////////////
void DensityAndThixotropyBCAlgorithm::addDistributionsH(SPtr<DistributionArray3D> distributions)
{
	this->distributionsH = distributions;
}
//////////////////////////////////////////////////////////////////////////
void DensityAndThixotropyBCAlgorithm::applyBC()
{
   using namespace D3Q27System;

	LBMReal f[D3Q27System::ENDF + 1];
	LBMReal feq[D3Q27System::ENDF + 1];
	LBMReal h[D3Q27System::ENDF + 1];
	LBMReal heq[D3Q27System::ENDF + 1];
	distributions->getDistributionInv(f, x1, x2, x3);
	distributionsH->getDistributionInv(h, x1, x2, x3);
	
	LBMReal rho, vx1, vx2, vx3;
	
	calcMacrosFct(f, rho, vx1, vx2, vx3);
	calcFeqFct(feq, rho, vx1, vx2, vx3);

	LBMReal lambda = D3Q27System::getDensity(h);
	D3Q27System::calcCompFeq(heq, lambda, vx1, vx2, vx3);


	int nx1 = x1;
	int nx2 = x2;
	int nx3 = x3;

	//flag points in direction of fluid
	if (bcPtr->hasDensityBoundaryFlag(D3Q27System::E)) { nx1 -= 1; }
	else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::W)) { nx1 += 1; }
	else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::N)) { nx2 -= 1; }
	else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::S)) { nx2 += 1; }
	else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::T)) { nx3 -= 1; }
	else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::B)) { nx3 += 1; }
	else	 UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on density boundary..."));

	LBMReal rhoBC = bcPtr->getBoundaryDensity();

	for (int fdir = D3Q27System::STARTF; fdir <= D3Q27System::ENDF; fdir++)
	{
		if (bcPtr->hasDensityBoundaryFlag(fdir))
		{
			LBMReal ftemp = calcFeqsForDirFct(fdir, rho, vx1, vx2, vx3);
			ftemp = calcFeqsForDirFct(fdir, rhoBC, vx1, vx2, vx3) + f[fdir] - ftemp;
			distributions->setDistributionForDirection(ftemp, nx1, nx2, nx3, fdir);

			LBMReal htemp = D3Q27System::getCompFeqForDirection(fdir, lambda, vx1, vx2, vx3);
			htemp = D3Q27System::getCompFeqForDirection(fdir,lambdaBC, vx1, vx2, vx3) + h[fdir] - htemp;
			distributionsH->setDistributionForDirection(htemp, nx1, nx2, nx3, fdir);
		}
	}
}
