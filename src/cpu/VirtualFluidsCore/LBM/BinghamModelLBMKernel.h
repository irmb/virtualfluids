#ifndef BinghamModelLBMKernel_H
#define BinghamModelLBMKernel_H

#include "ThixotropyModelLBMKernel2.h"
#include "Thixotropy.h"

//! \brief    Cumulant LBM kernel + Bingham plastic model 
//! \author K. Kutscher, M. Geier
class BinghamModelLBMKernel : public ThixotropyModelLBMKernel2
{
public:
	BinghamModelLBMKernel() {};
	~BinghamModelLBMKernel() {};
	SPtr<LBMKernel> clone() override
	{
		SPtr<LBMKernel> kernel(new BinghamModelLBMKernel());
		kernel->setNX(nx);
		kernel->setCollisionFactor(collFactor);
		dynamicPointerCast<BinghamModelLBMKernel>(kernel)->initDataSet();
		kernel->setBCProcessor(bcProcessor->clone(kernel));
		kernel->setWithForcing(withForcing);
		kernel->setForcingX1(muForcingX1);
		kernel->setForcingX2(muForcingX2);
		kernel->setForcingX3(muForcingX3);
		kernel->setIndex(ix1, ix2, ix3);
		kernel->setDeltaT(deltaT);

		return kernel;
	}
protected:	
	LBMReal getThyxotropyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho) const override
	{
		return Thixotropy::getBinghamCollFactor(omegaInf, shearRate, drho);
	}
};


#endif
