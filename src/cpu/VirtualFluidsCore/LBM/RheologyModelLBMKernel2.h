#ifndef RheologyModelLBMKernel2_H
#define RheologyModelLBMKernel2_H

#include "LBMKernel.h"
#include "BCProcessor.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

class RheologyModelLBMKernel2;

//! \brief Base class for model of thixotropy based on K16. Use Template Method design pattern for Implementation of different models. 
//! \author K. Kutscher, M. Geier
class RheologyModelLBMKernel2 : public LBMKernel
{
public:
	RheologyModelLBMKernel2();
	virtual ~RheologyModelLBMKernel2();
	void calculate(int step);
	virtual SPtr<LBMKernel> clone() { UB_THROW(UbException("SPtr<LBMKernel> clone() - belongs in the derived class")); };
	double getCalculationTime();

	void swapDistributions();

protected:
	void initDataSet();

	virtual real getRheologyCollFactor(real omegaInf, real shearRate, real drho) const { UB_THROW(UbException("real getRheologyCollFactor() - belongs in the derived class")); }

	real f[D3Q27System::ENDF + 1];

	UbTimer timer;

	real OxyyMxzz;
	
	CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF;
	CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF;
	CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsF;

	mu::value_type muX1, muX2, muX3;
	mu::value_type muDeltaT;
	mu::value_type muNu;
	real forcingX1;
	real forcingX2;
	real forcingX3;

	bool test;
};

#endif
