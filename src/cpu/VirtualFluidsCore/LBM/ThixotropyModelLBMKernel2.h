#ifndef ThixotropyModelLBMKernel2_H
#define ThixotropyModelLBMKernel2_H

#include "LBMKernel.h"
#include "BCProcessor.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

class ThixotropyModelLBMKernel2;

//! \brief Base class for model of thixotropy based on K16. Use Template Method design pattern for Implementation of different models. 
//! \author K. Kutscher, M. Geier
class ThixotropyModelLBMKernel2 : public LBMKernel
{
public:
	ThixotropyModelLBMKernel2();
	virtual ~ThixotropyModelLBMKernel2();
	void calculate(int step);
	virtual SPtr<LBMKernel> clone() { UB_THROW(UbException("SPtr<LBMKernel> clone() - belongs in the derived class")); };
	double getCalculationTime();

	void swapDistributions();

protected:
	void initDataSet();

	virtual LBMReal getThyxotropyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho) const { UB_THROW(UbException("LBMReal getThyxotropyCollFactor() - belongs in the derived class")); }

	LBMReal f[D3Q27System::ENDF + 1];

	UbTimer timer;

	LBMReal OxyyMxzz;
	
	CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF;
	CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF;
	CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsF;

	mu::value_type muX1, muX2, muX3;
	mu::value_type muDeltaT;
	mu::value_type muNu;
	LBMReal forcingX1;
	LBMReal forcingX2;
	LBMReal forcingX3;

	bool test;
};

#endif
