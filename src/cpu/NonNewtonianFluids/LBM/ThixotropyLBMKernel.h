#ifndef ThixotropyLBMKernel_H
#define ThixotropyLBMKernel_H

#include "LBMKernel.h"
#include "BCSet.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

class ThixotropyLBMKernel;

//! \brief    Cumulant + Fact. Central LBM kernel. 
//! \author Hussein
class ThixotropyLBMKernel : public LBMKernel
{
public:
	//! This option set relaxation parameter: NORMAL  
	enum Parameter { NORMAL, MAGIC };
public:
	ThixotropyLBMKernel();
	virtual ~ThixotropyLBMKernel(void);
	virtual void calculate(int step);
	virtual SPtr<LBMKernel> clone();
	real getCalculationTime();
 
	void setCollisionFactorF(real collFactor);
   void setCollisionFactorH(real collFactor);
   real getCollisionFactorF() const;
   real getCollisionFactorH() const;

	void setAlpha(real alpha);
	real getAlpha() const;

	void setTheta(real theta);
	real getTheta() const;

	void swapDistributions();

protected:
	virtual void initDataSet();
	real f[D3Q27System::ENDF + 1];

	UbTimer timer;

	real OxyyMxzz;
	Parameter parameter;

	CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF;
	CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF;
	CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsF;

	CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsH;
	CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsH;
	CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsH;

	mu::value_type muX1, muX2, muX3;
	mu::value_type muDeltaT;
	mu::value_type muNu;
	real forcingX1;
	real forcingX2;
	real forcingX3;

	real collFactorF;
   real collFactorH;

	real theta;
	real alpha;
};

#endif

