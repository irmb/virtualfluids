#ifndef ThixotropyExpLBMKernel_H
#define ThixotropyExpLBMKernel_H

#include "LBMKernel.h"
#include "BCProcessor.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

class ThixotropyExpLBMKernel;

//! \brief    Cumulant + Fact. Central LBM kernel. 
//! \author Hussein
class ThixotropyExpLBMKernel : public LBMKernel
{
public:
	//! This option set relaxation parameter: NORMAL  
	enum Parameter { NORMAL, MAGIC };
public:
	ThixotropyExpLBMKernel();
	virtual ~ThixotropyExpLBMKernel(void);
	virtual void calculate(int step);
	virtual SPtr<LBMKernel> clone();
	double getCalculationTime();
 
	void setCollisionFactorF(double collFactor);
   void setCollisionFactorH(double collFactor);
   double getCollisionFactorF() const;
   double getCollisionFactorH() const;

	void setAlpha(double alpha);
	double getAlpha() const;

	void setTheta(double theta);
	double getTheta() const;

	void swapDistributions();

protected:
	virtual void initDataSet();
	LBMReal f[D3Q27System::ENDF + 1];

	UbTimer timer;

	LBMReal OxyyMxzz;
	Parameter parameter;

	CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF;
	CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF;
	CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsF;

	CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsH;
	CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsH;
	CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsH;

	mu::value_type muX1, muX2, muX3;
	mu::value_type muDeltaT;
	mu::value_type muNu;
	LBMReal forcingX1;
	LBMReal forcingX2;
	LBMReal forcingX3;

	LBMReal collFactorF;
   LBMReal collFactorH;

	LBMReal theta;
	LBMReal alpha;
};

#endif

