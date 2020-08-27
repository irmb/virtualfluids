#ifndef DensityAndThixotropyBCAlgorithm_h__
#define DensityAndThixotropyBCAlgorithm_h__

#include "BCAlgorithm.h"


class DensityAndThixotropyBCAlgorithm : public BCAlgorithm
{
public:
	DensityAndThixotropyBCAlgorithm();
	virtual ~DensityAndThixotropyBCAlgorithm();
	SPtr<BCAlgorithm> clone();
	void addDistributions(SPtr<DistributionArray3D> distributions);
	//void addDistributionsF(SPtr<DistributionArray3D> distributions);
	void addDistributionsH(SPtr<DistributionArray3D> distributions);
	void applyBC();
	void setLambdaBC(LBMReal lambda) { this->lambdaBC = lambda; }
	LBMReal getLambdaBC() { return this->lambdaBC; }
protected:
	SPtr<DistributionArray3D> distributionsH;
private:
	LBMReal lambdaBC;
};
#endif // DensityAndThixotropyBCAlgorithm_h__

