#ifndef VelocityAndThixotropyBCAlgorithm_h__
#define VelocityAndThixotropyBCAlgorithm_h__

#include "BCAlgorithm.h"


class VelocityAndThixotropyBCAlgorithm : public BCAlgorithm
{
public:
	VelocityAndThixotropyBCAlgorithm();
	virtual ~VelocityAndThixotropyBCAlgorithm();
	SPtr<BCAlgorithm> clone();
	void addDistributions(SPtr<DistributionArray3D> distributions);
	void addDistributionsH(SPtr<DistributionArray3D> distributions);
	void applyBC();
	void setLambdaBC(LBMReal lambda) { this->lambdaBC = lambda; }
	LBMReal getLambdaBC() { return this->lambdaBC; }
protected:
	SPtr<DistributionArray3D> distributionsH;
private:
	LBMReal lambdaBC;
};
#endif // VelocityAndThixotropyBCAlgorithm_h__

