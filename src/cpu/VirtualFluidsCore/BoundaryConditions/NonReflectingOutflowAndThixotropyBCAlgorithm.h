#ifndef NonReflectingOutflowAndThixotropyBCAlgorithm_h__
#define NonReflectingOutflowAndThixotropyBCAlgorithm_h__

#include "BCAlgorithm.h"


class NonReflectingOutflowAndThixotropyBCAlgorithm : public BCAlgorithm
{
public:
	NonReflectingOutflowAndThixotropyBCAlgorithm();
	virtual ~NonReflectingOutflowAndThixotropyBCAlgorithm();
	SPtr<BCAlgorithm> clone();
	void addDistributions(SPtr<DistributionArray3D> distributions);
	//void addDistributionsF(SPtr<DistributionArray3D> distributions);
	void addDistributionsH(SPtr<DistributionArray3D> distributions);
	void applyBC();
	//void setLambdaBC(LBMReal lambda) { this->lambdaBC = lambda; }
	//LBMReal getLambdaBC() { return this->lambdaBC; }
protected:
	SPtr<DistributionArray3D> distributionsH;
private:
	//LBMReal lambdaBC;
};
#endif // NonReflectingOutflowAndThixotropyBCAlgorithm_h__

