#ifndef NoSlipAndThixotropyBCAlgorithm_h__
#define NoSlipAndThixotropyBCAlgorithm_h__

#include "BCAlgorithm.h"

class NoSlipAndThixotropyBCAlgorithm : public BCAlgorithm
{
public:
	NoSlipAndThixotropyBCAlgorithm();
	virtual ~NoSlipAndThixotropyBCAlgorithm();
	SPtr<BCAlgorithm> clone();
	void addDistributions(SPtr<DistributionArray3D> distributions);
	//void addDistributionsF(DistributionArray3DPtr distributions);
	void addDistributionsH(SPtr<DistributionArray3D> distributions);
	void applyBC();
protected:
	SPtr<DistributionArray3D> distributionsH;
private:

};
#endif // NoSlipAndThixotropyBCAlgorithm_h__

