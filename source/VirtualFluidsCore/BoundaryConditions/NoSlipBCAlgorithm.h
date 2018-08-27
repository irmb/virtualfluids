#ifndef NoSlipBCAlgorithm_h__
#define NoSlipBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class NoSlipBCAlgorithm : public BCAlgorithm
{
public:
   NoSlipBCAlgorithm();
   virtual ~NoSlipBCAlgorithm();
   SPtr<BCAlgorithm> clone();
   void addDistributions(SPtr<DistributionArray3D> distributions);
   void applyBC() override;
private:
};
#endif // NoSlipBCAlgorithm_h__
