#ifndef EqDensityBCAlgorithm_h__
#define EqDensityBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class EqDensityBCAlgorithm : public BCAlgorithm
{
public:
   EqDensityBCAlgorithm();
   ~EqDensityBCAlgorithm();
   SPtr<BCAlgorithm> clone();
   void addDistributions(SPtr<DistributionArray3D> distributions);
   void applyBC() override;
};
#endif // EqDensityBCAlgorithm_h__
