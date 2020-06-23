#ifndef NonEqDensityBCAlgorithm_h__
#define NonEqDensityBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class NonEqDensityBCAlgorithm : public BCAlgorithm
{
public:
   NonEqDensityBCAlgorithm();
   ~NonEqDensityBCAlgorithm();
   SPtr<BCAlgorithm> clone();
   void addDistributions(SPtr<DistributionArray3D> distributions);
   void applyBC() override;
};
#endif // NonEqDensityBCAlgorithm_h__
