#ifndef NonReflectingOutflowBCAlgorithm_h__
#define NonReflectingOutflowBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class NonReflectingOutflowBCAlgorithm : public BCAlgorithm
{
public:
   NonReflectingOutflowBCAlgorithm();
   ~NonReflectingOutflowBCAlgorithm();
   SPtr<BCAlgorithm> clone();
   void addDistributions(SPtr<DistributionArray3D> distributions);
   void applyBC() override;
};
#endif // NonReflectingDensityBCAlgorithm_h__
