#ifndef SlipBCAlgorithm_h__
#define SlipBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class SlipBCAlgorithm : public BCAlgorithm
{
public:
   SlipBCAlgorithm();
   virtual ~SlipBCAlgorithm();
   SPtr<BCAlgorithm> clone();
   void addDistributions(SPtr<DistributionArray3D> distributions);
   void applyBC() override;
};
#endif // SlipBCAlgorithm_h__
