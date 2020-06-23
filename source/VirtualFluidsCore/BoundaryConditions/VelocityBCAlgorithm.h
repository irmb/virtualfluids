#ifndef VelocityBoundaryCondition_h__
#define VelocityBoundaryCondition_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class VelocityBCAlgorithm : public BCAlgorithm
{
public:
   VelocityBCAlgorithm();
   ~VelocityBCAlgorithm();
   SPtr<BCAlgorithm> clone();
   void addDistributions(SPtr<DistributionArray3D> distributions);
   void applyBC() override;
};

#endif // VelocityBoundaryCondition_h__

