#ifndef RheologicalVelocityBCAlgorithm_h__
#define RheologicalVelocityBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class RheologicalVelocityBCAlgorithm : public BCAlgorithm
{
public:
   RheologicalVelocityBCAlgorithm();
   ~RheologicalVelocityBCAlgorithm();
   virtual SPtr<BCAlgorithm> clone() override { UB_THROW(UbException("LBMReal clone() - belongs in the derived class")); }
   void addDistributions(SPtr<DistributionArray3D> distributions) override;
   void applyBC() override;
protected:
   virtual LBMReal getThyxotropyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho) const = 0; // { UB_THROW(UbException("LBMReal getThyxotropyCollFactor() - belongs in the derived class")); }
};

#endif // VelocityBoundaryCondition_h__

