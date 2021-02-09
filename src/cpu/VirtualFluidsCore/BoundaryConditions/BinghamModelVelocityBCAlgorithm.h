#ifndef BinghamModelVelocityBCAlgorithm_h__
#define BinghamModelVelocityBCAlgorithm_h__

#include "RheologicalVelocityBCAlgorithm.h"
#include "Thixotropy.h"

class BinghamModelVelocityBCAlgorithm : public RheologicalVelocityBCAlgorithm
{
public:
   BinghamModelVelocityBCAlgorithm()
   {
      BCAlgorithm::type = BCAlgorithm::BinghamModelVelocityBCAlgorithm;
      BCAlgorithm::preCollision = false;
   }
   ~BinghamModelVelocityBCAlgorithm() {}
   SPtr<BCAlgorithm> clone() override
   {
      SPtr<BCAlgorithm> bc(new BinghamModelVelocityBCAlgorithm());
      return bc;
   }
protected:
   LBMReal getThyxotropyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho) const override 
   { 
      return Thixotropy::getBinghamCollFactor(omegaInf, shearRate, drho);
   }
};
#endif // BinghamModelVelocityBCAlgorithm_h__
