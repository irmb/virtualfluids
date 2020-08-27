#ifndef BinghamModelNoSlipBCAlgorithm_h__
#define BinghamModelNoSlipBCAlgorithm_h__

#include "ThixotropyNoSlipBCAlgorithm.h"
#include "Thixotropy.h"

class BinghamModelNoSlipBCAlgorithm : public ThixotropyNoSlipBCAlgorithm
{
public:
   BinghamModelNoSlipBCAlgorithm()
   {
      BCAlgorithm::type = BCAlgorithm::BinghamModelNoSlipBCAlgorithm;
      BCAlgorithm::preCollision = true;
   }
   ~BinghamModelNoSlipBCAlgorithm() {}
   SPtr<BCAlgorithm> clone() override
   {
      SPtr<BCAlgorithm> bc(new BinghamModelNoSlipBCAlgorithm());
      return bc;
   }
protected:
   LBMReal getThyxotropyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho) const override 
   { 
      return Thixotropy::getBinghamCollFactor(omegaInf, shearRate, drho);
   }
};
#endif // BinghamModelNoSlipBCAlgorithm_h__
