#ifndef PowellEyringModelNoSlipBCAlgorithm_h__
#define PowellEyringModelNoSlipBCAlgorithm_h__

#include "ThixotropyNoSlipBCAlgorithm.h"
#include "Thixotropy.h"

class PowellEyringModelNoSlipBCAlgorithm : public ThixotropyNoSlipBCAlgorithm
{
public:
   PowellEyringModelNoSlipBCAlgorithm() 
   {
      BCAlgorithm::type = BCAlgorithm::PowellEyringModelNoSlipBCAlgorithm;
      BCAlgorithm::preCollision = true;
   }
   ~PowellEyringModelNoSlipBCAlgorithm() {}
   SPtr<BCAlgorithm> clone() override
   {
      SPtr<BCAlgorithm> bc(new PowellEyringModelNoSlipBCAlgorithm());
      return bc;
   }
protected:
   LBMReal getThyxotropyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho) const override
   {
      return Thixotropy::getHerschelBulkleyCollFactor(omegaInf, shearRate, drho);
   }
};
#endif // PowellEyringModelNoSlipBCAlgorithm_h__