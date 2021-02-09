#ifndef HerschelBulkleyModelNoSlipBCAlgorithm_h__
#define HerschelBulkleyModelNoSlipBCAlgorithm_h__

#include "ThixotropyNoSlipBCAlgorithm.h"
#include "Thixotropy.h"

class HerschelBulkleyModelNoSlipBCAlgorithm : public ThixotropyNoSlipBCAlgorithm
{
public:
   HerschelBulkleyModelNoSlipBCAlgorithm() 
   {
      BCAlgorithm::type = BCAlgorithm::HerschelBulkleyModelNoSlipBCAlgorithm;
      BCAlgorithm::preCollision = true;
   }
   ~HerschelBulkleyModelNoSlipBCAlgorithm() {}
   SPtr<BCAlgorithm> clone() override
   {
      SPtr<BCAlgorithm> bc(new HerschelBulkleyModelNoSlipBCAlgorithm());
      return bc;
   }
protected:
   LBMReal getThyxotropyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho) const override
   {
      return Thixotropy::getHerschelBulkleyCollFactor(omegaInf, shearRate, drho);
   }
};
#endif // HerschelBulkleyModelNoSlipBCAlgorithm_h__