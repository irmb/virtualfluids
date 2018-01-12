#ifndef NoSlipBCAlgorithm_h__
#define NoSlipBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class NoSlipBCAlgorithm : public BCAlgorithm
{
public:
   NoSlipBCAlgorithm();
   virtual ~NoSlipBCAlgorithm();
   SPtr<BCAlgorithm> clone();
   void addDistributions(SPtr<DistributionArray3D> distributions);
   void applyBC() override;
private:
   //friend class boost::serialization::access;
   //template<class Archive>
   //void serialize(Archive & ar, const unsigned int version)
   //{
   //   ar & boost::serialization::base_object<BCAlgorithm>(*this);
   //}
};
#endif // NoSlipBCAlgorithm_h__
