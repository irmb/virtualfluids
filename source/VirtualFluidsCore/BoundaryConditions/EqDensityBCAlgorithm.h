#ifndef EqDensityBCAlgorithm_h__
#define EqDensityBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class EqDensityBCAlgorithm : public BCAlgorithm
{
public:
   EqDensityBCAlgorithm();
   ~EqDensityBCAlgorithm();
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
#endif // EqDensityBCAlgorithm_h__
