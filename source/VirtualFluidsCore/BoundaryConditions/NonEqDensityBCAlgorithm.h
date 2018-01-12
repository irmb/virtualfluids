#ifndef NonEqDensityBCAlgorithm_h__
#define NonEqDensityBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class NonEqDensityBCAlgorithm : public BCAlgorithm
{
public:
   NonEqDensityBCAlgorithm();
   ~NonEqDensityBCAlgorithm();
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
#endif // NonEqDensityBCAlgorithm_h__
