#ifndef EqDensityBCAlgorithm_h__
#define EqDensityBCAlgorithm_h__

#include "BCAlgorithm.h"

class DistributionArray3D;

class EqDensityBCAlgorithm;
typedef std::shared_ptr<EqDensityBCAlgorithm> EqDensityBCAlgorithmPtr;

class EqDensityBCAlgorithm : public BCAlgorithm
{
public:
   EqDensityBCAlgorithm();
   ~EqDensityBCAlgorithm();
   BCAlgorithmPtr clone();
   void addDistributions(std::shared_ptr<DistributionArray3D> distributions);
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
