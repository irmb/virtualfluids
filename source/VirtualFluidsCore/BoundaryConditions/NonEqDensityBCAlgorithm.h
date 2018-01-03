#ifndef NonEqDensityBCAlgorithm_h__
#define NonEqDensityBCAlgorithm_h__

#include "BCAlgorithm.h"

class DistributionArray3D;

class NonEqDensityBCAlgorithm;
typedef std::shared_ptr<NonEqDensityBCAlgorithm> NonEqDensityBCAlgorithmPtr;

class NonEqDensityBCAlgorithm : public BCAlgorithm
{
public:
   NonEqDensityBCAlgorithm();
   ~NonEqDensityBCAlgorithm();
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
#endif // NonEqDensityBCAlgorithm_h__
