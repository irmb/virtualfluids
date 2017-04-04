#ifndef NonReflectingDensityBCAlgorithm_h__
#define NonReflectingDensityBCAlgorithm_h__

#include "BCAlgorithm.h"

class NonReflectingDensityBCAlgorithm;
typedef boost::shared_ptr<NonReflectingDensityBCAlgorithm> NonReflectingDensityBCAlgorithmPtr;

class NonReflectingDensityBCAlgorithm : public BCAlgorithm
{
public:
   NonReflectingDensityBCAlgorithm();
   ~NonReflectingDensityBCAlgorithm();
   BCAlgorithmPtr clone();
   void addDistributions(DistributionArray3DPtr distributions);
protected:
   void applyBC();
private:
   //friend class boost::serialization::access;
   //template<class Archive>
   //void serialize(Archive & ar, const unsigned int version)
   //{
   //   ar & boost::serialization::base_object<BCAlgorithm>(*this);
   //}
};
#endif // NonReflectingDensityBCAlgorithm_h__
