#ifndef NonEqDensityBCAlgorithm_h__
#define NonEqDensityBCAlgorithm_h__

#include "BCAlgorithm.h"

class NonEqDensityBCAlgorithm;
typedef boost::shared_ptr<NonEqDensityBCAlgorithm> NonEqDensityBCAlgorithmPtr;

class NonEqDensityBCAlgorithm : public BCAlgorithm
{
public:
   NonEqDensityBCAlgorithm();
   ~NonEqDensityBCAlgorithm();
   BCAlgorithmPtr clone();
   void addDistributions(DistributionArray3DPtr distributions);
   void applyBC();
protected:
   
private:
   //friend class boost::serialization::access;
   //template<class Archive>
   //void serialize(Archive & ar, const unsigned int version)
   //{
   //   ar & boost::serialization::base_object<BCAlgorithm>(*this);
   //}
};
#endif // NonEqDensityBCAlgorithm_h__
