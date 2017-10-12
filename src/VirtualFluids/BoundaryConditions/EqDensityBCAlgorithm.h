#ifndef EqDensityBCAlgorithm_h__
#define EqDensityBCAlgorithm_h__

#include "BCAlgorithm.h"

class EqDensityBCAlgorithm;
typedef boost::shared_ptr<EqDensityBCAlgorithm> EqDensityBCAlgorithmPtr;

class EqDensityBCAlgorithm : public BCAlgorithm
{
public:
   EqDensityBCAlgorithm();
   ~EqDensityBCAlgorithm();
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
#endif // EqDensityBCAlgorithm_h__
