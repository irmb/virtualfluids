#ifndef NonReflectingSlipBCAlgorithm_h__
#define NonReflectingSlipBCAlgorithm_h__

#include "BCAlgorithm.h"

class NonReflectingSlipBCAlgorithm;
typedef boost::shared_ptr<NonReflectingSlipBCAlgorithm> NonReflectingSlipBCAlgorithmPtr;

class NonReflectingSlipBCAlgorithm : public BCAlgorithm
{
public:
   NonReflectingSlipBCAlgorithm();
   virtual ~NonReflectingSlipBCAlgorithm();
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
#endif // NonReflectingSlipBCAlgorithm_h__

