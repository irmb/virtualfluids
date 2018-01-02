#ifndef NoSlipBCAlgorithm_h__
#define NoSlipBCAlgorithm_h__

#include "BCAlgorithm.h"

class NoSlipBCAlgorithm;
typedef boost::shared_ptr<NoSlipBCAlgorithm> NoSlipBCAlgorithmPtr;

class NoSlipBCAlgorithm : public BCAlgorithm
{
public:
   NoSlipBCAlgorithm();
   virtual ~NoSlipBCAlgorithm();
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
#endif // NoSlipBCAlgorithm_h__
