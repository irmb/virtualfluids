#ifndef NoSlipBCAlgorithm_h__
#define NoSlipBCAlgorithm_h__

#include "BCAlgorithm.h"

class DistributionArray3D;

class NoSlipBCAlgorithm;
typedef std::shared_ptr<NoSlipBCAlgorithm> NoSlipBCAlgorithmPtr;

class NoSlipBCAlgorithm : public BCAlgorithm
{
public:
   NoSlipBCAlgorithm();
   virtual ~NoSlipBCAlgorithm();
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
#endif // NoSlipBCAlgorithm_h__
