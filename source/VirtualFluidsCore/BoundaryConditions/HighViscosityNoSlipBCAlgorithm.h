#ifndef HighViscosityNoSlipBCAlgorithm_h__
#define HighViscosityNoSlipBCAlgorithm_h__

#include "BCAlgorithm.h"

class HighViscosityNoSlipBCAlgorithm;
typedef boost::shared_ptr<HighViscosityNoSlipBCAlgorithm> HighViscosityNoSlipBCAlgorithmPtr;

class HighViscosityNoSlipBCAlgorithm : public BCAlgorithm
{
public:
   HighViscosityNoSlipBCAlgorithm();
   ~HighViscosityNoSlipBCAlgorithm();
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
#endif // HighViscosityNoSlipBCAlgorithm_h__

