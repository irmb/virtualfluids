#ifndef HighViscosityNoSlipBCAlgorithm_h__
#define HighViscosityNoSlipBCAlgorithm_h__

#include "BCAlgorithm.h"

class DistributionArray3D;

class HighViscosityNoSlipBCAlgorithm;
typedef std::shared_ptr<HighViscosityNoSlipBCAlgorithm> HighViscosityNoSlipBCAlgorithmPtr;

class HighViscosityNoSlipBCAlgorithm : public BCAlgorithm
{
public:
   HighViscosityNoSlipBCAlgorithm();
   ~HighViscosityNoSlipBCAlgorithm();
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
#endif // HighViscosityNoSlipBCAlgorithm_h__

