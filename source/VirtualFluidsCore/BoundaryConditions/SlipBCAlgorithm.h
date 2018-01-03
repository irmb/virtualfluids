#ifndef SlipBCAlgorithm_h__
#define SlipBCAlgorithm_h__

#include "BCAlgorithm.h"

class DistributionArray3D;

class SlipBCAlgorithm;
typedef std::shared_ptr<SlipBCAlgorithm> SlipBCAlgorithmPtr;

class SlipBCAlgorithm : public BCAlgorithm
{
public:
   SlipBCAlgorithm();
   virtual ~SlipBCAlgorithm();
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
#endif // SlipBCAlgorithm_h__
