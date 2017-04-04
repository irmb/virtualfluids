#ifndef SlipBCAlgorithm_h__
#define SlipBCAlgorithm_h__

#include "BCAlgorithm.h"

class SlipBCAlgorithm;
typedef boost::shared_ptr<SlipBCAlgorithm> SlipBCAlgorithmPtr;

class SlipBCAlgorithm : public BCAlgorithm
{
public:
   SlipBCAlgorithm();
   virtual ~SlipBCAlgorithm();
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
#endif // SlipBCAlgorithm_h__
