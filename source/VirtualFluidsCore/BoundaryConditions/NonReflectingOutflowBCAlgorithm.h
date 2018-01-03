#ifndef NonReflectingOutflowBCAlgorithm_h__
#define NonReflectingOutflowBCAlgorithm_h__

#include "BCAlgorithm.h"

class DistributionArray3D;

class NonReflectingOutflowBCAlgorithm;
typedef std::shared_ptr<NonReflectingOutflowBCAlgorithm> NonReflectingOutflowBCAlgorithmPtr;

class NonReflectingOutflowBCAlgorithm : public BCAlgorithm
{
public:
   NonReflectingOutflowBCAlgorithm();
   ~NonReflectingOutflowBCAlgorithm();
   BCAlgorithmPtr clone();
   void addDistributions(std::shared_ptr<DistributionArray3D> distributions);
   void applyBC() override;
private:
   int step;
   //friend class boost::serialization::access;
   //template<class Archive>
   //void serialize(Archive & ar, const unsigned int version)
   //{
   //   ar & boost::serialization::base_object<BCAlgorithm>(*this);
   //}
};
#endif // NonReflectingDensityBCAlgorithm_h__
