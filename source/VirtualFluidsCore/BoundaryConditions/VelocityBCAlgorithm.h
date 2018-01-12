#ifndef VelocityBoundaryCondition_h__
#define VelocityBoundaryCondition_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class VelocityBCAlgorithm : public BCAlgorithm
{
public:
   VelocityBCAlgorithm();
   ~VelocityBCAlgorithm();
   SPtr<BCAlgorithm> clone();
   void addDistributions(SPtr<DistributionArray3D> distributions);

   void applyBC() override;
private:
   //friend class boost::serialization::access;
   //template<class Archive>
   //void serialize(Archive & ar, const unsigned int version)
   //{
   //   ar & boost::serialization::base_object<BCAlgorithm>(*this);
   //}
};

#endif // VelocityBoundaryCondition_h__

