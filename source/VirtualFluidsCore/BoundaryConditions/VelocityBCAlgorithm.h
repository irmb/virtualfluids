#ifndef VelocityBoundaryCondition_h__
#define VelocityBoundaryCondition_h__

#include "BCAlgorithm.h"

class DistributionArray3D;

class VelocityBCAlgorithm;
typedef std::shared_ptr<VelocityBCAlgorithm> VelocityBCAlgorithmPtr;

class VelocityBCAlgorithm : public BCAlgorithm
{
public:
   VelocityBCAlgorithm();
   ~VelocityBCAlgorithm();
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

#endif // VelocityBoundaryCondition_h__

