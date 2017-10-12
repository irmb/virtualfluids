#ifndef VelocityBoundaryCondition_h__
#define VelocityBoundaryCondition_h__

#include "BCAlgorithm.h"

class VelocityBCAlgorithm;
typedef boost::shared_ptr<VelocityBCAlgorithm> VelocityBCAlgorithmPtr;

class VelocityBCAlgorithm : public BCAlgorithm
{
public:
   VelocityBCAlgorithm();
   ~VelocityBCAlgorithm();
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

#endif // VelocityBoundaryCondition_h__

