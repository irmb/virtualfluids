#ifndef DensityNEQBoundaryCondition_h__
#define DensityNEQBoundaryCondition_h__

#include "BoundaryCondition.h"

class NonEqDensityBoundaryCondition;
typedef boost::shared_ptr<NonEqDensityBoundaryCondition> NonEqDensityBoundaryConditionPtr;

class NonEqDensityBoundaryCondition : public BoundaryCondition
{
public:
   NonEqDensityBoundaryCondition();
   ~NonEqDensityBoundaryCondition();
   BoundaryConditionPtr clone();
   void addDistributions(DistributionArray3DPtr distributions);
protected:
   void applyBC();
private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<BoundaryCondition>(*this);
   }
};
#endif // DensityNEQBoundaryCondition_h__
