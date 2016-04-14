#ifndef SimpleDensityBoundaryCondition_h__
#define SimpleDensityBoundaryCondition_h__

#include "BoundaryCondition.h"

class EqDensityBoundaryCondition;
typedef boost::shared_ptr<EqDensityBoundaryCondition> SimpleDensityBoundaryConditionPtr;

class EqDensityBoundaryCondition : public BoundaryCondition
{
public:
   EqDensityBoundaryCondition();
   ~EqDensityBoundaryCondition();
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
#endif // SimpleDensityBoundaryCondition_h__
