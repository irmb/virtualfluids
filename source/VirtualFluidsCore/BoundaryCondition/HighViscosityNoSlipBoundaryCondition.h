#ifndef HighViscosityNoSlipBoundaryCondition_h__
#define HighViscosityNoSlipBoundaryCondition_h__

#include "BoundaryCondition.h"

class NoSlipBoundaryCondition;
typedef boost::shared_ptr<NoSlipBoundaryCondition> NoSlipBoundaryConditionPtr;

class HighViscosityNoSlipBoundaryCondition : public BoundaryCondition
{
public:
   HighViscosityNoSlipBoundaryCondition();
   ~HighViscosityNoSlipBoundaryCondition();
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
#endif // HighViscosityNoSlipBoundaryCondition_h__

