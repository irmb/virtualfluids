#ifndef SlipBoundaryCondition_h__
#define SlipBoundaryCondition_h__

#include "BoundaryCondition.h"

class SlipBoundaryCondition;
typedef boost::shared_ptr<SlipBoundaryCondition> SlipBoundaryConditionPtr;

class SlipBoundaryCondition : public BoundaryCondition
{
public:
   SlipBoundaryCondition();
   virtual ~SlipBoundaryCondition();
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
#endif // SlipBoundaryCondition_h__
