#ifndef NonReflectingDensityBoundaryCondition_h__
#define NonReflectingDensityBoundaryCondition_h__

#include "BoundaryCondition.h"

class NonReflectingDensityBoundaryCondition;
typedef boost::shared_ptr<NonReflectingDensityBoundaryCondition> NonReflectingDensityBoundaryConditionPtr;

class NonReflectingDensityBoundaryCondition : public BoundaryCondition
{
public:
   NonReflectingDensityBoundaryCondition();
   ~NonReflectingDensityBoundaryCondition();
   BoundaryConditionPtr clone();
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
#endif // NonReflectingDensityBoundaryCondition_h__
