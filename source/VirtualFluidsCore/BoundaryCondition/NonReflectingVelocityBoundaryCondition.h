//!  \file NonReflectingVelocityBoundaryCondition.h
//!  \brief Class implements velocity bc for non reflecting pressure bc.
//!  \author Konstantin Kutscher

#ifndef NonReflectingVelocityBoundaryCondition_h__
#define NonReflectingVelocityBoundaryCondition_h__

#include "BoundaryCondition.h"

class NonReflectingVelocityBoundaryCondition;
typedef boost::shared_ptr<NonReflectingVelocityBoundaryCondition> NonReflectingVelocityBoundaryConditionPtr;

//!  \brief Class implements velocity boundary condition for non reflecting pressure boundary condition

class NonReflectingVelocityBoundaryCondition : public BoundaryCondition
{
public:
   NonReflectingVelocityBoundaryCondition();
   ~NonReflectingVelocityBoundaryCondition();
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
#endif // NonReflectingVelocityBoundaryCondition_h__
