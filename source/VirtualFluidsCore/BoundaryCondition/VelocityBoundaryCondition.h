#ifndef VelocityBoundaryCondition_h__
#define VelocityBoundaryCondition_h__

#include "BoundaryCondition.h"

class VelocityBoundaryCondition;
typedef boost::shared_ptr<VelocityBoundaryCondition> VelocityBoundaryConditionPtr;

class VelocityBoundaryCondition : public BoundaryCondition
{
public:
   VelocityBoundaryCondition();
   ~VelocityBoundaryCondition();
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

#endif // VelocityBoundaryCondition_h__

