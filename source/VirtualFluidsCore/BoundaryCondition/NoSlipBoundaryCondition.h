#ifndef NoSlipBoundaryCondition_h__
#define NoSlipBoundaryCondition_h__

#include "BoundaryCondition.h"

class NoSlipBoundaryCondition;
typedef boost::shared_ptr<NoSlipBoundaryCondition> NoSlipBoundaryConditionPtr;

class NoSlipBoundaryCondition : public BoundaryCondition
{
public:
   NoSlipBoundaryCondition();
   virtual ~NoSlipBoundaryCondition();
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
#endif // NoSlipBoundaryCondition_h__
