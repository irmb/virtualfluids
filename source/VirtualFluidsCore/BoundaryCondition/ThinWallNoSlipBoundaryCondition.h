#ifndef ThinWallNoSlipBoundaryCondition_h__
#define ThinWallNoSlipBoundaryCondition_h__

#include "NoSlipBoundaryCondition.h"

class ThinWallNoSlipBoundaryCondition;
typedef boost::shared_ptr<ThinWallNoSlipBoundaryCondition> ThinWallNoSlipBoundaryConditionPtr;

class ThinWallNoSlipBoundaryCondition : public BoundaryCondition
{
public:
   ThinWallNoSlipBoundaryCondition();
   virtual ~ThinWallNoSlipBoundaryCondition();
   BoundaryConditionPtr clone();
   void addDistributions(DistributionArray3DPtr distributions);
   void setPass(int pass);
protected:
   void applyBC();
private:
   int pass;
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<BoundaryCondition>(*this);
   }
};
#endif // ThinWallNoSlipBoundaryCondition_h__
