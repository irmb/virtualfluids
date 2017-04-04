//!  \file NonReflectingVelocityBCAlgorithm.h
//!  \brief Class implements velocity bc for non reflecting pressure bc.
//!  \author Konstantin Kutscher

#ifndef NonReflectingVelocityBCAlgorithm_h__
#define NonReflectingVelocityBCAlgorithm_h__

#include "BCAlgorithm.h"

class NonReflectingVelocityBCAlgorithm;
typedef boost::shared_ptr<NonReflectingVelocityBCAlgorithm> NonReflectingVelocityBCAlgorithmPtr;

//!  \brief Class implements velocity boundary condition for non reflecting pressure boundary condition

class NonReflectingVelocityBCAlgorithm : public BCAlgorithm
{
public:
   NonReflectingVelocityBCAlgorithm();
   ~NonReflectingVelocityBCAlgorithm();
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
#endif // NonReflectingVelocityBCAlgorithm_h__
