//!  \file NonReflectingVelocityBCAlgorithm.h
//!  \brief Class implements velocity bc for non reflecting pressure bc.
//!  \author Konstantin Kutscher

#ifndef VelocityWithDensityBCAlgorithm_h__
#define VelocityWithDensityBCAlgorithm_h__

#include "BCAlgorithm.h"

class VelocityWithDensityBCAlgorithm;
typedef boost::shared_ptr<VelocityWithDensityBCAlgorithm> VelocityWithDensityBCAlgorithmPtr;

//!  \brief Class implements Dirichlet boundary condition for velocity. Set density in system. It is used together with non reflecting outflow.  

class VelocityWithDensityBCAlgorithm : public BCAlgorithm
{
public:
   VelocityWithDensityBCAlgorithm();
   ~VelocityWithDensityBCAlgorithm();
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
