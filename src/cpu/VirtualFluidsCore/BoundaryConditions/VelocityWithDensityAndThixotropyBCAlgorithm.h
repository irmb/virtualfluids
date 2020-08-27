//!  \file VelocityWithDensityAndThixotropyBCAlgorithm.h
//!  \brief Class implements velocity bc for non reflecting pressure bc.
//!  \author Konstantin Kutscher

#ifndef VelocityWithDensityAndThixotropyBCAlgorithm_h__
#define VelocityWithDensityAndThixotropyBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

//!  \brief Class implements Dirichlet boundary condition for velocity. Set density in system. It is used together with non reflecting outflow.  

class VelocityWithDensityAndThixotropyBCAlgorithm : public BCAlgorithm
{
public:
   VelocityWithDensityAndThixotropyBCAlgorithm();
   ~VelocityWithDensityAndThixotropyBCAlgorithm();
   SPtr<BCAlgorithm> clone();
   void addDistributions(SPtr<DistributionArray3D> distributions);
   void addDistributionsH(SPtr<DistributionArray3D> distributions);
   void applyBC();
   void setLambdaBC(LBMReal lambda) { this->lambdaBC = lambda; }
   LBMReal getLambdaBC() { return this->lambdaBC; }
protected:
   SPtr<DistributionArray3D> distributionsH;
private:
   LBMReal lambdaBC;
};
#endif // VelocityWithDensityAndThixotropyBCAlgorithm_h__
