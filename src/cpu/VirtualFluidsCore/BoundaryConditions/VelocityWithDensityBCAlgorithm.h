//!  \file VelocityWithDensityBCAlgorithm.h
//!  \brief Class implements velocity bc for non reflecting pressure bc.
//!  \author Konstantin Kutscher

#ifndef VelocityWithDensityBCAlgorithm_h__
#define VelocityWithDensityBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

//!  \brief Class implements Dirichlet boundary condition for velocity. Set density in system. It is used together with
//!  non reflecting outflow.

class VelocityWithDensityBCAlgorithm : public BCAlgorithm
{
public:
    VelocityWithDensityBCAlgorithm();
    ~VelocityWithDensityBCAlgorithm() override;
    SPtr<BCAlgorithm> clone() override;
    void addDistributions(SPtr<DistributionArray3D> distributions) override;
    void applyBC() override;
};
#endif // VelocityWithDensityBCAlgorithm_h__
