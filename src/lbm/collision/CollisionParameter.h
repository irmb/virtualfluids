#ifndef LBM_KERNEL_PARAMETER_H
#define LBM_KERNEL_PARAMETER_H

#include <basics/DataTypes.h>

namespace vf::lbm
{

struct CollisionParameter
{
    real distribution[27];
    real omega;
    real* quadricLimiter;
    real forceX;
    real forceY;
    real forceZ;
};

struct TurbulentViscosity
{
    real SGSconstant;
    real value { 0. };
};

struct MacroscopicValues
{
    real vx { 0. };
    real vy { 0. };
    real vz { 0. };
    real rho { 0. };
};

} // namespace vf::lbm

#endif
