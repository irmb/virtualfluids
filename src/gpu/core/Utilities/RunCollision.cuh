#ifndef GPU_RUNCOLLISIONKERNEL_CUH
#define GPU_RUNCOLLISIONKERNEL_CUH

#include <cuda_runtime.h>

#include <basics/DataTypes.h>

#include <lbm/collision/CollisionParameter.h>
#include <lbm/collision/TurbulentViscosity.h>

#include "Utilities/KernelUtilities.h"

namespace vf::gpu
{

struct GPUCollisionParameter
{
    real omega;
    unsigned int* neighborX;
    unsigned int* neighborY;
    unsigned int* neighborZ;
    real* distributions;
    real* rho;
    real* vx;
    real* vy;
    real* vz;
    real* turbulentViscosity;
    real SGSconstant;
    int numberOfLBnodes;
    real forceFactor;
    real* forces;
    real* bodyForceX;
    real* bodyForceY;
    real* bodyForceZ;
    real* quadricLimiters;
    bool isEvenTimestep;
    const uint* fluidNodeIndices;
    uint numberOfFluidNodes;
};

template <typename CollisionFunctor, vf::lbm::TurbulenceModel turbulenceModel, bool writeMacroscopicVariables, bool applyBodyForce>
__global__ void runCollision(CollisionFunctor collision, GPUCollisionParameter collisionParameter)
{
    const unsigned nodeIndex = getNodeIndex();

    if (nodeIndex >= collisionParameter.numberOfFluidNodes)
        return;

    const unsigned k_000 = collisionParameter.fluidNodeIndices[nodeIndex];

    vf::lbm::CollisionParameter para;
    para.omega = collisionParameter.omega;
    para.quadricLimiter = collisionParameter.quadricLimiters;

    if (applyBodyForce) {
        para.forceX = (collisionParameter.forces[0] + collisionParameter.bodyForceX[k_000]) * c1o2 * collisionParameter.forceFactor;
        para.forceY = (collisionParameter.forces[1] + collisionParameter.bodyForceY[k_000]) * c1o2 * collisionParameter.forceFactor;
        para.forceZ = (collisionParameter.forces[2] + collisionParameter.bodyForceZ[k_000]) * c1o2 * collisionParameter.forceFactor;

        // Reset body force. To be used when not using round-off correction.
        collisionParameter.bodyForceX[k_000] = c0o1;
        collisionParameter.bodyForceY[k_000] = c0o1;
        collisionParameter.bodyForceZ[k_000] = c0o1;

        ////////////////////////////////////////////////////////////////////////////////////
        //!> Round-off correction
        //!
        //!> Similar to Kahan summation algorithm (https://en.wikipedia.org/wiki/Kahan_summation_algorithm)
        //!> Essentially computes the round-off error of the applied force and adds it in the next time step as a compensation.
        //!> Seems to be necesseary at very high Re boundary layers, where the forcing and velocity can
        //!> differ by several orders of magnitude.
        //!> \note 16/05/2022: Testing, still ongoing!
        //!
        // bodyForceX[k_000] = (acc_x-(vvx-vx))*factor*c2o1;
        // bodyForceY[k_000] = (acc_y-(vvy-vy))*factor*c2o1;
        // bodyForceZ[k_000] = (acc_z-(vvz-vz))*factor*c2o1;
    } else {
        para.forceX = collisionParameter.forces[0] * c1o2 * collisionParameter.forceFactor;
        para.forceY = collisionParameter.forces[1] * c1o2 * collisionParameter.forceFactor;
        para.forceZ = collisionParameter.forces[2] * c1o2 * collisionParameter.forceFactor;
    }

    vf::lbm::TurbulentViscosity turbulentViscosity;
    if (turbulenceModel != vf::lbm::TurbulenceModel::None) {
        turbulentViscosity.value = collisionParameter.turbulentViscosity[k_000];
        turbulentViscosity.SGSconstant = collisionParameter.SGSconstant;
    }

    Distributions27 dist;
    getPointersToDistributions(dist, collisionParameter.distributions, collisionParameter.numberOfLBnodes,
                               collisionParameter.isEvenTimestep);

    ListIndices listIndices(k_000, collisionParameter.neighborX, collisionParameter.neighborY, collisionParameter.neighborZ);

    getPreCollisionDistribution(para.distribution, dist, listIndices);

    vf::lbm::MacroscopicValues macroscopicValues;
    collision(para, macroscopicValues, turbulentViscosity);

    if (writeMacroscopicVariables || turbulenceModel == vf::lbm::TurbulenceModel::AMD) {
        collisionParameter.vx[k_000] = macroscopicValues.vx;
        collisionParameter.vy[k_000] = macroscopicValues.vy;
        collisionParameter.vz[k_000] = macroscopicValues.vz;
        collisionParameter.rho[k_000] = macroscopicValues.rho;
    }
    if (turbulenceModel != vf::lbm::TurbulenceModel::None)
        collisionParameter.turbulentViscosity[k_000] = turbulentViscosity.value;

    setPostCollisionDistribution(dist, listIndices, para.distribution);
}

} // namespace vf::gpu

#endif
