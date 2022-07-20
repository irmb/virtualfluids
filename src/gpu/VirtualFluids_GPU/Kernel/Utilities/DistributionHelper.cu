#include "DistributionHelper.cuh"

#include <cuda_runtime.h>


#include <lbm/constants/NumericConstants.h>
#include "lbm/constants/D3Q27.h"
using namespace vf::lbm::dir;

namespace vf::gpu
{

__device__ DistributionWrapper::DistributionWrapper(real *distributions, unsigned int size_Mat, bool isEvenTimestep,
                                                    uint k, uint *neighborX, uint *neighborY, uint *neighborZ)
    : distribution_references(getDistributionReferences27(distributions, size_Mat, isEvenTimestep)), k(k), kw(neighborX[k]), ks(neighborY[k]),
      kb(neighborZ[k]), ksw(neighborY[kw]), kbw(neighborZ[kw]), kbs(neighborZ[ks]), kbsw(neighborZ[ksw])
{
    read();
}

__device__ void DistributionWrapper::read()
{
    distribution.f[vf::lbm::dir::PZZ] = (distribution_references.f[E])[k];
    distribution.f[vf::lbm::dir::MZZ] = (distribution_references.f[W])[kw];
    distribution.f[vf::lbm::dir::ZPZ] = (distribution_references.f[N])[k];
    distribution.f[vf::lbm::dir::ZMZ] = (distribution_references.f[S])[ks];
    distribution.f[vf::lbm::dir::ZZP] = (distribution_references.f[T])[k];
    distribution.f[vf::lbm::dir::ZZM] = (distribution_references.f[B])[kb];
    distribution.f[vf::lbm::dir::PPZ] = (distribution_references.f[NE])[k];
    distribution.f[vf::lbm::dir::MMZ] = (distribution_references.f[SW])[ksw];
    distribution.f[vf::lbm::dir::PMZ] = (distribution_references.f[SE])[ks];
    distribution.f[vf::lbm::dir::MPZ] = (distribution_references.f[NW])[kw];
    distribution.f[vf::lbm::dir::PZP] = (distribution_references.f[TE])[k];
    distribution.f[vf::lbm::dir::MZM] = (distribution_references.f[BW])[kbw];
    distribution.f[vf::lbm::dir::PZM] = (distribution_references.f[BE])[kb];
    distribution.f[vf::lbm::dir::MZP] = (distribution_references.f[TW])[kw];
    distribution.f[vf::lbm::dir::ZPP] = (distribution_references.f[TN])[k];
    distribution.f[vf::lbm::dir::ZMM] = (distribution_references.f[BS])[kbs];
    distribution.f[vf::lbm::dir::ZPM] = (distribution_references.f[BN])[kb];
    distribution.f[vf::lbm::dir::ZMP] = (distribution_references.f[TS])[ks];
    distribution.f[vf::lbm::dir::PPP] = (distribution_references.f[TNE])[k];
    distribution.f[vf::lbm::dir::MPP] = (distribution_references.f[TNW])[kw];
    distribution.f[vf::lbm::dir::PMP] = (distribution_references.f[TSE])[ks];
    distribution.f[vf::lbm::dir::MMP] = (distribution_references.f[TSW])[ksw];
    distribution.f[vf::lbm::dir::PPM] = (distribution_references.f[BNE])[kb];
    distribution.f[vf::lbm::dir::MPM] = (distribution_references.f[BNW])[kbw];
    distribution.f[vf::lbm::dir::PMM] = (distribution_references.f[BSE])[kbs];
    distribution.f[vf::lbm::dir::MMM] = (distribution_references.f[BSW])[kbsw];
    distribution.f[vf::lbm::dir::ZZZ] = (distribution_references.f[REST])[k];
}

__device__ void DistributionWrapper::write()
{
    (distribution_references.f[E])[k]      = distribution.f[vf::lbm::dir::PZZ];
    (distribution_references.f[W])[kw]     = distribution.f[vf::lbm::dir::MZZ];
    (distribution_references.f[N])[k]      = distribution.f[vf::lbm::dir::ZPZ];
    (distribution_references.f[S])[ks]     = distribution.f[vf::lbm::dir::ZMZ];
    (distribution_references.f[T])[k]      = distribution.f[vf::lbm::dir::ZZP];
    (distribution_references.f[B])[kb]     = distribution.f[vf::lbm::dir::ZZM];
    (distribution_references.f[NE])[k]     = distribution.f[vf::lbm::dir::PPZ];
    (distribution_references.f[SW])[ksw]   = distribution.f[vf::lbm::dir::MMZ];
    (distribution_references.f[SE])[ks]    = distribution.f[vf::lbm::dir::PMZ];
    (distribution_references.f[NW])[kw]    = distribution.f[vf::lbm::dir::MPZ];
    (distribution_references.f[TE])[k]     = distribution.f[vf::lbm::dir::PZP];
    (distribution_references.f[BW])[kbw]   = distribution.f[vf::lbm::dir::MZM];
    (distribution_references.f[BE])[kb]    = distribution.f[vf::lbm::dir::PZM];
    (distribution_references.f[TW])[kw]    = distribution.f[vf::lbm::dir::MZP];
    (distribution_references.f[TN])[k]     = distribution.f[vf::lbm::dir::ZPP];
    (distribution_references.f[BS])[kbs]   = distribution.f[vf::lbm::dir::ZMM];
    (distribution_references.f[BN])[kb]    = distribution.f[vf::lbm::dir::ZPM];
    (distribution_references.f[TS])[ks]    = distribution.f[vf::lbm::dir::ZMP];
    (distribution_references.f[TNE])[k]    = distribution.f[vf::lbm::dir::PPP];
    (distribution_references.f[TNW])[kw]   = distribution.f[vf::lbm::dir::MPP];
    (distribution_references.f[TSE])[ks]   = distribution.f[vf::lbm::dir::PMP];
    (distribution_references.f[TSW])[ksw]  = distribution.f[vf::lbm::dir::MMP];
    (distribution_references.f[BNE])[kb]   = distribution.f[vf::lbm::dir::PPM];
    (distribution_references.f[BNW])[kbw]  = distribution.f[vf::lbm::dir::MPM];
    (distribution_references.f[BSE])[kbs]  = distribution.f[vf::lbm::dir::PMM];
    (distribution_references.f[BSW])[kbsw] = distribution.f[vf::lbm::dir::MMM];
    (distribution_references.f[REST])[k]   = distribution.f[vf::lbm::dir::ZZZ];
}

__device__ bool isValidFluidNode(uint nodeType)
{
    return (nodeType == GEO_FLUID || nodeType == GEO_PM_0 || nodeType == GEO_PM_1 || nodeType == GEO_PM_2);
}


}