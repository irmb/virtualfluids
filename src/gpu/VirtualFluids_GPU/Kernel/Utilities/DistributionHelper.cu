#include "DistributionHelper.cuh"

#include <cuda_runtime.h>


#include <lbm/constants/NumericConstants.h>
#include <lbm/constants/D3Q27.h>
#include "lbm/constants/D3Q27.h"
using namespace vf::lbm::dir;

namespace vf::gpu
{

__device__ __host__ DistributionReferences27 getDistributionReferences27(real *distributions, unsigned int size_Mat, bool isEvenTimestep)
{
    DistributionReferences27 distribution_references;

    if (isEvenTimestep) {
        distribution_references.f[E]    = &distributions[E * size_Mat];
        distribution_references.f[W]    = &distributions[W * size_Mat];
        distribution_references.f[N]    = &distributions[N * size_Mat];
        distribution_references.f[S]    = &distributions[S * size_Mat];
        distribution_references.f[T]    = &distributions[T * size_Mat];
        distribution_references.f[B]    = &distributions[B * size_Mat];
        distribution_references.f[NE]   = &distributions[NE * size_Mat];
        distribution_references.f[SW]   = &distributions[SW * size_Mat];
        distribution_references.f[SE]   = &distributions[SE * size_Mat];
        distribution_references.f[NW]   = &distributions[NW * size_Mat];
        distribution_references.f[TE]   = &distributions[TE * size_Mat];
        distribution_references.f[BW]   = &distributions[BW * size_Mat];
        distribution_references.f[BE]   = &distributions[BE * size_Mat];
        distribution_references.f[TW]   = &distributions[TW * size_Mat];
        distribution_references.f[TN]   = &distributions[TN * size_Mat];
        distribution_references.f[BS]   = &distributions[BS * size_Mat];
        distribution_references.f[BN]   = &distributions[BN * size_Mat];
        distribution_references.f[TS]   = &distributions[TS * size_Mat];
        distribution_references.f[REST] = &distributions[REST * size_Mat];
        distribution_references.f[TNE]  = &distributions[TNE * size_Mat];
        distribution_references.f[TSW]  = &distributions[TSW * size_Mat];
        distribution_references.f[TSE]  = &distributions[TSE * size_Mat];
        distribution_references.f[TNW]  = &distributions[TNW * size_Mat];
        distribution_references.f[BNE]  = &distributions[BNE * size_Mat];
        distribution_references.f[BSW]  = &distributions[BSW * size_Mat];
        distribution_references.f[BSE]  = &distributions[BSE * size_Mat];
        distribution_references.f[BNW]  = &distributions[BNW * size_Mat];
    } else {
        distribution_references.f[W]    = &distributions[E * size_Mat];
        distribution_references.f[E]    = &distributions[W * size_Mat];
        distribution_references.f[S]    = &distributions[N * size_Mat];
        distribution_references.f[N]    = &distributions[S * size_Mat];
        distribution_references.f[B]    = &distributions[T * size_Mat];
        distribution_references.f[T]    = &distributions[B * size_Mat];
        distribution_references.f[SW]   = &distributions[NE * size_Mat];
        distribution_references.f[NE]   = &distributions[SW * size_Mat];
        distribution_references.f[NW]   = &distributions[SE * size_Mat];
        distribution_references.f[SE]   = &distributions[NW * size_Mat];
        distribution_references.f[BW]   = &distributions[TE * size_Mat];
        distribution_references.f[TE]   = &distributions[BW * size_Mat];
        distribution_references.f[TW]   = &distributions[BE * size_Mat];
        distribution_references.f[BE]   = &distributions[TW * size_Mat];
        distribution_references.f[BS]   = &distributions[TN * size_Mat];
        distribution_references.f[TN]   = &distributions[BS * size_Mat];
        distribution_references.f[TS]   = &distributions[BN * size_Mat];
        distribution_references.f[BN]   = &distributions[TS * size_Mat];
        distribution_references.f[REST] = &distributions[REST * size_Mat];
        distribution_references.f[BSW]  = &distributions[TNE * size_Mat];
        distribution_references.f[BNE]  = &distributions[TSW * size_Mat];
        distribution_references.f[BNW]  = &distributions[TSE * size_Mat];
        distribution_references.f[BSE]  = &distributions[TNW * size_Mat];
        distribution_references.f[TSW]  = &distributions[BNE * size_Mat];
        distribution_references.f[TNE]  = &distributions[BSW * size_Mat];
        distribution_references.f[TNW]  = &distributions[BSE * size_Mat];
        distribution_references.f[TSE]  = &distributions[BNW * size_Mat];
    }
    return distribution_references;
}

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

__device__ unsigned int getNodeIndex()
{
    const unsigned x = threadIdx.x;
    const unsigned y = blockIdx.x;
    const unsigned z = blockIdx.y;

    const unsigned nx = blockDim.x;
    const unsigned ny = gridDim.x;

    return nx * (ny * z + y) + x;
}

__device__ bool isValidFluidNode(uint nodeType)
{
    return (nodeType == GEO_FLUID || nodeType == GEO_PM_0 || nodeType == GEO_PM_1 || nodeType == GEO_PM_2);
}


}