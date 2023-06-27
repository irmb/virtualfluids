#include "DistributionHelper.cuh"

#include <cuda_runtime.h>

#include "basics/constants/NumericConstants.h"
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
    distribution.f[vf::lbm::dir::DIR_P00] = (distribution_references.f[DIR_P00])[k];
    distribution.f[vf::lbm::dir::DIR_M00] = (distribution_references.f[DIR_M00])[kw];
    distribution.f[vf::lbm::dir::DIR_0P0] = (distribution_references.f[DIR_0P0])[k];
    distribution.f[vf::lbm::dir::DIR_0M0] = (distribution_references.f[DIR_0M0])[ks];
    distribution.f[vf::lbm::dir::DIR_00P] = (distribution_references.f[DIR_00P])[k];
    distribution.f[vf::lbm::dir::DIR_00M] = (distribution_references.f[DIR_00M])[kb];
    distribution.f[vf::lbm::dir::DIR_PP0] = (distribution_references.f[DIR_PP0])[k];
    distribution.f[vf::lbm::dir::DIR_MM0] = (distribution_references.f[DIR_MM0])[ksw];
    distribution.f[vf::lbm::dir::DIR_PM0] = (distribution_references.f[DIR_PM0])[ks];
    distribution.f[vf::lbm::dir::DIR_MP0] = (distribution_references.f[DIR_MP0])[kw];
    distribution.f[vf::lbm::dir::DIR_P0P] = (distribution_references.f[DIR_P0P])[k];
    distribution.f[vf::lbm::dir::DIR_M0M] = (distribution_references.f[DIR_M0M])[kbw];
    distribution.f[vf::lbm::dir::DIR_P0M] = (distribution_references.f[DIR_P0M])[kb];
    distribution.f[vf::lbm::dir::DIR_M0P] = (distribution_references.f[DIR_M0P])[kw];
    distribution.f[vf::lbm::dir::DIR_0PP] = (distribution_references.f[DIR_0PP])[k];
    distribution.f[vf::lbm::dir::DIR_0MM] = (distribution_references.f[DIR_0MM])[kbs];
    distribution.f[vf::lbm::dir::DIR_0PM] = (distribution_references.f[DIR_0PM])[kb];
    distribution.f[vf::lbm::dir::DIR_0MP] = (distribution_references.f[DIR_0MP])[ks];
    distribution.f[vf::lbm::dir::DIR_PPP] = (distribution_references.f[DIR_PPP])[k];
    distribution.f[vf::lbm::dir::DIR_MPP] = (distribution_references.f[DIR_MPP])[kw];
    distribution.f[vf::lbm::dir::DIR_PMP] = (distribution_references.f[DIR_PMP])[ks];
    distribution.f[vf::lbm::dir::DIR_MMP] = (distribution_references.f[DIR_MMP])[ksw];
    distribution.f[vf::lbm::dir::DIR_PPM] = (distribution_references.f[DIR_PPM])[kb];
    distribution.f[vf::lbm::dir::DIR_MPM] = (distribution_references.f[DIR_MPM])[kbw];
    distribution.f[vf::lbm::dir::DIR_PMM] = (distribution_references.f[DIR_PMM])[kbs];
    distribution.f[vf::lbm::dir::DIR_MMM] = (distribution_references.f[DIR_MMM])[kbsw];
    distribution.f[vf::lbm::dir::DIR_000] = (distribution_references.f[DIR_000])[k];
}

__device__ void DistributionWrapper::write()
{
    (distribution_references.f[DIR_P00])[k]      = distribution.f[vf::lbm::dir::PZZ];
    (distribution_references.f[DIR_M00])[kw]     = distribution.f[vf::lbm::dir::MZZ];
    (distribution_references.f[DIR_0P0])[k]      = distribution.f[vf::lbm::dir::ZPZ];
    (distribution_references.f[DIR_0M0])[ks]     = distribution.f[vf::lbm::dir::ZMZ];
    (distribution_references.f[DIR_00P])[k]      = distribution.f[vf::lbm::dir::ZZP];
    (distribution_references.f[DIR_00M])[kb]     = distribution.f[vf::lbm::dir::ZZM];
    (distribution_references.f[DIR_PP0])[k]     = distribution.f[vf::lbm::dir::PPZ];
    (distribution_references.f[DIR_MM0])[ksw]   = distribution.f[vf::lbm::dir::MMZ];
    (distribution_references.f[DIR_PM0])[ks]    = distribution.f[vf::lbm::dir::PMZ];
    (distribution_references.f[DIR_MP0])[kw]    = distribution.f[vf::lbm::dir::MPZ];
    (distribution_references.f[DIR_P0P])[k]     = distribution.f[vf::lbm::dir::PZP];
    (distribution_references.f[DIR_M0M])[kbw]   = distribution.f[vf::lbm::dir::MZM];
    (distribution_references.f[DIR_P0M])[kb]    = distribution.f[vf::lbm::dir::PZM];
    (distribution_references.f[DIR_M0P])[kw]    = distribution.f[vf::lbm::dir::MZP];
    (distribution_references.f[DIR_0PP])[k]     = distribution.f[vf::lbm::dir::ZPP];
    (distribution_references.f[DIR_0MM])[kbs]   = distribution.f[vf::lbm::dir::ZMM];
    (distribution_references.f[DIR_0PM])[kb]    = distribution.f[vf::lbm::dir::ZPM];
    (distribution_references.f[DIR_0MP])[ks]    = distribution.f[vf::lbm::dir::ZMP];
    (distribution_references.f[DIR_PPP])[k]    = distribution.f[vf::lbm::dir::PPP];
    (distribution_references.f[DIR_MPP])[kw]   = distribution.f[vf::lbm::dir::MPP];
    (distribution_references.f[DIR_PMP])[ks]   = distribution.f[vf::lbm::dir::PMP];
    (distribution_references.f[DIR_MMP])[ksw]  = distribution.f[vf::lbm::dir::MMP];
    (distribution_references.f[DIR_PPM])[kb]   = distribution.f[vf::lbm::dir::PPM];
    (distribution_references.f[DIR_MPM])[kbw]  = distribution.f[vf::lbm::dir::MPM];
    (distribution_references.f[DIR_PMM])[kbs]  = distribution.f[vf::lbm::dir::PMM];
    (distribution_references.f[DIR_MMM])[kbsw] = distribution.f[vf::lbm::dir::MMM];
    (distribution_references.f[DIR_000])[k]   = distribution.f[vf::lbm::dir::ZZZ];
}

}