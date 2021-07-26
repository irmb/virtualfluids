#include "DistributionHelper.cuh"

#include <cuda_runtime.h>

#include "LBM/D3Q27.h"

#include <lbm/constants/NumericConstants.h>
#include <lbm/constants/D3Q27.h>

namespace vf
{
namespace gpu
{

__device__ __host__ DistributionReferences27 getDistributionReferences27(real *distributions, unsigned int size_Mat, bool isEvenTimestep)
{
    DistributionReferences27 distribution_references;

    if (isEvenTimestep) {
        distribution_references.f[dirE]    = &distributions[dirE * size_Mat];
        distribution_references.f[dirW]    = &distributions[dirW * size_Mat];
        distribution_references.f[dirN]    = &distributions[dirN * size_Mat];
        distribution_references.f[dirS]    = &distributions[dirS * size_Mat];
        distribution_references.f[dirT]    = &distributions[dirT * size_Mat];
        distribution_references.f[dirB]    = &distributions[dirB * size_Mat];
        distribution_references.f[dirNE]   = &distributions[dirNE * size_Mat];
        distribution_references.f[dirSW]   = &distributions[dirSW * size_Mat];
        distribution_references.f[dirSE]   = &distributions[dirSE * size_Mat];
        distribution_references.f[dirNW]   = &distributions[dirNW * size_Mat];
        distribution_references.f[dirTE]   = &distributions[dirTE * size_Mat];
        distribution_references.f[dirBW]   = &distributions[dirBW * size_Mat];
        distribution_references.f[dirBE]   = &distributions[dirBE * size_Mat];
        distribution_references.f[dirTW]   = &distributions[dirTW * size_Mat];
        distribution_references.f[dirTN]   = &distributions[dirTN * size_Mat];
        distribution_references.f[dirBS]   = &distributions[dirBS * size_Mat];
        distribution_references.f[dirBN]   = &distributions[dirBN * size_Mat];
        distribution_references.f[dirTS]   = &distributions[dirTS * size_Mat];
        distribution_references.f[dirREST] = &distributions[dirREST * size_Mat];
        distribution_references.f[dirTNE]  = &distributions[dirTNE * size_Mat];
        distribution_references.f[dirTSW]  = &distributions[dirTSW * size_Mat];
        distribution_references.f[dirTSE]  = &distributions[dirTSE * size_Mat];
        distribution_references.f[dirTNW]  = &distributions[dirTNW * size_Mat];
        distribution_references.f[dirBNE]  = &distributions[dirBNE * size_Mat];
        distribution_references.f[dirBSW]  = &distributions[dirBSW * size_Mat];
        distribution_references.f[dirBSE]  = &distributions[dirBSE * size_Mat];
        distribution_references.f[dirBNW]  = &distributions[dirBNW * size_Mat];
    } else {
        distribution_references.f[dirW]    = &distributions[dirE * size_Mat];
        distribution_references.f[dirE]    = &distributions[dirW * size_Mat];
        distribution_references.f[dirS]    = &distributions[dirN * size_Mat];
        distribution_references.f[dirN]    = &distributions[dirS * size_Mat];
        distribution_references.f[dirB]    = &distributions[dirT * size_Mat];
        distribution_references.f[dirT]    = &distributions[dirB * size_Mat];
        distribution_references.f[dirSW]   = &distributions[dirNE * size_Mat];
        distribution_references.f[dirNE]   = &distributions[dirSW * size_Mat];
        distribution_references.f[dirNW]   = &distributions[dirSE * size_Mat];
        distribution_references.f[dirSE]   = &distributions[dirNW * size_Mat];
        distribution_references.f[dirBW]   = &distributions[dirTE * size_Mat];
        distribution_references.f[dirTE]   = &distributions[dirBW * size_Mat];
        distribution_references.f[dirTW]   = &distributions[dirBE * size_Mat];
        distribution_references.f[dirBE]   = &distributions[dirTW * size_Mat];
        distribution_references.f[dirBS]   = &distributions[dirTN * size_Mat];
        distribution_references.f[dirTN]   = &distributions[dirBS * size_Mat];
        distribution_references.f[dirTS]   = &distributions[dirBN * size_Mat];
        distribution_references.f[dirBN]   = &distributions[dirTS * size_Mat];
        distribution_references.f[dirREST] = &distributions[dirREST * size_Mat];
        distribution_references.f[dirBSW]  = &distributions[dirTNE * size_Mat];
        distribution_references.f[dirBNE]  = &distributions[dirTSW * size_Mat];
        distribution_references.f[dirBNW]  = &distributions[dirTSE * size_Mat];
        distribution_references.f[dirBSE]  = &distributions[dirTNW * size_Mat];
        distribution_references.f[dirTSW]  = &distributions[dirBNE * size_Mat];
        distribution_references.f[dirTNE]  = &distributions[dirBSW * size_Mat];
        distribution_references.f[dirTNW]  = &distributions[dirBSE * size_Mat];
        distribution_references.f[dirTSE]  = &distributions[dirBNW * size_Mat];
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
    distribution.f[vf::lbm::dir::PZZ] = (distribution_references.f[dirE])[k];
    distribution.f[vf::lbm::dir::MZZ] = (distribution_references.f[dirW])[kw];
    distribution.f[vf::lbm::dir::ZPZ] = (distribution_references.f[dirN])[k];
    distribution.f[vf::lbm::dir::ZMZ] = (distribution_references.f[dirS])[ks];
    distribution.f[vf::lbm::dir::ZZP] = (distribution_references.f[dirT])[k];
    distribution.f[vf::lbm::dir::ZZM] = (distribution_references.f[dirB])[kb];
    distribution.f[vf::lbm::dir::PPZ] = (distribution_references.f[dirNE])[k];
    distribution.f[vf::lbm::dir::MMZ] = (distribution_references.f[dirSW])[ksw];
    distribution.f[vf::lbm::dir::PMZ] = (distribution_references.f[dirSE])[ks];
    distribution.f[vf::lbm::dir::MPZ] = (distribution_references.f[dirNW])[kw];
    distribution.f[vf::lbm::dir::PZP] = (distribution_references.f[dirTE])[k];
    distribution.f[vf::lbm::dir::MZM] = (distribution_references.f[dirBW])[kbw];
    distribution.f[vf::lbm::dir::PZM] = (distribution_references.f[dirBE])[kb];
    distribution.f[vf::lbm::dir::MZP] = (distribution_references.f[dirTW])[kw];
    distribution.f[vf::lbm::dir::ZPP] = (distribution_references.f[dirTN])[k];
    distribution.f[vf::lbm::dir::ZMM] = (distribution_references.f[dirBS])[kbs];
    distribution.f[vf::lbm::dir::ZPM] = (distribution_references.f[dirBN])[kb];
    distribution.f[vf::lbm::dir::ZMP] = (distribution_references.f[dirTS])[ks];
    distribution.f[vf::lbm::dir::PPP] = (distribution_references.f[dirTNE])[k];
    distribution.f[vf::lbm::dir::MPP] = (distribution_references.f[dirTNW])[kw];
    distribution.f[vf::lbm::dir::PMP] = (distribution_references.f[dirTSE])[ks];
    distribution.f[vf::lbm::dir::MMP] = (distribution_references.f[dirTSW])[ksw];
    distribution.f[vf::lbm::dir::PPM] = (distribution_references.f[dirBNE])[kb];
    distribution.f[vf::lbm::dir::MPM] = (distribution_references.f[dirBNW])[kbw];
    distribution.f[vf::lbm::dir::PMM] = (distribution_references.f[dirBSE])[kbs];
    distribution.f[vf::lbm::dir::MMM] = (distribution_references.f[dirBSW])[kbsw];
    distribution.f[vf::lbm::dir::ZZZ] = (distribution_references.f[dirREST])[k];
}

__device__ void DistributionWrapper::write()
{
    (distribution_references.f[dirE])[k]      = distribution.f[vf::lbm::dir::PZZ];
    (distribution_references.f[dirW])[kw]     = distribution.f[vf::lbm::dir::MZZ];
    (distribution_references.f[dirN])[k]      = distribution.f[vf::lbm::dir::ZPZ];
    (distribution_references.f[dirS])[ks]     = distribution.f[vf::lbm::dir::ZMZ];
    (distribution_references.f[dirT])[k]      = distribution.f[vf::lbm::dir::ZZP];
    (distribution_references.f[dirB])[kb]     = distribution.f[vf::lbm::dir::ZZM];
    (distribution_references.f[dirNE])[k]     = distribution.f[vf::lbm::dir::PPZ];
    (distribution_references.f[dirSW])[ksw]   = distribution.f[vf::lbm::dir::MMZ];
    (distribution_references.f[dirSE])[ks]    = distribution.f[vf::lbm::dir::PMZ];
    (distribution_references.f[dirNW])[kw]    = distribution.f[vf::lbm::dir::MPZ];
    (distribution_references.f[dirTE])[k]     = distribution.f[vf::lbm::dir::PZP];
    (distribution_references.f[dirBW])[kbw]   = distribution.f[vf::lbm::dir::MZM];
    (distribution_references.f[dirBE])[kb]    = distribution.f[vf::lbm::dir::PZM];
    (distribution_references.f[dirTW])[kw]    = distribution.f[vf::lbm::dir::MZP];
    (distribution_references.f[dirTN])[k]     = distribution.f[vf::lbm::dir::ZPP];
    (distribution_references.f[dirBS])[kbs]   = distribution.f[vf::lbm::dir::ZMM];
    (distribution_references.f[dirBN])[kb]    = distribution.f[vf::lbm::dir::ZPM];
    (distribution_references.f[dirTS])[ks]    = distribution.f[vf::lbm::dir::ZMP];
    (distribution_references.f[dirTNE])[k]    = distribution.f[vf::lbm::dir::PPP];
    (distribution_references.f[dirTNW])[kw]   = distribution.f[vf::lbm::dir::MPP];
    (distribution_references.f[dirTSE])[ks]   = distribution.f[vf::lbm::dir::PMP];
    (distribution_references.f[dirTSW])[ksw]  = distribution.f[vf::lbm::dir::MMP];
    (distribution_references.f[dirBNE])[kb]   = distribution.f[vf::lbm::dir::PPM];
    (distribution_references.f[dirBNW])[kbw]  = distribution.f[vf::lbm::dir::MPM];
    (distribution_references.f[dirBSE])[kbs]  = distribution.f[vf::lbm::dir::PMM];
    (distribution_references.f[dirBSW])[kbsw] = distribution.f[vf::lbm::dir::MMM];
    (distribution_references.f[dirREST])[k]   = distribution.f[vf::lbm::dir::ZZZ];
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


} // namespace gpu
} // namespace vf
