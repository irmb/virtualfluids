#include "DistributionHelper.cuh"

#include "LBM/D3Q27.h"

namespace vf
{
namespace gpu
{

__device__ Distributions27 getDistributions27(real* distributions, unsigned int size_Mat, bool isEvenTimestep)
{
    Distributions27 dist;
    if (isEvenTimestep)
    {
        dist.f[dirE   ] = &distributions[dirE   *size_Mat];
        dist.f[dirW   ] = &distributions[dirW   *size_Mat];
        dist.f[dirN   ] = &distributions[dirN   *size_Mat];
        dist.f[dirS   ] = &distributions[dirS   *size_Mat];
        dist.f[dirT   ] = &distributions[dirT   *size_Mat];
        dist.f[dirB   ] = &distributions[dirB   *size_Mat];
        dist.f[dirNE  ] = &distributions[dirNE  *size_Mat];
        dist.f[dirSW  ] = &distributions[dirSW  *size_Mat];
        dist.f[dirSE  ] = &distributions[dirSE  *size_Mat];
        dist.f[dirNW  ] = &distributions[dirNW  *size_Mat];
        dist.f[dirTE  ] = &distributions[dirTE  *size_Mat];
        dist.f[dirBW  ] = &distributions[dirBW  *size_Mat];
        dist.f[dirBE  ] = &distributions[dirBE  *size_Mat];
        dist.f[dirTW  ] = &distributions[dirTW  *size_Mat];
        dist.f[dirTN  ] = &distributions[dirTN  *size_Mat];
        dist.f[dirBS  ] = &distributions[dirBS  *size_Mat];
        dist.f[dirBN  ] = &distributions[dirBN  *size_Mat];
        dist.f[dirTS  ] = &distributions[dirTS  *size_Mat];
        dist.f[dirREST] = &distributions[dirREST*size_Mat];
        dist.f[dirTNE ] = &distributions[dirTNE *size_Mat];
        dist.f[dirTSW ] = &distributions[dirTSW *size_Mat];
        dist.f[dirTSE ] = &distributions[dirTSE *size_Mat];
        dist.f[dirTNW ] = &distributions[dirTNW *size_Mat];
        dist.f[dirBNE ] = &distributions[dirBNE *size_Mat];
        dist.f[dirBSW ] = &distributions[dirBSW *size_Mat];
        dist.f[dirBSE ] = &distributions[dirBSE *size_Mat];
        dist.f[dirBNW ] = &distributions[dirBNW *size_Mat];
    }
    else
    {
        dist.f[dirW   ] = &distributions[dirE   *size_Mat];
        dist.f[dirE   ] = &distributions[dirW   *size_Mat];
        dist.f[dirS   ] = &distributions[dirN   *size_Mat];
        dist.f[dirN   ] = &distributions[dirS   *size_Mat];
        dist.f[dirB   ] = &distributions[dirT   *size_Mat];
        dist.f[dirT   ] = &distributions[dirB   *size_Mat];
        dist.f[dirSW  ] = &distributions[dirNE  *size_Mat];
        dist.f[dirNE  ] = &distributions[dirSW  *size_Mat];
        dist.f[dirNW  ] = &distributions[dirSE  *size_Mat];
        dist.f[dirSE  ] = &distributions[dirNW  *size_Mat];
        dist.f[dirBW  ] = &distributions[dirTE  *size_Mat];
        dist.f[dirTE  ] = &distributions[dirBW  *size_Mat];
        dist.f[dirTW  ] = &distributions[dirBE  *size_Mat];
        dist.f[dirBE  ] = &distributions[dirTW  *size_Mat];
        dist.f[dirBS  ] = &distributions[dirTN  *size_Mat];
        dist.f[dirTN  ] = &distributions[dirBS  *size_Mat];
        dist.f[dirTS  ] = &distributions[dirBN  *size_Mat];
        dist.f[dirBN  ] = &distributions[dirTS  *size_Mat];
        dist.f[dirREST] = &distributions[dirREST*size_Mat];
        dist.f[dirBSW ] = &distributions[dirTNE *size_Mat];
        dist.f[dirBNE ] = &distributions[dirTSW *size_Mat];
        dist.f[dirBNW ] = &distributions[dirTSE *size_Mat];
        dist.f[dirBSE ] = &distributions[dirTNW *size_Mat];
        dist.f[dirTSW ] = &distributions[dirBNE *size_Mat];
        dist.f[dirTNE ] = &distributions[dirBSW *size_Mat];
        dist.f[dirTNW ] = &distributions[dirBSE *size_Mat];
        dist.f[dirTSE ] = &distributions[dirBNW *size_Mat];
    }
    return dist;
}


__device__ DistributionWrapper::DistributionWrapper(
    real* distributions,
    unsigned int size_Mat,
    bool isEvenTimestep,
    uint k,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ) :
    dist(getDistributions27(distributions, size_Mat, isEvenTimestep)),
    k(k),
    kw  (neighborX[k]),
    ks  (neighborY[k]),
    kb  (neighborZ[k]),
    ksw (neighborY[kw]),
    kbw (neighborZ[kw]),
    kbs (neighborZ[ks]),
    kbsw(neighborZ[ksw])
{ 
    read();
}

__device__ void DistributionWrapper::read()
{
    distribution.f[vf::lbm::dir::PZZ] = (dist.f[dirE   ])[k];
    distribution.f[vf::lbm::dir::MZZ] = (dist.f[dirW   ])[kw];
    distribution.f[vf::lbm::dir::ZPZ] = (dist.f[dirN   ])[k];
    distribution.f[vf::lbm::dir::ZMZ] = (dist.f[dirS   ])[ks];
    distribution.f[vf::lbm::dir::ZZP] = (dist.f[dirT   ])[k];
    distribution.f[vf::lbm::dir::ZZM] = (dist.f[dirB   ])[kb];
    distribution.f[vf::lbm::dir::PPZ] = (dist.f[dirNE  ])[k];
    distribution.f[vf::lbm::dir::MMZ] = (dist.f[dirSW  ])[ksw];
    distribution.f[vf::lbm::dir::PMZ] = (dist.f[dirSE  ])[ks];
    distribution.f[vf::lbm::dir::MPZ] = (dist.f[dirNW  ])[kw];
    distribution.f[vf::lbm::dir::PZP] = (dist.f[dirTE  ])[k];
    distribution.f[vf::lbm::dir::MZM] = (dist.f[dirBW  ])[kbw];
    distribution.f[vf::lbm::dir::PZM] = (dist.f[dirBE  ])[kb];
    distribution.f[vf::lbm::dir::MZP] = (dist.f[dirTW  ])[kw];
    distribution.f[vf::lbm::dir::ZPP] = (dist.f[dirTN  ])[k];
    distribution.f[vf::lbm::dir::ZMM] = (dist.f[dirBS  ])[kbs];
    distribution.f[vf::lbm::dir::ZPM] = (dist.f[dirBN  ])[kb];
    distribution.f[vf::lbm::dir::ZMP] = (dist.f[dirTS  ])[ks];
    distribution.f[vf::lbm::dir::PPP] = (dist.f[dirTNE ])[k];
    distribution.f[vf::lbm::dir::MPP] = (dist.f[dirTNW ])[kw];
    distribution.f[vf::lbm::dir::PMP] = (dist.f[dirTSE ])[ks];
    distribution.f[vf::lbm::dir::MMP] = (dist.f[dirTSW ])[ksw];
    distribution.f[vf::lbm::dir::PPM] = (dist.f[dirBNE ])[kb];
    distribution.f[vf::lbm::dir::MPM] = (dist.f[dirBNW ])[kbw];
    distribution.f[vf::lbm::dir::PMM] = (dist.f[dirBSE ])[kbs];
    distribution.f[vf::lbm::dir::MMM] = (dist.f[dirBSW ])[kbsw];
    distribution.f[vf::lbm::dir::ZZZ] = (dist.f[dirREST])[k];
}

__device__ void DistributionWrapper::write()
{
    (dist.f[dirE   ])[k]    = distribution.f[vf::lbm::dir::PZZ];
    (dist.f[dirW   ])[kw]   = distribution.f[vf::lbm::dir::MZZ];
    (dist.f[dirN   ])[k]    = distribution.f[vf::lbm::dir::ZPZ];
    (dist.f[dirS   ])[ks]   = distribution.f[vf::lbm::dir::ZMZ];
    (dist.f[dirT   ])[k]    = distribution.f[vf::lbm::dir::ZZP];
    (dist.f[dirB   ])[kb]   = distribution.f[vf::lbm::dir::ZZM];
    (dist.f[dirNE  ])[k]    = distribution.f[vf::lbm::dir::PPZ];
    (dist.f[dirSW  ])[ksw]  = distribution.f[vf::lbm::dir::MMZ];
    (dist.f[dirSE  ])[ks]   = distribution.f[vf::lbm::dir::PMZ];
    (dist.f[dirNW  ])[kw]   = distribution.f[vf::lbm::dir::MPZ];
    (dist.f[dirTE  ])[k]    = distribution.f[vf::lbm::dir::PZP];
    (dist.f[dirBW  ])[kbw]  = distribution.f[vf::lbm::dir::MZM];
    (dist.f[dirBE  ])[kb]   = distribution.f[vf::lbm::dir::PZM];
    (dist.f[dirTW  ])[kw]   = distribution.f[vf::lbm::dir::MZP];
    (dist.f[dirTN  ])[k]    = distribution.f[vf::lbm::dir::ZPP];
    (dist.f[dirBS  ])[kbs]  = distribution.f[vf::lbm::dir::ZMM];
    (dist.f[dirBN  ])[kb]   = distribution.f[vf::lbm::dir::ZPM];
    (dist.f[dirTS  ])[ks]   = distribution.f[vf::lbm::dir::ZMP];
    (dist.f[dirTNE ])[k]    = distribution.f[vf::lbm::dir::PPP];
    (dist.f[dirTNW ])[kw]   = distribution.f[vf::lbm::dir::MPP];
    (dist.f[dirTSE ])[ks]   = distribution.f[vf::lbm::dir::PMP];
    (dist.f[dirTSW ])[ksw]  = distribution.f[vf::lbm::dir::MMP];
    (dist.f[dirBNE ])[kb]   = distribution.f[vf::lbm::dir::PPM];
    (dist.f[dirBNW ])[kbw]  = distribution.f[vf::lbm::dir::MPM];
    (dist.f[dirBSE ])[kbs]  = distribution.f[vf::lbm::dir::PMM];
    (dist.f[dirBSW ])[kbsw] = distribution.f[vf::lbm::dir::MMM];
    (dist.f[dirREST])[k]    = distribution.f[vf::lbm::dir::ZZZ];
}

__device__ unsigned int getNodeIndex()
{
    const unsigned  x = threadIdx.x;
    const unsigned  y = blockIdx.x;
    const unsigned  z = blockIdx.y;

    const unsigned nx = blockDim.x;
    const unsigned ny = gridDim.x;

    return nx*(ny*z + y) + x;
}

__device__ bool isValidFluidNode(uint k, int size_Mat, uint nodeType)
{
    return (k < size_Mat) && (nodeType == GEO_FLUID);
}

__device__ void getLevelForce(real fx, real fy, real fz, int level, real* forces)
{
    real fx_t {1.}, fy_t {1.}, fz_t {1.};
    for (int i = 0; i < level; i++)
    {
        fx_t *= vf::lbm::constant::c2o1;
        fy_t *= vf::lbm::constant::c2o1;
        fz_t *= vf::lbm::constant::c2o1;
    }

    forces[0] = fx / fx_t;
    forces[1] = fy / fy_t;
    forces[2] = fz / fz_t;
}


}
}
