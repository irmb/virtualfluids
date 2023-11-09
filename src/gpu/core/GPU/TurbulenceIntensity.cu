//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////

/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include "basics/constants/NumericConstants.h"

#include "lbm/MacroscopicQuantities.h"
#include "LBM/GPUHelperFunctions/KernelUtilities.h"


using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

//////////////////////////////////////////////////////////////////////////////
__global__ void CalcTurbulenceIntensity(
   real* vxx,
   real* vyy,
   real* vzz,
   real* vxy,
   real* vxz,
   real* vyz,
   real* vx_mean,
   real* vy_mean,
   real* vz_mean, 
   real *distributions, 
   uint* typeOfGridNode, 
   unsigned int* neighborX,
   unsigned int* neighborY,
   unsigned int* neighborZ,
   unsigned long long numberOfLBnodes, 
   bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    if (nodeIndex >= numberOfLBnodes)
        return;

    if (!isValidFluidNode(typeOfGridNode[nodeIndex]))
        return;

    Distributions27 dist;
    vf::gpu::getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
    vf::gpu::ListIndices listIndices(nodeIndex, neighborX, neighborY, neighborZ);

    real distribution[27];
    vf::gpu::getPreCollisionDistribution(distribution, dist, listIndices);

    // analogue to LBCalcMacCompSP27
    real rho = vf::lbm::getDensity(distribution);
    real vx = vf::lbm::getCompressibleVelocityX1(distribution, rho);
    real vy = vf::lbm::getCompressibleVelocityX2(distribution, rho);
    real vz = vf::lbm::getCompressibleVelocityX3(distribution, rho);

    // compute subtotals:
    // fluctuations
    vxx[nodeIndex] = vxx[nodeIndex] + vx * vx;
    vyy[nodeIndex] = vyy[nodeIndex] + vy * vy;
    vzz[nodeIndex] = vzz[nodeIndex] + vz * vz;
    vxy[nodeIndex] = vxy[nodeIndex] + vx * vy;
    vxz[nodeIndex] = vxz[nodeIndex] + vx * vz;
    vyz[nodeIndex] = vyz[nodeIndex] + vy * vz;

    // velocity (for mean velocity)
    vx_mean[nodeIndex] = vx_mean[nodeIndex] + vx;
    vy_mean[nodeIndex] = vy_mean[nodeIndex] + vy;
    vz_mean[nodeIndex] = vz_mean[nodeIndex] + vz;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
