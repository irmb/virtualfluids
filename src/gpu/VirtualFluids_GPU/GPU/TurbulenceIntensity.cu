//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////

/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
#include "lbm/MacroscopicQuantities.h"
#include "../Kernel/Utilities/DistributionHelper.cuh"


using namespace vf::lbm::constant;

//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void CalcTurbulenceIntensity(
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
   unsigned int size_Mat, 
   bool isEvenTimestep)
{
   const unsigned k = vf::gpu::getNodeIndex();

   if (k >= size_Mat)
       return;

   if (!vf::gpu::isValidFluidNode(typeOfGridNode[k]))
       return;

   vf::gpu::DistributionWrapper distr_wrapper(distributions, size_Mat, isEvenTimestep, k, neighborX, neighborY,
                                              neighborZ);
   const auto &distribution = distr_wrapper.distribution;

   // analogue to LBCalcMacCompSP27
   real rho   = vf::lbm::getDensity(distribution.f);
   real vx    = vf::lbm::getCompressibleVelocityX1(distribution.f, rho);
   real vy    = vf::lbm::getCompressibleVelocityX2(distribution.f, rho);
   real vz    = vf::lbm::getCompressibleVelocityX3(distribution.f, rho);   


   // compute subtotals:
   // fluctuations
   vxx[k] = vxx[k] + vx * vx;
   vyy[k] = vyy[k] + vy * vy;
   vzz[k] = vzz[k] + vz * vz;
   vxy[k] = vxy[k] + vx * vy;
   vxz[k] = vxz[k] + vx * vz;
   vyz[k] = vyz[k] + vy * vz;

   // velocity (for mean velocity)
   vx_mean[k] = vx_mean[k] + vx;
   vy_mean[k] = vy_mean[k] + vy;
   vz_mean[k] = vz_mean[k] + vz; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
