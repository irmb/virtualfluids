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
#include "lbm/constants/NumericConstants.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

#include "lbm/MacroscopicQuantities.h"

#include "../Kernel/Utilities/DistributionHelper.cuh"


////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMac27( real* vxD,
                                        real* vyD,
                                        real* vzD,
                                        real* rhoD,
                                        unsigned int* geoD,
                                        unsigned int* neighborX,
                                        unsigned int* neighborY,
                                        unsigned int* neighborZ,
                                        unsigned int size_Mat,
                                        real* distributions,
                                        bool isEvenTimestep)
{
   const unsigned int tx = threadIdx.x;    // Thread index = lokaler i index
   const unsigned int by = blockIdx.x;     // Block index x
   const unsigned int bz = blockIdx.y;     // Block index y
   const unsigned int x = tx + STARTOFFX;  // Globaler x-Index 
   const unsigned int y = by + STARTOFFY;  // Globaler y-Index 
   const unsigned int z = bz + STARTOFFZ;  // Globaler z-Index 

   const unsigned nx = blockDim.x + 2 * STARTOFFX;
   const unsigned ny = gridDim.x + 2 * STARTOFFY;

   const unsigned int k = nx*(ny*z + y) + x; // Zugriff auf arrays im device


   if(k >= size_Mat)
      return;

   if(!vf::gpu::isValidFluidNode(geoD[k]))
      return;

   rhoD[k] = c0o1;
   vxD[k]  = c0o1;
   vyD[k]  = c0o1;
   vzD[k]  = c0o1;

   vf::gpu::DistributionWrapper distr_wrapper(distributions, size_Mat, isEvenTimestep, k, neighborX, neighborY, neighborZ);
   const auto& distribution = distr_wrapper.distribution;

   rhoD[k] = vf::lbm::getDensity(distribution.f);
   vxD[k] = vf::lbm::getIncompressibleVelocityX1(distribution.f);
   vyD[k] = vf::lbm::getIncompressibleVelocityX2(distribution.f);
   vzD[k] = vf::lbm::getIncompressibleVelocityX3(distribution.f);

}





////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMacSP27( real* vxD,
                                          real* vyD,
                                          real* vzD,
                                          real* rhoD,
                                          real* pressD,
                                          unsigned int* geoD,
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
                                          unsigned int size_Mat,
                                          real* DD,
                                          bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00   ] = &DD[DIR_P00   *size_Mat];
      D.f[DIR_M00   ] = &DD[DIR_M00   *size_Mat];
      D.f[DIR_0P0   ] = &DD[DIR_0P0   *size_Mat];
      D.f[DIR_0M0   ] = &DD[DIR_0M0   *size_Mat];
      D.f[DIR_00P   ] = &DD[DIR_00P   *size_Mat];
      D.f[DIR_00M   ] = &DD[DIR_00M   *size_Mat];
      D.f[DIR_PP0  ] = &DD[DIR_PP0  *size_Mat];
      D.f[DIR_MM0  ] = &DD[DIR_MM0  *size_Mat];
      D.f[DIR_PM0  ] = &DD[DIR_PM0  *size_Mat];
      D.f[DIR_MP0  ] = &DD[DIR_MP0  *size_Mat];
      D.f[DIR_P0P  ] = &DD[DIR_P0P  *size_Mat];
      D.f[DIR_M0M  ] = &DD[DIR_M0M  *size_Mat];
      D.f[DIR_P0M  ] = &DD[DIR_P0M  *size_Mat];
      D.f[DIR_M0P  ] = &DD[DIR_M0P  *size_Mat];
      D.f[DIR_0PP  ] = &DD[DIR_0PP  *size_Mat];
      D.f[DIR_0MM  ] = &DD[DIR_0MM  *size_Mat];
      D.f[DIR_0PM  ] = &DD[DIR_0PM  *size_Mat];
      D.f[DIR_0MP  ] = &DD[DIR_0MP  *size_Mat];
      D.f[DIR_000] = &DD[DIR_000*size_Mat];
      D.f[DIR_PPP ] = &DD[DIR_PPP *size_Mat];
      D.f[DIR_MMP ] = &DD[DIR_MMP *size_Mat];
      D.f[DIR_PMP ] = &DD[DIR_PMP *size_Mat];
      D.f[DIR_MPP ] = &DD[DIR_MPP *size_Mat];
      D.f[DIR_PPM ] = &DD[DIR_PPM *size_Mat];
      D.f[DIR_MMM ] = &DD[DIR_MMM *size_Mat];
      D.f[DIR_PMM ] = &DD[DIR_PMM *size_Mat];
      D.f[DIR_MPM ] = &DD[DIR_MPM *size_Mat];
   } 
   else
   {
      D.f[DIR_M00   ] = &DD[DIR_P00   *size_Mat];
      D.f[DIR_P00   ] = &DD[DIR_M00   *size_Mat];
      D.f[DIR_0M0   ] = &DD[DIR_0P0   *size_Mat];
      D.f[DIR_0P0   ] = &DD[DIR_0M0   *size_Mat];
      D.f[DIR_00M   ] = &DD[DIR_00P   *size_Mat];
      D.f[DIR_00P   ] = &DD[DIR_00M   *size_Mat];
      D.f[DIR_MM0  ] = &DD[DIR_PP0  *size_Mat];
      D.f[DIR_PP0  ] = &DD[DIR_MM0  *size_Mat];
      D.f[DIR_MP0  ] = &DD[DIR_PM0  *size_Mat];
      D.f[DIR_PM0  ] = &DD[DIR_MP0  *size_Mat];
      D.f[DIR_M0M  ] = &DD[DIR_P0P  *size_Mat];
      D.f[DIR_P0P  ] = &DD[DIR_M0M  *size_Mat];
      D.f[DIR_M0P  ] = &DD[DIR_P0M  *size_Mat];
      D.f[DIR_P0M  ] = &DD[DIR_M0P  *size_Mat];
      D.f[DIR_0MM  ] = &DD[DIR_0PP  *size_Mat];
      D.f[DIR_0PP  ] = &DD[DIR_0MM  *size_Mat];
      D.f[DIR_0MP  ] = &DD[DIR_0PM  *size_Mat];
      D.f[DIR_0PM  ] = &DD[DIR_0MP  *size_Mat];
      D.f[DIR_000] = &DD[DIR_000*size_Mat];
      D.f[DIR_PPP ] = &DD[DIR_MMM *size_Mat];
      D.f[DIR_MMP ] = &DD[DIR_PPM *size_Mat];
      D.f[DIR_PMP ] = &DD[DIR_MPM *size_Mat];
      D.f[DIR_MPP ] = &DD[DIR_PMM *size_Mat];
      D.f[DIR_PPM ] = &DD[DIR_MMP *size_Mat];
      D.f[DIR_MMM ] = &DD[DIR_PPP *size_Mat];
      D.f[DIR_PMM ] = &DD[DIR_MPP *size_Mat];
      D.f[DIR_MPM ] = &DD[DIR_PMP *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= k;
      unsigned int ke   = k;
      unsigned int kw   = neighborX[k];
      unsigned int kn   = k;
      unsigned int ks   = neighborY[k];
      unsigned int kt   = k;
      unsigned int kb   = neighborZ[k];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = k;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = k;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = k;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = k;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
      pressD[k] = c0o1;
	  rhoD[k]   = c0o1;
	  vxD[k]    = c0o1;
	  vyD[k]    = c0o1;
	  vzD[k]    = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         rhoD[k]    =   (D.f[DIR_P00   ])[ke  ]+ (D.f[DIR_M00   ])[kw  ]+ 
                        (D.f[DIR_0P0   ])[kn  ]+ (D.f[DIR_0M0   ])[ks  ]+
                        (D.f[DIR_00P   ])[kt  ]+ (D.f[DIR_00M   ])[kb  ]+
                        (D.f[DIR_PP0  ])[kne ]+ (D.f[DIR_MM0  ])[ksw ]+
                        (D.f[DIR_PM0  ])[kse ]+ (D.f[DIR_MP0  ])[knw ]+
                        (D.f[DIR_P0P  ])[kte ]+ (D.f[DIR_M0M  ])[kbw ]+
                        (D.f[DIR_P0M  ])[kbe ]+ (D.f[DIR_M0P  ])[ktw ]+
                        (D.f[DIR_0PP  ])[ktn ]+ (D.f[DIR_0MM  ])[kbs ]+
                        (D.f[DIR_0PM  ])[kbn ]+ (D.f[DIR_0MP  ])[kts ]+
                        (D.f[DIR_000])[kzero]+ 
                        (D.f[DIR_PPP ])[ktne]+ (D.f[DIR_MMP ])[ktsw]+ 
                        (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]+ 
                        (D.f[DIR_PPM ])[kbne]+ (D.f[DIR_MMM ])[kbsw]+ 
                        (D.f[DIR_PMM ])[kbse]+ (D.f[DIR_MPM ])[kbnw];

         vxD[k]     =   (D.f[DIR_P00   ])[ke  ]- (D.f[DIR_M00   ])[kw  ]+ 
                        (D.f[DIR_PP0  ])[kne ]- (D.f[DIR_MM0  ])[ksw ]+
                        (D.f[DIR_PM0  ])[kse ]- (D.f[DIR_MP0  ])[knw ]+
                        (D.f[DIR_P0P  ])[kte ]- (D.f[DIR_M0M  ])[kbw ]+
                        (D.f[DIR_P0M  ])[kbe ]- (D.f[DIR_M0P  ])[ktw ]+
                        (D.f[DIR_PPP ])[ktne]- (D.f[DIR_MMP ])[ktsw]+ 
                        (D.f[DIR_PMP ])[ktse]- (D.f[DIR_MPP ])[ktnw]+ 
                        (D.f[DIR_PPM ])[kbne]- (D.f[DIR_MMM ])[kbsw]+ 
                        (D.f[DIR_PMM ])[kbse]- (D.f[DIR_MPM ])[kbnw];

         vyD[k]     =   (D.f[DIR_0P0   ])[kn  ]- (D.f[DIR_0M0   ])[ks  ]+
                        (D.f[DIR_PP0  ])[kne ]- (D.f[DIR_MM0  ])[ksw ]-
                        (D.f[DIR_PM0  ])[kse ]+ (D.f[DIR_MP0  ])[knw ]+
                        (D.f[DIR_0PP  ])[ktn ]- (D.f[DIR_0MM  ])[kbs ]+
                        (D.f[DIR_0PM  ])[kbn ]- (D.f[DIR_0MP  ])[kts ]+
                        (D.f[DIR_PPP ])[ktne]- (D.f[DIR_MMP ])[ktsw]- 
                        (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]+ 
                        (D.f[DIR_PPM ])[kbne]- (D.f[DIR_MMM ])[kbsw]- 
                        (D.f[DIR_PMM ])[kbse]+ (D.f[DIR_MPM ])[kbnw];

         vzD[k]     =   (D.f[DIR_00P   ])[kt  ]- (D.f[DIR_00M   ])[kb  ]+
                        (D.f[DIR_P0P  ])[kte ]- (D.f[DIR_M0M  ])[kbw ]-
                        (D.f[DIR_P0M  ])[kbe ]+ (D.f[DIR_M0P  ])[ktw ]+
                        (D.f[DIR_0PP  ])[ktn ]- (D.f[DIR_0MM  ])[kbs ]-
                        (D.f[DIR_0PM  ])[kbn ]+ (D.f[DIR_0MP  ])[kts ]+
                        (D.f[DIR_PPP ])[ktne]+ (D.f[DIR_MMP ])[ktsw]+ 
                        (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]- 
                        (D.f[DIR_PPM ])[kbne]- (D.f[DIR_MMM ])[kbsw]- 
                        (D.f[DIR_PMM ])[kbse]- (D.f[DIR_MPM ])[kbnw];

         pressD[k]  =  ((D.f[DIR_P00   ])[ke  ]+ (D.f[DIR_M00   ])[kw  ]+ 
                        (D.f[DIR_0P0   ])[kn  ]+ (D.f[DIR_0M0   ])[ks  ]+
                        (D.f[DIR_00P   ])[kt  ]+ (D.f[DIR_00M   ])[kb  ]+
                        2.f*(
                        (D.f[DIR_PP0  ])[kne ]+ (D.f[DIR_MM0  ])[ksw ]+
                        (D.f[DIR_PM0  ])[kse ]+ (D.f[DIR_MP0  ])[knw ]+
                        (D.f[DIR_P0P  ])[kte ]+ (D.f[DIR_M0M  ])[kbw ]+
                        (D.f[DIR_P0M  ])[kbe ]+ (D.f[DIR_M0P  ])[ktw ]+
                        (D.f[DIR_0PP  ])[ktn ]+ (D.f[DIR_0MM  ])[kbs ]+
                        (D.f[DIR_0PM  ])[kbn ]+ (D.f[DIR_0MP  ])[kts ])+
                        3.f*(
                        (D.f[DIR_PPP ])[ktne]+ (D.f[DIR_MMP ])[ktsw]+ 
                        (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]+ 
                        (D.f[DIR_PPM ])[kbne]+ (D.f[DIR_MMM ])[kbsw]+ 
                        (D.f[DIR_PMM ])[kbse]+ (D.f[DIR_MPM ])[kbnw])-
                        rhoD[k]-(vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1+c0o1*rhoD[k])) * c1o2+rhoD[k]; // times zero for incompressible case   
         //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5

      }
   }
}


////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMacCompSP27(
   real *vxD,
   real *vyD,
   real *vzD,
   real *rhoD,
   real *pressD,
   unsigned int *geoD,
   unsigned int *neighborX,
   unsigned int *neighborY,
   unsigned int *neighborZ,
   unsigned int size_Mat,
   real *distributions,
   bool isEvenTimestep)
{
    const unsigned k = vf::gpu::getNodeIndex();

    if(k >= size_Mat)
        return;

    pressD[k] = c0o1;
    rhoD[k]   = c0o1;
    vxD[k]    = c0o1;
    vyD[k]    = c0o1;
    vzD[k]    = c0o1;

    if (!vf::gpu::isValidFluidNode(geoD[k]))
        return;

    vf::gpu::DistributionWrapper distr_wrapper(distributions, size_Mat, isEvenTimestep, k, neighborX, neighborY,
                                               neighborZ);
    const auto &distribution = distr_wrapper.distribution;

    rhoD[k]   = vf::lbm::getDensity(distribution.f);
    vxD[k]    = vf::lbm::getCompressibleVelocityX1(distribution.f, rhoD[k]);
    vyD[k]    = vf::lbm::getCompressibleVelocityX2(distribution.f, rhoD[k]);
    vzD[k]    = vf::lbm::getCompressibleVelocityX3(distribution.f, rhoD[k]);
    pressD[k] = vf::lbm::getPressure(distribution.f, rhoD[k], vxD[k], vyD[k], vzD[k]); 
}




































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMedSP27( real* vxD,
                                          real* vyD,
                                          real* vzD,
                                          real* rhoD,
                                          real* pressD,
                                          unsigned int* geoD,
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
                                          unsigned int size_Mat,
                                          real* DD,
                                          bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00   ] = &DD[DIR_P00   *size_Mat];
      D.f[DIR_M00   ] = &DD[DIR_M00   *size_Mat];
      D.f[DIR_0P0   ] = &DD[DIR_0P0   *size_Mat];
      D.f[DIR_0M0   ] = &DD[DIR_0M0   *size_Mat];
      D.f[DIR_00P   ] = &DD[DIR_00P   *size_Mat];
      D.f[DIR_00M   ] = &DD[DIR_00M   *size_Mat];
      D.f[DIR_PP0  ] = &DD[DIR_PP0  *size_Mat];
      D.f[DIR_MM0  ] = &DD[DIR_MM0  *size_Mat];
      D.f[DIR_PM0  ] = &DD[DIR_PM0  *size_Mat];
      D.f[DIR_MP0  ] = &DD[DIR_MP0  *size_Mat];
      D.f[DIR_P0P  ] = &DD[DIR_P0P  *size_Mat];
      D.f[DIR_M0M  ] = &DD[DIR_M0M  *size_Mat];
      D.f[DIR_P0M  ] = &DD[DIR_P0M  *size_Mat];
      D.f[DIR_M0P  ] = &DD[DIR_M0P  *size_Mat];
      D.f[DIR_0PP  ] = &DD[DIR_0PP  *size_Mat];
      D.f[DIR_0MM  ] = &DD[DIR_0MM  *size_Mat];
      D.f[DIR_0PM  ] = &DD[DIR_0PM  *size_Mat];
      D.f[DIR_0MP  ] = &DD[DIR_0MP  *size_Mat];
      D.f[DIR_000] = &DD[DIR_000*size_Mat];
      D.f[DIR_PPP ] = &DD[DIR_PPP *size_Mat];
      D.f[DIR_MMP ] = &DD[DIR_MMP *size_Mat];
      D.f[DIR_PMP ] = &DD[DIR_PMP *size_Mat];
      D.f[DIR_MPP ] = &DD[DIR_MPP *size_Mat];
      D.f[DIR_PPM ] = &DD[DIR_PPM *size_Mat];
      D.f[DIR_MMM ] = &DD[DIR_MMM *size_Mat];
      D.f[DIR_PMM ] = &DD[DIR_PMM *size_Mat];
      D.f[DIR_MPM ] = &DD[DIR_MPM *size_Mat];
   } 
   else
   {
      D.f[DIR_M00   ] = &DD[DIR_P00   *size_Mat];
      D.f[DIR_P00   ] = &DD[DIR_M00   *size_Mat];
      D.f[DIR_0M0   ] = &DD[DIR_0P0   *size_Mat];
      D.f[DIR_0P0   ] = &DD[DIR_0M0   *size_Mat];
      D.f[DIR_00M   ] = &DD[DIR_00P   *size_Mat];
      D.f[DIR_00P   ] = &DD[DIR_00M   *size_Mat];
      D.f[DIR_MM0  ] = &DD[DIR_PP0  *size_Mat];
      D.f[DIR_PP0  ] = &DD[DIR_MM0  *size_Mat];
      D.f[DIR_MP0  ] = &DD[DIR_PM0  *size_Mat];
      D.f[DIR_PM0  ] = &DD[DIR_MP0  *size_Mat];
      D.f[DIR_M0M  ] = &DD[DIR_P0P  *size_Mat];
      D.f[DIR_P0P  ] = &DD[DIR_M0M  *size_Mat];
      D.f[DIR_M0P  ] = &DD[DIR_P0M  *size_Mat];
      D.f[DIR_P0M  ] = &DD[DIR_M0P  *size_Mat];
      D.f[DIR_0MM  ] = &DD[DIR_0PP  *size_Mat];
      D.f[DIR_0PP  ] = &DD[DIR_0MM  *size_Mat];
      D.f[DIR_0MP  ] = &DD[DIR_0PM  *size_Mat];
      D.f[DIR_0PM  ] = &DD[DIR_0MP  *size_Mat];
      D.f[DIR_000] = &DD[DIR_000*size_Mat];
      D.f[DIR_PPP ] = &DD[DIR_MMM *size_Mat];
      D.f[DIR_MMP ] = &DD[DIR_PPM *size_Mat];
      D.f[DIR_PMP ] = &DD[DIR_MPM *size_Mat];
      D.f[DIR_MPP ] = &DD[DIR_PMM *size_Mat];
      D.f[DIR_PPM ] = &DD[DIR_MMP *size_Mat];
      D.f[DIR_MMM ] = &DD[DIR_PPP *size_Mat];
      D.f[DIR_PMM ] = &DD[DIR_MPP *size_Mat];
      D.f[DIR_MPM ] = &DD[DIR_PMP *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= k;
      unsigned int ke   = k;
      unsigned int kw   = neighborX[k];
      unsigned int kn   = k;
      unsigned int ks   = neighborY[k];
      unsigned int kt   = k;
      unsigned int kb   = neighborZ[k];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = k;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = k;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = k;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = k;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
      real PRESS = pressD[k];
      real RHO   = rhoD[k];
      real VX    = vxD[k];
      real VY    = vyD[k];
      real VZ    = vzD[k];
      //////////////////////////////////////////////////////////////////////////
      pressD[k] = c0o1;
	  rhoD[k]   = c0o1;
	  vxD[k]    = c0o1;
	  vyD[k]    = c0o1;
	  vzD[k]    = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         rhoD[k]    =   (D.f[DIR_P00   ])[ke  ]+ (D.f[DIR_M00   ])[kw  ]+ 
                        (D.f[DIR_0P0   ])[kn  ]+ (D.f[DIR_0M0   ])[ks  ]+
                        (D.f[DIR_00P   ])[kt  ]+ (D.f[DIR_00M   ])[kb  ]+
                        (D.f[DIR_PP0  ])[kne ]+ (D.f[DIR_MM0  ])[ksw ]+
                        (D.f[DIR_PM0  ])[kse ]+ (D.f[DIR_MP0  ])[knw ]+
                        (D.f[DIR_P0P  ])[kte ]+ (D.f[DIR_M0M  ])[kbw ]+
                        (D.f[DIR_P0M  ])[kbe ]+ (D.f[DIR_M0P  ])[ktw ]+
                        (D.f[DIR_0PP  ])[ktn ]+ (D.f[DIR_0MM  ])[kbs ]+
                        (D.f[DIR_0PM  ])[kbn ]+ (D.f[DIR_0MP  ])[kts ]+
                        (D.f[DIR_000])[kzero]+ 
                        (D.f[DIR_PPP ])[ktne]+ (D.f[DIR_MMP ])[ktsw]+ 
                        (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]+ 
                        (D.f[DIR_PPM ])[kbne]+ (D.f[DIR_MMM ])[kbsw]+ 
                        (D.f[DIR_PMM ])[kbse]+ (D.f[DIR_MPM ])[kbnw]+
                        RHO;

         vxD[k]     =   (D.f[DIR_P00   ])[ke  ]- (D.f[DIR_M00   ])[kw  ]+ 
                        (D.f[DIR_PP0  ])[kne ]- (D.f[DIR_MM0  ])[ksw ]+
                        (D.f[DIR_PM0  ])[kse ]- (D.f[DIR_MP0  ])[knw ]+
                        (D.f[DIR_P0P  ])[kte ]- (D.f[DIR_M0M  ])[kbw ]+
                        (D.f[DIR_P0M  ])[kbe ]- (D.f[DIR_M0P  ])[ktw ]+
                        (D.f[DIR_PPP ])[ktne]- (D.f[DIR_MMP ])[ktsw]+ 
                        (D.f[DIR_PMP ])[ktse]- (D.f[DIR_MPP ])[ktnw]+ 
                        (D.f[DIR_PPM ])[kbne]- (D.f[DIR_MMM ])[kbsw]+ 
                        (D.f[DIR_PMM ])[kbse]- (D.f[DIR_MPM ])[kbnw]+
                        VX;

         vyD[k]     =   (D.f[DIR_0P0   ])[kn  ]- (D.f[DIR_0M0   ])[ks  ]+
                        (D.f[DIR_PP0  ])[kne ]- (D.f[DIR_MM0  ])[ksw ]-
                        (D.f[DIR_PM0  ])[kse ]+ (D.f[DIR_MP0  ])[knw ]+
                        (D.f[DIR_0PP  ])[ktn ]- (D.f[DIR_0MM  ])[kbs ]+
                        (D.f[DIR_0PM  ])[kbn ]- (D.f[DIR_0MP  ])[kts ]+
                        (D.f[DIR_PPP ])[ktne]- (D.f[DIR_MMP ])[ktsw]- 
                        (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]+ 
                        (D.f[DIR_PPM ])[kbne]- (D.f[DIR_MMM ])[kbsw]- 
                        (D.f[DIR_PMM ])[kbse]+ (D.f[DIR_MPM ])[kbnw]+
                        VY;

         vzD[k]     =   (D.f[DIR_00P   ])[kt  ]- (D.f[DIR_00M   ])[kb  ]+
                        (D.f[DIR_P0P  ])[kte ]- (D.f[DIR_M0M  ])[kbw ]-
                        (D.f[DIR_P0M  ])[kbe ]+ (D.f[DIR_M0P  ])[ktw ]+
                        (D.f[DIR_0PP  ])[ktn ]- (D.f[DIR_0MM  ])[kbs ]-
                        (D.f[DIR_0PM  ])[kbn ]+ (D.f[DIR_0MP  ])[kts ]+
                        (D.f[DIR_PPP ])[ktne]+ (D.f[DIR_MMP ])[ktsw]+ 
                        (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]- 
                        (D.f[DIR_PPM ])[kbne]- (D.f[DIR_MMM ])[kbsw]- 
                        (D.f[DIR_PMM ])[kbse]- (D.f[DIR_MPM ])[kbnw]+
                        VZ;

         pressD[k]  =   ((D.f[DIR_P00   ])[ke  ]+ (D.f[DIR_M00   ])[kw  ]+ 
                        (D.f[DIR_0P0   ])[kn  ]+ (D.f[DIR_0M0   ])[ks  ]+
                        (D.f[DIR_00P   ])[kt  ]+ (D.f[DIR_00M   ])[kb  ]+
                        c2o1*(
                        (D.f[DIR_PP0  ])[kne ]+ (D.f[DIR_MM0  ])[ksw ]+
                        (D.f[DIR_PM0  ])[kse ]+ (D.f[DIR_MP0  ])[knw ]+
                        (D.f[DIR_P0P  ])[kte ]+ (D.f[DIR_M0M  ])[kbw ]+
                        (D.f[DIR_P0M  ])[kbe ]+ (D.f[DIR_M0P  ])[ktw ]+
                        (D.f[DIR_0PP  ])[ktn ]+ (D.f[DIR_0MM  ])[kbs ]+
                        (D.f[DIR_0PM  ])[kbn ]+ (D.f[DIR_0MP  ])[kts ])+
                        c3o1*(
                        (D.f[DIR_PPP ])[ktne]+ (D.f[DIR_MMP ])[ktsw]+ 
                        (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]+ 
                        (D.f[DIR_PPM ])[kbne]+ (D.f[DIR_MMM ])[kbsw]+ 
                        (D.f[DIR_PMM ])[kbse]+ (D.f[DIR_MPM ])[kbnw])-
                        rhoD[k]-(vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1+rhoD[k])) * c1o2+rhoD[k]+
                        PRESS;    
         //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
      }
   }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMedCompSP27( real* vxD,
											  real* vyD,
											  real* vzD,
											  real* rhoD,
											  real* pressD,
											  unsigned int* geoD,
											  unsigned int* neighborX,
											  unsigned int* neighborY,
											  unsigned int* neighborZ,
											  unsigned int size_Mat,
											  real* DD,
											  bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00   ] = &DD[DIR_P00   *size_Mat];
      D.f[DIR_M00   ] = &DD[DIR_M00   *size_Mat];
      D.f[DIR_0P0   ] = &DD[DIR_0P0   *size_Mat];
      D.f[DIR_0M0   ] = &DD[DIR_0M0   *size_Mat];
      D.f[DIR_00P   ] = &DD[DIR_00P   *size_Mat];
      D.f[DIR_00M   ] = &DD[DIR_00M   *size_Mat];
      D.f[DIR_PP0  ] = &DD[DIR_PP0  *size_Mat];
      D.f[DIR_MM0  ] = &DD[DIR_MM0  *size_Mat];
      D.f[DIR_PM0  ] = &DD[DIR_PM0  *size_Mat];
      D.f[DIR_MP0  ] = &DD[DIR_MP0  *size_Mat];
      D.f[DIR_P0P  ] = &DD[DIR_P0P  *size_Mat];
      D.f[DIR_M0M  ] = &DD[DIR_M0M  *size_Mat];
      D.f[DIR_P0M  ] = &DD[DIR_P0M  *size_Mat];
      D.f[DIR_M0P  ] = &DD[DIR_M0P  *size_Mat];
      D.f[DIR_0PP  ] = &DD[DIR_0PP  *size_Mat];
      D.f[DIR_0MM  ] = &DD[DIR_0MM  *size_Mat];
      D.f[DIR_0PM  ] = &DD[DIR_0PM  *size_Mat];
      D.f[DIR_0MP  ] = &DD[DIR_0MP  *size_Mat];
      D.f[DIR_000] = &DD[DIR_000*size_Mat];
      D.f[DIR_PPP ] = &DD[DIR_PPP *size_Mat];
      D.f[DIR_MMP ] = &DD[DIR_MMP *size_Mat];
      D.f[DIR_PMP ] = &DD[DIR_PMP *size_Mat];
      D.f[DIR_MPP ] = &DD[DIR_MPP *size_Mat];
      D.f[DIR_PPM ] = &DD[DIR_PPM *size_Mat];
      D.f[DIR_MMM ] = &DD[DIR_MMM *size_Mat];
      D.f[DIR_PMM ] = &DD[DIR_PMM *size_Mat];
      D.f[DIR_MPM ] = &DD[DIR_MPM *size_Mat];
   } 
   else
   {
      D.f[DIR_M00   ] = &DD[DIR_P00   *size_Mat];
      D.f[DIR_P00   ] = &DD[DIR_M00   *size_Mat];
      D.f[DIR_0M0   ] = &DD[DIR_0P0   *size_Mat];
      D.f[DIR_0P0   ] = &DD[DIR_0M0   *size_Mat];
      D.f[DIR_00M   ] = &DD[DIR_00P   *size_Mat];
      D.f[DIR_00P   ] = &DD[DIR_00M   *size_Mat];
      D.f[DIR_MM0  ] = &DD[DIR_PP0  *size_Mat];
      D.f[DIR_PP0  ] = &DD[DIR_MM0  *size_Mat];
      D.f[DIR_MP0  ] = &DD[DIR_PM0  *size_Mat];
      D.f[DIR_PM0  ] = &DD[DIR_MP0  *size_Mat];
      D.f[DIR_M0M  ] = &DD[DIR_P0P  *size_Mat];
      D.f[DIR_P0P  ] = &DD[DIR_M0M  *size_Mat];
      D.f[DIR_M0P  ] = &DD[DIR_P0M  *size_Mat];
      D.f[DIR_P0M  ] = &DD[DIR_M0P  *size_Mat];
      D.f[DIR_0MM  ] = &DD[DIR_0PP  *size_Mat];
      D.f[DIR_0PP  ] = &DD[DIR_0MM  *size_Mat];
      D.f[DIR_0MP  ] = &DD[DIR_0PM  *size_Mat];
      D.f[DIR_0PM  ] = &DD[DIR_0MP  *size_Mat];
      D.f[DIR_000] = &DD[DIR_000*size_Mat];
      D.f[DIR_PPP ] = &DD[DIR_MMM *size_Mat];
      D.f[DIR_MMP ] = &DD[DIR_PPM *size_Mat];
      D.f[DIR_PMP ] = &DD[DIR_MPM *size_Mat];
      D.f[DIR_MPP ] = &DD[DIR_PMM *size_Mat];
      D.f[DIR_PPM ] = &DD[DIR_MMP *size_Mat];
      D.f[DIR_MMM ] = &DD[DIR_PPP *size_Mat];
      D.f[DIR_PMM ] = &DD[DIR_MPP *size_Mat];
      D.f[DIR_MPM ] = &DD[DIR_PMP *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      //unsigned int kzero= k;
      unsigned int ke   = k;
      unsigned int kw   = neighborX[k];
      unsigned int kn   = k;
      unsigned int ks   = neighborY[k];
      unsigned int kt   = k;
      unsigned int kb   = neighborZ[k];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = k;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = k;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = k;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = k;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
      real PRESS = pressD[k];
      real RHO   = rhoD[k];
      real VX    = vxD[k];
      real VY    = vyD[k];
      real VZ    = vzD[k];
      //////////////////////////////////////////////////////////////////////////
      pressD[k] = c0o1;
	  rhoD[k]   = c0o1;
	  vxD[k]    = c0o1;
	  vyD[k]    = c0o1;
	  vzD[k]    = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
		  real mfcbb = (D.f[DIR_P00])[k];//[ke   ];
		  real mfabb = (D.f[DIR_M00])[kw];//[kw   ];  
		  real mfbcb = (D.f[DIR_0P0])[k];//[kn   ];
		  real mfbab = (D.f[DIR_0M0])[ks];//[ks   ];  
		  real mfbbc = (D.f[DIR_00P])[k];//[kt   ];
		  real mfbba = (D.f[DIR_00M])[kb];//[kb   ];  
		  real mfccb = (D.f[DIR_PP0])[k];//[kne  ];  
		  real mfaab = (D.f[DIR_MM0])[ksw];//[ksw  ];
		  real mfcab = (D.f[DIR_PM0])[ks];//[kse  ]; 
		  real mfacb = (D.f[DIR_MP0])[kw];//[knw  ]; 
		  real mfcbc = (D.f[DIR_P0P])[k];//[kte  ];  
		  real mfaba = (D.f[DIR_M0M])[kbw];//[kbw  ];
		  real mfcba = (D.f[DIR_P0M])[kb];//[kbe  ]; 
		  real mfabc = (D.f[DIR_M0P])[kw];//[ktw  ]; 
		  real mfbcc = (D.f[DIR_0PP])[k];//[ktn  ];  
		  real mfbaa = (D.f[DIR_0MM])[kbs];//[kbs  ];
		  real mfbca = (D.f[DIR_0PM])[kb];//[kbn  ]; 
		  real mfbac = (D.f[DIR_0MP])[ks];//[kts  ]; 
		  real mfbbb = (D.f[DIR_000])[k];//[kzero];
		  real mfccc = (D.f[DIR_PPP])[k];//[ktne ]; 
		  real mfaac = (D.f[DIR_MMP])[ksw];//[ktsw ]; 
		  real mfcac = (D.f[DIR_PMP])[ks];//[ktse ];
		  real mfacc = (D.f[DIR_MPP])[kw];//[ktnw ];
		  real mfcca = (D.f[DIR_PPM])[kb];//[kbne ];
		  real mfaaa = (D.f[DIR_MMM])[kbsw];//[kbsw ];
		  real mfcaa = (D.f[DIR_PMM])[kbs];//[kbse ]; 
		  real mfaca = (D.f[DIR_MPM])[kbw];//[kbnw ]; 
		  ////////////////////////////////////////////////////////////////////////////////////
		  real drho = 
			  ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
			  (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
			  ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

		  real rho = c1o1 + drho;
		  
		  rhoD[k] = drho + RHO;

		  vxD[k] = 
			  (((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
			  (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
			  (mfcbb - mfabb)) / rho) + VX;
		  vyD[k] = 
			  (((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
			  (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
			  (mfbcb - mfbab)) / rho) + VY;
		  vzD[k] = 
			  (((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
			  (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
			  (mfbbc - mfbba)) / rho) + VZ;

		  //rhoD[k] =
			 // (D.f[DIR_P00])[ke] + (D.f[DIR_M00])[kw] +
			 // (D.f[DIR_0P0])[kn] + (D.f[DIR_0M0])[ks] +
			 // (D.f[DIR_00P])[kt] + (D.f[DIR_00M])[kb] +
			 // (D.f[DIR_PP0])[kne] + (D.f[DIR_MM0])[ksw] +
			 // (D.f[DIR_PM0])[kse] + (D.f[DIR_MP0])[knw] +
			 // (D.f[DIR_P0P])[kte] + (D.f[DIR_M0M])[kbw] +
			 // (D.f[DIR_P0M])[kbe] + (D.f[DIR_M0P])[ktw] +
			 // (D.f[DIR_0PP])[ktn] + (D.f[DIR_0MM])[kbs] +
			 // (D.f[DIR_0PM])[kbn] + (D.f[DIR_0MP])[kts] +
			 // (D.f[DIR_000])[kzero] +
			 // (D.f[DIR_PPP])[ktne] + (D.f[DIR_MMP])[ktsw] +
			 // (D.f[DIR_PMP])[ktse] + (D.f[DIR_MPP])[ktnw] +
			 // (D.f[DIR_PPM])[kbne] + (D.f[DIR_MMM])[kbsw] +
			 // (D.f[DIR_PMM])[kbse] + (D.f[DIR_MPM])[kbnw];// +RHO;

    //     vxD[k] =  
			 //((D.f[DIR_P00  ])[ke  ]- (D.f[DIR_M00   ])[kw  ]+ 
    //         (D.f[DIR_PP0  ])[kne ]- (D.f[DIR_MM0  ])[ksw ]+
    //         (D.f[DIR_PM0  ])[kse ]- (D.f[DIR_MP0  ])[knw ]+
    //         (D.f[DIR_P0P  ])[kte ]- (D.f[DIR_M0M  ])[kbw ]+
    //         (D.f[DIR_P0M  ])[kbe ]- (D.f[DIR_M0P  ])[ktw ]+
    //         (D.f[DIR_PPP ])[ktne]- (D.f[DIR_MMP ])[ktsw]+ 
    //         (D.f[DIR_PMP ])[ktse]- (D.f[DIR_MPP ])[ktnw]+ 
    //         (D.f[DIR_PPM ])[kbne]- (D.f[DIR_MMM ])[kbsw]+ 
    //         (D.f[DIR_PMM ])[kbse]- (D.f[DIR_MPM ])[kbnw]) / (one + rhoD[k])+
    //         VX;

    //     vyD[k] =  
			 //((D.f[DIR_0P0  ])[kn  ]- (D.f[DIR_0M0   ])[ks  ]+
    //         (D.f[DIR_PP0  ])[kne ]- (D.f[DIR_MM0  ])[ksw ]-
    //         (D.f[DIR_PM0  ])[kse ]+ (D.f[DIR_MP0  ])[knw ]+
    //         (D.f[DIR_0PP  ])[ktn ]- (D.f[DIR_0MM  ])[kbs ]+
    //         (D.f[DIR_0PM  ])[kbn ]- (D.f[DIR_0MP  ])[kts ]+
    //         (D.f[DIR_PPP ])[ktne]- (D.f[DIR_MMP ])[ktsw]- 
    //         (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]+ 
    //         (D.f[DIR_PPM ])[kbne]- (D.f[DIR_MMM ])[kbsw]- 
    //         (D.f[DIR_PMM ])[kbse]+ (D.f[DIR_MPM ])[kbnw]) / (one + rhoD[k])+
    //         VY;

    //     vzD[k] =  
			 //((D.f[DIR_00P  ])[kt  ]- (D.f[DIR_00M   ])[kb  ]+
    //         (D.f[DIR_P0P  ])[kte ]- (D.f[DIR_M0M  ])[kbw ]-
    //         (D.f[DIR_P0M  ])[kbe ]+ (D.f[DIR_M0P  ])[ktw ]+
    //         (D.f[DIR_0PP  ])[ktn ]- (D.f[DIR_0MM  ])[kbs ]-
    //         (D.f[DIR_0PM  ])[kbn ]+ (D.f[DIR_0MP  ])[kts ]+
    //         (D.f[DIR_PPP ])[ktne]+ (D.f[DIR_MMP ])[ktsw]+ 
    //         (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]- 
    //         (D.f[DIR_PPM ])[kbne]- (D.f[DIR_MMM ])[kbsw]- 
    //         (D.f[DIR_PMM ])[kbse]- (D.f[DIR_MPM ])[kbnw]) / (one + rhoD[k])+
    //         VZ;

         pressD[k]  =  ((D.f[DIR_P00   ])[ke  ]+ (D.f[DIR_M00   ])[kw  ]+ 
                        (D.f[DIR_0P0   ])[kn  ]+ (D.f[DIR_0M0   ])[ks  ]+
                        (D.f[DIR_00P   ])[kt  ]+ (D.f[DIR_00M   ])[kb  ]+
                        c2o1*(
                        (D.f[DIR_PP0  ])[kne ]+ (D.f[DIR_MM0  ])[ksw ]+
                        (D.f[DIR_PM0  ])[kse ]+ (D.f[DIR_MP0  ])[knw ]+
                        (D.f[DIR_P0P  ])[kte ]+ (D.f[DIR_M0M  ])[kbw ]+
                        (D.f[DIR_P0M  ])[kbe ]+ (D.f[DIR_M0P  ])[ktw ]+
                        (D.f[DIR_0PP  ])[ktn ]+ (D.f[DIR_0MM  ])[kbs ]+
                        (D.f[DIR_0PM  ])[kbn ]+ (D.f[DIR_0MP  ])[kts ])+
                        c3o1*(
                        (D.f[DIR_PPP ])[ktne]+ (D.f[DIR_MMP ])[ktsw]+ 
                        (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]+ 
                        (D.f[DIR_PPM ])[kbne]+ (D.f[DIR_MMM ])[kbsw]+ 
                        (D.f[DIR_PMM ])[kbse]+ (D.f[DIR_MPM ])[kbnw])-
                        rhoD[k]-(vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1+rhoD[k])) * c1o2+rhoD[k]+
                        PRESS;    
         //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
      }
   }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMedCompAD27(
	real* vxD,
	real* vyD,
	real* vzD,
	real* rhoD,
	real* pressD,
	real* concD,
	unsigned int* geoD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int size_Mat,
	real* DD,
	real* DD_AD,
	bool isEvenTimestep)
{
	Distributions27 D;
	if (isEvenTimestep == true)
	{
		D.f[DIR_P00] = &DD[DIR_P00   *size_Mat];
		D.f[DIR_M00] = &DD[DIR_M00   *size_Mat];
		D.f[DIR_0P0] = &DD[DIR_0P0   *size_Mat];
		D.f[DIR_0M0] = &DD[DIR_0M0   *size_Mat];
		D.f[DIR_00P] = &DD[DIR_00P   *size_Mat];
		D.f[DIR_00M] = &DD[DIR_00M   *size_Mat];
		D.f[DIR_PP0] = &DD[DIR_PP0  *size_Mat];
		D.f[DIR_MM0] = &DD[DIR_MM0  *size_Mat];
		D.f[DIR_PM0] = &DD[DIR_PM0  *size_Mat];
		D.f[DIR_MP0] = &DD[DIR_MP0  *size_Mat];
		D.f[DIR_P0P] = &DD[DIR_P0P  *size_Mat];
		D.f[DIR_M0M] = &DD[DIR_M0M  *size_Mat];
		D.f[DIR_P0M] = &DD[DIR_P0M  *size_Mat];
		D.f[DIR_M0P] = &DD[DIR_M0P  *size_Mat];
		D.f[DIR_0PP] = &DD[DIR_0PP  *size_Mat];
		D.f[DIR_0MM] = &DD[DIR_0MM  *size_Mat];
		D.f[DIR_0PM] = &DD[DIR_0PM  *size_Mat];
		D.f[DIR_0MP] = &DD[DIR_0MP  *size_Mat];
		D.f[DIR_000] = &DD[DIR_000*size_Mat];
		D.f[DIR_PPP] = &DD[DIR_PPP *size_Mat];
		D.f[DIR_MMP] = &DD[DIR_MMP *size_Mat];
		D.f[DIR_PMP] = &DD[DIR_PMP *size_Mat];
		D.f[DIR_MPP] = &DD[DIR_MPP *size_Mat];
		D.f[DIR_PPM] = &DD[DIR_PPM *size_Mat];
		D.f[DIR_MMM] = &DD[DIR_MMM *size_Mat];
		D.f[DIR_PMM] = &DD[DIR_PMM *size_Mat];
		D.f[DIR_MPM] = &DD[DIR_MPM *size_Mat];
	}
	else
	{
		D.f[DIR_M00] = &DD[DIR_P00   *size_Mat];
		D.f[DIR_P00] = &DD[DIR_M00   *size_Mat];
		D.f[DIR_0M0] = &DD[DIR_0P0   *size_Mat];
		D.f[DIR_0P0] = &DD[DIR_0M0   *size_Mat];
		D.f[DIR_00M] = &DD[DIR_00P   *size_Mat];
		D.f[DIR_00P] = &DD[DIR_00M   *size_Mat];
		D.f[DIR_MM0] = &DD[DIR_PP0  *size_Mat];
		D.f[DIR_PP0] = &DD[DIR_MM0  *size_Mat];
		D.f[DIR_MP0] = &DD[DIR_PM0  *size_Mat];
		D.f[DIR_PM0] = &DD[DIR_MP0  *size_Mat];
		D.f[DIR_M0M] = &DD[DIR_P0P  *size_Mat];
		D.f[DIR_P0P] = &DD[DIR_M0M  *size_Mat];
		D.f[DIR_M0P] = &DD[DIR_P0M  *size_Mat];
		D.f[DIR_P0M] = &DD[DIR_M0P  *size_Mat];
		D.f[DIR_0MM] = &DD[DIR_0PP  *size_Mat];
		D.f[DIR_0PP] = &DD[DIR_0MM  *size_Mat];
		D.f[DIR_0MP] = &DD[DIR_0PM  *size_Mat];
		D.f[DIR_0PM] = &DD[DIR_0MP  *size_Mat];
		D.f[DIR_000] = &DD[DIR_000*size_Mat];
		D.f[DIR_PPP] = &DD[DIR_MMM *size_Mat];
		D.f[DIR_MMP] = &DD[DIR_PPM *size_Mat];
		D.f[DIR_PMP] = &DD[DIR_MPM *size_Mat];
		D.f[DIR_MPP] = &DD[DIR_PMM *size_Mat];
		D.f[DIR_PPM] = &DD[DIR_MMP *size_Mat];
		D.f[DIR_MMM] = &DD[DIR_PPP *size_Mat];
		D.f[DIR_PMM] = &DD[DIR_MPP *size_Mat];
		D.f[DIR_MPM] = &DD[DIR_PMP *size_Mat];
	}
	////////////////////////////////////////////////////////////////////////////////
	Distributions27 Dad;
	if (isEvenTimestep == true)
	{
		Dad.f[DIR_P00]    = &DD_AD[DIR_P00   *size_Mat];
		Dad.f[DIR_M00]    = &DD_AD[DIR_M00   *size_Mat];
		Dad.f[DIR_0P0]    = &DD_AD[DIR_0P0   *size_Mat];
		Dad.f[DIR_0M0]    = &DD_AD[DIR_0M0   *size_Mat];
		Dad.f[DIR_00P]    = &DD_AD[DIR_00P   *size_Mat];
		Dad.f[DIR_00M]    = &DD_AD[DIR_00M   *size_Mat];
		Dad.f[DIR_PP0]   = &DD_AD[DIR_PP0  *size_Mat];
		Dad.f[DIR_MM0]   = &DD_AD[DIR_MM0  *size_Mat];
		Dad.f[DIR_PM0]   = &DD_AD[DIR_PM0  *size_Mat];
		Dad.f[DIR_MP0]   = &DD_AD[DIR_MP0  *size_Mat];
		Dad.f[DIR_P0P]   = &DD_AD[DIR_P0P  *size_Mat];
		Dad.f[DIR_M0M]   = &DD_AD[DIR_M0M  *size_Mat];
		Dad.f[DIR_P0M]   = &DD_AD[DIR_P0M  *size_Mat];
		Dad.f[DIR_M0P]   = &DD_AD[DIR_M0P  *size_Mat];
		Dad.f[DIR_0PP]   = &DD_AD[DIR_0PP  *size_Mat];
		Dad.f[DIR_0MM]   = &DD_AD[DIR_0MM  *size_Mat];
		Dad.f[DIR_0PM]   = &DD_AD[DIR_0PM  *size_Mat];
		Dad.f[DIR_0MP]   = &DD_AD[DIR_0MP  *size_Mat];
		Dad.f[DIR_000] = &DD_AD[DIR_000*size_Mat];
		Dad.f[DIR_PPP]  = &DD_AD[DIR_PPP *size_Mat];
		Dad.f[DIR_MMP]  = &DD_AD[DIR_MMP *size_Mat];
		Dad.f[DIR_PMP]  = &DD_AD[DIR_PMP *size_Mat];
		Dad.f[DIR_MPP]  = &DD_AD[DIR_MPP *size_Mat];
		Dad.f[DIR_PPM]  = &DD_AD[DIR_PPM *size_Mat];
		Dad.f[DIR_MMM]  = &DD_AD[DIR_MMM *size_Mat];
		Dad.f[DIR_PMM]  = &DD_AD[DIR_PMM *size_Mat];
		Dad.f[DIR_MPM]  = &DD_AD[DIR_MPM *size_Mat];
	}						
	else					
	{						
		Dad.f[DIR_M00]    = &DD_AD[DIR_P00   *size_Mat];
		Dad.f[DIR_P00]    = &DD_AD[DIR_M00   *size_Mat];
		Dad.f[DIR_0M0]    = &DD_AD[DIR_0P0   *size_Mat];
		Dad.f[DIR_0P0]    = &DD_AD[DIR_0M0   *size_Mat];
		Dad.f[DIR_00M]    = &DD_AD[DIR_00P   *size_Mat];
		Dad.f[DIR_00P]    = &DD_AD[DIR_00M   *size_Mat];
		Dad.f[DIR_MM0]   = &DD_AD[DIR_PP0  *size_Mat];
		Dad.f[DIR_PP0]   = &DD_AD[DIR_MM0  *size_Mat];
		Dad.f[DIR_MP0]   = &DD_AD[DIR_PM0  *size_Mat];
		Dad.f[DIR_PM0]   = &DD_AD[DIR_MP0  *size_Mat];
		Dad.f[DIR_M0M]   = &DD_AD[DIR_P0P  *size_Mat];
		Dad.f[DIR_P0P]   = &DD_AD[DIR_M0M  *size_Mat];
		Dad.f[DIR_M0P]   = &DD_AD[DIR_P0M  *size_Mat];
		Dad.f[DIR_P0M]   = &DD_AD[DIR_M0P  *size_Mat];
		Dad.f[DIR_0MM]   = &DD_AD[DIR_0PP  *size_Mat];
		Dad.f[DIR_0PP]   = &DD_AD[DIR_0MM  *size_Mat];
		Dad.f[DIR_0MP]   = &DD_AD[DIR_0PM  *size_Mat];
		Dad.f[DIR_0PM]   = &DD_AD[DIR_0MP  *size_Mat];
		Dad.f[DIR_000] = &DD_AD[DIR_000*size_Mat];
		Dad.f[DIR_PPP]  = &DD_AD[DIR_MMM *size_Mat];
		Dad.f[DIR_MMP]  = &DD_AD[DIR_PPM *size_Mat];
		Dad.f[DIR_PMP]  = &DD_AD[DIR_MPM *size_Mat];
		Dad.f[DIR_MPP]  = &DD_AD[DIR_PMM *size_Mat];
		Dad.f[DIR_PPM]  = &DD_AD[DIR_MMP *size_Mat];
		Dad.f[DIR_MMM]  = &DD_AD[DIR_PPP *size_Mat];
		Dad.f[DIR_PMM]  = &DD_AD[DIR_MPP *size_Mat];
		Dad.f[DIR_MPM]  = &DD_AD[DIR_PMP *size_Mat];
	}
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if (k < size_Mat)
	{
		//////////////////////////////////////////////////////////////////////////
		//index
		//unsigned int kzero = k;
		unsigned int ke = k;
		unsigned int kw = neighborX[k];
		unsigned int kn = k;
		unsigned int ks = neighborY[k];
		unsigned int kt = k;
		unsigned int kb = neighborZ[k];
		unsigned int ksw = neighborY[kw];
		unsigned int kne = k;
		unsigned int kse = ks;
		unsigned int knw = kw;
		unsigned int kbw = neighborZ[kw];
		unsigned int kte = k;
		unsigned int kbe = kb;
		unsigned int ktw = kw;
		unsigned int kbs = neighborZ[ks];
		unsigned int ktn = k;
		unsigned int kbn = kb;
		unsigned int kts = ks;
		unsigned int ktse = ks;
		unsigned int kbnw = kbw;
		unsigned int ktnw = kw;
		unsigned int kbse = kbs;
		unsigned int ktsw = ksw;
		unsigned int kbne = kb;
		unsigned int ktne = k;
		unsigned int kbsw = neighborZ[ksw];
		//////////////////////////////////////////////////////////////////////////
		real CONC  = concD[k];
		real PRESS = pressD[k];
		real RHO   = rhoD[k];
		real VX    = vxD[k];
		real VY    = vyD[k];
		real VZ    = vzD[k];
		//////////////////////////////////////////////////////////////////////////
		concD[k] = c0o1;
		pressD[k] = c0o1;
		rhoD[k] = c0o1;
		vxD[k] = c0o1;
		vyD[k] = c0o1;
		vzD[k] = c0o1;

		if (geoD[k] == GEO_FLUID)
		{
			real mfcbb = (D.f[DIR_P00])[k];//[ke   ];
			real mfabb = (D.f[DIR_M00])[kw];//[kw   ];  
			real mfbcb = (D.f[DIR_0P0])[k];//[kn   ];
			real mfbab = (D.f[DIR_0M0])[ks];//[ks   ];  
			real mfbbc = (D.f[DIR_00P])[k];//[kt   ];
			real mfbba = (D.f[DIR_00M])[kb];//[kb   ];  
			real mfccb = (D.f[DIR_PP0])[k];//[kne  ];  
			real mfaab = (D.f[DIR_MM0])[ksw];//[ksw  ];
			real mfcab = (D.f[DIR_PM0])[ks];//[kse  ]; 
			real mfacb = (D.f[DIR_MP0])[kw];//[knw  ]; 
			real mfcbc = (D.f[DIR_P0P])[k];//[kte  ];  
			real mfaba = (D.f[DIR_M0M])[kbw];//[kbw  ];
			real mfcba = (D.f[DIR_P0M])[kb];//[kbe  ]; 
			real mfabc = (D.f[DIR_M0P])[kw];//[ktw  ]; 
			real mfbcc = (D.f[DIR_0PP])[k];//[ktn  ];  
			real mfbaa = (D.f[DIR_0MM])[kbs];//[kbs  ];
			real mfbca = (D.f[DIR_0PM])[kb];//[kbn  ]; 
			real mfbac = (D.f[DIR_0MP])[ks];//[kts  ]; 
			real mfbbb = (D.f[DIR_000])[k];//[kzero];
			real mfccc = (D.f[DIR_PPP])[k];//[ktne ]; 
			real mfaac = (D.f[DIR_MMP])[ksw];//[ktsw ]; 
			real mfcac = (D.f[DIR_PMP])[ks];//[ktse ];
			real mfacc = (D.f[DIR_MPP])[kw];//[ktnw ];
			real mfcca = (D.f[DIR_PPM])[kb];//[kbne ];
			real mfaaa = (D.f[DIR_MMM])[kbsw];//[kbsw ];
			real mfcaa = (D.f[DIR_PMM])[kbs];//[kbse ]; 
			real mfaca = (D.f[DIR_MPM])[kbw];//[kbnw ]; 
			////////////////////////////////////////////////////////////////////////////////////
			real drho =
				((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				 (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
				  ((mfabb + mfcbb) + (mfbab + mfbcb)  +  (mfbba + mfbbc))) + mfbbb;
			real rho = c1o1 + drho;
			////////////////////////////////////////////////////////////////////////////////////

			rhoD[k] = drho + RHO;

			vxD[k] =
				(((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
				(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
					(mfcbb - mfabb)) / rho) + VX;
			
			vyD[k] =
				(((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
				(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
					(mfbcb - mfbab)) / rho) + VY;
			
			vzD[k] =
				(((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
				(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
					(mfbbc - mfbba)) / rho) + VZ;

			pressD[k] = 
				((D.f[DIR_P00])[ke] + (D.f[DIR_M00])[kw] +
				 (D.f[DIR_0P0])[kn] + (D.f[DIR_0M0])[ks] +
				 (D.f[DIR_00P])[kt] + (D.f[DIR_00M])[kb] +
				 c2o1*(
				 (D.f[DIR_PP0])[kne] + (D.f[DIR_MM0])[ksw] +
				 (D.f[DIR_PM0])[kse] + (D.f[DIR_MP0])[knw] +
				 (D.f[DIR_P0P])[kte] + (D.f[DIR_M0M])[kbw] +
				 (D.f[DIR_P0M])[kbe] + (D.f[DIR_M0P])[ktw] +
				 (D.f[DIR_0PP])[ktn] + (D.f[DIR_0MM])[kbs] +
				 (D.f[DIR_0PM])[kbn] + (D.f[DIR_0MP])[kts]) +
				 c3o1*(
				 (D.f[DIR_PPP])[ktne] + (D.f[DIR_MMP])[ktsw] +
				 (D.f[DIR_PMP])[ktse] + (D.f[DIR_MPP])[ktnw] +
				 (D.f[DIR_PPM])[kbne] + (D.f[DIR_MMM])[kbsw] +
				 (D.f[DIR_PMM])[kbse] + (D.f[DIR_MPM])[kbnw]) -
				 rhoD[k] - (vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1 + rhoD[k])) * c1o2 + rhoD[k] +
				 PRESS;
				 //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
			//////////////////////////////////////////////////////////////////////////
			mfcbb = (Dad.f[DIR_P00   ])[k   ];
			mfabb = (Dad.f[DIR_M00   ])[kw  ];
			mfbcb = (Dad.f[DIR_0P0   ])[k   ];
			mfbab = (Dad.f[DIR_0M0   ])[ks  ];
			mfbbc = (Dad.f[DIR_00P   ])[k   ];
			mfbba = (Dad.f[DIR_00M   ])[kb  ];
			mfccb = (Dad.f[DIR_PP0  ])[k   ];
			mfaab = (Dad.f[DIR_MM0  ])[ksw ];
			mfcab = (Dad.f[DIR_PM0  ])[ks  ];
			mfacb = (Dad.f[DIR_MP0  ])[kw  ];
			mfcbc = (Dad.f[DIR_P0P  ])[k   ];
			mfaba = (Dad.f[DIR_M0M  ])[kbw ];
			mfcba = (Dad.f[DIR_P0M  ])[kb  ];
			mfabc = (Dad.f[DIR_M0P  ])[kw  ];
			mfbcc = (Dad.f[DIR_0PP  ])[k   ];
			mfbaa = (Dad.f[DIR_0MM  ])[kbs ];
			mfbca = (Dad.f[DIR_0PM  ])[kb  ];
			mfbac = (Dad.f[DIR_0MP  ])[ks  ];
			mfbbb = (Dad.f[DIR_000])[k   ];
			mfccc = (Dad.f[DIR_PPP ])[k   ];
			mfaac = (Dad.f[DIR_MMP ])[ksw ];
			mfcac = (Dad.f[DIR_PMP ])[ks  ];
			mfacc = (Dad.f[DIR_MPP ])[kw  ];
			mfcca = (Dad.f[DIR_PPM ])[kb  ];
			mfaaa = (Dad.f[DIR_MMM ])[kbsw];
			mfcaa = (Dad.f[DIR_PMM ])[kbs ];
			mfaca = (Dad.f[DIR_MPM ])[kbw ];
			//////////////////////////////////////////////////////////////////////////
			concD[k] = 
				((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa)   + (mfaac + mfcca))) +
				 (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba)   + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
				  ((mfabb + mfcbb) + (mfbab + mfbcb)  +  (mfbba + mfbbc))) +  mfbbb + CONC;
		}
	}
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMacMedSP27( real* vxD,
                                             real* vyD,
                                             real* vzD,
                                             real* rhoD,
                                             real* pressD,
                                             unsigned int* geoD,
                                             unsigned int* neighborX,
                                             unsigned int* neighborY,
                                             unsigned int* neighborZ,
                                             unsigned int tdiff,
                                             unsigned int size_Mat,
                                             bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
      real PRESS = pressD[k];
      real RHO   = rhoD[k];
      real VX    = vxD[k];
      real VY    = vyD[k];
      real VZ    = vzD[k];
      //////////////////////////////////////////////////////////////////////////
      pressD[k] = c0o1;
      rhoD[k]   = c0o1;
      vxD[k]    = c0o1;
      vyD[k]    = c0o1;
      vzD[k]    = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         rhoD[k]    =   RHO   / tdiff;
         vxD[k]     =   VX    / tdiff;
         vyD[k]     =   VY    / tdiff;
         vzD[k]     =   VZ    / tdiff;
         pressD[k]  =   PRESS / tdiff;    
      }
   }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
__global__ void LBResetMedianValuesSP27(
	real* vxD,
	real* vyD,
	real* vzD,
	real* rhoD,
	real* pressD,
	unsigned int size_Mat,
	bool isEvenTimestep)
{
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if (k<size_Mat)
	{
		//////////////////////////////////////////////////////////////////////////
		pressD[k] = c0o1;
		rhoD[k] = c0o1;
		vxD[k] = c0o1;
		vyD[k] = c0o1;
		vzD[k] = c0o1;
	}
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
__global__ void LBResetMedianValuesAD27(
	real* vxD,
	real* vyD,
	real* vzD,
	real* rhoD,
	real* pressD,
	real* concD,
	unsigned int size_Mat,
	bool isEvenTimestep)
{
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if (k < size_Mat)
	{
		//////////////////////////////////////////////////////////////////////////
		concD[k]  = c0o1;
		pressD[k] = c0o1;
		rhoD[k]   = c0o1;
		vxD[k]    = c0o1;
		vyD[k]    = c0o1;
		vzD[k]    = c0o1;
	}
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMeasurePoints( real* vxMP,
												real* vyMP,
												real* vzMP,
												real* rhoMP,
												unsigned int* kMP,
												unsigned int numberOfPointskMP,
												unsigned int MPClockCycle,
												unsigned int t,
												unsigned int* geoD,
												unsigned int* neighborX,
												unsigned int* neighborY,
												unsigned int* neighborZ,
												unsigned int size_Mat,
												real* DD,
												bool isEvenTimestep)
{
	Distributions27 D;
	if (isEvenTimestep==true)
	{
		D.f[DIR_P00   ] = &DD[DIR_P00   *size_Mat];
		D.f[DIR_M00   ] = &DD[DIR_M00   *size_Mat];
		D.f[DIR_0P0   ] = &DD[DIR_0P0   *size_Mat];
		D.f[DIR_0M0   ] = &DD[DIR_0M0   *size_Mat];
		D.f[DIR_00P   ] = &DD[DIR_00P   *size_Mat];
		D.f[DIR_00M   ] = &DD[DIR_00M   *size_Mat];
		D.f[DIR_PP0  ] = &DD[DIR_PP0  *size_Mat];
		D.f[DIR_MM0  ] = &DD[DIR_MM0  *size_Mat];
		D.f[DIR_PM0  ] = &DD[DIR_PM0  *size_Mat];
		D.f[DIR_MP0  ] = &DD[DIR_MP0  *size_Mat];
		D.f[DIR_P0P  ] = &DD[DIR_P0P  *size_Mat];
		D.f[DIR_M0M  ] = &DD[DIR_M0M  *size_Mat];
		D.f[DIR_P0M  ] = &DD[DIR_P0M  *size_Mat];
		D.f[DIR_M0P  ] = &DD[DIR_M0P  *size_Mat];
		D.f[DIR_0PP  ] = &DD[DIR_0PP  *size_Mat];
		D.f[DIR_0MM  ] = &DD[DIR_0MM  *size_Mat];
		D.f[DIR_0PM  ] = &DD[DIR_0PM  *size_Mat];
		D.f[DIR_0MP  ] = &DD[DIR_0MP  *size_Mat];
		D.f[DIR_000] = &DD[DIR_000*size_Mat];
		D.f[DIR_PPP ] = &DD[DIR_PPP *size_Mat];
		D.f[DIR_MMP ] = &DD[DIR_MMP *size_Mat];
		D.f[DIR_PMP ] = &DD[DIR_PMP *size_Mat];
		D.f[DIR_MPP ] = &DD[DIR_MPP *size_Mat];
		D.f[DIR_PPM ] = &DD[DIR_PPM *size_Mat];
		D.f[DIR_MMM ] = &DD[DIR_MMM *size_Mat];
		D.f[DIR_PMM ] = &DD[DIR_PMM *size_Mat];
		D.f[DIR_MPM ] = &DD[DIR_MPM *size_Mat];
	} 
	else
	{
		D.f[DIR_M00   ] = &DD[DIR_P00   *size_Mat];
		D.f[DIR_P00   ] = &DD[DIR_M00   *size_Mat];
		D.f[DIR_0M0   ] = &DD[DIR_0P0   *size_Mat];
		D.f[DIR_0P0   ] = &DD[DIR_0M0   *size_Mat];
		D.f[DIR_00M   ] = &DD[DIR_00P   *size_Mat];
		D.f[DIR_00P   ] = &DD[DIR_00M   *size_Mat];
		D.f[DIR_MM0  ] = &DD[DIR_PP0  *size_Mat];
		D.f[DIR_PP0  ] = &DD[DIR_MM0  *size_Mat];
		D.f[DIR_MP0  ] = &DD[DIR_PM0  *size_Mat];
		D.f[DIR_PM0  ] = &DD[DIR_MP0  *size_Mat];
		D.f[DIR_M0M  ] = &DD[DIR_P0P  *size_Mat];
		D.f[DIR_P0P  ] = &DD[DIR_M0M  *size_Mat];
		D.f[DIR_M0P  ] = &DD[DIR_P0M  *size_Mat];
		D.f[DIR_P0M  ] = &DD[DIR_M0P  *size_Mat];
		D.f[DIR_0MM  ] = &DD[DIR_0PP  *size_Mat];
		D.f[DIR_0PP  ] = &DD[DIR_0MM  *size_Mat];
		D.f[DIR_0MP  ] = &DD[DIR_0PM  *size_Mat];
		D.f[DIR_0PM  ] = &DD[DIR_0MP  *size_Mat];
		D.f[DIR_000] = &DD[DIR_000*size_Mat];
		D.f[DIR_PPP ] = &DD[DIR_MMM *size_Mat];
		D.f[DIR_MMP ] = &DD[DIR_PPM *size_Mat];
		D.f[DIR_PMP ] = &DD[DIR_MPM *size_Mat];
		D.f[DIR_MPP ] = &DD[DIR_PMM *size_Mat];
		D.f[DIR_PPM ] = &DD[DIR_MMP *size_Mat];
		D.f[DIR_MMM ] = &DD[DIR_PPP *size_Mat];
		D.f[DIR_PMM ] = &DD[DIR_MPP *size_Mat];
		D.f[DIR_MPM ] = &DD[DIR_PMP *size_Mat];
	}
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if(k<numberOfPointskMP)
	{
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= kMP[k];//k;
      unsigned int ke   = kzero;
      unsigned int kw   = neighborX[kzero];
      unsigned int kn   = kzero;
      unsigned int ks   = neighborY[kzero];
      unsigned int kt   = kzero;
      unsigned int kb   = neighborZ[kzero];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = kzero;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = kzero;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = kzero;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = kzero;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
	  unsigned int kMac = k*MPClockCycle + t;
	  //////////////////////////////////////////////////////////////////////////

      if(geoD[kzero] == GEO_FLUID)
      {
         rhoMP[kMac]=   (D.f[DIR_P00   ])[ke  ]+ (D.f[DIR_M00   ])[kw  ]+ 
                        (D.f[DIR_0P0   ])[kn  ]+ (D.f[DIR_0M0   ])[ks  ]+
                        (D.f[DIR_00P   ])[kt  ]+ (D.f[DIR_00M   ])[kb  ]+
                        (D.f[DIR_PP0  ])[kne ]+ (D.f[DIR_MM0  ])[ksw ]+
                        (D.f[DIR_PM0  ])[kse ]+ (D.f[DIR_MP0  ])[knw ]+
                        (D.f[DIR_P0P  ])[kte ]+ (D.f[DIR_M0M  ])[kbw ]+
                        (D.f[DIR_P0M  ])[kbe ]+ (D.f[DIR_M0P  ])[ktw ]+
                        (D.f[DIR_0PP  ])[ktn ]+ (D.f[DIR_0MM  ])[kbs ]+
                        (D.f[DIR_0PM  ])[kbn ]+ (D.f[DIR_0MP  ])[kts ]+
                        (D.f[DIR_000])[kzero]+ 
                        (D.f[DIR_PPP ])[ktne]+ (D.f[DIR_MMP ])[ktsw]+ 
                        (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]+ 
                        (D.f[DIR_PPM ])[kbne]+ (D.f[DIR_MMM ])[kbsw]+ 
                        (D.f[DIR_PMM ])[kbse]+ (D.f[DIR_MPM ])[kbnw];

         vxMP[kMac] =   (D.f[DIR_P00   ])[ke  ]- (D.f[DIR_M00   ])[kw  ]+ 
                        (D.f[DIR_PP0  ])[kne ]- (D.f[DIR_MM0  ])[ksw ]+
                        (D.f[DIR_PM0  ])[kse ]- (D.f[DIR_MP0  ])[knw ]+
                        (D.f[DIR_P0P  ])[kte ]- (D.f[DIR_M0M  ])[kbw ]+
                        (D.f[DIR_P0M  ])[kbe ]- (D.f[DIR_M0P  ])[ktw ]+
                        (D.f[DIR_PPP ])[ktne]- (D.f[DIR_MMP ])[ktsw]+ 
                        (D.f[DIR_PMP ])[ktse]- (D.f[DIR_MPP ])[ktnw]+ 
                        (D.f[DIR_PPM ])[kbne]- (D.f[DIR_MMM ])[kbsw]+ 
                        (D.f[DIR_PMM ])[kbse]- (D.f[DIR_MPM ])[kbnw];

         vyMP[kMac] =   (D.f[DIR_0P0   ])[kn  ]- (D.f[DIR_0M0   ])[ks  ]+
                        (D.f[DIR_PP0  ])[kne ]- (D.f[DIR_MM0  ])[ksw ]-
                        (D.f[DIR_PM0  ])[kse ]+ (D.f[DIR_MP0  ])[knw ]+
                        (D.f[DIR_0PP  ])[ktn ]- (D.f[DIR_0MM  ])[kbs ]+
                        (D.f[DIR_0PM  ])[kbn ]- (D.f[DIR_0MP  ])[kts ]+
                        (D.f[DIR_PPP ])[ktne]- (D.f[DIR_MMP ])[ktsw]- 
                        (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]+ 
                        (D.f[DIR_PPM ])[kbne]- (D.f[DIR_MMM ])[kbsw]- 
                        (D.f[DIR_PMM ])[kbse]+ (D.f[DIR_MPM ])[kbnw];

         vzMP[kMac] =   (D.f[DIR_00P   ])[kt  ]- (D.f[DIR_00M   ])[kb  ]+
                        (D.f[DIR_P0P  ])[kte ]- (D.f[DIR_M0M  ])[kbw ]-
                        (D.f[DIR_P0M  ])[kbe ]+ (D.f[DIR_M0P  ])[ktw ]+
                        (D.f[DIR_0PP  ])[ktn ]- (D.f[DIR_0MM  ])[kbs ]-
                        (D.f[DIR_0PM  ])[kbn ]+ (D.f[DIR_0MP  ])[kts ]+
                        (D.f[DIR_PPP ])[ktne]+ (D.f[DIR_MMP ])[ktsw]+ 
                        (D.f[DIR_PMP ])[ktse]+ (D.f[DIR_MPP ])[ktnw]- 
                        (D.f[DIR_PPM ])[kbne]- (D.f[DIR_MMM ])[kbsw]- 
                        (D.f[DIR_PMM ])[kbse]- (D.f[DIR_MPM ])[kbnw];
      }
   }
}
////////////////////////////////////////////////////////////////////////////////





































////////////////////////////////////////////////////////////////////////////////
__global__ void LBSetOutputWallVelocitySP27( real* vxD,
														real* vyD,
														real* vzD,
														real* vxWall,
														real* vyWall,
														real* vzWall,
														int numberOfWallNodes, 
														int* kWallNodes, 
														real* rhoD,
														real* pressD,
														unsigned int* geoD,
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														unsigned int size_Mat,
														real* DD,
														bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfWallNodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KWN  = kWallNodes[k];
      //////////////////////////////////////////////////////////////////////////
      vxD[KWN] = 0.0;//vxWall[k];
      vyD[KWN] = 0.0;//vyWall[k];
      vzD[KWN] = 0.0;//vzWall[k];
   }
}





























