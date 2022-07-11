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
extern "C" __global__ void LBCalcMac27( real* vxD,
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
extern "C" __global__ void LBCalcMacSP27( real* vxD,
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
      D.f[E   ] = &DD[E   *size_Mat];
      D.f[W   ] = &DD[W   *size_Mat];
      D.f[N   ] = &DD[N   *size_Mat];
      D.f[S   ] = &DD[S   *size_Mat];
      D.f[T   ] = &DD[T   *size_Mat];
      D.f[B   ] = &DD[B   *size_Mat];
      D.f[NE  ] = &DD[NE  *size_Mat];
      D.f[SW  ] = &DD[SW  *size_Mat];
      D.f[SE  ] = &DD[SE  *size_Mat];
      D.f[NW  ] = &DD[NW  *size_Mat];
      D.f[TE  ] = &DD[TE  *size_Mat];
      D.f[BW  ] = &DD[BW  *size_Mat];
      D.f[BE  ] = &DD[BE  *size_Mat];
      D.f[TW  ] = &DD[TW  *size_Mat];
      D.f[TN  ] = &DD[TN  *size_Mat];
      D.f[BS  ] = &DD[BS  *size_Mat];
      D.f[BN  ] = &DD[BN  *size_Mat];
      D.f[TS  ] = &DD[TS  *size_Mat];
      D.f[REST] = &DD[REST*size_Mat];
      D.f[TNE ] = &DD[TNE *size_Mat];
      D.f[TSW ] = &DD[TSW *size_Mat];
      D.f[TSE ] = &DD[TSE *size_Mat];
      D.f[TNW ] = &DD[TNW *size_Mat];
      D.f[BNE ] = &DD[BNE *size_Mat];
      D.f[BSW ] = &DD[BSW *size_Mat];
      D.f[BSE ] = &DD[BSE *size_Mat];
      D.f[BNW ] = &DD[BNW *size_Mat];
   } 
   else
   {
      D.f[W   ] = &DD[E   *size_Mat];
      D.f[E   ] = &DD[W   *size_Mat];
      D.f[S   ] = &DD[N   *size_Mat];
      D.f[N   ] = &DD[S   *size_Mat];
      D.f[B   ] = &DD[T   *size_Mat];
      D.f[T   ] = &DD[B   *size_Mat];
      D.f[SW  ] = &DD[NE  *size_Mat];
      D.f[NE  ] = &DD[SW  *size_Mat];
      D.f[NW  ] = &DD[SE  *size_Mat];
      D.f[SE  ] = &DD[NW  *size_Mat];
      D.f[BW  ] = &DD[TE  *size_Mat];
      D.f[TE  ] = &DD[BW  *size_Mat];
      D.f[TW  ] = &DD[BE  *size_Mat];
      D.f[BE  ] = &DD[TW  *size_Mat];
      D.f[BS  ] = &DD[TN  *size_Mat];
      D.f[TN  ] = &DD[BS  *size_Mat];
      D.f[TS  ] = &DD[BN  *size_Mat];
      D.f[BN  ] = &DD[TS  *size_Mat];
      D.f[REST] = &DD[REST*size_Mat];
      D.f[TNE ] = &DD[BSW *size_Mat];
      D.f[TSW ] = &DD[BNE *size_Mat];
      D.f[TSE ] = &DD[BNW *size_Mat];
      D.f[TNW ] = &DD[BSE *size_Mat];
      D.f[BNE ] = &DD[TSW *size_Mat];
      D.f[BSW ] = &DD[TNE *size_Mat];
      D.f[BSE ] = &DD[TNW *size_Mat];
      D.f[BNW ] = &DD[TSE *size_Mat];
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
         rhoD[k]    =   (D.f[E   ])[ke  ]+ (D.f[W   ])[kw  ]+ 
                        (D.f[N   ])[kn  ]+ (D.f[S   ])[ks  ]+
                        (D.f[T   ])[kt  ]+ (D.f[B   ])[kb  ]+
                        (D.f[NE  ])[kne ]+ (D.f[SW  ])[ksw ]+
                        (D.f[SE  ])[kse ]+ (D.f[NW  ])[knw ]+
                        (D.f[TE  ])[kte ]+ (D.f[BW  ])[kbw ]+
                        (D.f[BE  ])[kbe ]+ (D.f[TW  ])[ktw ]+
                        (D.f[TN  ])[ktn ]+ (D.f[BS  ])[kbs ]+
                        (D.f[BN  ])[kbn ]+ (D.f[TS  ])[kts ]+
                        (D.f[REST])[kzero]+ 
                        (D.f[TNE ])[ktne]+ (D.f[TSW ])[ktsw]+ 
                        (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]+ 
                        (D.f[BNE ])[kbne]+ (D.f[BSW ])[kbsw]+ 
                        (D.f[BSE ])[kbse]+ (D.f[BNW ])[kbnw];

         vxD[k]     =   (D.f[E   ])[ke  ]- (D.f[W   ])[kw  ]+ 
                        (D.f[NE  ])[kne ]- (D.f[SW  ])[ksw ]+
                        (D.f[SE  ])[kse ]- (D.f[NW  ])[knw ]+
                        (D.f[TE  ])[kte ]- (D.f[BW  ])[kbw ]+
                        (D.f[BE  ])[kbe ]- (D.f[TW  ])[ktw ]+
                        (D.f[TNE ])[ktne]- (D.f[TSW ])[ktsw]+ 
                        (D.f[TSE ])[ktse]- (D.f[TNW ])[ktnw]+ 
                        (D.f[BNE ])[kbne]- (D.f[BSW ])[kbsw]+ 
                        (D.f[BSE ])[kbse]- (D.f[BNW ])[kbnw];

         vyD[k]     =   (D.f[N   ])[kn  ]- (D.f[S   ])[ks  ]+
                        (D.f[NE  ])[kne ]- (D.f[SW  ])[ksw ]-
                        (D.f[SE  ])[kse ]+ (D.f[NW  ])[knw ]+
                        (D.f[TN  ])[ktn ]- (D.f[BS  ])[kbs ]+
                        (D.f[BN  ])[kbn ]- (D.f[TS  ])[kts ]+
                        (D.f[TNE ])[ktne]- (D.f[TSW ])[ktsw]- 
                        (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]+ 
                        (D.f[BNE ])[kbne]- (D.f[BSW ])[kbsw]- 
                        (D.f[BSE ])[kbse]+ (D.f[BNW ])[kbnw];

         vzD[k]     =   (D.f[T   ])[kt  ]- (D.f[B   ])[kb  ]+
                        (D.f[TE  ])[kte ]- (D.f[BW  ])[kbw ]-
                        (D.f[BE  ])[kbe ]+ (D.f[TW  ])[ktw ]+
                        (D.f[TN  ])[ktn ]- (D.f[BS  ])[kbs ]-
                        (D.f[BN  ])[kbn ]+ (D.f[TS  ])[kts ]+
                        (D.f[TNE ])[ktne]+ (D.f[TSW ])[ktsw]+ 
                        (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]- 
                        (D.f[BNE ])[kbne]- (D.f[BSW ])[kbsw]- 
                        (D.f[BSE ])[kbse]- (D.f[BNW ])[kbnw];

         pressD[k]  =  ((D.f[E   ])[ke  ]+ (D.f[W   ])[kw  ]+ 
                        (D.f[N   ])[kn  ]+ (D.f[S   ])[ks  ]+
                        (D.f[T   ])[kt  ]+ (D.f[B   ])[kb  ]+
                        2.f*(
                        (D.f[NE  ])[kne ]+ (D.f[SW  ])[ksw ]+
                        (D.f[SE  ])[kse ]+ (D.f[NW  ])[knw ]+
                        (D.f[TE  ])[kte ]+ (D.f[BW  ])[kbw ]+
                        (D.f[BE  ])[kbe ]+ (D.f[TW  ])[ktw ]+
                        (D.f[TN  ])[ktn ]+ (D.f[BS  ])[kbs ]+
                        (D.f[BN  ])[kbn ]+ (D.f[TS  ])[kts ])+
                        3.f*(
                        (D.f[TNE ])[ktne]+ (D.f[TSW ])[ktsw]+ 
                        (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]+ 
                        (D.f[BNE ])[kbne]+ (D.f[BSW ])[kbsw]+ 
                        (D.f[BSE ])[kbse]+ (D.f[BNW ])[kbnw])-
                        rhoD[k]-(vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1+c0o1*rhoD[k])) * c1o2+rhoD[k]; // times zero for incompressible case   
         //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5

      }
   }
}


////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMacCompSP27(real *vxD, real *vyD, real *vzD, real *rhoD, real *pressD,
                                             unsigned int *geoD, unsigned int *neighborX, unsigned int *neighborY,
                                             unsigned int *neighborZ, unsigned int size_Mat, real *distributions,
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
extern "C" __global__ void LBCalcMedSP27( real* vxD,
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
      D.f[E   ] = &DD[E   *size_Mat];
      D.f[W   ] = &DD[W   *size_Mat];
      D.f[N   ] = &DD[N   *size_Mat];
      D.f[S   ] = &DD[S   *size_Mat];
      D.f[T   ] = &DD[T   *size_Mat];
      D.f[B   ] = &DD[B   *size_Mat];
      D.f[NE  ] = &DD[NE  *size_Mat];
      D.f[SW  ] = &DD[SW  *size_Mat];
      D.f[SE  ] = &DD[SE  *size_Mat];
      D.f[NW  ] = &DD[NW  *size_Mat];
      D.f[TE  ] = &DD[TE  *size_Mat];
      D.f[BW  ] = &DD[BW  *size_Mat];
      D.f[BE  ] = &DD[BE  *size_Mat];
      D.f[TW  ] = &DD[TW  *size_Mat];
      D.f[TN  ] = &DD[TN  *size_Mat];
      D.f[BS  ] = &DD[BS  *size_Mat];
      D.f[BN  ] = &DD[BN  *size_Mat];
      D.f[TS  ] = &DD[TS  *size_Mat];
      D.f[REST] = &DD[REST*size_Mat];
      D.f[TNE ] = &DD[TNE *size_Mat];
      D.f[TSW ] = &DD[TSW *size_Mat];
      D.f[TSE ] = &DD[TSE *size_Mat];
      D.f[TNW ] = &DD[TNW *size_Mat];
      D.f[BNE ] = &DD[BNE *size_Mat];
      D.f[BSW ] = &DD[BSW *size_Mat];
      D.f[BSE ] = &DD[BSE *size_Mat];
      D.f[BNW ] = &DD[BNW *size_Mat];
   } 
   else
   {
      D.f[W   ] = &DD[E   *size_Mat];
      D.f[E   ] = &DD[W   *size_Mat];
      D.f[S   ] = &DD[N   *size_Mat];
      D.f[N   ] = &DD[S   *size_Mat];
      D.f[B   ] = &DD[T   *size_Mat];
      D.f[T   ] = &DD[B   *size_Mat];
      D.f[SW  ] = &DD[NE  *size_Mat];
      D.f[NE  ] = &DD[SW  *size_Mat];
      D.f[NW  ] = &DD[SE  *size_Mat];
      D.f[SE  ] = &DD[NW  *size_Mat];
      D.f[BW  ] = &DD[TE  *size_Mat];
      D.f[TE  ] = &DD[BW  *size_Mat];
      D.f[TW  ] = &DD[BE  *size_Mat];
      D.f[BE  ] = &DD[TW  *size_Mat];
      D.f[BS  ] = &DD[TN  *size_Mat];
      D.f[TN  ] = &DD[BS  *size_Mat];
      D.f[TS  ] = &DD[BN  *size_Mat];
      D.f[BN  ] = &DD[TS  *size_Mat];
      D.f[REST] = &DD[REST*size_Mat];
      D.f[TNE ] = &DD[BSW *size_Mat];
      D.f[TSW ] = &DD[BNE *size_Mat];
      D.f[TSE ] = &DD[BNW *size_Mat];
      D.f[TNW ] = &DD[BSE *size_Mat];
      D.f[BNE ] = &DD[TSW *size_Mat];
      D.f[BSW ] = &DD[TNE *size_Mat];
      D.f[BSE ] = &DD[TNW *size_Mat];
      D.f[BNW ] = &DD[TSE *size_Mat];
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
         rhoD[k]    =   (D.f[E   ])[ke  ]+ (D.f[W   ])[kw  ]+ 
                        (D.f[N   ])[kn  ]+ (D.f[S   ])[ks  ]+
                        (D.f[T   ])[kt  ]+ (D.f[B   ])[kb  ]+
                        (D.f[NE  ])[kne ]+ (D.f[SW  ])[ksw ]+
                        (D.f[SE  ])[kse ]+ (D.f[NW  ])[knw ]+
                        (D.f[TE  ])[kte ]+ (D.f[BW  ])[kbw ]+
                        (D.f[BE  ])[kbe ]+ (D.f[TW  ])[ktw ]+
                        (D.f[TN  ])[ktn ]+ (D.f[BS  ])[kbs ]+
                        (D.f[BN  ])[kbn ]+ (D.f[TS  ])[kts ]+
                        (D.f[REST])[kzero]+ 
                        (D.f[TNE ])[ktne]+ (D.f[TSW ])[ktsw]+ 
                        (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]+ 
                        (D.f[BNE ])[kbne]+ (D.f[BSW ])[kbsw]+ 
                        (D.f[BSE ])[kbse]+ (D.f[BNW ])[kbnw]+
                        RHO;

         vxD[k]     =   (D.f[E   ])[ke  ]- (D.f[W   ])[kw  ]+ 
                        (D.f[NE  ])[kne ]- (D.f[SW  ])[ksw ]+
                        (D.f[SE  ])[kse ]- (D.f[NW  ])[knw ]+
                        (D.f[TE  ])[kte ]- (D.f[BW  ])[kbw ]+
                        (D.f[BE  ])[kbe ]- (D.f[TW  ])[ktw ]+
                        (D.f[TNE ])[ktne]- (D.f[TSW ])[ktsw]+ 
                        (D.f[TSE ])[ktse]- (D.f[TNW ])[ktnw]+ 
                        (D.f[BNE ])[kbne]- (D.f[BSW ])[kbsw]+ 
                        (D.f[BSE ])[kbse]- (D.f[BNW ])[kbnw]+
                        VX;

         vyD[k]     =   (D.f[N   ])[kn  ]- (D.f[S   ])[ks  ]+
                        (D.f[NE  ])[kne ]- (D.f[SW  ])[ksw ]-
                        (D.f[SE  ])[kse ]+ (D.f[NW  ])[knw ]+
                        (D.f[TN  ])[ktn ]- (D.f[BS  ])[kbs ]+
                        (D.f[BN  ])[kbn ]- (D.f[TS  ])[kts ]+
                        (D.f[TNE ])[ktne]- (D.f[TSW ])[ktsw]- 
                        (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]+ 
                        (D.f[BNE ])[kbne]- (D.f[BSW ])[kbsw]- 
                        (D.f[BSE ])[kbse]+ (D.f[BNW ])[kbnw]+
                        VY;

         vzD[k]     =   (D.f[T   ])[kt  ]- (D.f[B   ])[kb  ]+
                        (D.f[TE  ])[kte ]- (D.f[BW  ])[kbw ]-
                        (D.f[BE  ])[kbe ]+ (D.f[TW  ])[ktw ]+
                        (D.f[TN  ])[ktn ]- (D.f[BS  ])[kbs ]-
                        (D.f[BN  ])[kbn ]+ (D.f[TS  ])[kts ]+
                        (D.f[TNE ])[ktne]+ (D.f[TSW ])[ktsw]+ 
                        (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]- 
                        (D.f[BNE ])[kbne]- (D.f[BSW ])[kbsw]- 
                        (D.f[BSE ])[kbse]- (D.f[BNW ])[kbnw]+
                        VZ;

         pressD[k]  =   ((D.f[E   ])[ke  ]+ (D.f[W   ])[kw  ]+ 
                        (D.f[N   ])[kn  ]+ (D.f[S   ])[ks  ]+
                        (D.f[T   ])[kt  ]+ (D.f[B   ])[kb  ]+
                        c2o1*(
                        (D.f[NE  ])[kne ]+ (D.f[SW  ])[ksw ]+
                        (D.f[SE  ])[kse ]+ (D.f[NW  ])[knw ]+
                        (D.f[TE  ])[kte ]+ (D.f[BW  ])[kbw ]+
                        (D.f[BE  ])[kbe ]+ (D.f[TW  ])[ktw ]+
                        (D.f[TN  ])[ktn ]+ (D.f[BS  ])[kbs ]+
                        (D.f[BN  ])[kbn ]+ (D.f[TS  ])[kts ])+
                        c3o1*(
                        (D.f[TNE ])[ktne]+ (D.f[TSW ])[ktsw]+ 
                        (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]+ 
                        (D.f[BNE ])[kbne]+ (D.f[BSW ])[kbsw]+ 
                        (D.f[BSE ])[kbse]+ (D.f[BNW ])[kbnw])-
                        rhoD[k]-(vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1+rhoD[k])) * c1o2+rhoD[k]+
                        PRESS;    
         //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
      }
   }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMedCompSP27( real* vxD,
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
      D.f[E   ] = &DD[E   *size_Mat];
      D.f[W   ] = &DD[W   *size_Mat];
      D.f[N   ] = &DD[N   *size_Mat];
      D.f[S   ] = &DD[S   *size_Mat];
      D.f[T   ] = &DD[T   *size_Mat];
      D.f[B   ] = &DD[B   *size_Mat];
      D.f[NE  ] = &DD[NE  *size_Mat];
      D.f[SW  ] = &DD[SW  *size_Mat];
      D.f[SE  ] = &DD[SE  *size_Mat];
      D.f[NW  ] = &DD[NW  *size_Mat];
      D.f[TE  ] = &DD[TE  *size_Mat];
      D.f[BW  ] = &DD[BW  *size_Mat];
      D.f[BE  ] = &DD[BE  *size_Mat];
      D.f[TW  ] = &DD[TW  *size_Mat];
      D.f[TN  ] = &DD[TN  *size_Mat];
      D.f[BS  ] = &DD[BS  *size_Mat];
      D.f[BN  ] = &DD[BN  *size_Mat];
      D.f[TS  ] = &DD[TS  *size_Mat];
      D.f[REST] = &DD[REST*size_Mat];
      D.f[TNE ] = &DD[TNE *size_Mat];
      D.f[TSW ] = &DD[TSW *size_Mat];
      D.f[TSE ] = &DD[TSE *size_Mat];
      D.f[TNW ] = &DD[TNW *size_Mat];
      D.f[BNE ] = &DD[BNE *size_Mat];
      D.f[BSW ] = &DD[BSW *size_Mat];
      D.f[BSE ] = &DD[BSE *size_Mat];
      D.f[BNW ] = &DD[BNW *size_Mat];
   } 
   else
   {
      D.f[W   ] = &DD[E   *size_Mat];
      D.f[E   ] = &DD[W   *size_Mat];
      D.f[S   ] = &DD[N   *size_Mat];
      D.f[N   ] = &DD[S   *size_Mat];
      D.f[B   ] = &DD[T   *size_Mat];
      D.f[T   ] = &DD[B   *size_Mat];
      D.f[SW  ] = &DD[NE  *size_Mat];
      D.f[NE  ] = &DD[SW  *size_Mat];
      D.f[NW  ] = &DD[SE  *size_Mat];
      D.f[SE  ] = &DD[NW  *size_Mat];
      D.f[BW  ] = &DD[TE  *size_Mat];
      D.f[TE  ] = &DD[BW  *size_Mat];
      D.f[TW  ] = &DD[BE  *size_Mat];
      D.f[BE  ] = &DD[TW  *size_Mat];
      D.f[BS  ] = &DD[TN  *size_Mat];
      D.f[TN  ] = &DD[BS  *size_Mat];
      D.f[TS  ] = &DD[BN  *size_Mat];
      D.f[BN  ] = &DD[TS  *size_Mat];
      D.f[REST] = &DD[REST*size_Mat];
      D.f[TNE ] = &DD[BSW *size_Mat];
      D.f[TSW ] = &DD[BNE *size_Mat];
      D.f[TSE ] = &DD[BNW *size_Mat];
      D.f[TNW ] = &DD[BSE *size_Mat];
      D.f[BNE ] = &DD[TSW *size_Mat];
      D.f[BSW ] = &DD[TNE *size_Mat];
      D.f[BSE ] = &DD[TNW *size_Mat];
      D.f[BNW ] = &DD[TSE *size_Mat];
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
		  real mfcbb = (D.f[E])[k];//[ke   ];
		  real mfabb = (D.f[W])[kw];//[kw   ];  
		  real mfbcb = (D.f[N])[k];//[kn   ];
		  real mfbab = (D.f[S])[ks];//[ks   ];  
		  real mfbbc = (D.f[T])[k];//[kt   ];
		  real mfbba = (D.f[B])[kb];//[kb   ];  
		  real mfccb = (D.f[NE])[k];//[kne  ];  
		  real mfaab = (D.f[SW])[ksw];//[ksw  ];
		  real mfcab = (D.f[SE])[ks];//[kse  ]; 
		  real mfacb = (D.f[NW])[kw];//[knw  ]; 
		  real mfcbc = (D.f[TE])[k];//[kte  ];  
		  real mfaba = (D.f[BW])[kbw];//[kbw  ];
		  real mfcba = (D.f[BE])[kb];//[kbe  ]; 
		  real mfabc = (D.f[TW])[kw];//[ktw  ]; 
		  real mfbcc = (D.f[TN])[k];//[ktn  ];  
		  real mfbaa = (D.f[BS])[kbs];//[kbs  ];
		  real mfbca = (D.f[BN])[kb];//[kbn  ]; 
		  real mfbac = (D.f[TS])[ks];//[kts  ]; 
		  real mfbbb = (D.f[REST])[k];//[kzero];
		  real mfccc = (D.f[TNE])[k];//[ktne ]; 
		  real mfaac = (D.f[TSW])[ksw];//[ktsw ]; 
		  real mfcac = (D.f[TSE])[ks];//[ktse ];
		  real mfacc = (D.f[TNW])[kw];//[ktnw ];
		  real mfcca = (D.f[BNE])[kb];//[kbne ];
		  real mfaaa = (D.f[BSW])[kbsw];//[kbsw ];
		  real mfcaa = (D.f[BSE])[kbs];//[kbse ]; 
		  real mfaca = (D.f[BNW])[kbw];//[kbnw ]; 
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
			 // (D.f[E])[ke] + (D.f[W])[kw] +
			 // (D.f[N])[kn] + (D.f[S])[ks] +
			 // (D.f[T])[kt] + (D.f[B])[kb] +
			 // (D.f[NE])[kne] + (D.f[SW])[ksw] +
			 // (D.f[SE])[kse] + (D.f[NW])[knw] +
			 // (D.f[TE])[kte] + (D.f[BW])[kbw] +
			 // (D.f[BE])[kbe] + (D.f[TW])[ktw] +
			 // (D.f[TN])[ktn] + (D.f[BS])[kbs] +
			 // (D.f[BN])[kbn] + (D.f[TS])[kts] +
			 // (D.f[REST])[kzero] +
			 // (D.f[TNE])[ktne] + (D.f[TSW])[ktsw] +
			 // (D.f[TSE])[ktse] + (D.f[TNW])[ktnw] +
			 // (D.f[BNE])[kbne] + (D.f[BSW])[kbsw] +
			 // (D.f[BSE])[kbse] + (D.f[BNW])[kbnw];// +RHO;

    //     vxD[k] =  
			 //((D.f[E  ])[ke  ]- (D.f[W   ])[kw  ]+ 
    //         (D.f[NE  ])[kne ]- (D.f[SW  ])[ksw ]+
    //         (D.f[SE  ])[kse ]- (D.f[NW  ])[knw ]+
    //         (D.f[TE  ])[kte ]- (D.f[BW  ])[kbw ]+
    //         (D.f[BE  ])[kbe ]- (D.f[TW  ])[ktw ]+
    //         (D.f[TNE ])[ktne]- (D.f[TSW ])[ktsw]+ 
    //         (D.f[TSE ])[ktse]- (D.f[TNW ])[ktnw]+ 
    //         (D.f[BNE ])[kbne]- (D.f[BSW ])[kbsw]+ 
    //         (D.f[BSE ])[kbse]- (D.f[BNW ])[kbnw]) / (one + rhoD[k])+
    //         VX;

    //     vyD[k] =  
			 //((D.f[N  ])[kn  ]- (D.f[S   ])[ks  ]+
    //         (D.f[NE  ])[kne ]- (D.f[SW  ])[ksw ]-
    //         (D.f[SE  ])[kse ]+ (D.f[NW  ])[knw ]+
    //         (D.f[TN  ])[ktn ]- (D.f[BS  ])[kbs ]+
    //         (D.f[BN  ])[kbn ]- (D.f[TS  ])[kts ]+
    //         (D.f[TNE ])[ktne]- (D.f[TSW ])[ktsw]- 
    //         (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]+ 
    //         (D.f[BNE ])[kbne]- (D.f[BSW ])[kbsw]- 
    //         (D.f[BSE ])[kbse]+ (D.f[BNW ])[kbnw]) / (one + rhoD[k])+
    //         VY;

    //     vzD[k] =  
			 //((D.f[T  ])[kt  ]- (D.f[B   ])[kb  ]+
    //         (D.f[TE  ])[kte ]- (D.f[BW  ])[kbw ]-
    //         (D.f[BE  ])[kbe ]+ (D.f[TW  ])[ktw ]+
    //         (D.f[TN  ])[ktn ]- (D.f[BS  ])[kbs ]-
    //         (D.f[BN  ])[kbn ]+ (D.f[TS  ])[kts ]+
    //         (D.f[TNE ])[ktne]+ (D.f[TSW ])[ktsw]+ 
    //         (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]- 
    //         (D.f[BNE ])[kbne]- (D.f[BSW ])[kbsw]- 
    //         (D.f[BSE ])[kbse]- (D.f[BNW ])[kbnw]) / (one + rhoD[k])+
    //         VZ;

         pressD[k]  =  ((D.f[E   ])[ke  ]+ (D.f[W   ])[kw  ]+ 
                        (D.f[N   ])[kn  ]+ (D.f[S   ])[ks  ]+
                        (D.f[T   ])[kt  ]+ (D.f[B   ])[kb  ]+
                        c2o1*(
                        (D.f[NE  ])[kne ]+ (D.f[SW  ])[ksw ]+
                        (D.f[SE  ])[kse ]+ (D.f[NW  ])[knw ]+
                        (D.f[TE  ])[kte ]+ (D.f[BW  ])[kbw ]+
                        (D.f[BE  ])[kbe ]+ (D.f[TW  ])[ktw ]+
                        (D.f[TN  ])[ktn ]+ (D.f[BS  ])[kbs ]+
                        (D.f[BN  ])[kbn ]+ (D.f[TS  ])[kts ])+
                        c3o1*(
                        (D.f[TNE ])[ktne]+ (D.f[TSW ])[ktsw]+ 
                        (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]+ 
                        (D.f[BNE ])[kbne]+ (D.f[BSW ])[kbsw]+ 
                        (D.f[BSE ])[kbse]+ (D.f[BNW ])[kbnw])-
                        rhoD[k]-(vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1+rhoD[k])) * c1o2+rhoD[k]+
                        PRESS;    
         //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
      }
   }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMedCompAD27(
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
		D.f[E] = &DD[E   *size_Mat];
		D.f[W] = &DD[W   *size_Mat];
		D.f[N] = &DD[N   *size_Mat];
		D.f[S] = &DD[S   *size_Mat];
		D.f[T] = &DD[T   *size_Mat];
		D.f[B] = &DD[B   *size_Mat];
		D.f[NE] = &DD[NE  *size_Mat];
		D.f[SW] = &DD[SW  *size_Mat];
		D.f[SE] = &DD[SE  *size_Mat];
		D.f[NW] = &DD[NW  *size_Mat];
		D.f[TE] = &DD[TE  *size_Mat];
		D.f[BW] = &DD[BW  *size_Mat];
		D.f[BE] = &DD[BE  *size_Mat];
		D.f[TW] = &DD[TW  *size_Mat];
		D.f[TN] = &DD[TN  *size_Mat];
		D.f[BS] = &DD[BS  *size_Mat];
		D.f[BN] = &DD[BN  *size_Mat];
		D.f[TS] = &DD[TS  *size_Mat];
		D.f[REST] = &DD[REST*size_Mat];
		D.f[TNE] = &DD[TNE *size_Mat];
		D.f[TSW] = &DD[TSW *size_Mat];
		D.f[TSE] = &DD[TSE *size_Mat];
		D.f[TNW] = &DD[TNW *size_Mat];
		D.f[BNE] = &DD[BNE *size_Mat];
		D.f[BSW] = &DD[BSW *size_Mat];
		D.f[BSE] = &DD[BSE *size_Mat];
		D.f[BNW] = &DD[BNW *size_Mat];
	}
	else
	{
		D.f[W] = &DD[E   *size_Mat];
		D.f[E] = &DD[W   *size_Mat];
		D.f[S] = &DD[N   *size_Mat];
		D.f[N] = &DD[S   *size_Mat];
		D.f[B] = &DD[T   *size_Mat];
		D.f[T] = &DD[B   *size_Mat];
		D.f[SW] = &DD[NE  *size_Mat];
		D.f[NE] = &DD[SW  *size_Mat];
		D.f[NW] = &DD[SE  *size_Mat];
		D.f[SE] = &DD[NW  *size_Mat];
		D.f[BW] = &DD[TE  *size_Mat];
		D.f[TE] = &DD[BW  *size_Mat];
		D.f[TW] = &DD[BE  *size_Mat];
		D.f[BE] = &DD[TW  *size_Mat];
		D.f[BS] = &DD[TN  *size_Mat];
		D.f[TN] = &DD[BS  *size_Mat];
		D.f[TS] = &DD[BN  *size_Mat];
		D.f[BN] = &DD[TS  *size_Mat];
		D.f[REST] = &DD[REST*size_Mat];
		D.f[TNE] = &DD[BSW *size_Mat];
		D.f[TSW] = &DD[BNE *size_Mat];
		D.f[TSE] = &DD[BNW *size_Mat];
		D.f[TNW] = &DD[BSE *size_Mat];
		D.f[BNE] = &DD[TSW *size_Mat];
		D.f[BSW] = &DD[TNE *size_Mat];
		D.f[BSE] = &DD[TNW *size_Mat];
		D.f[BNW] = &DD[TSE *size_Mat];
	}
	////////////////////////////////////////////////////////////////////////////////
	Distributions27 Dad;
	if (isEvenTimestep == true)
	{
		Dad.f[E]    = &DD_AD[E   *size_Mat];
		Dad.f[W]    = &DD_AD[W   *size_Mat];
		Dad.f[N]    = &DD_AD[N   *size_Mat];
		Dad.f[S]    = &DD_AD[S   *size_Mat];
		Dad.f[T]    = &DD_AD[T   *size_Mat];
		Dad.f[B]    = &DD_AD[B   *size_Mat];
		Dad.f[NE]   = &DD_AD[NE  *size_Mat];
		Dad.f[SW]   = &DD_AD[SW  *size_Mat];
		Dad.f[SE]   = &DD_AD[SE  *size_Mat];
		Dad.f[NW]   = &DD_AD[NW  *size_Mat];
		Dad.f[TE]   = &DD_AD[TE  *size_Mat];
		Dad.f[BW]   = &DD_AD[BW  *size_Mat];
		Dad.f[BE]   = &DD_AD[BE  *size_Mat];
		Dad.f[TW]   = &DD_AD[TW  *size_Mat];
		Dad.f[TN]   = &DD_AD[TN  *size_Mat];
		Dad.f[BS]   = &DD_AD[BS  *size_Mat];
		Dad.f[BN]   = &DD_AD[BN  *size_Mat];
		Dad.f[TS]   = &DD_AD[TS  *size_Mat];
		Dad.f[REST] = &DD_AD[REST*size_Mat];
		Dad.f[TNE]  = &DD_AD[TNE *size_Mat];
		Dad.f[TSW]  = &DD_AD[TSW *size_Mat];
		Dad.f[TSE]  = &DD_AD[TSE *size_Mat];
		Dad.f[TNW]  = &DD_AD[TNW *size_Mat];
		Dad.f[BNE]  = &DD_AD[BNE *size_Mat];
		Dad.f[BSW]  = &DD_AD[BSW *size_Mat];
		Dad.f[BSE]  = &DD_AD[BSE *size_Mat];
		Dad.f[BNW]  = &DD_AD[BNW *size_Mat];
	}						
	else					
	{						
		Dad.f[W]    = &DD_AD[E   *size_Mat];
		Dad.f[E]    = &DD_AD[W   *size_Mat];
		Dad.f[S]    = &DD_AD[N   *size_Mat];
		Dad.f[N]    = &DD_AD[S   *size_Mat];
		Dad.f[B]    = &DD_AD[T   *size_Mat];
		Dad.f[T]    = &DD_AD[B   *size_Mat];
		Dad.f[SW]   = &DD_AD[NE  *size_Mat];
		Dad.f[NE]   = &DD_AD[SW  *size_Mat];
		Dad.f[NW]   = &DD_AD[SE  *size_Mat];
		Dad.f[SE]   = &DD_AD[NW  *size_Mat];
		Dad.f[BW]   = &DD_AD[TE  *size_Mat];
		Dad.f[TE]   = &DD_AD[BW  *size_Mat];
		Dad.f[TW]   = &DD_AD[BE  *size_Mat];
		Dad.f[BE]   = &DD_AD[TW  *size_Mat];
		Dad.f[BS]   = &DD_AD[TN  *size_Mat];
		Dad.f[TN]   = &DD_AD[BS  *size_Mat];
		Dad.f[TS]   = &DD_AD[BN  *size_Mat];
		Dad.f[BN]   = &DD_AD[TS  *size_Mat];
		Dad.f[REST] = &DD_AD[REST*size_Mat];
		Dad.f[TNE]  = &DD_AD[BSW *size_Mat];
		Dad.f[TSW]  = &DD_AD[BNE *size_Mat];
		Dad.f[TSE]  = &DD_AD[BNW *size_Mat];
		Dad.f[TNW]  = &DD_AD[BSE *size_Mat];
		Dad.f[BNE]  = &DD_AD[TSW *size_Mat];
		Dad.f[BSW]  = &DD_AD[TNE *size_Mat];
		Dad.f[BSE]  = &DD_AD[TNW *size_Mat];
		Dad.f[BNW]  = &DD_AD[TSE *size_Mat];
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
			real mfcbb = (D.f[E])[k];//[ke   ];
			real mfabb = (D.f[W])[kw];//[kw   ];  
			real mfbcb = (D.f[N])[k];//[kn   ];
			real mfbab = (D.f[S])[ks];//[ks   ];  
			real mfbbc = (D.f[T])[k];//[kt   ];
			real mfbba = (D.f[B])[kb];//[kb   ];  
			real mfccb = (D.f[NE])[k];//[kne  ];  
			real mfaab = (D.f[SW])[ksw];//[ksw  ];
			real mfcab = (D.f[SE])[ks];//[kse  ]; 
			real mfacb = (D.f[NW])[kw];//[knw  ]; 
			real mfcbc = (D.f[TE])[k];//[kte  ];  
			real mfaba = (D.f[BW])[kbw];//[kbw  ];
			real mfcba = (D.f[BE])[kb];//[kbe  ]; 
			real mfabc = (D.f[TW])[kw];//[ktw  ]; 
			real mfbcc = (D.f[TN])[k];//[ktn  ];  
			real mfbaa = (D.f[BS])[kbs];//[kbs  ];
			real mfbca = (D.f[BN])[kb];//[kbn  ]; 
			real mfbac = (D.f[TS])[ks];//[kts  ]; 
			real mfbbb = (D.f[REST])[k];//[kzero];
			real mfccc = (D.f[TNE])[k];//[ktne ]; 
			real mfaac = (D.f[TSW])[ksw];//[ktsw ]; 
			real mfcac = (D.f[TSE])[ks];//[ktse ];
			real mfacc = (D.f[TNW])[kw];//[ktnw ];
			real mfcca = (D.f[BNE])[kb];//[kbne ];
			real mfaaa = (D.f[BSW])[kbsw];//[kbsw ];
			real mfcaa = (D.f[BSE])[kbs];//[kbse ]; 
			real mfaca = (D.f[BNW])[kbw];//[kbnw ]; 
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
				((D.f[E])[ke] + (D.f[W])[kw] +
				 (D.f[N])[kn] + (D.f[S])[ks] +
				 (D.f[T])[kt] + (D.f[B])[kb] +
				 c2o1*(
				 (D.f[NE])[kne] + (D.f[SW])[ksw] +
				 (D.f[SE])[kse] + (D.f[NW])[knw] +
				 (D.f[TE])[kte] + (D.f[BW])[kbw] +
				 (D.f[BE])[kbe] + (D.f[TW])[ktw] +
				 (D.f[TN])[ktn] + (D.f[BS])[kbs] +
				 (D.f[BN])[kbn] + (D.f[TS])[kts]) +
				 c3o1*(
				 (D.f[TNE])[ktne] + (D.f[TSW])[ktsw] +
				 (D.f[TSE])[ktse] + (D.f[TNW])[ktnw] +
				 (D.f[BNE])[kbne] + (D.f[BSW])[kbsw] +
				 (D.f[BSE])[kbse] + (D.f[BNW])[kbnw]) -
				 rhoD[k] - (vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1 + rhoD[k])) * c1o2 + rhoD[k] +
				 PRESS;
				 //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
			//////////////////////////////////////////////////////////////////////////
			mfcbb = (Dad.f[E   ])[k   ];
			mfabb = (Dad.f[W   ])[kw  ];
			mfbcb = (Dad.f[N   ])[k   ];
			mfbab = (Dad.f[S   ])[ks  ];
			mfbbc = (Dad.f[T   ])[k   ];
			mfbba = (Dad.f[B   ])[kb  ];
			mfccb = (Dad.f[NE  ])[k   ];
			mfaab = (Dad.f[SW  ])[ksw ];
			mfcab = (Dad.f[SE  ])[ks  ];
			mfacb = (Dad.f[NW  ])[kw  ];
			mfcbc = (Dad.f[TE  ])[k   ];
			mfaba = (Dad.f[BW  ])[kbw ];
			mfcba = (Dad.f[BE  ])[kb  ];
			mfabc = (Dad.f[TW  ])[kw  ];
			mfbcc = (Dad.f[TN  ])[k   ];
			mfbaa = (Dad.f[BS  ])[kbs ];
			mfbca = (Dad.f[BN  ])[kb  ];
			mfbac = (Dad.f[TS  ])[ks  ];
			mfbbb = (Dad.f[REST])[k   ];
			mfccc = (Dad.f[TNE ])[k   ];
			mfaac = (Dad.f[TSW ])[ksw ];
			mfcac = (Dad.f[TSE ])[ks  ];
			mfacc = (Dad.f[TNW ])[kw  ];
			mfcca = (Dad.f[BNE ])[kb  ];
			mfaaa = (Dad.f[BSW ])[kbsw];
			mfcaa = (Dad.f[BSE ])[kbs ];
			mfaca = (Dad.f[BNW ])[kbw ];
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
extern "C" __global__ void LBCalcMacMedSP27( real* vxD,
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
extern "C" __global__ void LBResetMedianValuesSP27(
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
extern "C" __global__ void LBResetMedianValuesAD27(
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
extern "C" __global__ void LBCalcMeasurePoints( real* vxMP,
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
		D.f[E   ] = &DD[E   *size_Mat];
		D.f[W   ] = &DD[W   *size_Mat];
		D.f[N   ] = &DD[N   *size_Mat];
		D.f[S   ] = &DD[S   *size_Mat];
		D.f[T   ] = &DD[T   *size_Mat];
		D.f[B   ] = &DD[B   *size_Mat];
		D.f[NE  ] = &DD[NE  *size_Mat];
		D.f[SW  ] = &DD[SW  *size_Mat];
		D.f[SE  ] = &DD[SE  *size_Mat];
		D.f[NW  ] = &DD[NW  *size_Mat];
		D.f[TE  ] = &DD[TE  *size_Mat];
		D.f[BW  ] = &DD[BW  *size_Mat];
		D.f[BE  ] = &DD[BE  *size_Mat];
		D.f[TW  ] = &DD[TW  *size_Mat];
		D.f[TN  ] = &DD[TN  *size_Mat];
		D.f[BS  ] = &DD[BS  *size_Mat];
		D.f[BN  ] = &DD[BN  *size_Mat];
		D.f[TS  ] = &DD[TS  *size_Mat];
		D.f[REST] = &DD[REST*size_Mat];
		D.f[TNE ] = &DD[TNE *size_Mat];
		D.f[TSW ] = &DD[TSW *size_Mat];
		D.f[TSE ] = &DD[TSE *size_Mat];
		D.f[TNW ] = &DD[TNW *size_Mat];
		D.f[BNE ] = &DD[BNE *size_Mat];
		D.f[BSW ] = &DD[BSW *size_Mat];
		D.f[BSE ] = &DD[BSE *size_Mat];
		D.f[BNW ] = &DD[BNW *size_Mat];
	} 
	else
	{
		D.f[W   ] = &DD[E   *size_Mat];
		D.f[E   ] = &DD[W   *size_Mat];
		D.f[S   ] = &DD[N   *size_Mat];
		D.f[N   ] = &DD[S   *size_Mat];
		D.f[B   ] = &DD[T   *size_Mat];
		D.f[T   ] = &DD[B   *size_Mat];
		D.f[SW  ] = &DD[NE  *size_Mat];
		D.f[NE  ] = &DD[SW  *size_Mat];
		D.f[NW  ] = &DD[SE  *size_Mat];
		D.f[SE  ] = &DD[NW  *size_Mat];
		D.f[BW  ] = &DD[TE  *size_Mat];
		D.f[TE  ] = &DD[BW  *size_Mat];
		D.f[TW  ] = &DD[BE  *size_Mat];
		D.f[BE  ] = &DD[TW  *size_Mat];
		D.f[BS  ] = &DD[TN  *size_Mat];
		D.f[TN  ] = &DD[BS  *size_Mat];
		D.f[TS  ] = &DD[BN  *size_Mat];
		D.f[BN  ] = &DD[TS  *size_Mat];
		D.f[REST] = &DD[REST*size_Mat];
		D.f[TNE ] = &DD[BSW *size_Mat];
		D.f[TSW ] = &DD[BNE *size_Mat];
		D.f[TSE ] = &DD[BNW *size_Mat];
		D.f[TNW ] = &DD[BSE *size_Mat];
		D.f[BNE ] = &DD[TSW *size_Mat];
		D.f[BSW ] = &DD[TNE *size_Mat];
		D.f[BSE ] = &DD[TNW *size_Mat];
		D.f[BNW ] = &DD[TSE *size_Mat];
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
         rhoMP[kMac]=   (D.f[E   ])[ke  ]+ (D.f[W   ])[kw  ]+ 
                        (D.f[N   ])[kn  ]+ (D.f[S   ])[ks  ]+
                        (D.f[T   ])[kt  ]+ (D.f[B   ])[kb  ]+
                        (D.f[NE  ])[kne ]+ (D.f[SW  ])[ksw ]+
                        (D.f[SE  ])[kse ]+ (D.f[NW  ])[knw ]+
                        (D.f[TE  ])[kte ]+ (D.f[BW  ])[kbw ]+
                        (D.f[BE  ])[kbe ]+ (D.f[TW  ])[ktw ]+
                        (D.f[TN  ])[ktn ]+ (D.f[BS  ])[kbs ]+
                        (D.f[BN  ])[kbn ]+ (D.f[TS  ])[kts ]+
                        (D.f[REST])[kzero]+ 
                        (D.f[TNE ])[ktne]+ (D.f[TSW ])[ktsw]+ 
                        (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]+ 
                        (D.f[BNE ])[kbne]+ (D.f[BSW ])[kbsw]+ 
                        (D.f[BSE ])[kbse]+ (D.f[BNW ])[kbnw];

         vxMP[kMac] =   (D.f[E   ])[ke  ]- (D.f[W   ])[kw  ]+ 
                        (D.f[NE  ])[kne ]- (D.f[SW  ])[ksw ]+
                        (D.f[SE  ])[kse ]- (D.f[NW  ])[knw ]+
                        (D.f[TE  ])[kte ]- (D.f[BW  ])[kbw ]+
                        (D.f[BE  ])[kbe ]- (D.f[TW  ])[ktw ]+
                        (D.f[TNE ])[ktne]- (D.f[TSW ])[ktsw]+ 
                        (D.f[TSE ])[ktse]- (D.f[TNW ])[ktnw]+ 
                        (D.f[BNE ])[kbne]- (D.f[BSW ])[kbsw]+ 
                        (D.f[BSE ])[kbse]- (D.f[BNW ])[kbnw];

         vyMP[kMac] =   (D.f[N   ])[kn  ]- (D.f[S   ])[ks  ]+
                        (D.f[NE  ])[kne ]- (D.f[SW  ])[ksw ]-
                        (D.f[SE  ])[kse ]+ (D.f[NW  ])[knw ]+
                        (D.f[TN  ])[ktn ]- (D.f[BS  ])[kbs ]+
                        (D.f[BN  ])[kbn ]- (D.f[TS  ])[kts ]+
                        (D.f[TNE ])[ktne]- (D.f[TSW ])[ktsw]- 
                        (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]+ 
                        (D.f[BNE ])[kbne]- (D.f[BSW ])[kbsw]- 
                        (D.f[BSE ])[kbse]+ (D.f[BNW ])[kbnw];

         vzMP[kMac] =   (D.f[T   ])[kt  ]- (D.f[B   ])[kb  ]+
                        (D.f[TE  ])[kte ]- (D.f[BW  ])[kbw ]-
                        (D.f[BE  ])[kbe ]+ (D.f[TW  ])[ktw ]+
                        (D.f[TN  ])[ktn ]- (D.f[BS  ])[kbs ]-
                        (D.f[BN  ])[kbn ]+ (D.f[TS  ])[kts ]+
                        (D.f[TNE ])[ktne]+ (D.f[TSW ])[ktsw]+ 
                        (D.f[TSE ])[ktse]+ (D.f[TNW ])[ktnw]- 
                        (D.f[BNE ])[kbne]- (D.f[BSW ])[kbsw]- 
                        (D.f[BSE ])[kbse]- (D.f[BNW ])[kbnw];
      }
   }
}
////////////////////////////////////////////////////////////////////////////////





































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBSetOutputWallVelocitySP27( real* vxD,
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





























