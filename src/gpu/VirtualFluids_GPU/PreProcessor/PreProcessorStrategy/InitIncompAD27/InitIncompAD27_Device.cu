#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void LB_Init_Incomp_AD_27(unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int* geoD,
	real* Conc,
	real* ux,
	real* uy,
	real* uz,
	unsigned int size_Mat,
	real* DD27,
	bool EvenOrOdd)
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
      ////////////////////////////////////////////////////////////////////////////////
      unsigned int BC;
      BC        =   geoD[k];

      if( BC != GEO_SOLID && BC != GEO_VOID)
      {
         Distributions27 D27;
         if (EvenOrOdd==true)
         {
            D27.f[dP00   ] = &DD27[dP00   *size_Mat];
            D27.f[dM00   ] = &DD27[dM00   *size_Mat];
            D27.f[DIR_0P0   ] = &DD27[DIR_0P0   *size_Mat];
            D27.f[DIR_0M0   ] = &DD27[DIR_0M0   *size_Mat];
            D27.f[DIR_00P   ] = &DD27[DIR_00P   *size_Mat];
            D27.f[DIR_00M   ] = &DD27[DIR_00M   *size_Mat];
            D27.f[DIR_PP0  ] = &DD27[DIR_PP0  *size_Mat];
            D27.f[DIR_MM0  ] = &DD27[DIR_MM0  *size_Mat];
            D27.f[DIR_PM0  ] = &DD27[DIR_PM0  *size_Mat];
            D27.f[DIR_MP0  ] = &DD27[DIR_MP0  *size_Mat];
            D27.f[DIR_P0P  ] = &DD27[DIR_P0P  *size_Mat];
            D27.f[DIR_M0M  ] = &DD27[DIR_M0M  *size_Mat];
            D27.f[DIR_P0M  ] = &DD27[DIR_P0M  *size_Mat];
            D27.f[DIR_M0P  ] = &DD27[DIR_M0P  *size_Mat];
            D27.f[DIR_0PP  ] = &DD27[DIR_0PP  *size_Mat];
            D27.f[DIR_0MM  ] = &DD27[DIR_0MM  *size_Mat];
            D27.f[DIR_0PM  ] = &DD27[DIR_0PM  *size_Mat];
            D27.f[DIR_0MP  ] = &DD27[DIR_0MP  *size_Mat];
            D27.f[d000] = &DD27[d000*size_Mat];
            D27.f[DIR_PPP ] = &DD27[DIR_PPP *size_Mat];
            D27.f[DIR_MMP ] = &DD27[DIR_MMP *size_Mat];
            D27.f[DIR_PMP ] = &DD27[DIR_PMP *size_Mat];
            D27.f[DIR_MPP ] = &DD27[DIR_MPP *size_Mat];
            D27.f[DIR_PPM ] = &DD27[DIR_PPM *size_Mat];
            D27.f[DIR_MMM ] = &DD27[DIR_MMM *size_Mat];
            D27.f[DIR_PMM ] = &DD27[DIR_PMM *size_Mat];
            D27.f[DIR_MPM ] = &DD27[DIR_MPM *size_Mat];
         }
         else
         {
            D27.f[dM00   ] = &DD27[dP00   *size_Mat];
            D27.f[dP00   ] = &DD27[dM00   *size_Mat];
            D27.f[DIR_0M0   ] = &DD27[DIR_0P0   *size_Mat];
            D27.f[DIR_0P0   ] = &DD27[DIR_0M0   *size_Mat];
            D27.f[DIR_00M   ] = &DD27[DIR_00P   *size_Mat];
            D27.f[DIR_00P   ] = &DD27[DIR_00M   *size_Mat];
            D27.f[DIR_MM0  ] = &DD27[DIR_PP0  *size_Mat];
            D27.f[DIR_PP0  ] = &DD27[DIR_MM0  *size_Mat];
            D27.f[DIR_MP0  ] = &DD27[DIR_PM0  *size_Mat];
            D27.f[DIR_PM0  ] = &DD27[DIR_MP0  *size_Mat];
            D27.f[DIR_M0M  ] = &DD27[DIR_P0P  *size_Mat];
            D27.f[DIR_P0P  ] = &DD27[DIR_M0M  *size_Mat];
            D27.f[DIR_M0P  ] = &DD27[DIR_P0M  *size_Mat];
            D27.f[DIR_P0M  ] = &DD27[DIR_M0P  *size_Mat];
            D27.f[DIR_0MM  ] = &DD27[DIR_0PP  *size_Mat];
            D27.f[DIR_0PP  ] = &DD27[DIR_0MM  *size_Mat];
            D27.f[DIR_0MP  ] = &DD27[DIR_0PM  *size_Mat];
            D27.f[DIR_0PM  ] = &DD27[DIR_0MP  *size_Mat];
            D27.f[d000] = &DD27[d000*size_Mat];
            D27.f[DIR_MMM ] = &DD27[DIR_PPP *size_Mat];
            D27.f[DIR_PPM ] = &DD27[DIR_MMP *size_Mat];
            D27.f[DIR_MPM ] = &DD27[DIR_PMP *size_Mat];
            D27.f[DIR_PMM ] = &DD27[DIR_MPP *size_Mat];
            D27.f[DIR_MMP ] = &DD27[DIR_PPM *size_Mat];
            D27.f[DIR_PPP ] = &DD27[DIR_MMM *size_Mat];
            D27.f[DIR_MPP ] = &DD27[DIR_PMM *size_Mat];
            D27.f[DIR_PMP ] = &DD27[DIR_MPM *size_Mat];
         }
         //////////////////////////////////////////////////////////////////////////
         real ConcD = Conc[k];
         real   vx1 = ux[k];
         real   vx2 = uy[k];
         real   vx3 = uz[k];
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //D3Q27
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         (D27.f[d000])[kzero] =   c8o27* ConcD*(c1o1-cu_sq);
         (D27.f[dP00   ])[ke   ] =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
         (D27.f[dM00   ])[kw   ] =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
         (D27.f[DIR_0P0   ])[kn   ] =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
         (D27.f[DIR_0M0   ])[ks   ] =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
         (D27.f[DIR_00P   ])[kt   ] =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
         (D27.f[DIR_00M   ])[kb   ] =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
         (D27.f[DIR_PP0  ])[kne  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         (D27.f[DIR_MM0  ])[ksw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         (D27.f[DIR_PM0  ])[kse  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         (D27.f[DIR_MP0  ])[knw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         (D27.f[DIR_P0P  ])[kte  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         (D27.f[DIR_M0M  ])[kbw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         (D27.f[DIR_P0M  ])[kbe  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         (D27.f[DIR_M0P  ])[ktw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         (D27.f[DIR_0PP  ])[ktn  ] =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         (D27.f[DIR_0MM  ])[kbs  ] =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         (D27.f[DIR_0PM  ])[kbn  ] =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         (D27.f[DIR_0MP  ])[kts  ] =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         (D27.f[DIR_PPP ])[ktne ] =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         (D27.f[DIR_MMM ])[kbsw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         (D27.f[DIR_PPM ])[kbne ] =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         (D27.f[DIR_MMP ])[ktsw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         (D27.f[DIR_PMP ])[ktse ] =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         (D27.f[DIR_MPM ])[kbnw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         (D27.f[DIR_PMM ])[kbse ] =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         (D27.f[DIR_MPP ])[ktnw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      }
   }
}