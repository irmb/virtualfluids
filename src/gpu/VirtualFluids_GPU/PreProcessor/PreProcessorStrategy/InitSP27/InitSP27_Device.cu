#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void LB_Init_SP_27(unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int* geoD,
	real* rho,
	real* ux,
	real* uy,
	real* uz,
	unsigned int size_Mat,
	real* DD,
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

      if( BC != GEO_SOLID &&  BC != GEO_VOID)
      {
         Distributions27 D;
         if (EvenOrOdd==true)
         {
            D.f[dP00   ] = &DD[dP00   *size_Mat];
            D.f[dM00   ] = &DD[dM00   *size_Mat];
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
            D.f[d000] = &DD[d000*size_Mat];
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
            D.f[dM00   ] = &DD[dP00   *size_Mat];
            D.f[dP00   ] = &DD[dM00   *size_Mat];
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
            D.f[d000] = &DD[d000*size_Mat];
            D.f[DIR_MMM ] = &DD[DIR_PPP *size_Mat];
            D.f[DIR_PPM ] = &DD[DIR_MMP *size_Mat];
            D.f[DIR_MPM ] = &DD[DIR_PMP *size_Mat];
            D.f[DIR_PMM ] = &DD[DIR_MPP *size_Mat];
            D.f[DIR_MMP ] = &DD[DIR_PPM *size_Mat];
            D.f[DIR_PPP ] = &DD[DIR_MMM *size_Mat];
            D.f[DIR_MPP ] = &DD[DIR_PMM *size_Mat];
            D.f[DIR_PMP ] = &DD[DIR_MPM *size_Mat];
         }
         //////////////////////////////////////////////////////////////////////////
         real drho = rho[k];//0.0f;//
         real  vx1 = ux[k]; //0.0f;//
         real  vx2 = uy[k]; //0.0f;//
         real  vx3 = uz[k]; //0.0f;//
         //////////////////////////////////////////////////////////////////////////
         //index
         //////////////////////////////////////////////////////////////////////////
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
         real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         (D.f[d000])[kzero] =   c8o27* (drho-cu_sq);
         (D.f[dP00   ])[ke   ] =   c2o27* (drho+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
         (D.f[dM00   ])[kw   ] =   c2o27* (drho+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
         (D.f[DIR_0P0   ])[kn   ] =   c2o27* (drho+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
         (D.f[DIR_0M0   ])[ks   ] =   c2o27* (drho+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
         (D.f[DIR_00P   ])[kt   ] =   c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
         (D.f[DIR_00M   ])[kb   ] =   c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
         (D.f[DIR_PP0  ])[kne  ] =   c1o54* (drho+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         (D.f[DIR_MM0  ])[ksw  ] =   c1o54* (drho+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         (D.f[DIR_PM0  ])[kse  ] =   c1o54* (drho+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         (D.f[DIR_MP0  ])[knw  ] =   c1o54* (drho+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         (D.f[DIR_P0P  ])[kte  ] =   c1o54* (drho+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         (D.f[DIR_M0M  ])[kbw  ] =   c1o54* (drho+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         (D.f[DIR_P0M  ])[kbe  ] =   c1o54* (drho+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         (D.f[DIR_M0P  ])[ktw  ] =   c1o54* (drho+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         (D.f[DIR_0PP  ])[ktn  ] =   c1o54* (drho+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         (D.f[DIR_0MM  ])[kbs  ] =   c1o54* (drho+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         (D.f[DIR_0PM  ])[kbn  ] =   c1o54* (drho+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         (D.f[DIR_0MP  ])[kts  ] =   c1o54* (drho+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         (D.f[DIR_PPP ])[ktne ] =   c1o216*(drho+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         (D.f[DIR_MMM ])[kbsw ] =   c1o216*(drho+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         (D.f[DIR_PPM ])[kbne ] =   c1o216*(drho+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         (D.f[DIR_MMP ])[ktsw ] =   c1o216*(drho+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         (D.f[DIR_PMP ])[ktse ] =   c1o216*(drho+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         (D.f[DIR_MPM ])[kbnw ] =   c1o216*(drho+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         (D.f[DIR_PMM ])[kbse ] =   c1o216*(drho+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         (D.f[DIR_MPP ])[ktnw ] =   c1o216*(drho+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      }
	  else
	  {
		  //////////////////////////////////////////////////////////////////////////
		  Distributions27 D;
		  D.f[d000] = &DD[d000*size_Mat];
		  //////////////////////////////////////////////////////////////////////////
		  (D.f[d000])[k] = c96o1;
		  //////////////////////////////////////////////////////////////////////////
	  }
   }
}