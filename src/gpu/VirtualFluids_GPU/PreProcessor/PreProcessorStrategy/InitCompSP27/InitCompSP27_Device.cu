#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
#include "math.h"

#include <stdio.h>

__global__ void LB_Init_Comp_SP_27(unsigned int* neighborX,
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

         (D.f[DIR_000])[kzero] =   c8o27* (drho-cu_sq*(c1o1+drho));
         (D.f[DIR_P00   ])[ke   ] =   c2o27* (drho+ (c1o1+drho) * (c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq));
         (D.f[DIR_M00   ])[kw   ] =   c2o27* (drho+ (c1o1+drho) * (c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq));
         (D.f[DIR_0P0   ])[kn   ] =   c2o27* (drho+ (c1o1+drho) * (c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq));
         (D.f[DIR_0M0   ])[ks   ] =   c2o27* (drho+ (c1o1+drho) * (c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq));
         (D.f[DIR_00P   ])[kt   ] =   c2o27* (drho+ (c1o1+drho) * (c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq));
         (D.f[DIR_00M   ])[kb   ] =   c2o27* (drho+ (c1o1+drho) * (c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq));
         (D.f[DIR_PP0  ])[kne  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
         (D.f[DIR_MM0  ])[ksw  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
         (D.f[DIR_PM0  ])[kse  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
         (D.f[DIR_MP0  ])[knw  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
         (D.f[DIR_P0P  ])[kte  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
         (D.f[DIR_M0M  ])[kbw  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
         (D.f[DIR_P0M  ])[kbe  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
         (D.f[DIR_M0P  ])[ktw  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
         (D.f[DIR_0PP  ])[ktn  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
         (D.f[DIR_0MM  ])[kbs  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
         (D.f[DIR_0PM  ])[kbn  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
         (D.f[DIR_0MP  ])[kts  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
         (D.f[DIR_PPP ])[ktne ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
         (D.f[DIR_MMM ])[kbsw ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
         (D.f[DIR_PPM ])[kbne ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
         (D.f[DIR_MMP ])[ktsw ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
         (D.f[DIR_PMP ])[ktse ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
         (D.f[DIR_MPM ])[kbnw ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
         (D.f[DIR_PMM ])[kbse ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
         (D.f[DIR_MPP ])[ktnw ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));
      }
   }
}










////////////////////////////////////////////////////////////////////////////////
__global__ void LB_Init_Comp_Neq_SP_27( unsigned int* neighborX,
                                                   unsigned int* neighborY,
                                                   unsigned int* neighborZ,
                                                   unsigned int* neighborWSB,
                                                   unsigned int* geoD,
                                                   real* rho,
                                                   real* ux,
                                                   real* uy,
                                                   real* uz,
                                                   unsigned int size_Mat,
                                                   real* DD,
                                                   real omega,
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
        BC = geoD[k];

        if( BC != GEO_SOLID &&  BC != GEO_VOID)
        {
            Distributions27 D;
            if (EvenOrOdd==true)
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
	        //////////////////////////////////////////////////////////////////////////////
	        //neighbor index
	        uint kPx   = neighborX[k];
	        uint kPy   = neighborY[k];
	        uint kPz   = neighborZ[k];
	        uint kMxyz = neighborWSB[k];
	        uint kMx   = neighborZ[neighborY[kMxyz]];
	        uint kMy   = neighborZ[neighborX[kMxyz]];
	        uint kMz   = neighborY[neighborX[kMxyz]];
            //////////////////////////////////////////////////////////////////////////
	        //getVeloX//
	        real vx1NeighborPx = ux[kPx];
	        real vx1NeighborMx = ux[kMx];
	        real vx1NeighborPy = ux[kPy];
	        real vx1NeighborMy = ux[kMy];
	        real vx1NeighborPz = ux[kPz];
	        real vx1NeighborMz = ux[kMz];
	        //getVeloY//
	        real vx2NeighborPx = uy[kPx];
	        real vx2NeighborMx = uy[kMx];
	        real vx2NeighborPy = uy[kPy];
	        real vx2NeighborMy = uy[kMy];
	        real vx2NeighborPz = uy[kPz];
	        real vx2NeighborMz = uy[kMz];
	        //getVeloZ//
	        real vx3NeighborPx = uz[kPx];
	        real vx3NeighborMx = uz[kMx];
	        real vx3NeighborPy = uz[kPy];
	        real vx3NeighborMy = uz[kMy];
	        real vx3NeighborPz = uz[kPz];
	        real vx3NeighborMz = uz[kMz];
            //////////////////////////////////////////////////////////////////////////

	        real dvx1dx = (vx1NeighborPx - vx1NeighborMx) / c2o1;
	        real dvx1dy = (vx1NeighborPy - vx1NeighborMy) / c2o1;
	        real dvx1dz = (vx1NeighborPz - vx1NeighborMz) / c2o1;

	        real dvx2dx = (vx2NeighborPx - vx2NeighborMx) / c2o1;
	        real dvx2dy = (vx2NeighborPy - vx2NeighborMy) / c2o1;
	        real dvx2dz = (vx2NeighborPz - vx2NeighborMz) / c2o1;

	        real dvx3dx = (vx3NeighborPx - vx3NeighborMx) / c2o1;
	        real dvx3dy = (vx3NeighborPy - vx3NeighborMy) / c2o1;
	        real dvx3dz = (vx3NeighborPz - vx3NeighborMz) / c2o1;

            //////////////////////////////////////////////////////////////////////////

            // the following code is copy and pasted from VirtualFluidsCore/Visitors/InitDistributionsBlockVisitor.cpp
            // i.e. Konstantins code

            real ax = dvx1dx;
            real ay = dvx1dy;
            real az = dvx1dz;

            real bx = dvx2dx;
            real by = dvx2dy;
            real bz = dvx2dz;

            real cx = dvx3dx;
            real cy = dvx3dy;
            real cz = dvx3dz;

            real eps_new = c1o1;
            real op      = c1o1;
            real o       = omega;

            real f_E    =            eps_new *((5.*ax*o + 5.*by*o + 5.*cz*o - 8.*ax*op + 4.*by*op + 4.*cz*op)/(54.*o*op));

            real f_N    =    f_E   + eps_new *((2.*(ax - by))/(9.*o));
            real f_T    =    f_E   + eps_new *((2.*(ax - cz))/(9.*o));
            real f_NE   =            eps_new *(-(5.*cz*o + 3.*(ay + bx)*op - 2.*cz*op + ax*(5.*o + op) + by*(5.*o + op))/(54.*o*op));
            real f_SE   =    f_NE  + eps_new *((  ay + bx )/(9.*o));
            real f_TE   =            eps_new *(-(5.*cz*o + by*(5.*o - 2.*op) + 3.*(az + cx)*op + cz*op + ax*(5.*o + op))/(54.*o*op));
            real f_BE   =    f_TE  + eps_new *((  az + cx )/(9.*o));
            real f_TN   =            eps_new *(-(5.*ax*o + 5.*by*o + 5.*cz*o - 2.*ax*op + by*op + 3.*bz*op + 3.*cy*op + cz*op)/(54.*o*op));
            real f_BN   =    f_TN  + eps_new *((  bz + cy )/(9.*o));
            real f_ZERO =            eps_new *((5.*(ax + by + cz))/(9.*op));
            real f_TNE  =            eps_new *(-(ay + az + bx + bz + cx + cy)/(72.*o));
            real f_TSW  =  - f_TNE - eps_new *((ay + bx)/(36.*o));
            real f_TSE  =  - f_TNE - eps_new *((az + cx)/(36.*o));
            real f_TNW  =  - f_TNE - eps_new *((bz + cy)/(36.*o));

            //////////////////////////////////////////////////////////////////////////
            real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

            (D.f[DIR_000])[kzero] =   c8o27* (drho-cu_sq*(c1o1+drho));
            (D.f[DIR_P00   ])[ke   ] =   c2o27* (drho+ (c1o1+drho) * (c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq));
            (D.f[DIR_M00   ])[kw   ] =   c2o27* (drho+ (c1o1+drho) * (c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq));
            (D.f[DIR_0P0   ])[kn   ] =   c2o27* (drho+ (c1o1+drho) * (c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq));
            (D.f[DIR_0M0   ])[ks   ] =   c2o27* (drho+ (c1o1+drho) * (c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq));
            (D.f[DIR_00P   ])[kt   ] =   c2o27* (drho+ (c1o1+drho) * (c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq));
            (D.f[DIR_00M   ])[kb   ] =   c2o27* (drho+ (c1o1+drho) * (c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq));
            (D.f[DIR_PP0  ])[kne  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
            (D.f[DIR_MM0  ])[ksw  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
            (D.f[DIR_PM0  ])[kse  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
            (D.f[DIR_MP0  ])[knw  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
            (D.f[DIR_P0P  ])[kte  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
            (D.f[DIR_M0M  ])[kbw  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
            (D.f[DIR_P0M  ])[kbe  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
            (D.f[DIR_M0P  ])[ktw  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
            (D.f[DIR_0PP  ])[ktn  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
            (D.f[DIR_0MM  ])[kbs  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
            (D.f[DIR_0PM  ])[kbn  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
            (D.f[DIR_0MP  ])[kts  ] =   c1o54* (drho+ (c1o1+drho) * (c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
            (D.f[DIR_PPP ])[ktne ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
            (D.f[DIR_MMM ])[kbsw ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
            (D.f[DIR_PPM ])[kbne ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
            (D.f[DIR_MMP ])[ktsw ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
            (D.f[DIR_PMP ])[ktse ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
            (D.f[DIR_MPM ])[kbnw ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
            (D.f[DIR_PMM ])[kbse ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
            (D.f[DIR_MPP ])[ktnw ] =   c1o216*(drho+ (c1o1+drho) * (c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));

            //////////////////////////////////////////////////////////////////////////

            (D.f[DIR_000])[kzero] += (c1o1+drho) * f_ZERO;
            (D.f[DIR_P00   ])[ke   ] += (c1o1+drho) * f_E   ;
            (D.f[DIR_M00   ])[kw   ] += (c1o1+drho) * f_E   ;
            (D.f[DIR_0P0   ])[kn   ] += (c1o1+drho) * f_N   ;
            (D.f[DIR_0M0   ])[ks   ] += (c1o1+drho) * f_N   ;
            (D.f[DIR_00P   ])[kt   ] += (c1o1+drho) * f_T   ;
            (D.f[DIR_00M   ])[kb   ] += (c1o1+drho) * f_T   ;
            (D.f[DIR_PP0  ])[kne  ] += (c1o1+drho) * f_NE  ;
            (D.f[DIR_MM0  ])[ksw  ] += (c1o1+drho) * f_NE  ;
            (D.f[DIR_PM0  ])[kse  ] += (c1o1+drho) * f_SE  ;
            (D.f[DIR_MP0  ])[knw  ] += (c1o1+drho) * f_SE  ;
            (D.f[DIR_P0P  ])[kte  ] += (c1o1+drho) * f_TE  ;
            (D.f[DIR_M0M  ])[kbw  ] += (c1o1+drho) * f_TE  ;
            (D.f[DIR_P0M  ])[kbe  ] += (c1o1+drho) * f_BE  ;
            (D.f[DIR_M0P  ])[ktw  ] += (c1o1+drho) * f_BE  ;
            (D.f[DIR_0PP  ])[ktn  ] += (c1o1+drho) * f_TN  ;
            (D.f[DIR_0MM  ])[kbs  ] += (c1o1+drho) * f_TN  ;
            (D.f[DIR_0PM  ])[kbn  ] += (c1o1+drho) * f_BN  ;
            (D.f[DIR_0MP  ])[kts  ] += (c1o1+drho) * f_BN  ;
            (D.f[DIR_PPP ])[ktne ] += (c1o1+drho) * f_TNE ;
            (D.f[DIR_MMM ])[kbsw ] += (c1o1+drho) * f_TNE ;
            (D.f[DIR_PPM ])[kbne ] += (c1o1+drho) * f_TSW ;
            (D.f[DIR_MMP ])[ktsw ] += (c1o1+drho) * f_TSW ;
            (D.f[DIR_PMP ])[ktse ] += (c1o1+drho) * f_TSE ;
            (D.f[DIR_MPM ])[kbnw ] += (c1o1+drho) * f_TSE ;
            (D.f[DIR_PMM ])[kbse ] += (c1o1+drho) * f_TNW ;
            (D.f[DIR_MPP ])[ktnw ] += (c1o1+drho) * f_TNW ;

            //////////////////////////////////////////////////////////////////////////
        }
	    else
	    {
		    //////////////////////////////////////////////////////////////////////////
		    Distributions27 D;
		    D.f[DIR_000] = &DD[DIR_000*size_Mat];
		    //////////////////////////////////////////////////////////////////////////
		    (D.f[DIR_000])[k] = c96o1;
		    //////////////////////////////////////////////////////////////////////////
	    }
   }
}