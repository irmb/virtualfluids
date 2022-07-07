#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
#include "math.h"

extern "C" __global__ void LB_Init_Incomp_AD_27(unsigned int* neighborX,
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
            D27.f[E   ] = &DD27[E   *size_Mat];
            D27.f[W   ] = &DD27[W   *size_Mat];
            D27.f[N   ] = &DD27[N   *size_Mat];
            D27.f[S   ] = &DD27[S   *size_Mat];
            D27.f[T   ] = &DD27[T   *size_Mat];
            D27.f[B   ] = &DD27[B   *size_Mat];
            D27.f[NE  ] = &DD27[NE  *size_Mat];
            D27.f[SW  ] = &DD27[SW  *size_Mat];
            D27.f[SE  ] = &DD27[SE  *size_Mat];
            D27.f[NW  ] = &DD27[NW  *size_Mat];
            D27.f[TE  ] = &DD27[TE  *size_Mat];
            D27.f[BW  ] = &DD27[BW  *size_Mat];
            D27.f[BE  ] = &DD27[BE  *size_Mat];
            D27.f[TW  ] = &DD27[TW  *size_Mat];
            D27.f[TN  ] = &DD27[TN  *size_Mat];
            D27.f[BS  ] = &DD27[BS  *size_Mat];
            D27.f[BN  ] = &DD27[BN  *size_Mat];
            D27.f[TS  ] = &DD27[TS  *size_Mat];
            D27.f[dirREST] = &DD27[dirREST*size_Mat];
            D27.f[TNE ] = &DD27[TNE *size_Mat];
            D27.f[TSW ] = &DD27[TSW *size_Mat];
            D27.f[TSE ] = &DD27[TSE *size_Mat];
            D27.f[TNW ] = &DD27[TNW *size_Mat];
            D27.f[BNE ] = &DD27[BNE *size_Mat];
            D27.f[BSW ] = &DD27[BSW *size_Mat];
            D27.f[BSE ] = &DD27[BSE *size_Mat];
            D27.f[BNW ] = &DD27[BNW *size_Mat];
         }
         else
         {
            D27.f[W   ] = &DD27[E   *size_Mat];
            D27.f[E   ] = &DD27[W   *size_Mat];
            D27.f[S   ] = &DD27[N   *size_Mat];
            D27.f[N   ] = &DD27[S   *size_Mat];
            D27.f[B   ] = &DD27[T   *size_Mat];
            D27.f[T   ] = &DD27[B   *size_Mat];
            D27.f[SW  ] = &DD27[NE  *size_Mat];
            D27.f[NE  ] = &DD27[SW  *size_Mat];
            D27.f[NW  ] = &DD27[SE  *size_Mat];
            D27.f[SE  ] = &DD27[NW  *size_Mat];
            D27.f[BW  ] = &DD27[TE  *size_Mat];
            D27.f[TE  ] = &DD27[BW  *size_Mat];
            D27.f[TW  ] = &DD27[BE  *size_Mat];
            D27.f[BE  ] = &DD27[TW  *size_Mat];
            D27.f[BS  ] = &DD27[TN  *size_Mat];
            D27.f[TN  ] = &DD27[BS  *size_Mat];
            D27.f[TS  ] = &DD27[BN  *size_Mat];
            D27.f[BN  ] = &DD27[TS  *size_Mat];
            D27.f[dirREST] = &DD27[dirREST*size_Mat];
            D27.f[BSW ] = &DD27[TNE *size_Mat];
            D27.f[BNE ] = &DD27[TSW *size_Mat];
            D27.f[BNW ] = &DD27[TSE *size_Mat];
            D27.f[BSE ] = &DD27[TNW *size_Mat];
            D27.f[TSW ] = &DD27[BNE *size_Mat];
            D27.f[TNE ] = &DD27[BSW *size_Mat];
            D27.f[TNW ] = &DD27[BSE *size_Mat];
            D27.f[TSE ] = &DD27[BNW *size_Mat];
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

         (D27.f[dirREST])[kzero] =   c8o27* ConcD*(c1o1-cu_sq);
         (D27.f[E   ])[ke   ] =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
         (D27.f[W   ])[kw   ] =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
         (D27.f[N   ])[kn   ] =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
         (D27.f[S   ])[ks   ] =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
         (D27.f[T   ])[kt   ] =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
         (D27.f[B   ])[kb   ] =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
         (D27.f[NE  ])[kne  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         (D27.f[SW  ])[ksw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         (D27.f[SE  ])[kse  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         (D27.f[NW  ])[knw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         (D27.f[TE  ])[kte  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         (D27.f[BW  ])[kbw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         (D27.f[BE  ])[kbe  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         (D27.f[TW  ])[ktw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         (D27.f[TN  ])[ktn  ] =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         (D27.f[BS  ])[kbs  ] =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         (D27.f[BN  ])[kbn  ] =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         (D27.f[TS  ])[kts  ] =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         (D27.f[TNE ])[ktne ] =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         (D27.f[BSW ])[kbsw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         (D27.f[BNE ])[kbne ] =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         (D27.f[TSW ])[ktsw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         (D27.f[TSE ])[ktse ] =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         (D27.f[BNW ])[kbnw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         (D27.f[BSE ])[kbse ] =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         (D27.f[TNW ])[ktnw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      }
   }
}