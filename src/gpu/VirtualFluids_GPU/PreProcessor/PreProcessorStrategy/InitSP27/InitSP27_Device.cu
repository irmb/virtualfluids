#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
#include "math.h"

extern "C" __global__ void LB_Init_SP_27(unsigned int* neighborX,
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
            D.f[BSW ] = &DD[TNE *size_Mat];
            D.f[BNE ] = &DD[TSW *size_Mat];
            D.f[BNW ] = &DD[TSE *size_Mat];
            D.f[BSE ] = &DD[TNW *size_Mat];
            D.f[TSW ] = &DD[BNE *size_Mat];
            D.f[TNE ] = &DD[BSW *size_Mat];
            D.f[TNW ] = &DD[BSE *size_Mat];
            D.f[TSE ] = &DD[BNW *size_Mat];
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

         (D.f[REST])[kzero] =   c8o27* (drho-cu_sq);
         (D.f[E   ])[ke   ] =   c2o27* (drho+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
         (D.f[W   ])[kw   ] =   c2o27* (drho+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
         (D.f[N   ])[kn   ] =   c2o27* (drho+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
         (D.f[S   ])[ks   ] =   c2o27* (drho+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
         (D.f[T   ])[kt   ] =   c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
         (D.f[B   ])[kb   ] =   c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
         (D.f[NE  ])[kne  ] =   c1o54* (drho+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         (D.f[SW  ])[ksw  ] =   c1o54* (drho+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         (D.f[SE  ])[kse  ] =   c1o54* (drho+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         (D.f[NW  ])[knw  ] =   c1o54* (drho+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         (D.f[TE  ])[kte  ] =   c1o54* (drho+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         (D.f[BW  ])[kbw  ] =   c1o54* (drho+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         (D.f[BE  ])[kbe  ] =   c1o54* (drho+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         (D.f[TW  ])[ktw  ] =   c1o54* (drho+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         (D.f[TN  ])[ktn  ] =   c1o54* (drho+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         (D.f[BS  ])[kbs  ] =   c1o54* (drho+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         (D.f[BN  ])[kbn  ] =   c1o54* (drho+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         (D.f[TS  ])[kts  ] =   c1o54* (drho+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         (D.f[TNE ])[ktne ] =   c1o216*(drho+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         (D.f[BSW ])[kbsw ] =   c1o216*(drho+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         (D.f[BNE ])[kbne ] =   c1o216*(drho+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         (D.f[TSW ])[ktsw ] =   c1o216*(drho+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         (D.f[TSE ])[ktse ] =   c1o216*(drho+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         (D.f[BNW ])[kbnw ] =   c1o216*(drho+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         (D.f[BSE ])[kbse ] =   c1o216*(drho+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         (D.f[TNW ])[ktnw ] =   c1o216*(drho+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      }
	  else
	  {
		  //////////////////////////////////////////////////////////////////////////
		  Distributions27 D;
		  D.f[REST] = &DD[REST*size_Mat];
		  //////////////////////////////////////////////////////////////////////////
		  (D.f[REST])[k] = c96o1;
		  //////////////////////////////////////////////////////////////////////////
	  }
   }
}