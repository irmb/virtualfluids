/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include "lbm/constants/NumericConstants.h"
#include "KernelUtilities.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

//////////////////////////////////////////////////////////////////////////////
__global__ void QSlipDevice27(real* DD, 
                                         int* k_Q, 
                                         real* QQ,
                                         unsigned int numberOfBCnodes,
                                         real om1, 
                                         unsigned int* neighborX,
                                         unsigned int* neighborY,
                                         unsigned int* neighborZ,
                                         unsigned int size_Mat, 
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

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[E   * numberOfBCnodes];
      q_dirW   = &QQ[W   * numberOfBCnodes];
      q_dirN   = &QQ[N   * numberOfBCnodes];
      q_dirS   = &QQ[S   * numberOfBCnodes];
      q_dirT   = &QQ[T   * numberOfBCnodes];
      q_dirB   = &QQ[B   * numberOfBCnodes];
      q_dirNE  = &QQ[NE  * numberOfBCnodes];
      q_dirSW  = &QQ[SW  * numberOfBCnodes];
      q_dirSE  = &QQ[SE  * numberOfBCnodes];
      q_dirNW  = &QQ[NW  * numberOfBCnodes];
      q_dirTE  = &QQ[TE  * numberOfBCnodes];
      q_dirBW  = &QQ[BW  * numberOfBCnodes];
      q_dirBE  = &QQ[BE  * numberOfBCnodes];
      q_dirTW  = &QQ[TW  * numberOfBCnodes];
      q_dirTN  = &QQ[TN  * numberOfBCnodes];
      q_dirBS  = &QQ[BS  * numberOfBCnodes];
      q_dirBN  = &QQ[BN  * numberOfBCnodes];
      q_dirTS  = &QQ[TS  * numberOfBCnodes];
      q_dirTNE = &QQ[TNE * numberOfBCnodes];
      q_dirTSW = &QQ[TSW * numberOfBCnodes];
      q_dirTSE = &QQ[TSE * numberOfBCnodes];
      q_dirTNW = &QQ[TNW * numberOfBCnodes];
      q_dirBNE = &QQ[BNE * numberOfBCnodes];
      q_dirBSW = &QQ[BSW * numberOfBCnodes];
      q_dirBSE = &QQ[BSE * numberOfBCnodes];
      q_dirBNW = &QQ[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real f_W    = (D.f[E   ])[ke   ];
      real f_E    = (D.f[W   ])[kw   ];
      real f_S    = (D.f[N   ])[kn   ];
      real f_N    = (D.f[S   ])[ks   ];
      real f_B    = (D.f[T   ])[kt   ];
      real f_T    = (D.f[B   ])[kb   ];
      real f_SW   = (D.f[NE  ])[kne  ];
      real f_NE   = (D.f[SW  ])[ksw  ];
      real f_NW   = (D.f[SE  ])[kse  ];
      real f_SE   = (D.f[NW  ])[knw  ];
      real f_BW   = (D.f[TE  ])[kte  ];
      real f_TE   = (D.f[BW  ])[kbw  ];
      real f_TW   = (D.f[BE  ])[kbe  ];
      real f_BE   = (D.f[TW  ])[ktw  ];
      real f_BS   = (D.f[TN  ])[ktn  ];
      real f_TN   = (D.f[BS  ])[kbs  ];
      real f_TS   = (D.f[BN  ])[kbn  ];
      real f_BN   = (D.f[TS  ])[kts  ];
      real f_BSW  = (D.f[TNE ])[ktne ];
      real f_BNE  = (D.f[TSW ])[ktsw ];
      real f_BNW  = (D.f[TSE ])[ktse ];
      real f_BSE  = (D.f[TNW ])[ktnw ];
      real f_TSW  = (D.f[BNE ])[kbne ];
      real f_TNE  = (D.f[BSW ])[kbsw ];
      real f_TNW  = (D.f[BSE ])[kbse ];
      real f_TSE  = (D.f[BNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[REST])[kzero]); 

      vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W); 
         

      vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S); 

      vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[REST])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real fac = c1o1;//c99o100;
	  real VeloX = fac*vx1;
	  real VeloY = fac*vx2;
	  real VeloZ = fac*vx3;
	  bool x = false;
	  bool y = false;
	  bool z = false;

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = c0o1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 x = true;
         feq=c2o27* (drho+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq); 
         (D.f[W])[kw]=(c1o1-q)/(c1o1+q)*(f_E-feq*om1)/(c1o1-om1)+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q);
         //(D.f[W])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = c0o1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 x = true;
         feq=c2o27* (drho+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq); 
         (D.f[E])[ke]=(c1o1-q)/(c1o1+q)*(f_W-feq*om1)/(c1o1-om1)+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q);
         //(D.f[E])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
		 VeloY = c0o1;
	     VeloZ = fac*vx3;
		 y = true;
         feq=c2o27* (drho+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq); 
         (D.f[S])[ks]=(c1o1-q)/(c1o1+q)*(f_N-feq*om1)/(c1o1-om1)+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q);
         //(D.f[S])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
		 VeloY = c0o1;
	     VeloZ = fac*vx3;
		 y = true;
         feq=c2o27* (drho+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq); 
         (D.f[N])[kn]=(c1o1-q)/(c1o1+q)*(f_S-feq*om1)/(c1o1-om1)+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q);
         //(D.f[N])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
		 VeloZ = c0o1;
		 z = true;
         feq=c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq); 
         (D.f[B])[kb]=(c1o1-q)/(c1o1+q)*(f_T-feq*om1)/(c1o1-om1)+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q);
         //(D.f[B])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
		 VeloZ = c0o1;
		 z = true;
         feq=c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq); 
         (D.f[T])[kt]=(c1o1-q)/(c1o1+q)*(f_B-feq*om1)/(c1o1-om1)+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q);
         //(D.f[T])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
         feq=c1o54* (drho+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         (D.f[SW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-feq*om1)/(c1o1-om1)+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q);
         //(D.f[SW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
         feq=c1o54* (drho+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         (D.f[NE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-feq*om1)/(c1o1-om1)+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
         //(D.f[NE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
         feq=c1o54* (drho+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         (D.f[NW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-feq*om1)/(c1o1-om1)+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
         //(D.f[NW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
         feq=c1o54* (drho+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         (D.f[SE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-feq*om1)/(c1o1-om1)+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
         //(D.f[SE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         (D.f[BW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-feq*om1)/(c1o1-om1)+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
         //(D.f[BW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         (D.f[TE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-feq*om1)/(c1o1-om1)+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
         //(D.f[TE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         (D.f[TW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-feq*om1)/(c1o1-om1)+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
         //(D.f[TW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         (D.f[BE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-feq*om1)/(c1o1-om1)+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
         //(D.f[BE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         (D.f[BS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-feq*om1)/(c1o1-om1)+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
         //(D.f[BS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         (D.f[TN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-feq*om1)/(c1o1-om1)+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q);
         //(D.f[TN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         (D.f[TS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-feq*om1)/(c1o1-om1)+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
         //(D.f[TS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         (D.f[BN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-feq*om1)/(c1o1-om1)+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q);
         //(D.f[BN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         (D.f[BSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-feq*om1)/(c1o1-om1)+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
         //(D.f[BSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         (D.f[TNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-feq*om1)/(c1o1-om1)+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
         //(D.f[TNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         (D.f[TSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-feq*om1)/(c1o1-om1)+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
         //(D.f[TSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         (D.f[BNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-feq*om1)/(c1o1-om1)+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
         //(D.f[BNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         (D.f[BNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-feq*om1)/(c1o1-om1)+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
         //(D.f[BNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         (D.f[TSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-feq*om1)/(c1o1-om1)+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
         //(D.f[TSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         (D.f[TNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-feq*om1)/(c1o1-om1)+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
         //(D.f[TNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         (D.f[BSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-feq*om1)/(c1o1-om1)+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);
         //(D.f[BSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QSlipDeviceComp27(
                                    real* distributions, 
                                    int* subgridDistanceIndices, 
                                    real* subgridDistances,
                                    unsigned int numberOfBCnodes,
                                    real omega, 
                                    unsigned int* neighborX,
                                    unsigned int* neighborY,
                                    unsigned int* neighborZ,
                                    unsigned int numberOfLBnodes, 
                                    bool isEvenTimestep)
{
   //! The slip boundary condition is executed in the following steps
   //!
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned  x = threadIdx.x;  // global x-index 
   const unsigned  y = blockIdx.x;   // global y-index 
   const unsigned  z = blockIdx.y;   // global z-index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;

   if(k < numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist;
      getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local subgrid distances (q's)
      //!
      SubgridDistances27 subgridD;
      getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);
      
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int indexOfBCnode  = subgridDistanceIndices[k];
      unsigned int kzero= indexOfBCnode;
      unsigned int ke   = indexOfBCnode;
      unsigned int kw   = neighborX[indexOfBCnode];
      unsigned int kn   = indexOfBCnode;
      unsigned int ks   = neighborY[indexOfBCnode];
      unsigned int kt   = indexOfBCnode;
      unsigned int kb   = neighborZ[indexOfBCnode];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = indexOfBCnode;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = indexOfBCnode;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = indexOfBCnode;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = indexOfBCnode;
      unsigned int kbsw = neighborZ[ksw];
      
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local distributions
      //!
      real f_W    = (dist.f[E   ])[ke   ];
      real f_E    = (dist.f[W   ])[kw   ];
      real f_S    = (dist.f[N   ])[kn   ];
      real f_N    = (dist.f[S   ])[ks   ];
      real f_B    = (dist.f[T   ])[kt   ];
      real f_T    = (dist.f[B   ])[kb   ];
      real f_SW   = (dist.f[NE  ])[kne  ];
      real f_NE   = (dist.f[SW  ])[ksw  ];
      real f_NW   = (dist.f[SE  ])[kse  ];
      real f_SE   = (dist.f[NW  ])[knw  ];
      real f_BW   = (dist.f[TE  ])[kte  ];
      real f_TE   = (dist.f[BW  ])[kbw  ];
      real f_TW   = (dist.f[BE  ])[kbe  ];
      real f_BE   = (dist.f[TW  ])[ktw  ];
      real f_BS   = (dist.f[TN  ])[ktn  ];
      real f_TN   = (dist.f[BS  ])[kbs  ];
      real f_TS   = (dist.f[BN  ])[kbn  ];
      real f_BN   = (dist.f[TS  ])[kts  ];
      real f_BSW  = (dist.f[TNE ])[ktne ];
      real f_BNE  = (dist.f[TSW ])[ktsw ];
      real f_BNW  = (dist.f[TSE ])[ktse ];
      real f_BSE  = (dist.f[TNW ])[ktnw ];
      real f_TSW  = (dist.f[BNE ])[kbne ];
      real f_TNE  = (dist.f[BSW ])[kbsw ];
      real f_TNW  = (dist.f[BSE ])[kbse ];
      real f_TSE  = (dist.f[BNW ])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Calculate macroscopic quantities
      //!
      real drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                  f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[REST])[kzero]); 

      real vx1  = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                   (f_E - f_W)) / (c1o1 + drho);

      real vx2  = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                   (f_N - f_S)) / (c1o1 + drho);

      real vx3  = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                   (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                   (f_T - f_B)) / (c1o1 + drho);

      real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3) * (c1o1 + drho);

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Multiply the local velocities by the slipLength
      //!
      real slipLength = c1o1;
      real VeloX = slipLength*vx1;
      real VeloY = slipLength*vx2;
      real VeloZ = slipLength*vx3;

      ////////////////////////////////////////////////////////////////////////////////
      //! - Update distributions with subgrid distance (q) between zero and one
      //!
      real feq, q, velocityLB, velocityBC;

      bool x = false;
      bool y = false;
      bool z = false;

      q = (subgridD.q[E])[k];
      if (q>=c0o1 && q<=c1o1)  // only update distribution for q between zero and one
      {
         VeloX = c0o1;
         x = true;

         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloX;
         (dist.f[W])[kw] = getInterpolatedDistributionForVeloBC(q, f_E, f_W, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[W])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = c0o1;
         x = true;

         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloX;
         (dist.f[E])[ke] = getInterpolatedDistributionForVeloBC(q, f_W, f_E, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[N])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = c0o1;
         y = true;

         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloY;
         (dist.f[S])[ks] = getInterpolatedDistributionForVeloBC(q, f_N, f_S, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[S])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = c0o1;
         y = true;

         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloY;
         (dist.f[N])[kn] = getInterpolatedDistributionForVeloBC(q, f_S, f_N, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[T])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = c0o1;
         z = true;

         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloZ;
         (dist.f[B])[kb] = getInterpolatedDistributionForVeloBC(q, f_T, f_B, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[B])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = c0o1;
         z = true;

         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloZ;
         (dist.f[T])[kt] = getInterpolatedDistributionForVeloBC(q, f_B, f_T, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[NE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloY;
         (dist.f[SW])[ksw] = getInterpolatedDistributionForVeloBC(q, f_NE, f_SW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[SW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloY;
         (dist.f[NE])[kne] = getInterpolatedDistributionForVeloBC(q, f_SW, f_NE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[SE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloY;
         (dist.f[NW])[knw] = getInterpolatedDistributionForVeloBC(q, f_SE, f_NW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[NW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloY;
         (dist.f[SE])[kse] = getInterpolatedDistributionForVeloBC(q, f_NW, f_SE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[TE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloZ;
         (dist.f[BW])[kbw] = getInterpolatedDistributionForVeloBC(q, f_TE, f_BW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[BW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
        VeloX = slipLength*vx1;
        VeloZ = slipLength*vx3;
        if (x == true) VeloX = c0o1;
        if (z == true) VeloZ = c0o1;

         velocityLB = -vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloZ;
         (dist.f[TE])[kte] = getInterpolatedDistributionForVeloBC(q, f_BW, f_TE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[BE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloZ;
         (dist.f[TW])[ktw] = getInterpolatedDistributionForVeloBC(q, f_BE, f_TW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[TW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloZ;
         (dist.f[BE])[kbe] = getInterpolatedDistributionForVeloBC(q, f_TW, f_BE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[TN])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY + VeloZ;
         (dist.f[BS])[kbs] = getInterpolatedDistributionForVeloBC(q, f_TN, f_BS, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[BS])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY - VeloZ;
         (dist.f[TN])[ktn] = getInterpolatedDistributionForVeloBC(q, f_BS, f_TN, feq, omega, velocityBC, c1o54);
      }


      q = (subgridD.q[BN])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY - VeloZ;
         (dist.f[TS])[kts] = getInterpolatedDistributionForVeloBC(q, f_BN, f_TS, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[TS])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY + VeloZ;
         (dist.f[BN])[kbn] = getInterpolatedDistributionForVeloBC(q, f_TS, f_BN, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[TNE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY + VeloZ;
         (dist.f[BSW])[kbsw] = getInterpolatedDistributionForVeloBC(q, f_TNE, f_BSW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[BSW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY - VeloZ;
         (dist.f[TNE])[ktne] = getInterpolatedDistributionForVeloBC(q, f_BSW, f_TNE, feq, omega, velocityBC, c1o216);
      }


      q = (subgridD.q[BNE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY - VeloZ;
         (dist.f[TSW])[ktsw] = getInterpolatedDistributionForVeloBC(q, f_BNE, f_TSW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[TSW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY + VeloZ;
         (dist.f[BNE])[kbne] = getInterpolatedDistributionForVeloBC(q, f_TSW, f_BNE, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[TSE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY + VeloZ;
         (dist.f[BNW])[kbnw] = getInterpolatedDistributionForVeloBC(q, f_TSE, f_BNW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[BNW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY - VeloZ;
         (dist.f[TSE])[ktse] = getInterpolatedDistributionForVeloBC(q, f_BNW, f_TSE, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[BSE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY - VeloZ;
         (dist.f[TNW])[ktnw] = getInterpolatedDistributionForVeloBC(q, f_BSE, f_TNW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[TNW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY + VeloZ;
         (dist.f[BSE])[kbse] = getInterpolatedDistributionForVeloBC(q, f_TNW, f_BSE, feq, omega, velocityBC, c1o216);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
__global__ void QSlipDeviceComp27TurbViscosity(real* DD, 
											 int* k_Q, 
											 real* QQ,
											 unsigned int numberOfBCnodes,
											 real om1, 
											 unsigned int* neighborX,
											 unsigned int* neighborY,
											 unsigned int* neighborZ,
                                  real* turbViscosity,
											 unsigned int size_Mat, 
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

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[E   * numberOfBCnodes];
      q_dirW   = &QQ[W   * numberOfBCnodes];
      q_dirN   = &QQ[N   * numberOfBCnodes];
      q_dirS   = &QQ[S   * numberOfBCnodes];
      q_dirT   = &QQ[T   * numberOfBCnodes];
      q_dirB   = &QQ[B   * numberOfBCnodes];
      q_dirNE  = &QQ[NE  * numberOfBCnodes];
      q_dirSW  = &QQ[SW  * numberOfBCnodes];
      q_dirSE  = &QQ[SE  * numberOfBCnodes];
      q_dirNW  = &QQ[NW  * numberOfBCnodes];
      q_dirTE  = &QQ[TE  * numberOfBCnodes];
      q_dirBW  = &QQ[BW  * numberOfBCnodes];
      q_dirBE  = &QQ[BE  * numberOfBCnodes];
      q_dirTW  = &QQ[TW  * numberOfBCnodes];
      q_dirTN  = &QQ[TN  * numberOfBCnodes];
      q_dirBS  = &QQ[BS  * numberOfBCnodes];
      q_dirBN  = &QQ[BN  * numberOfBCnodes];
      q_dirTS  = &QQ[TS  * numberOfBCnodes];
      q_dirTNE = &QQ[TNE * numberOfBCnodes];
      q_dirTSW = &QQ[TSW * numberOfBCnodes];
      q_dirTSE = &QQ[TSE * numberOfBCnodes];
      q_dirTNW = &QQ[TNW * numberOfBCnodes];
      q_dirBNE = &QQ[BNE * numberOfBCnodes];
      q_dirBSW = &QQ[BSW * numberOfBCnodes];
      q_dirBSE = &QQ[BSE * numberOfBCnodes];
      q_dirBNW = &QQ[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      
      ////////////////////////////////////////////////////////////////////////////////
      real f_W    = (D.f[E   ])[ke   ];
      real f_E    = (D.f[W   ])[kw   ];
      real f_S    = (D.f[N   ])[kn   ];
      real f_N    = (D.f[S   ])[ks   ];
      real f_B    = (D.f[T   ])[kt   ];
      real f_T    = (D.f[B   ])[kb   ];
      real f_SW   = (D.f[NE  ])[kne  ];
      real f_NE   = (D.f[SW  ])[ksw  ];
      real f_NW   = (D.f[SE  ])[kse  ];
      real f_SE   = (D.f[NW  ])[knw  ];
      real f_BW   = (D.f[TE  ])[kte  ];
      real f_TE   = (D.f[BW  ])[kbw  ];
      real f_TW   = (D.f[BE  ])[kbe  ];
      real f_BE   = (D.f[TW  ])[ktw  ];
      real f_BS   = (D.f[TN  ])[ktn  ];
      real f_TN   = (D.f[BS  ])[kbs  ];
      real f_TS   = (D.f[BN  ])[kbn  ];
      real f_BN   = (D.f[TS  ])[kts  ];
      real f_BSW  = (D.f[TNE ])[ktne ];
      real f_BNE  = (D.f[TSW ])[ktsw ];
      real f_BNW  = (D.f[TSE ])[ktse ];
      real f_BSE  = (D.f[TNW ])[ktnw ];
      real f_TSW  = (D.f[BNE ])[kbne ];
      real f_TNE  = (D.f[BSW ])[kbsw ];
      real f_TNW  = (D.f[BSE ])[kbse ];
      real f_TSE  = (D.f[BNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[REST])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (c1o1 + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (c1o1 + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (c1o1 + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[REST])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real om_turb = om1 / (c1o1 + c3o1*om1*max(c0o1, turbViscosity[k_Q[k]]));
     
     real fac = c1o1;//c99o100;
	  real VeloX = fac*vx1;
	  real VeloY = fac*vx2;
	  real VeloZ = fac*vx3;
	  bool x = false;
	  bool y = false;
	  bool z = false;

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = c0o1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 x = true;
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[W])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q) - c2o27 * drho;
         //feq=c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq); 
         //(D.f[W])[kw]=(one-q)/(one+q)*(f_E-feq*om1)/(one-om1)+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);
         //(D.f[W])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = c0o1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 x = true;
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[E])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q) - c2o27 * drho;
         //feq=c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
         //(D.f[E])[ke]=(one-q)/(one+q)*(f_W-feq*om_turb)/(one-om_turb)+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);
         //(D.f[E])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
		 VeloY = c0o1;
	     VeloZ = fac*vx3;
		 y = true;
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[S])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q) - c2o27 * drho;
         //feq=c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
         //(D.f[S])[ks]=(one-q)/(one+q)*(f_N-feq*om_turb)/(one-om_turb)+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);
         //(D.f[S])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
		 VeloY = c0o1;
	     VeloZ = fac*vx3;
		 y = true;
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[N])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q) - c2o27 * drho;
         //feq=c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
         //(D.f[N])[kn]=(one-q)/(one+q)*(f_S-feq*om_turb)/(one-om_turb)+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);
         //(D.f[N])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
		 VeloZ = c0o1;
		 z = true;
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[B])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q) - c2o27 * drho;
         //feq=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq); 
         //(D.f[B])[kb]=(one-q)/(one+q)*(f_T-feq*om_turb)/(one-om_turb)+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);
         //(D.f[B])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
		 VeloZ = c0o1;
		 z = true;
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[T])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q) - c2o27 * drho;
         //feq=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
         //(D.f[T])[kt]=(one-q)/(one+q)*(f_B-feq*om_turb)/(one-om_turb)+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);
         //(D.f[T])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         //(D.f[SW])[ksw]=(one-q)/(one+q)*(f_NE-feq*om_turb)/(one-om_turb)+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);
         //(D.f[SW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         //(D.f[NE])[kne]=(one-q)/(one+q)*(f_SW-feq*om_turb)/(one-om_turb)+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);
         //(D.f[NE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         //(D.f[NW])[knw]=(one-q)/(one+q)*(f_SE-feq*om_turb)/(one-om_turb)+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);
         //(D.f[NW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         //(D.f[SE])[kse]=(one-q)/(one+q)*(f_NW-feq*om_turb)/(one-om_turb)+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);
         //(D.f[SE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (z == true) VeloZ = c0o1;
      //  if (k==10000) printf("AFTER x: %u \t  y: %u \t z: %u \n  VeloX: %f \t VeloY: %f \t VeloZ: %f \n\n", x,y,z, VeloX,VeloY,VeloZ);
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         //(D.f[BW])[kbw]=(one-q)/(one+q)*(f_TE-feq*om_turb)/(one-om_turb)+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);
         //(D.f[BW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         //(D.f[TE])[kte]=(one-q)/(one+q)*(f_BW-feq*om_turb)/(one-om_turb)+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);
         //(D.f[TE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         //(D.f[TW])[ktw]=(one-q)/(one+q)*(f_BE-feq*om_turb)/(one-om_turb)+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);
         //(D.f[TW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         //(D.f[BE])[kbe]=(one-q)/(one+q)*(f_TW-feq*om_turb)/(one-om_turb)+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);
         //(D.f[BE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         //(D.f[BS])[kbs]=(one-q)/(one+q)*(f_TN-feq*om_turb)/(one-om_turb)+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);
         //(D.f[BS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         //(D.f[TN])[ktn]=(one-q)/(one+q)*(f_BS-feq*om_turb)/(one-om_turb)+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);
         //(D.f[TN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         //(D.f[TS])[kts]=(one-q)/(one+q)*(f_BN-feq*om_turb)/(one-om_turb)+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);
         //(D.f[TS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         //(D.f[BN])[kbn]=(one-q)/(one+q)*(f_TS-feq*om_turb)/(one-om_turb)+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);
         //(D.f[BN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         //(D.f[BSW])[kbsw]=(one-q)/(one+q)*(f_TNE-feq*om_turb)/(one-om_turb)+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);
         //(D.f[BSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         //(D.f[TNE])[ktne]=(one-q)/(one+q)*(f_BSW-feq*om_turb)/(one-om_turb)+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);
         //(D.f[TNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         //(D.f[TSW])[ktsw]=(one-q)/(one+q)*(f_BNE-feq*om_turb)/(one-om_turb)+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);
         //(D.f[TSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         //(D.f[BNE])[kbne]=(one-q)/(one+q)*(f_TSW-feq*om_turb)/(one-om_turb)+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);
         //(D.f[BNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         //(D.f[BNW])[kbnw]=(one-q)/(one+q)*(f_TSE-feq*om_turb)/(one-om_turb)+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);
         //(D.f[BNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         //(D.f[TSE])[ktse]=(one-q)/(one+q)*(f_BNW-feq*om_turb)/(one-om_turb)+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);
         //(D.f[TSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         //(D.f[TNW])[ktnw]=(one-q)/(one+q)*(f_BSE-feq*om_turb)/(one-om_turb)+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);
         //(D.f[TNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (y == true) VeloY = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         //(D.f[BSE])[kbse]=(one-q)/(one+q)*(f_TNW-feq*om_turb)/(one-om_turb)+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);
         //(D.f[BSE])[kbse]=zero;
      }
   }
}






































//////////////////////////////////////////////////////////////////////////////
__global__ void QSlipGeomDeviceComp27(real* DD, 
												 int* k_Q, 
												 real* QQ,
												 unsigned int  numberOfBCnodes,
												 real om1, 
												 real* NormalX,
												 real* NormalY,
												 real* NormalZ,
												 unsigned int* neighborX,
												 unsigned int* neighborY,
												 unsigned int* neighborZ,
												 unsigned int size_Mat, 
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

   if(k< numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[E   * numberOfBCnodes];
      q_dirW   = &QQ[W   * numberOfBCnodes];
      q_dirN   = &QQ[N   * numberOfBCnodes];
      q_dirS   = &QQ[S   * numberOfBCnodes];
      q_dirT   = &QQ[T   * numberOfBCnodes];
      q_dirB   = &QQ[B   * numberOfBCnodes];
      q_dirNE  = &QQ[NE  * numberOfBCnodes];
      q_dirSW  = &QQ[SW  * numberOfBCnodes];
      q_dirSE  = &QQ[SE  * numberOfBCnodes];
      q_dirNW  = &QQ[NW  * numberOfBCnodes];
      q_dirTE  = &QQ[TE  * numberOfBCnodes];
      q_dirBW  = &QQ[BW  * numberOfBCnodes];
      q_dirBE  = &QQ[BE  * numberOfBCnodes];
      q_dirTW  = &QQ[TW  * numberOfBCnodes];
      q_dirTN  = &QQ[TN  * numberOfBCnodes];
      q_dirBS  = &QQ[BS  * numberOfBCnodes];
      q_dirBN  = &QQ[BN  * numberOfBCnodes];
      q_dirTS  = &QQ[TS  * numberOfBCnodes];
      q_dirTNE = &QQ[TNE * numberOfBCnodes];
      q_dirTSW = &QQ[TSW * numberOfBCnodes];
      q_dirTSE = &QQ[TSE * numberOfBCnodes];
      q_dirTNW = &QQ[TNW * numberOfBCnodes];
      q_dirBNE = &QQ[BNE * numberOfBCnodes];
      q_dirBSW = &QQ[BSW * numberOfBCnodes];
      q_dirBSE = &QQ[BSE * numberOfBCnodes];
      q_dirBNW = &QQ[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      real *nx_dirE,   *nx_dirW,   *nx_dirN,   *nx_dirS,   *nx_dirT,   *nx_dirB, 
              *nx_dirNE,  *nx_dirSW,  *nx_dirSE,  *nx_dirNW,  *nx_dirTE,  *nx_dirBW,
              *nx_dirBE,  *nx_dirTW,  *nx_dirTN,  *nx_dirBS,  *nx_dirBN,  *nx_dirTS,
              *nx_dirTNE, *nx_dirTSW, *nx_dirTSE, *nx_dirTNW, *nx_dirBNE, *nx_dirBSW,
              *nx_dirBSE, *nx_dirBNW; 
      nx_dirE   = &NormalX[E   * numberOfBCnodes];
      nx_dirW   = &NormalX[W   * numberOfBCnodes];
      nx_dirN   = &NormalX[N   * numberOfBCnodes];
      nx_dirS   = &NormalX[S   * numberOfBCnodes];
      nx_dirT   = &NormalX[T   * numberOfBCnodes];
      nx_dirB   = &NormalX[B   * numberOfBCnodes];
      nx_dirNE  = &NormalX[NE  * numberOfBCnodes];
      nx_dirSW  = &NormalX[SW  * numberOfBCnodes];
      nx_dirSE  = &NormalX[SE  * numberOfBCnodes];
      nx_dirNW  = &NormalX[NW  * numberOfBCnodes];
      nx_dirTE  = &NormalX[TE  * numberOfBCnodes];
      nx_dirBW  = &NormalX[BW  * numberOfBCnodes];
      nx_dirBE  = &NormalX[BE  * numberOfBCnodes];
      nx_dirTW  = &NormalX[TW  * numberOfBCnodes];
      nx_dirTN  = &NormalX[TN  * numberOfBCnodes];
      nx_dirBS  = &NormalX[BS  * numberOfBCnodes];
      nx_dirBN  = &NormalX[BN  * numberOfBCnodes];
      nx_dirTS  = &NormalX[TS  * numberOfBCnodes];
      nx_dirTNE = &NormalX[TNE * numberOfBCnodes];
      nx_dirTSW = &NormalX[TSW * numberOfBCnodes];
      nx_dirTSE = &NormalX[TSE * numberOfBCnodes];
      nx_dirTNW = &NormalX[TNW * numberOfBCnodes];
      nx_dirBNE = &NormalX[BNE * numberOfBCnodes];
      nx_dirBSW = &NormalX[BSW * numberOfBCnodes];
      nx_dirBSE = &NormalX[BSE * numberOfBCnodes];
      nx_dirBNW = &NormalX[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      real *ny_dirE,   *ny_dirW,   *ny_dirN,   *ny_dirS,   *ny_dirT,   *ny_dirB, 
              *ny_dirNE,  *ny_dirSW,  *ny_dirSE,  *ny_dirNW,  *ny_dirTE,  *ny_dirBW,
              *ny_dirBE,  *ny_dirTW,  *ny_dirTN,  *ny_dirBS,  *ny_dirBN,  *ny_dirTS,
              *ny_dirTNE, *ny_dirTSW, *ny_dirTSE, *ny_dirTNW, *ny_dirBNE, *ny_dirBSW,
              *ny_dirBSE, *ny_dirBNW; 
      ny_dirE   = &NormalY[E   * numberOfBCnodes];
      ny_dirW   = &NormalY[W   * numberOfBCnodes];
      ny_dirN   = &NormalY[N   * numberOfBCnodes];
      ny_dirS   = &NormalY[S   * numberOfBCnodes];
      ny_dirT   = &NormalY[T   * numberOfBCnodes];
      ny_dirB   = &NormalY[B   * numberOfBCnodes];
      ny_dirNE  = &NormalY[NE  * numberOfBCnodes];
      ny_dirSW  = &NormalY[SW  * numberOfBCnodes];
      ny_dirSE  = &NormalY[SE  * numberOfBCnodes];
      ny_dirNW  = &NormalY[NW  * numberOfBCnodes];
      ny_dirTE  = &NormalY[TE  * numberOfBCnodes];
      ny_dirBW  = &NormalY[BW  * numberOfBCnodes];
      ny_dirBE  = &NormalY[BE  * numberOfBCnodes];
      ny_dirTW  = &NormalY[TW  * numberOfBCnodes];
      ny_dirTN  = &NormalY[TN  * numberOfBCnodes];
      ny_dirBS  = &NormalY[BS  * numberOfBCnodes];
      ny_dirBN  = &NormalY[BN  * numberOfBCnodes];
      ny_dirTS  = &NormalY[TS  * numberOfBCnodes];
      ny_dirTNE = &NormalY[TNE * numberOfBCnodes];
      ny_dirTSW = &NormalY[TSW * numberOfBCnodes];
      ny_dirTSE = &NormalY[TSE * numberOfBCnodes];
      ny_dirTNW = &NormalY[TNW * numberOfBCnodes];
      ny_dirBNE = &NormalY[BNE * numberOfBCnodes];
      ny_dirBSW = &NormalY[BSW * numberOfBCnodes];
      ny_dirBSE = &NormalY[BSE * numberOfBCnodes];
      ny_dirBNW = &NormalY[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      real *nz_dirE,   *nz_dirW,   *nz_dirN,   *nz_dirS,   *nz_dirT,   *nz_dirB, 
              *nz_dirNE,  *nz_dirSW,  *nz_dirSE,  *nz_dirNW,  *nz_dirTE,  *nz_dirBW,
              *nz_dirBE,  *nz_dirTW,  *nz_dirTN,  *nz_dirBS,  *nz_dirBN,  *nz_dirTS,
              *nz_dirTNE, *nz_dirTSW, *nz_dirTSE, *nz_dirTNW, *nz_dirBNE, *nz_dirBSW,
              *nz_dirBSE, *nz_dirBNW; 
      nz_dirE   = &NormalZ[E   * numberOfBCnodes];
      nz_dirW   = &NormalZ[W   * numberOfBCnodes];
      nz_dirN   = &NormalZ[N   * numberOfBCnodes];
      nz_dirS   = &NormalZ[S   * numberOfBCnodes];
      nz_dirT   = &NormalZ[T   * numberOfBCnodes];
      nz_dirB   = &NormalZ[B   * numberOfBCnodes];
      nz_dirNE  = &NormalZ[NE  * numberOfBCnodes];
      nz_dirSW  = &NormalZ[SW  * numberOfBCnodes];
      nz_dirSE  = &NormalZ[SE  * numberOfBCnodes];
      nz_dirNW  = &NormalZ[NW  * numberOfBCnodes];
      nz_dirTE  = &NormalZ[TE  * numberOfBCnodes];
      nz_dirBW  = &NormalZ[BW  * numberOfBCnodes];
      nz_dirBE  = &NormalZ[BE  * numberOfBCnodes];
      nz_dirTW  = &NormalZ[TW  * numberOfBCnodes];
      nz_dirTN  = &NormalZ[TN  * numberOfBCnodes];
      nz_dirBS  = &NormalZ[BS  * numberOfBCnodes];
      nz_dirBN  = &NormalZ[BN  * numberOfBCnodes];
      nz_dirTS  = &NormalZ[TS  * numberOfBCnodes];
      nz_dirTNE = &NormalZ[TNE * numberOfBCnodes];
      nz_dirTSW = &NormalZ[TSW * numberOfBCnodes];
      nz_dirTSE = &NormalZ[TSE * numberOfBCnodes];
      nz_dirTNW = &NormalZ[TNW * numberOfBCnodes];
      nz_dirBNE = &NormalZ[BNE * numberOfBCnodes];
      nz_dirBSW = &NormalZ[BSW * numberOfBCnodes];
      nz_dirBSE = &NormalZ[BSE * numberOfBCnodes];
      nz_dirBNW = &NormalZ[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real f_W    = (D.f[E   ])[ke   ];
      real f_E    = (D.f[W   ])[kw   ];
      real f_S    = (D.f[N   ])[kn   ];
      real f_N    = (D.f[S   ])[ks   ];
      real f_B    = (D.f[T   ])[kt   ];
      real f_T    = (D.f[B   ])[kb   ];
      real f_SW   = (D.f[NE  ])[kne  ];
      real f_NE   = (D.f[SW  ])[ksw  ];
      real f_NW   = (D.f[SE  ])[kse  ];
      real f_SE   = (D.f[NW  ])[knw  ];
      real f_BW   = (D.f[TE  ])[kte  ];
      real f_TE   = (D.f[BW  ])[kbw  ];
      real f_TW   = (D.f[BE  ])[kbe  ];
      real f_BE   = (D.f[TW  ])[ktw  ];
      real f_BS   = (D.f[TN  ])[ktn  ];
      real f_TN   = (D.f[BS  ])[kbs  ];
      real f_TS   = (D.f[BN  ])[kbn  ];
      real f_BN   = (D.f[TS  ])[kts  ];
      real f_BSW  = (D.f[TNE ])[ktne ];
      real f_BNE  = (D.f[TSW ])[ktsw ];
      real f_BNW  = (D.f[TSE ])[ktse ];
      real f_BSE  = (D.f[TNW ])[ktnw ];
      real f_TSW  = (D.f[BNE ])[kbne ];
      real f_TNE  = (D.f[BSW ])[kbsw ];
      real f_TNW  = (D.f[BSE ])[kbse ];
      real f_TSE  = (D.f[BNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[REST])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (c1o1 + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (c1o1 + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (c1o1 + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real VeloX = vx1;
	  real VeloY = vx2;
	  real VeloZ = vx3;
	  real fac = c0o1;//0.5;
 	  real phi = c0o1;
	  //real alpha = c1o100;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real kxyFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_SW+f_BSW+f_TSW-f_NW-f_BNW-f_TNW-f_SE-f_BSE-f_TSE+f_NE+f_BNE+f_TNE ) / (c1o1 + drho) - ((vx1*vx2)));
      real kyzFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_BS+f_BSE+f_BSW-f_TS-f_TSE-f_TSW-f_BN-f_BNE-f_BNW+f_TN+f_TNE+f_TNW ) / (c1o1 + drho) - ((vx2*vx3)));
      real kxzFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_BW+f_BSW+f_BNW-f_TW-f_TSW-f_TNW-f_BE-f_BSE-f_BNE+f_TE+f_TSE+f_TNE ) / (c1o1 + drho) - ((vx1*vx3)));

	  real kxxFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_E+f_NE+f_SE+f_TE+f_BE+f_W+f_NW+f_SW+f_TW+f_BW+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (c1o1 + drho) - ((c1o3*drho + vx1*vx1)));
	  real kyyFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_N+f_NE+f_NW+f_TN+f_BN+f_S+f_SE+f_SW+f_TS+f_BS+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (c1o1 + drho) - ((c1o3*drho + vx2*vx2)));
	  real kzzFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_T+f_TE+f_TW+f_TN+f_BS+f_B+f_BE+f_BW+f_BN+f_BS+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (c1o1 + drho) - ((c1o3*drho + vx3*vx3)));

	  real magS = sqrtf(kxyFromfcNEQ*kxyFromfcNEQ + kyzFromfcNEQ*kyzFromfcNEQ + kxzFromfcNEQ*kxzFromfcNEQ + kxxFromfcNEQ*kxxFromfcNEQ + kyyFromfcNEQ*kyyFromfcNEQ + kzzFromfcNEQ*kzzFromfcNEQ);

	  //fac = fac * magS / (c1o3 * (one / om1 - c1o2));
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //real *facAst = &QQ[REST * numberOfBCnodes];

	  //fac = fac * alpha + facAst[k] * (one - alpha);
	  //facAst[k] = fac;
	  //(&QQ[REST * numberOfBCnodes])[KQK] = fac;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ////real uk = sqrtf(vx1*vx1 + vx2*vx2 + vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //real phi = expf(magS/0.01f) - one;
	  //phi = (phi > one) ? one:phi;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //real C = five;
	  //real kappa = 0.41f;
	  //real phi = (C * kappa * c1o2 * logf(magS / (c1o3 * (one / om1 - c1o2))) - one) / (C * kappa * c1o2 * logf(magS / (c1o3 * (one / om1 - c1o2))));
	  //phi = (phi < zero) ? zero:phi;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //real sum = zero, count = zero;
   //   q = q_dirE   [k]; if (q>=zero && q<=one) sum += (q *   nx_dirE[k] ); count += one;
   //   q = q_dirW   [k]; if (q>=zero && q<=one) sum += (q * (-nx_dirW[k])); count += one;
   //   q = q_dirN   [k]; if (q>=zero && q<=one) sum += (q *   nx_dirN[k] ); count += one;
   //   q = q_dirS   [k]; if (q>=zero && q<=one) sum += (q * (-nx_dirS[k])); count += one;
   //   q = q_dirT   [k]; if (q>=zero && q<=one) sum += (q *   nx_dirT[k] ); count += one;
   //   q = q_dirB   [k]; if (q>=zero && q<=one) sum += (q * (-nx_dirB[k])); count += one;
   //   q = q_dirNE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirNE[k]  + ny_dirNE[k])/(sqrtf(two))); count += one;
   //   q = q_dirSW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirSW[k]) - ny_dirSW[k])/(sqrtf(two))); count += one;
   //   q = q_dirSE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirSE[k]  - ny_dirSE[k])/(sqrtf(two))); count += one;
   //   q = q_dirNW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirNW[k]) + ny_dirNW[k])/(sqrtf(two))); count += one;
   //   q = q_dirTE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirTE[k]  + nz_dirTE[k])/(sqrtf(two))); count += one;
   //   q = q_dirBW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirBW[k]) - nz_dirBW[k])/(sqrtf(two))); count += one;
   //   q = q_dirBE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirBE[k]  - nz_dirBE[k])/(sqrtf(two))); count += one;
   //   q = q_dirTW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirTW[k]) + nz_dirTW[k])/(sqrtf(two))); count += one;
   //   q = q_dirTN  [k]; if (q>=zero && q<=one) sum += (q * (  ny_dirTN[k]  + nz_dirTN[k])/(sqrtf(two))); count += one;
   //   q = q_dirBS  [k]; if (q>=zero && q<=one) sum += (q * ((-ny_dirBS[k]) - nz_dirBS[k])/(sqrtf(two))); count += one;
   //   q = q_dirBN  [k]; if (q>=zero && q<=one) sum += (q * (  ny_dirBN[k]  - nz_dirBN[k])/(sqrtf(two))); count += one;
   //   q = q_dirTS  [k]; if (q>=zero && q<=one) sum += (q * ((-ny_dirTS[k]) + nz_dirTS[k])/(sqrtf(two))); count += one;
   //   q = q_dirTNE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirTNE[k] + ny_dirTNE[k] + nz_dirTNE[k])/(sqrtf(three))); count += one;
   //   q = q_dirTSW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirTSW[k])- ny_dirTSW[k] + nz_dirTSW[k])/(sqrtf(three))); count += one;
   //   q = q_dirTSE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirTSE[k] - ny_dirTSE[k] + nz_dirTSE[k])/(sqrtf(three))); count += one;
   //   q = q_dirTNW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirTNW[k])+ ny_dirTNW[k] + nz_dirTNW[k])/(sqrtf(three))); count += one;
   //   q = q_dirBNE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirBNE[k] + ny_dirBNE[k] - nz_dirBNE[k])/(sqrtf(three))); count += one;
   //   q = q_dirBSW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirBSW[k])- ny_dirBSW[k] - nz_dirBSW[k])/(sqrtf(three))); count += one;
   //   q = q_dirBSE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirBSE[k] - ny_dirBSE[k] - nz_dirBSE[k])/(sqrtf(three))); count += one;
   //   q = q_dirBNW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirBNW[k])+ ny_dirBNW[k] - nz_dirBNW[k])/(sqrtf(three))); count += one;
	  //real qMed = sum/count;
	  //real phi = fac / (qMed + fac);
	  //phi = (phi > one) ? one:one;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real testQ = c2o1;

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirE[k] + vx2 * ny_dirE[k] + vx3 * nz_dirE[k]) * nx_dirE[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( nx_dirE[k]) + fac);
		 VeloX *= phi;
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[W])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirW[k] + vx2 * ny_dirW[k] + vx3 * nz_dirW[k]) * nx_dirW[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-nx_dirW[k]) + fac);
		 VeloX *= phi;
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[E])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloY = vx2 - (vx1 * nx_dirN[k] + vx2 * ny_dirN[k] + vx3 * nz_dirN[k]) * ny_dirN[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( ny_dirN[k]) + fac);
		 VeloY *= phi;
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[S])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloY = vx2 - (vx1 * nx_dirS[k] + vx2 * ny_dirS[k] + vx3 * nz_dirS[k]) * ny_dirS[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-ny_dirS[k]) + fac);
		 VeloY *= phi;
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[N])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloZ = vx3 - (vx1 * nx_dirT[k] + vx2 * ny_dirT[k] + vx3 * nz_dirT[k]) * nz_dirT[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( nz_dirT[k]) + fac);
		 VeloZ *= phi;
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[B])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloZ = vx3 - (vx1 * nx_dirB[k] + vx2 * ny_dirB[k] + vx3 * nz_dirB[k]) * nz_dirB[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-nz_dirB[k]) + fac);
		 VeloZ *= phi;
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[T])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * nx_dirNE[k];
		 VeloY = vx2 - (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * ny_dirNE[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( nx_dirNE[k] + ny_dirNE[k]) + fac);
		 VeloX *= phi;
		 VeloY *= phi;
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * nx_dirSW[k];
		 VeloY = vx2 - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * ny_dirSW[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-nx_dirSW[k] - ny_dirSW[k]) + fac);
		 VeloX *= phi;
		 VeloY *= phi;
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * nx_dirSE[k];
		 VeloY = vx2 - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * ny_dirSE[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( nx_dirSE[k] - ny_dirSE[k]) + fac);
		 VeloX *= phi;
		 VeloY *= phi;
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * nx_dirNW[k];
		 VeloY = vx2 - (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * ny_dirNW[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-nx_dirNW[k] + ny_dirNW[k]) + fac);
		 VeloX *= phi;
		 VeloY *= phi;
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nx_dirTE[k];
		 VeloZ = vx3 - (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nz_dirTE[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( nx_dirTE[k] + nz_dirTE[k]) + fac);
		 VeloX *= phi;
		 VeloZ *= phi;
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nx_dirBW[k];
		 VeloZ = vx3 - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nz_dirBW[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-nx_dirBW[k] - nz_dirBW[k]) + fac);
		 VeloX *= phi;
		 VeloZ *= phi;
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nx_dirBE[k];
		 VeloZ = vx3 - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nz_dirBE[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( nx_dirBE[k] - nz_dirBE[k]) + fac);
		 VeloX *= phi;
		 VeloZ *= phi;
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nx_dirTW[k];
		 VeloZ = vx3 - (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nz_dirTW[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-nx_dirTW[k] + nz_dirTW[k]) + fac);
		 VeloX *= phi;
		 VeloZ *= phi;
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloY = vx2 - (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * ny_dirTN[k];
		 VeloZ = vx3 - (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * nz_dirTN[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( ny_dirTN[k] + nz_dirTN[k]) + fac);
		 VeloY *= phi;
		 VeloZ *= phi;
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloY = vx2 - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * ny_dirBS[k];
		 VeloZ = vx3 - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * nz_dirBS[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-ny_dirBS[k] - nz_dirBS[k]) + fac);
		 VeloY *= phi;
		 VeloZ *= phi;
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloY = vx2 - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * ny_dirBN[k];
		 VeloZ = vx3 - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * nz_dirBN[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( ny_dirBN[k] - nz_dirBN[k]) + fac);
		 VeloY *= phi;
		 VeloZ *= phi;
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloY = vx2 - (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * ny_dirTS[k];
		 VeloZ = vx3 - (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * nz_dirTS[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-ny_dirTS[k] + nz_dirTS[k]) + fac);
		 VeloY *= phi;
		 VeloZ *= phi;
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * nx_dirTNE[k];
		 VeloY = vx2 - (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * ny_dirTNE[k];
		 VeloZ = vx3 - (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * nz_dirTNE[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( nx_dirTNE[k] + ny_dirTNE[k] + nz_dirTNE[k]) + fac);
		 VeloX *= phi;
		 VeloY *= phi;
		 VeloZ *= phi;
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * nx_dirBSW[k];
		 VeloY = vx2 - (vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * ny_dirBSW[k];
		 VeloZ = vx3 - (vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * nz_dirBSW[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-nx_dirBSW[k] - ny_dirBSW[k] - nz_dirBSW[k]) + fac);
		 VeloX *= phi;
		 VeloY *= phi;
		 VeloZ *= phi;
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * nx_dirBNE[k];
		 VeloY = vx2 - (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * ny_dirBNE[k];
		 VeloZ = vx3 - (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * nz_dirBNE[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( nx_dirBNE[k] + ny_dirBNE[k] - nz_dirBNE[k]) + fac);
		 VeloX *= phi;
		 VeloY *= phi;
		 VeloZ *= phi;
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * nx_dirTSW[k];
		 VeloY = vx2 - (vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * ny_dirTSW[k];
		 VeloZ = vx3 - (vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * nz_dirTSW[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-nx_dirTSW[k] - ny_dirTSW[k] + nz_dirTSW[k]) + fac);
		 VeloX *= phi;
		 VeloY *= phi;
		 VeloZ *= phi;
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * nx_dirTSE[k];
		 VeloY = vx2 - (vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * ny_dirTSE[k];
		 VeloZ = vx3 - (vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * nz_dirTSE[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( nx_dirTSE[k] - ny_dirTSE[k] + nz_dirTSE[k]) + fac);
		 VeloX *= phi;
		 VeloY *= phi;
		 VeloZ *= phi;
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * nx_dirBNW[k];
		 VeloY = vx2 - (vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * ny_dirBNW[k];
		 VeloZ = vx3 - (vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * nz_dirBNW[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-nx_dirBNW[k] + ny_dirBNW[k] - nz_dirBNW[k]) + fac);
		 VeloX *= phi;
		 VeloY *= phi;
		 VeloZ *= phi;
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * nx_dirBSE[k];
		 VeloY = vx2 - (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * ny_dirBSE[k];
		 VeloZ = vx3 - (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * nz_dirBSE[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = (phi > one) ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( nx_dirBSE[k] - ny_dirBSE[k] - nz_dirBSE[k]) + fac);
		 VeloX *= phi;
		 VeloY *= phi;
		 VeloZ *= phi;
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * nx_dirTNW[k];
		 VeloY = vx2 - (vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * ny_dirTNW[k];
		 VeloZ = vx3 - (vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * nz_dirTNW[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-nx_dirTNW[k] + ny_dirTNW[k] + nz_dirTNW[k]) + fac);
		 VeloX *= phi;
		 VeloY *= phi;
		 VeloZ *= phi;
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QSlipNormDeviceComp27(real* DD, 
												 int* k_Q, 
												 real* QQ,
												 unsigned int  numberOfBCnodes,
												 real om1, 
												 real* NormalX,
												 real* NormalY,
												 real* NormalZ,
												 unsigned int* neighborX,
												 unsigned int* neighborY,
												 unsigned int* neighborZ,
												 unsigned int size_Mat, 
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

   if(k< numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[E   * numberOfBCnodes];
      q_dirW   = &QQ[W   * numberOfBCnodes];
      q_dirN   = &QQ[N   * numberOfBCnodes];
      q_dirS   = &QQ[S   * numberOfBCnodes];
      q_dirT   = &QQ[T   * numberOfBCnodes];
      q_dirB   = &QQ[B   * numberOfBCnodes];
      q_dirNE  = &QQ[NE  * numberOfBCnodes];
      q_dirSW  = &QQ[SW  * numberOfBCnodes];
      q_dirSE  = &QQ[SE  * numberOfBCnodes];
      q_dirNW  = &QQ[NW  * numberOfBCnodes];
      q_dirTE  = &QQ[TE  * numberOfBCnodes];
      q_dirBW  = &QQ[BW  * numberOfBCnodes];
      q_dirBE  = &QQ[BE  * numberOfBCnodes];
      q_dirTW  = &QQ[TW  * numberOfBCnodes];
      q_dirTN  = &QQ[TN  * numberOfBCnodes];
      q_dirBS  = &QQ[BS  * numberOfBCnodes];
      q_dirBN  = &QQ[BN  * numberOfBCnodes];
      q_dirTS  = &QQ[TS  * numberOfBCnodes];
      q_dirTNE = &QQ[TNE * numberOfBCnodes];
      q_dirTSW = &QQ[TSW * numberOfBCnodes];
      q_dirTSE = &QQ[TSE * numberOfBCnodes];
      q_dirTNW = &QQ[TNW * numberOfBCnodes];
      q_dirBNE = &QQ[BNE * numberOfBCnodes];
      q_dirBSW = &QQ[BSW * numberOfBCnodes];
      q_dirBSE = &QQ[BSE * numberOfBCnodes];
      q_dirBNW = &QQ[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      real *nx_dirE,   *nx_dirW,   *nx_dirN,   *nx_dirS,   *nx_dirT,   *nx_dirB, 
              *nx_dirNE,  *nx_dirSW,  *nx_dirSE,  *nx_dirNW,  *nx_dirTE,  *nx_dirBW,
              *nx_dirBE,  *nx_dirTW,  *nx_dirTN,  *nx_dirBS,  *nx_dirBN,  *nx_dirTS,
              *nx_dirTNE, *nx_dirTSW, *nx_dirTSE, *nx_dirTNW, *nx_dirBNE, *nx_dirBSW,
              *nx_dirBSE, *nx_dirBNW; 
      nx_dirE   = &NormalX[E   * numberOfBCnodes];
      nx_dirW   = &NormalX[W   * numberOfBCnodes];
      nx_dirN   = &NormalX[N   * numberOfBCnodes];
      nx_dirS   = &NormalX[S   * numberOfBCnodes];
      nx_dirT   = &NormalX[T   * numberOfBCnodes];
      nx_dirB   = &NormalX[B   * numberOfBCnodes];
      nx_dirNE  = &NormalX[NE  * numberOfBCnodes];
      nx_dirSW  = &NormalX[SW  * numberOfBCnodes];
      nx_dirSE  = &NormalX[SE  * numberOfBCnodes];
      nx_dirNW  = &NormalX[NW  * numberOfBCnodes];
      nx_dirTE  = &NormalX[TE  * numberOfBCnodes];
      nx_dirBW  = &NormalX[BW  * numberOfBCnodes];
      nx_dirBE  = &NormalX[BE  * numberOfBCnodes];
      nx_dirTW  = &NormalX[TW  * numberOfBCnodes];
      nx_dirTN  = &NormalX[TN  * numberOfBCnodes];
      nx_dirBS  = &NormalX[BS  * numberOfBCnodes];
      nx_dirBN  = &NormalX[BN  * numberOfBCnodes];
      nx_dirTS  = &NormalX[TS  * numberOfBCnodes];
      nx_dirTNE = &NormalX[TNE * numberOfBCnodes];
      nx_dirTSW = &NormalX[TSW * numberOfBCnodes];
      nx_dirTSE = &NormalX[TSE * numberOfBCnodes];
      nx_dirTNW = &NormalX[TNW * numberOfBCnodes];
      nx_dirBNE = &NormalX[BNE * numberOfBCnodes];
      nx_dirBSW = &NormalX[BSW * numberOfBCnodes];
      nx_dirBSE = &NormalX[BSE * numberOfBCnodes];
      nx_dirBNW = &NormalX[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      real *ny_dirE,   *ny_dirW,   *ny_dirN,   *ny_dirS,   *ny_dirT,   *ny_dirB, 
              *ny_dirNE,  *ny_dirSW,  *ny_dirSE,  *ny_dirNW,  *ny_dirTE,  *ny_dirBW,
              *ny_dirBE,  *ny_dirTW,  *ny_dirTN,  *ny_dirBS,  *ny_dirBN,  *ny_dirTS,
              *ny_dirTNE, *ny_dirTSW, *ny_dirTSE, *ny_dirTNW, *ny_dirBNE, *ny_dirBSW,
              *ny_dirBSE, *ny_dirBNW; 
      ny_dirE   = &NormalY[E   * numberOfBCnodes];
      ny_dirW   = &NormalY[W   * numberOfBCnodes];
      ny_dirN   = &NormalY[N   * numberOfBCnodes];
      ny_dirS   = &NormalY[S   * numberOfBCnodes];
      ny_dirT   = &NormalY[T   * numberOfBCnodes];
      ny_dirB   = &NormalY[B   * numberOfBCnodes];
      ny_dirNE  = &NormalY[NE  * numberOfBCnodes];
      ny_dirSW  = &NormalY[SW  * numberOfBCnodes];
      ny_dirSE  = &NormalY[SE  * numberOfBCnodes];
      ny_dirNW  = &NormalY[NW  * numberOfBCnodes];
      ny_dirTE  = &NormalY[TE  * numberOfBCnodes];
      ny_dirBW  = &NormalY[BW  * numberOfBCnodes];
      ny_dirBE  = &NormalY[BE  * numberOfBCnodes];
      ny_dirTW  = &NormalY[TW  * numberOfBCnodes];
      ny_dirTN  = &NormalY[TN  * numberOfBCnodes];
      ny_dirBS  = &NormalY[BS  * numberOfBCnodes];
      ny_dirBN  = &NormalY[BN  * numberOfBCnodes];
      ny_dirTS  = &NormalY[TS  * numberOfBCnodes];
      ny_dirTNE = &NormalY[TNE * numberOfBCnodes];
      ny_dirTSW = &NormalY[TSW * numberOfBCnodes];
      ny_dirTSE = &NormalY[TSE * numberOfBCnodes];
      ny_dirTNW = &NormalY[TNW * numberOfBCnodes];
      ny_dirBNE = &NormalY[BNE * numberOfBCnodes];
      ny_dirBSW = &NormalY[BSW * numberOfBCnodes];
      ny_dirBSE = &NormalY[BSE * numberOfBCnodes];
      ny_dirBNW = &NormalY[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      real *nz_dirE,   *nz_dirW,   *nz_dirN,   *nz_dirS,   *nz_dirT,   *nz_dirB, 
              *nz_dirNE,  *nz_dirSW,  *nz_dirSE,  *nz_dirNW,  *nz_dirTE,  *nz_dirBW,
              *nz_dirBE,  *nz_dirTW,  *nz_dirTN,  *nz_dirBS,  *nz_dirBN,  *nz_dirTS,
              *nz_dirTNE, *nz_dirTSW, *nz_dirTSE, *nz_dirTNW, *nz_dirBNE, *nz_dirBSW,
              *nz_dirBSE, *nz_dirBNW; 
      nz_dirE   = &NormalZ[E   * numberOfBCnodes];
      nz_dirW   = &NormalZ[W   * numberOfBCnodes];
      nz_dirN   = &NormalZ[N   * numberOfBCnodes];
      nz_dirS   = &NormalZ[S   * numberOfBCnodes];
      nz_dirT   = &NormalZ[T   * numberOfBCnodes];
      nz_dirB   = &NormalZ[B   * numberOfBCnodes];
      nz_dirNE  = &NormalZ[NE  * numberOfBCnodes];
      nz_dirSW  = &NormalZ[SW  * numberOfBCnodes];
      nz_dirSE  = &NormalZ[SE  * numberOfBCnodes];
      nz_dirNW  = &NormalZ[NW  * numberOfBCnodes];
      nz_dirTE  = &NormalZ[TE  * numberOfBCnodes];
      nz_dirBW  = &NormalZ[BW  * numberOfBCnodes];
      nz_dirBE  = &NormalZ[BE  * numberOfBCnodes];
      nz_dirTW  = &NormalZ[TW  * numberOfBCnodes];
      nz_dirTN  = &NormalZ[TN  * numberOfBCnodes];
      nz_dirBS  = &NormalZ[BS  * numberOfBCnodes];
      nz_dirBN  = &NormalZ[BN  * numberOfBCnodes];
      nz_dirTS  = &NormalZ[TS  * numberOfBCnodes];
      nz_dirTNE = &NormalZ[TNE * numberOfBCnodes];
      nz_dirTSW = &NormalZ[TSW * numberOfBCnodes];
      nz_dirTSE = &NormalZ[TSE * numberOfBCnodes];
      nz_dirTNW = &NormalZ[TNW * numberOfBCnodes];
      nz_dirBNE = &NormalZ[BNE * numberOfBCnodes];
      nz_dirBSW = &NormalZ[BSW * numberOfBCnodes];
      nz_dirBSE = &NormalZ[BSE * numberOfBCnodes];
      nz_dirBNW = &NormalZ[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real f_W    = (D.f[E   ])[ke   ];
      real f_E    = (D.f[W   ])[kw   ];
      real f_S    = (D.f[N   ])[kn   ];
      real f_N    = (D.f[S   ])[ks   ];
      real f_B    = (D.f[T   ])[kt   ];
      real f_T    = (D.f[B   ])[kb   ];
      real f_SW   = (D.f[NE  ])[kne  ];
      real f_NE   = (D.f[SW  ])[ksw  ];
      real f_NW   = (D.f[SE  ])[kse  ];
      real f_SE   = (D.f[NW  ])[knw  ];
      real f_BW   = (D.f[TE  ])[kte  ];
      real f_TE   = (D.f[BW  ])[kbw  ];
      real f_TW   = (D.f[BE  ])[kbe  ];
      real f_BE   = (D.f[TW  ])[ktw  ];
      real f_BS   = (D.f[TN  ])[ktn  ];
      real f_TN   = (D.f[BS  ])[kbs  ];
      real f_TS   = (D.f[BN  ])[kbn  ];
      real f_BN   = (D.f[TS  ])[kts  ];
      real f_BSW  = (D.f[TNE ])[ktne ];
      real f_BNE  = (D.f[TSW ])[ktsw ];
      real f_BNW  = (D.f[TSE ])[ktse ];
      real f_BSE  = (D.f[TNW ])[ktnw ];
      real f_TSW  = (D.f[BNE ])[kbne ];
      real f_TNE  = (D.f[BSW ])[kbsw ];
      real f_TNW  = (D.f[BSE ])[kbse ];
      real f_TSE  = (D.f[BNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[REST])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (c1o1 + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (c1o1 + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (c1o1 + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real VeloX = vx1;
	  real VeloY = vx2;
	  real VeloZ = vx3;
	  real fac = c1o100;//0.5;
 	  //real phi = c0o1;
	  real alpha = c1o100;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real kxyFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_SW+f_BSW+f_TSW-f_NW-f_BNW-f_TNW-f_SE-f_BSE-f_TSE+f_NE+f_BNE+f_TNE ) / (c1o1 + drho) - ((vx1*vx2)));
      real kyzFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_BS+f_BSE+f_BSW-f_TS-f_TSE-f_TSW-f_BN-f_BNE-f_BNW+f_TN+f_TNE+f_TNW ) / (c1o1 + drho) - ((vx2*vx3)));
      real kxzFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_BW+f_BSW+f_BNW-f_TW-f_TSW-f_TNW-f_BE-f_BSE-f_BNE+f_TE+f_TSE+f_TNE ) / (c1o1 + drho) - ((vx1*vx3)));

	  real kxxFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_E+f_NE+f_SE+f_TE+f_BE+f_W+f_NW+f_SW+f_TW+f_BW+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (c1o1 + drho) - ((c1o3*drho + vx1*vx1)));
	  real kyyFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_N+f_NE+f_NW+f_TN+f_BN+f_S+f_SE+f_SW+f_TS+f_BS+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (c1o1 + drho) - ((c1o3*drho + vx2*vx2)));
	  real kzzFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_T+f_TE+f_TW+f_TN+f_BS+f_B+f_BE+f_BW+f_BN+f_BS+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (c1o1 + drho) - ((c1o3*drho + vx3*vx3)));

	  real magS = sqrtf(kxyFromfcNEQ*kxyFromfcNEQ + kyzFromfcNEQ*kyzFromfcNEQ + kxzFromfcNEQ*kxzFromfcNEQ + kxxFromfcNEQ*kxxFromfcNEQ + kyyFromfcNEQ*kyyFromfcNEQ + kzzFromfcNEQ*kzzFromfcNEQ);

	  fac = fac * magS / (c1o3 * (c1o1 / om1 - c1o2));
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real *facAst = &QQ[REST * numberOfBCnodes];

	  fac = fac * alpha + facAst[k] * (c1o1 - alpha);
	  facAst[k] = fac;
	  //(&QQ[REST * numberOfBCnodes])[KQK] = fac;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ////real uk = sqrtf(vx1*vx1 + vx2*vx2 + vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //real phi = expf(magS/0.01f) - one;
	  //phi = (phi > one) ? one:phi;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //real C = five;
	  //real kappa = 0.41f;
	  //real phi = (C * kappa * c1o2 * logf(magS / (c1o3 * (one / om1 - c1o2))) - one) / (C * kappa * c1o2 * logf(magS / (c1o3 * (one / om1 - c1o2))));
	  //phi = (phi < zero) ? zero:phi;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //real sum = zero, count = zero;
   //   q = q_dirE   [k]; if (q>=zero && q<=one) sum += (q *   nx_dirE[k] ); count += one;
   //   q = q_dirW   [k]; if (q>=zero && q<=one) sum += (q * (-nx_dirW[k])); count += one;
   //   q = q_dirN   [k]; if (q>=zero && q<=one) sum += (q *   nx_dirN[k] ); count += one;
   //   q = q_dirS   [k]; if (q>=zero && q<=one) sum += (q * (-nx_dirS[k])); count += one;
   //   q = q_dirT   [k]; if (q>=zero && q<=one) sum += (q *   nx_dirT[k] ); count += one;
   //   q = q_dirB   [k]; if (q>=zero && q<=one) sum += (q * (-nx_dirB[k])); count += one;
   //   q = q_dirNE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirNE[k]  + ny_dirNE[k])/(sqrtf(two))); count += one;
   //   q = q_dirSW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirSW[k]) - ny_dirSW[k])/(sqrtf(two))); count += one;
   //   q = q_dirSE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirSE[k]  - ny_dirSE[k])/(sqrtf(two))); count += one;
   //   q = q_dirNW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirNW[k]) + ny_dirNW[k])/(sqrtf(two))); count += one;
   //   q = q_dirTE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirTE[k]  + nz_dirTE[k])/(sqrtf(two))); count += one;
   //   q = q_dirBW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirBW[k]) - nz_dirBW[k])/(sqrtf(two))); count += one;
   //   q = q_dirBE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirBE[k]  - nz_dirBE[k])/(sqrtf(two))); count += one;
   //   q = q_dirTW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirTW[k]) + nz_dirTW[k])/(sqrtf(two))); count += one;
   //   q = q_dirTN  [k]; if (q>=zero && q<=one) sum += (q * (  ny_dirTN[k]  + nz_dirTN[k])/(sqrtf(two))); count += one;
   //   q = q_dirBS  [k]; if (q>=zero && q<=one) sum += (q * ((-ny_dirBS[k]) - nz_dirBS[k])/(sqrtf(two))); count += one;
   //   q = q_dirBN  [k]; if (q>=zero && q<=one) sum += (q * (  ny_dirBN[k]  - nz_dirBN[k])/(sqrtf(two))); count += one;
   //   q = q_dirTS  [k]; if (q>=zero && q<=one) sum += (q * ((-ny_dirTS[k]) + nz_dirTS[k])/(sqrtf(two))); count += one;
   //   q = q_dirTNE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirTNE[k] + ny_dirTNE[k] + nz_dirTNE[k])/(sqrtf(three))); count += one;
   //   q = q_dirTSW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirTSW[k])- ny_dirTSW[k] + nz_dirTSW[k])/(sqrtf(three))); count += one;
   //   q = q_dirTSE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirTSE[k] - ny_dirTSE[k] + nz_dirTSE[k])/(sqrtf(three))); count += one;
   //   q = q_dirTNW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirTNW[k])+ ny_dirTNW[k] + nz_dirTNW[k])/(sqrtf(three))); count += one;
   //   q = q_dirBNE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirBNE[k] + ny_dirBNE[k] - nz_dirBNE[k])/(sqrtf(three))); count += one;
   //   q = q_dirBSW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirBSW[k])- ny_dirBSW[k] - nz_dirBSW[k])/(sqrtf(three))); count += one;
   //   q = q_dirBSE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirBSE[k] - ny_dirBSE[k] - nz_dirBSE[k])/(sqrtf(three))); count += one;
   //   q = q_dirBNW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirBNW[k])+ ny_dirBNW[k] - nz_dirBNW[k])/(sqrtf(three))); count += one;
	  //real qMed = sum/count;
	  //real phi = fac / (qMed + fac);
	  //phi = (phi > one) ? one:one;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real sliplength = 0.9f;//c1o2;
	  real qSlip = c0o1;
	  real un = c0o1;
	  real ut = c0o1;
	  real tangential = c0o1;
	  //real smallSingle = Op0000002;

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirE[k] + vx2 * ny_dirE[k] + vx3 * nz_dirE[k]) * nx_dirE[k];
		 un = fabs((vx1 * nx_dirE[k] + vx2 * ny_dirE[k] + vx3 * nz_dirE[k]) * nx_dirE[k]);
		 ut = fabs(VeloX);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[W])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirW[k] + vx2 * ny_dirW[k] + vx3 * nz_dirW[k]) * nx_dirW[k];
		 un = fabs(-(vx1 * nx_dirW[k] + vx2 * ny_dirW[k] + vx3 * nz_dirW[k]) * nx_dirW[k]);
		 ut = fabs(-VeloX);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[E])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloY = vx2 - (vx1 * nx_dirN[k] + vx2 * ny_dirN[k] + vx3 * nz_dirN[k]) * ny_dirN[k];
		 un = fabs( (vx1 * nx_dirN[k] + vx2 * ny_dirN[k] + vx3 * nz_dirN[k]) * ny_dirN[k]);
		 ut = fabs( VeloY);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( ny_dirN[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[S])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloY = vx2 - (vx1 * nx_dirS[k] + vx2 * ny_dirS[k] + vx3 * nz_dirS[k]) * ny_dirS[k];
		 un = fabs(-(vx1 * nx_dirS[k] + vx2 * ny_dirS[k] + vx3 * nz_dirS[k]) * ny_dirS[k]);
		 ut = fabs(-VeloY);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-ny_dirS[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[N])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloZ = vx3 - (vx1 * nx_dirT[k] + vx2 * ny_dirT[k] + vx3 * nz_dirT[k]) * nz_dirT[k];
		 un = fabs( (vx1 * nx_dirT[k] + vx2 * ny_dirT[k] + vx3 * nz_dirT[k]) * nz_dirT[k]);
		 ut = fabs( VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nz_dirT[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[B])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloZ = vx3 - (vx1 * nx_dirB[k] + vx2 * ny_dirB[k] + vx3 * nz_dirB[k]) * nz_dirB[k];
		 un = fabs(-(vx1 * nx_dirB[k] + vx2 * ny_dirB[k] + vx3 * nz_dirB[k]) * nz_dirB[k]);
		 ut = fabs(-VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nz_dirB[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[T])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * nx_dirNE[k];
		 VeloY = vx2 - (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * ny_dirNE[k];
		 un = fabs( (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * nx_dirNE[k] + (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * ny_dirNE[k]);
		 ut = fabs( VeloX + VeloY);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirNE[k]+ny_dirNE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * nx_dirSW[k];
		 VeloY = vx2 - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * ny_dirSW[k];
		 un = fabs(-(vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * nx_dirSW[k] - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * ny_dirSW[k]);
		 ut = fabs(-VeloX - VeloY);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirSW[k]-ny_dirSW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * nx_dirSE[k];
		 VeloY = vx2 - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * ny_dirSE[k];
		 un = fabs( (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * nx_dirSE[k] - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * ny_dirSE[k]);
		 ut = fabs( VeloX - VeloY);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirSE[k]-ny_dirSE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * nx_dirNW[k];
		 VeloY = vx2 - (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * ny_dirNW[k];
		 un = fabs(-(vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * nx_dirNW[k] + (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * ny_dirNW[k]);
		 ut = fabs(-VeloX + VeloY);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirNW[k]+ny_dirNW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nx_dirTE[k];
		 VeloZ = vx3 - (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nz_dirTE[k];
		 un = fabs( (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nx_dirTE[k] + (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nz_dirTE[k]);
		 ut = fabs( VeloX + VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirTE[k]+nz_dirTE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nx_dirBW[k];
		 VeloZ = vx3 - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nz_dirBW[k];
		 un = fabs(-(vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nx_dirBW[k] - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nz_dirBW[k]);
		 ut = fabs(-VeloX - VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirBW[k]-nz_dirBW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nx_dirBE[k];
		 VeloZ = vx3 - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nz_dirBE[k];
		 un = fabs( (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nx_dirBE[k] - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nz_dirBE[k]);
		 ut = fabs( VeloX - VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirBE[k]-nz_dirBE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nx_dirTW[k];
		 VeloZ = vx3 - (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nz_dirTW[k];
		 un = fabs(-(vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nx_dirTW[k] + (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nz_dirTW[k]);
		 ut = fabs(-VeloX + VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirTW[k]+nz_dirTW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloY = vx2 - (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * ny_dirTN[k];
		 VeloZ = vx3 - (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * nz_dirTN[k];
		 un = fabs( (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * ny_dirTN[k] + (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * nz_dirTN[k]);
		 ut = fabs( VeloY + VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( ny_dirTN[k]+nz_dirTN[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloY = vx2 - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * ny_dirBS[k];
		 VeloZ = vx3 - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * nz_dirBS[k];
		 un = fabs(-(vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * ny_dirBS[k] - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * nz_dirBS[k]);
		 ut = fabs(-VeloY - VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-ny_dirBS[k]-nz_dirBS[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloY = vx2 - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * ny_dirBN[k];
		 VeloZ = vx3 - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * nz_dirBN[k];
		 un = fabs( (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * ny_dirBN[k] - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * nz_dirBN[k]);
		 ut = fabs( VeloY - VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( ny_dirBN[k]-nz_dirBN[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloY = vx2 - (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * ny_dirTS[k];
		 VeloZ = vx3 - (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * nz_dirTS[k];
		 un = fabs(-(vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * ny_dirTS[k] + (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * nz_dirTS[k]);
		 ut = fabs(-VeloY + VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-ny_dirTS[k]+nz_dirTS[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * nx_dirTNE[k];
		 VeloY = vx2 - (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * ny_dirTNE[k];
		 VeloZ = vx3 - (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * nz_dirTNE[k];
		 un = fabs( (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * nx_dirTNE[k] 
				   +(vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * ny_dirTNE[k] 
				   +(vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * nz_dirTNE[k]);
		 ut = fabs( VeloX + VeloY + VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirTNE[k] + ny_dirTNE[k] + nz_dirTNE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * nx_dirBSW[k];
		 VeloY = vx2 - (vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * ny_dirBSW[k];
		 VeloZ = vx3 - (vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * nz_dirBSW[k];
		 un = fabs(-(vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * nx_dirBSW[k] 
				   -(vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * ny_dirBSW[k] 
				   -(vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * nz_dirBSW[k]);
		 ut = fabs(-VeloX - VeloY - VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirBSW[k] - ny_dirBSW[k] - nz_dirBSW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * nx_dirBNE[k];
		 VeloY = vx2 - (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * ny_dirBNE[k];
		 VeloZ = vx3 - (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * nz_dirBNE[k];
		 un = fabs( (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * nx_dirBNE[k] 
				   +(vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * ny_dirBNE[k] 
				   -(vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * nz_dirBNE[k]);
		 ut = fabs( VeloX + VeloY - VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirBNE[k] + ny_dirBNE[k] - nz_dirBNE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * nx_dirTSW[k];
		 VeloY = vx2 - (vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * ny_dirTSW[k];
		 VeloZ = vx3 - (vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * nz_dirTSW[k];
		 un = fabs(-(vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * nx_dirTSW[k] 
				   -(vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * ny_dirTSW[k] 
				   +(vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * nz_dirTSW[k]);
		 ut = fabs(-VeloX - VeloY + VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirTSW[k] - ny_dirTSW[k] + nz_dirTSW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * nx_dirTSE[k];
		 VeloY = vx2 - (vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * ny_dirTSE[k];
		 VeloZ = vx3 - (vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * nz_dirTSE[k];
		 un = fabs(+(vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * nx_dirTSE[k] 
				   -(vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * ny_dirTSE[k] 
				   +(vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * nz_dirTSE[k]);
		 ut = fabs(+VeloX - VeloY + VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirTSE[k] - ny_dirTSE[k] + nz_dirTSE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * nx_dirBNW[k];
		 VeloY = vx2 - (vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * ny_dirBNW[k];
		 VeloZ = vx3 - (vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * nz_dirBNW[k];
		 un = fabs(-(vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * nx_dirBNW[k] 
				   +(vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * ny_dirBNW[k] 
				   -(vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * nz_dirBNW[k]);
		 ut = fabs(-VeloX + VeloY - VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirBNW[k] + ny_dirBNW[k] - nz_dirBNW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * nx_dirBSE[k];
		 VeloY = vx2 - (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * ny_dirBSE[k];
		 VeloZ = vx3 - (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * nz_dirBSE[k];
		 un = fabs( (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * nx_dirBSE[k] 
				   -(vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * ny_dirBSE[k] 
				   -(vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * nz_dirBSE[k]);
		 ut = fabs( VeloX - VeloY - VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirBSE[k] - ny_dirBSE[k] - nz_dirBSE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = vx1 - (vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * nx_dirTNW[k];
		 VeloY = vx2 - (vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * ny_dirTNW[k];
		 VeloZ = vx3 - (vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * nz_dirTNW[k];
		 un = fabs(-(vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * nx_dirTNW[k] 
				   +(vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * ny_dirTNW[k] 
				   +(vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * nz_dirTNW[k]);
		 ut = fabs(-VeloX + VeloY + VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirTNW[k] + ny_dirTNW[k] + nz_dirTNW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (smallSingle + q));
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE))/(c1o1+q) - c1o216 * drho;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
