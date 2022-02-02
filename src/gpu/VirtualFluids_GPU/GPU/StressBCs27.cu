/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;

//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QStressDeviceComp27(real* DD, 
											 int* k_Q, 
											 real* QQ,
											 unsigned int sizeQ,
											 real om1, 
											 unsigned int* neighborX,
											 unsigned int* neighborY,
											 unsigned int* neighborZ,
											 unsigned int size_Mat, 
											 bool evenOrOdd)
{
   Distributions27 D;
   if (evenOrOdd==true)
   {
      D.f[dirE   ] = &DD[dirE   *size_Mat];
      D.f[dirW   ] = &DD[dirW   *size_Mat];
      D.f[dirN   ] = &DD[dirN   *size_Mat];
      D.f[dirS   ] = &DD[dirS   *size_Mat];
      D.f[dirT   ] = &DD[dirT   *size_Mat];
      D.f[dirB   ] = &DD[dirB   *size_Mat];
      D.f[dirNE  ] = &DD[dirNE  *size_Mat];
      D.f[dirSW  ] = &DD[dirSW  *size_Mat];
      D.f[dirSE  ] = &DD[dirSE  *size_Mat];
      D.f[dirNW  ] = &DD[dirNW  *size_Mat];
      D.f[dirTE  ] = &DD[dirTE  *size_Mat];
      D.f[dirBW  ] = &DD[dirBW  *size_Mat];
      D.f[dirBE  ] = &DD[dirBE  *size_Mat];
      D.f[dirTW  ] = &DD[dirTW  *size_Mat];
      D.f[dirTN  ] = &DD[dirTN  *size_Mat];
      D.f[dirBS  ] = &DD[dirBS  *size_Mat];
      D.f[dirBN  ] = &DD[dirBN  *size_Mat];
      D.f[dirTS  ] = &DD[dirTS  *size_Mat];
      D.f[dirZERO] = &DD[dirZERO*size_Mat];
      D.f[dirTNE ] = &DD[dirTNE *size_Mat];
      D.f[dirTSW ] = &DD[dirTSW *size_Mat];
      D.f[dirTSE ] = &DD[dirTSE *size_Mat];
      D.f[dirTNW ] = &DD[dirTNW *size_Mat];
      D.f[dirBNE ] = &DD[dirBNE *size_Mat];
      D.f[dirBSW ] = &DD[dirBSW *size_Mat];
      D.f[dirBSE ] = &DD[dirBSE *size_Mat];
      D.f[dirBNW ] = &DD[dirBNW *size_Mat];
   } 
   else
   {
      D.f[dirW   ] = &DD[dirE   *size_Mat];
      D.f[dirE   ] = &DD[dirW   *size_Mat];
      D.f[dirS   ] = &DD[dirN   *size_Mat];
      D.f[dirN   ] = &DD[dirS   *size_Mat];
      D.f[dirB   ] = &DD[dirT   *size_Mat];
      D.f[dirT   ] = &DD[dirB   *size_Mat];
      D.f[dirSW  ] = &DD[dirNE  *size_Mat];
      D.f[dirNE  ] = &DD[dirSW  *size_Mat];
      D.f[dirNW  ] = &DD[dirSE  *size_Mat];
      D.f[dirSE  ] = &DD[dirNW  *size_Mat];
      D.f[dirBW  ] = &DD[dirTE  *size_Mat];
      D.f[dirTE  ] = &DD[dirBW  *size_Mat];
      D.f[dirTW  ] = &DD[dirBE  *size_Mat];
      D.f[dirBE  ] = &DD[dirTW  *size_Mat];
      D.f[dirBS  ] = &DD[dirTN  *size_Mat];
      D.f[dirTN  ] = &DD[dirBS  *size_Mat];
      D.f[dirTS  ] = &DD[dirBN  *size_Mat];
      D.f[dirBN  ] = &DD[dirTS  *size_Mat];
      D.f[dirZERO] = &DD[dirZERO*size_Mat];
      D.f[dirTNE ] = &DD[dirBSW *size_Mat];
      D.f[dirTSW ] = &DD[dirBNE *size_Mat];
      D.f[dirTSE ] = &DD[dirBNW *size_Mat];
      D.f[dirTNW ] = &DD[dirBSE *size_Mat];
      D.f[dirBNE ] = &DD[dirTSW *size_Mat];
      D.f[dirBSW ] = &DD[dirTNE *size_Mat];
      D.f[dirBSE ] = &DD[dirTNW *size_Mat];
      D.f[dirBNW ] = &DD[dirTSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<sizeQ)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
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
      real f_W    = (D.f[dirE   ])[ke   ];
      real f_E    = (D.f[dirW   ])[kw   ];
      real f_S    = (D.f[dirN   ])[kn   ];
      real f_N    = (D.f[dirS   ])[ks   ];
      real f_B    = (D.f[dirT   ])[kt   ];
      real f_T    = (D.f[dirB   ])[kb   ];
      real f_SW   = (D.f[dirNE  ])[kne  ];
      real f_NE   = (D.f[dirSW  ])[ksw  ];
      real f_NW   = (D.f[dirSE  ])[kse  ];
      real f_SE   = (D.f[dirNW  ])[knw  ];
      real f_BW   = (D.f[dirTE  ])[kte  ];
      real f_TE   = (D.f[dirBW  ])[kbw  ];
      real f_TW   = (D.f[dirBE  ])[kbe  ];
      real f_BE   = (D.f[dirTW  ])[ktw  ];
      real f_BS   = (D.f[dirTN  ])[ktn  ];
      real f_TN   = (D.f[dirBS  ])[kbs  ];
      real f_TS   = (D.f[dirBN  ])[kbn  ];
      real f_BN   = (D.f[dirTS  ])[kts  ];
      real f_BSW  = (D.f[dirTNE ])[ktne ];
      real f_BNE  = (D.f[dirTSW ])[ktsw ];
      real f_BNW  = (D.f[dirTSE ])[ktse ];
      real f_BSE  = (D.f[dirTNW ])[ktnw ];
      real f_TSW  = (D.f[dirBNE ])[kbne ];
      real f_TNE  = (D.f[dirBSW ])[kbsw ];
      real f_TNW  = (D.f[dirBSE ])[kbse ];
      real f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

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
      if (evenOrOdd==false)
      {
         D.f[dirE   ] = &DD[dirE   *size_Mat];
         D.f[dirW   ] = &DD[dirW   *size_Mat];
         D.f[dirN   ] = &DD[dirN   *size_Mat];
         D.f[dirS   ] = &DD[dirS   *size_Mat];
         D.f[dirT   ] = &DD[dirT   *size_Mat];
         D.f[dirB   ] = &DD[dirB   *size_Mat];
         D.f[dirNE  ] = &DD[dirNE  *size_Mat];
         D.f[dirSW  ] = &DD[dirSW  *size_Mat];
         D.f[dirSE  ] = &DD[dirSE  *size_Mat];
         D.f[dirNW  ] = &DD[dirNW  *size_Mat];
         D.f[dirTE  ] = &DD[dirTE  *size_Mat];
         D.f[dirBW  ] = &DD[dirBW  *size_Mat];
         D.f[dirBE  ] = &DD[dirBE  *size_Mat];
         D.f[dirTW  ] = &DD[dirTW  *size_Mat];
         D.f[dirTN  ] = &DD[dirTN  *size_Mat];
         D.f[dirBS  ] = &DD[dirBS  *size_Mat];
         D.f[dirBN  ] = &DD[dirBN  *size_Mat];
         D.f[dirTS  ] = &DD[dirTS  *size_Mat];
         D.f[dirZERO] = &DD[dirZERO*size_Mat];
         D.f[dirTNE ] = &DD[dirTNE *size_Mat];
         D.f[dirTSW ] = &DD[dirTSW *size_Mat];
         D.f[dirTSE ] = &DD[dirTSE *size_Mat];
         D.f[dirTNW ] = &DD[dirTNW *size_Mat];
         D.f[dirBNE ] = &DD[dirBNE *size_Mat];
         D.f[dirBSW ] = &DD[dirBSW *size_Mat];
         D.f[dirBSE ] = &DD[dirBSE *size_Mat];
         D.f[dirBNW ] = &DD[dirBNW *size_Mat];
      } 
      else
      {
         D.f[dirW   ] = &DD[dirE   *size_Mat];
         D.f[dirE   ] = &DD[dirW   *size_Mat];
         D.f[dirS   ] = &DD[dirN   *size_Mat];
         D.f[dirN   ] = &DD[dirS   *size_Mat];
         D.f[dirB   ] = &DD[dirT   *size_Mat];
         D.f[dirT   ] = &DD[dirB   *size_Mat];
         D.f[dirSW  ] = &DD[dirNE  *size_Mat];
         D.f[dirNE  ] = &DD[dirSW  *size_Mat];
         D.f[dirNW  ] = &DD[dirSE  *size_Mat];
         D.f[dirSE  ] = &DD[dirNW  *size_Mat];
         D.f[dirBW  ] = &DD[dirTE  *size_Mat];
         D.f[dirTE  ] = &DD[dirBW  *size_Mat];
         D.f[dirTW  ] = &DD[dirBE  *size_Mat];
         D.f[dirBE  ] = &DD[dirTW  *size_Mat];
         D.f[dirBS  ] = &DD[dirTN  *size_Mat];
         D.f[dirTN  ] = &DD[dirBS  *size_Mat];
         D.f[dirTS  ] = &DD[dirBN  *size_Mat];
         D.f[dirBN  ] = &DD[dirTS  *size_Mat];
         D.f[dirZERO] = &DD[dirZERO*size_Mat];
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real fac = c1o1;//c99o100;
	  real VeloX = fac*vx1;
	  real VeloY = fac*vx2;
	  real VeloZ = fac*vx3;
	  bool x = false;
	  bool y = false;
	  bool z = false;
     
   //   printf("qs: \t %f \t %f \n", q_dirE[k],q_dirT[k] );

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = c0o1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 x = true;
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirW])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q) - c2o27 * drho;
         //feq=c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq); 
         //(D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-feq*om1)/(one-om1)+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = c0o1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 x = true;
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirE])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q) - c2o27 * drho;
         //feq=c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
         //(D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-feq*om1)/(one-om1)+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
		 VeloY = c0o1;
	     VeloZ = fac*vx3;
		 y = true;
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirS])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q) - c2o27 * drho;
         //feq=c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
         //(D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-feq*om1)/(one-om1)+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
		 VeloY = c0o1;
	     VeloZ = fac*vx3;
		 y = true;
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirN])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q) - c2o27 * drho;
         //feq=c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
         //(D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-feq*om1)/(one-om1)+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
		 VeloZ = c0o1;
		 z = true;
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirB])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q) - c2o27 * drho;
         //feq=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq); 
         //(D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-feq*om1)/(one-om1)+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
		 VeloZ = c0o1;
		 z = true;
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirT])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q) - c2o27 * drho;
         //feq=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
         //(D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-feq*om1)/(one-om1)+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);
         //(D.f[dirT])[kt]=zero;
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
         (D.f[dirSW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         //(D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-feq*om1)/(one-om1)+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);
         //(D.f[dirSW])[ksw]=zero;
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
         (D.f[dirNE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         //(D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-feq*om1)/(one-om1)+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);
         //(D.f[dirNE])[kne]=zero;
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
         (D.f[dirNW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         //(D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-feq*om1)/(one-om1)+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);
         //(D.f[dirNW])[knw]=zero;
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
         (D.f[dirSE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         //(D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-feq*om1)/(one-om1)+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = c0o1;
		 if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         //(D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-feq*om1)/(one-om1)+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);
         //(D.f[dirBW])[kbw]=zero;
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
         (D.f[dirTE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         //(D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-feq*om1)/(one-om1)+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);
         //(D.f[dirTE])[kte]=zero;
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
         (D.f[dirTW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         //(D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-feq*om1)/(one-om1)+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);
         //(D.f[dirTW])[ktw]=zero;
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
         (D.f[dirBE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         //(D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-feq*om1)/(one-om1)+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);
         //(D.f[dirBE])[kbe]=zero;
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
         (D.f[dirBS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         //(D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-feq*om1)/(one-om1)+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);
         //(D.f[dirBS])[kbs]=zero;
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
         (D.f[dirTN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         //(D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-feq*om1)/(one-om1)+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);
         //(D.f[dirTN])[ktn]=zero;
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
         (D.f[dirTS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         //(D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-feq*om1)/(one-om1)+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);
         //(D.f[dirTS])[kts]=zero;
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
         (D.f[dirBN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
         //feq=c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         //(D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-feq*om1)/(one-om1)+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);
         //(D.f[dirBN])[kbn]=zero;
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
         (D.f[dirBSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         //(D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-feq*om1)/(one-om1)+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);
         //(D.f[dirBSW])[kbsw]=zero;
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
         (D.f[dirTNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         //(D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-feq*om1)/(one-om1)+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);
         //(D.f[dirTNE])[ktne]=zero;
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
         (D.f[dirTSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         //(D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-feq*om1)/(one-om1)+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);
         //(D.f[dirTSW])[ktsw]=zero;
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
         (D.f[dirBNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         //(D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-feq*om1)/(one-om1)+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);
         //(D.f[dirBNE])[kbne]=zero;
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
         (D.f[dirBNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         //(D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-feq*om1)/(one-om1)+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);
         //(D.f[dirBNW])[kbnw]=zero;
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
         (D.f[dirTSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         //(D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-feq*om1)/(one-om1)+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);
         //(D.f[dirTSE])[ktse]=zero;
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
         (D.f[dirTNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         //(D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-feq*om1)/(one-om1)+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);
         //(D.f[dirTNW])[ktnw]=zero;
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
         (D.f[dirBSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
         //feq=c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         //(D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-feq*om1)/(one-om1)+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}