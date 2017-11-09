/* Device code */
#include "LBM/D3Q27.h"
#include "GPU/constant.h"

//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QSlipDevice27(doubflo* DD, 
                                         int* k_Q, 
                                         doubflo* QQ,
                                         unsigned int sizeQ,
                                         doubflo om1, 
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
      doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
      doubflo f_W    = (D.f[dirE   ])[ke   ];
      doubflo f_E    = (D.f[dirW   ])[kw   ];
      doubflo f_S    = (D.f[dirN   ])[kn   ];
      doubflo f_N    = (D.f[dirS   ])[ks   ];
      doubflo f_B    = (D.f[dirT   ])[kt   ];
      doubflo f_T    = (D.f[dirB   ])[kb   ];
      doubflo f_SW   = (D.f[dirNE  ])[kne  ];
      doubflo f_NE   = (D.f[dirSW  ])[ksw  ];
      doubflo f_NW   = (D.f[dirSE  ])[kse  ];
      doubflo f_SE   = (D.f[dirNW  ])[knw  ];
      doubflo f_BW   = (D.f[dirTE  ])[kte  ];
      doubflo f_TE   = (D.f[dirBW  ])[kbw  ];
      doubflo f_TW   = (D.f[dirBE  ])[kbe  ];
      doubflo f_BE   = (D.f[dirTW  ])[ktw  ];
      doubflo f_BS   = (D.f[dirTN  ])[ktn  ];
      doubflo f_TN   = (D.f[dirBS  ])[kbs  ];
      doubflo f_TS   = (D.f[dirBN  ])[kbn  ];
      doubflo f_BN   = (D.f[dirTS  ])[kts  ];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W); 
         

      vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S); 

      vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

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
	  doubflo fac = one;//c99o100;
	  doubflo VeloX = fac*vx1;
	  doubflo VeloY = fac*vx2;
	  doubflo VeloZ = fac*vx3;
	  bool x = false;
	  bool y = false;
	  bool z = false;

      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = zero;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 x = true;
         feq=c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq); 
         (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-feq*om1)/(one-om1)+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = zero;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 x = true;
         feq=c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
         (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-feq*om1)/(one-om1)+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
		 VeloY = zero;
	     VeloZ = fac*vx3;
		 y = true;
         feq=c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
         (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-feq*om1)/(one-om1)+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
		 VeloY = zero;
	     VeloZ = fac*vx3;
		 y = true;
         feq=c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
         (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-feq*om1)/(one-om1)+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
		 VeloZ = zero;
		 z = true;
         feq=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq); 
         (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-feq*om1)/(one-om1)+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
		 VeloZ = zero;
		 z = true;
         feq=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
         (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-feq*om1)/(one-om1)+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
         feq=c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-feq*om1)/(one-om1)+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
         feq=c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-feq*om1)/(one-om1)+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
         feq=c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-feq*om1)/(one-om1)+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
         feq=c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-feq*om1)/(one-om1)+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-feq*om1)/(one-om1)+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-feq*om1)/(one-om1)+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-feq*om1)/(one-om1)+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-feq*om1)/(one-om1)+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-feq*om1)/(one-om1)+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-feq*om1)/(one-om1)+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-feq*om1)/(one-om1)+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-feq*om1)/(one-om1)+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-feq*om1)/(one-om1)+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-feq*om1)/(one-om1)+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-feq*om1)/(one-om1)+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-feq*om1)/(one-om1)+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-feq*om1)/(one-om1)+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-feq*om1)/(one-om1)+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-feq*om1)/(one-om1)+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-feq*om1)/(one-om1)+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QSlipDeviceComp27(doubflo* DD, 
											 int* k_Q, 
											 doubflo* QQ,
											 unsigned int sizeQ,
											 doubflo om1, 
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
      doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
      doubflo f_W    = (D.f[dirE   ])[ke   ];
      doubflo f_E    = (D.f[dirW   ])[kw   ];
      doubflo f_S    = (D.f[dirN   ])[kn   ];
      doubflo f_N    = (D.f[dirS   ])[ks   ];
      doubflo f_B    = (D.f[dirT   ])[kt   ];
      doubflo f_T    = (D.f[dirB   ])[kb   ];
      doubflo f_SW   = (D.f[dirNE  ])[kne  ];
      doubflo f_NE   = (D.f[dirSW  ])[ksw  ];
      doubflo f_NW   = (D.f[dirSE  ])[kse  ];
      doubflo f_SE   = (D.f[dirNW  ])[knw  ];
      doubflo f_BW   = (D.f[dirTE  ])[kte  ];
      doubflo f_TE   = (D.f[dirBW  ])[kbw  ];
      doubflo f_TW   = (D.f[dirBE  ])[kbe  ];
      doubflo f_BE   = (D.f[dirTW  ])[ktw  ];
      doubflo f_BS   = (D.f[dirTN  ])[ktn  ];
      doubflo f_TN   = (D.f[dirBS  ])[kbs  ];
      doubflo f_TS   = (D.f[dirBN  ])[kbn  ];
      doubflo f_BN   = (D.f[dirTS  ])[kts  ];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (one + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (one + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (one + drho); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (one + drho);

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
	  doubflo fac = one;//c99o100;
	  doubflo VeloX = fac*vx1;
	  doubflo VeloY = fac*vx2;
	  doubflo VeloZ = fac*vx3;
	  bool x = false;
	  bool y = false;
	  bool z = false;

      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = zero;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 x = true;
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q) - c2over27 * drho;
         //feq=c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq); 
         //(D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-feq*om1)/(one-om1)+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = zero;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 x = true;
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q) - c2over27 * drho;
         //feq=c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
         //(D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-feq*om1)/(one-om1)+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
		 VeloY = zero;
	     VeloZ = fac*vx3;
		 y = true;
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q) - c2over27 * drho;
         //feq=c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
         //(D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-feq*om1)/(one-om1)+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
		 VeloY = zero;
	     VeloZ = fac*vx3;
		 y = true;
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q) - c2over27 * drho;
         //feq=c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
         //(D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-feq*om1)/(one-om1)+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
		 VeloZ = zero;
		 z = true;
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
         (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q) - c2over27 * drho;
         //feq=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq); 
         //(D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-feq*om1)/(one-om1)+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
		 VeloZ = zero;
		 z = true;
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
         (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q) - c2over27 * drho;
         //feq=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
         //(D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-feq*om1)/(one-om1)+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q) - c1over54 * drho;
         //feq=c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         //(D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-feq*om1)/(one-om1)+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q) - c1over54 * drho;
         //feq=c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         //(D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-feq*om1)/(one-om1)+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q) - c1over54 * drho;
         //feq=c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         //(D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-feq*om1)/(one-om1)+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q) - c1over54 * drho;
         //feq=c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         //(D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-feq*om1)/(one-om1)+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //feq=c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         //(D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-feq*om1)/(one-om1)+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //feq=c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         //(D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-feq*om1)/(one-om1)+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //feq=c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         //(D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-feq*om1)/(one-om1)+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //feq=c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         //(D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-feq*om1)/(one-om1)+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //feq=c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         //(D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-feq*om1)/(one-om1)+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //feq=c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         //(D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-feq*om1)/(one-om1)+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //feq=c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         //(D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-feq*om1)/(one-om1)+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //feq=c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         //(D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-feq*om1)/(one-om1)+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //feq=c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         //(D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-feq*om1)/(one-om1)+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //feq=c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         //(D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-feq*om1)/(one-om1)+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //feq=c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         //(D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-feq*om1)/(one-om1)+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //feq=c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         //(D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-feq*om1)/(one-om1)+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //feq=c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         //(D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-feq*om1)/(one-om1)+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //feq=c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         //(D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-feq*om1)/(one-om1)+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //feq=c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         //(D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-feq*om1)/(one-om1)+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = fac*vx1;
	     VeloY = fac*vx2;
	     VeloZ = fac*vx3;
		 if (x == true) VeloX = zero;
		 if (y == true) VeloY = zero;
		 if (z == true) VeloZ = zero;
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //feq=c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         //(D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-feq*om1)/(one-om1)+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QSlipGeomDeviceComp27(doubflo* DD, 
												 int* k_Q, 
												 doubflo* QQ,
												 unsigned int sizeQ,
												 doubflo om1, 
												 doubflo* NormalX,
												 doubflo* NormalY,
												 doubflo* NormalZ,
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
      doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
      doubflo *nx_dirE,   *nx_dirW,   *nx_dirN,   *nx_dirS,   *nx_dirT,   *nx_dirB, 
              *nx_dirNE,  *nx_dirSW,  *nx_dirSE,  *nx_dirNW,  *nx_dirTE,  *nx_dirBW,
              *nx_dirBE,  *nx_dirTW,  *nx_dirTN,  *nx_dirBS,  *nx_dirBN,  *nx_dirTS,
              *nx_dirTNE, *nx_dirTSW, *nx_dirTSE, *nx_dirTNW, *nx_dirBNE, *nx_dirBSW,
              *nx_dirBSE, *nx_dirBNW; 
      nx_dirE   = &NormalX[dirE   *sizeQ];
      nx_dirW   = &NormalX[dirW   *sizeQ];
      nx_dirN   = &NormalX[dirN   *sizeQ];
      nx_dirS   = &NormalX[dirS   *sizeQ];
      nx_dirT   = &NormalX[dirT   *sizeQ];
      nx_dirB   = &NormalX[dirB   *sizeQ];
      nx_dirNE  = &NormalX[dirNE  *sizeQ];
      nx_dirSW  = &NormalX[dirSW  *sizeQ];
      nx_dirSE  = &NormalX[dirSE  *sizeQ];
      nx_dirNW  = &NormalX[dirNW  *sizeQ];
      nx_dirTE  = &NormalX[dirTE  *sizeQ];
      nx_dirBW  = &NormalX[dirBW  *sizeQ];
      nx_dirBE  = &NormalX[dirBE  *sizeQ];
      nx_dirTW  = &NormalX[dirTW  *sizeQ];
      nx_dirTN  = &NormalX[dirTN  *sizeQ];
      nx_dirBS  = &NormalX[dirBS  *sizeQ];
      nx_dirBN  = &NormalX[dirBN  *sizeQ];
      nx_dirTS  = &NormalX[dirTS  *sizeQ];
      nx_dirTNE = &NormalX[dirTNE *sizeQ];
      nx_dirTSW = &NormalX[dirTSW *sizeQ];
      nx_dirTSE = &NormalX[dirTSE *sizeQ];
      nx_dirTNW = &NormalX[dirTNW *sizeQ];
      nx_dirBNE = &NormalX[dirBNE *sizeQ];
      nx_dirBSW = &NormalX[dirBSW *sizeQ];
      nx_dirBSE = &NormalX[dirBSE *sizeQ];
      nx_dirBNW = &NormalX[dirBNW *sizeQ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo *ny_dirE,   *ny_dirW,   *ny_dirN,   *ny_dirS,   *ny_dirT,   *ny_dirB, 
              *ny_dirNE,  *ny_dirSW,  *ny_dirSE,  *ny_dirNW,  *ny_dirTE,  *ny_dirBW,
              *ny_dirBE,  *ny_dirTW,  *ny_dirTN,  *ny_dirBS,  *ny_dirBN,  *ny_dirTS,
              *ny_dirTNE, *ny_dirTSW, *ny_dirTSE, *ny_dirTNW, *ny_dirBNE, *ny_dirBSW,
              *ny_dirBSE, *ny_dirBNW; 
      ny_dirE   = &NormalY[dirE   *sizeQ];
      ny_dirW   = &NormalY[dirW   *sizeQ];
      ny_dirN   = &NormalY[dirN   *sizeQ];
      ny_dirS   = &NormalY[dirS   *sizeQ];
      ny_dirT   = &NormalY[dirT   *sizeQ];
      ny_dirB   = &NormalY[dirB   *sizeQ];
      ny_dirNE  = &NormalY[dirNE  *sizeQ];
      ny_dirSW  = &NormalY[dirSW  *sizeQ];
      ny_dirSE  = &NormalY[dirSE  *sizeQ];
      ny_dirNW  = &NormalY[dirNW  *sizeQ];
      ny_dirTE  = &NormalY[dirTE  *sizeQ];
      ny_dirBW  = &NormalY[dirBW  *sizeQ];
      ny_dirBE  = &NormalY[dirBE  *sizeQ];
      ny_dirTW  = &NormalY[dirTW  *sizeQ];
      ny_dirTN  = &NormalY[dirTN  *sizeQ];
      ny_dirBS  = &NormalY[dirBS  *sizeQ];
      ny_dirBN  = &NormalY[dirBN  *sizeQ];
      ny_dirTS  = &NormalY[dirTS  *sizeQ];
      ny_dirTNE = &NormalY[dirTNE *sizeQ];
      ny_dirTSW = &NormalY[dirTSW *sizeQ];
      ny_dirTSE = &NormalY[dirTSE *sizeQ];
      ny_dirTNW = &NormalY[dirTNW *sizeQ];
      ny_dirBNE = &NormalY[dirBNE *sizeQ];
      ny_dirBSW = &NormalY[dirBSW *sizeQ];
      ny_dirBSE = &NormalY[dirBSE *sizeQ];
      ny_dirBNW = &NormalY[dirBNW *sizeQ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo *nz_dirE,   *nz_dirW,   *nz_dirN,   *nz_dirS,   *nz_dirT,   *nz_dirB, 
              *nz_dirNE,  *nz_dirSW,  *nz_dirSE,  *nz_dirNW,  *nz_dirTE,  *nz_dirBW,
              *nz_dirBE,  *nz_dirTW,  *nz_dirTN,  *nz_dirBS,  *nz_dirBN,  *nz_dirTS,
              *nz_dirTNE, *nz_dirTSW, *nz_dirTSE, *nz_dirTNW, *nz_dirBNE, *nz_dirBSW,
              *nz_dirBSE, *nz_dirBNW; 
      nz_dirE   = &NormalZ[dirE   *sizeQ];
      nz_dirW   = &NormalZ[dirW   *sizeQ];
      nz_dirN   = &NormalZ[dirN   *sizeQ];
      nz_dirS   = &NormalZ[dirS   *sizeQ];
      nz_dirT   = &NormalZ[dirT   *sizeQ];
      nz_dirB   = &NormalZ[dirB   *sizeQ];
      nz_dirNE  = &NormalZ[dirNE  *sizeQ];
      nz_dirSW  = &NormalZ[dirSW  *sizeQ];
      nz_dirSE  = &NormalZ[dirSE  *sizeQ];
      nz_dirNW  = &NormalZ[dirNW  *sizeQ];
      nz_dirTE  = &NormalZ[dirTE  *sizeQ];
      nz_dirBW  = &NormalZ[dirBW  *sizeQ];
      nz_dirBE  = &NormalZ[dirBE  *sizeQ];
      nz_dirTW  = &NormalZ[dirTW  *sizeQ];
      nz_dirTN  = &NormalZ[dirTN  *sizeQ];
      nz_dirBS  = &NormalZ[dirBS  *sizeQ];
      nz_dirBN  = &NormalZ[dirBN  *sizeQ];
      nz_dirTS  = &NormalZ[dirTS  *sizeQ];
      nz_dirTNE = &NormalZ[dirTNE *sizeQ];
      nz_dirTSW = &NormalZ[dirTSW *sizeQ];
      nz_dirTSE = &NormalZ[dirTSE *sizeQ];
      nz_dirTNW = &NormalZ[dirTNW *sizeQ];
      nz_dirBNE = &NormalZ[dirBNE *sizeQ];
      nz_dirBSW = &NormalZ[dirBSW *sizeQ];
      nz_dirBSE = &NormalZ[dirBSE *sizeQ];
      nz_dirBNW = &NormalZ[dirBNW *sizeQ];
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
      doubflo f_W    = (D.f[dirE   ])[ke   ];
      doubflo f_E    = (D.f[dirW   ])[kw   ];
      doubflo f_S    = (D.f[dirN   ])[kn   ];
      doubflo f_N    = (D.f[dirS   ])[ks   ];
      doubflo f_B    = (D.f[dirT   ])[kt   ];
      doubflo f_T    = (D.f[dirB   ])[kb   ];
      doubflo f_SW   = (D.f[dirNE  ])[kne  ];
      doubflo f_NE   = (D.f[dirSW  ])[ksw  ];
      doubflo f_NW   = (D.f[dirSE  ])[kse  ];
      doubflo f_SE   = (D.f[dirNW  ])[knw  ];
      doubflo f_BW   = (D.f[dirTE  ])[kte  ];
      doubflo f_TE   = (D.f[dirBW  ])[kbw  ];
      doubflo f_TW   = (D.f[dirBE  ])[kbe  ];
      doubflo f_BE   = (D.f[dirTW  ])[ktw  ];
      doubflo f_BS   = (D.f[dirTN  ])[ktn  ];
      doubflo f_TN   = (D.f[dirBS  ])[kbs  ];
      doubflo f_TS   = (D.f[dirBN  ])[kbn  ];
      doubflo f_BN   = (D.f[dirTS  ])[kts  ];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (one + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (one + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (one + drho); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (one + drho);

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
	  doubflo VeloX = vx1;
	  doubflo VeloY = vx2;
	  doubflo VeloZ = vx3;
	  doubflo fac = zero;//0.5;
 	  doubflo phi = zero;
	  doubflo alpha = c1o100;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo kxyFromfcNEQ = -(three * om1 / (one-om1))*((f_SW+f_BSW+f_TSW-f_NW-f_BNW-f_TNW-f_SE-f_BSE-f_TSE+f_NE+f_BNE+f_TNE ) / (one + drho) - ((vx1*vx2)));
      doubflo kyzFromfcNEQ = -(three * om1 / (one-om1))*((f_BS+f_BSE+f_BSW-f_TS-f_TSE-f_TSW-f_BN-f_BNE-f_BNW+f_TN+f_TNE+f_TNW ) / (one + drho) - ((vx2*vx3)));
      doubflo kxzFromfcNEQ = -(three * om1 / (one-om1))*((f_BW+f_BSW+f_BNW-f_TW-f_TSW-f_TNW-f_BE-f_BSE-f_BNE+f_TE+f_TSE+f_TNE ) / (one + drho) - ((vx1*vx3)));

	  doubflo kxxFromfcNEQ = -(three * om1 / (one-om1))*((f_E+f_NE+f_SE+f_TE+f_BE+f_W+f_NW+f_SW+f_TW+f_BW+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (one + drho) - ((c1o3*drho + vx1*vx1)));
	  doubflo kyyFromfcNEQ = -(three * om1 / (one-om1))*((f_N+f_NE+f_NW+f_TN+f_BN+f_S+f_SE+f_SW+f_TS+f_BS+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (one + drho) - ((c1o3*drho + vx2*vx2)));
	  doubflo kzzFromfcNEQ = -(three * om1 / (one-om1))*((f_T+f_TE+f_TW+f_TN+f_BS+f_B+f_BE+f_BW+f_BN+f_BS+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (one + drho) - ((c1o3*drho + vx3*vx3)));

	  doubflo magS = sqrtf(kxyFromfcNEQ*kxyFromfcNEQ + kyzFromfcNEQ*kyzFromfcNEQ + kxzFromfcNEQ*kxzFromfcNEQ + kxxFromfcNEQ*kxxFromfcNEQ + kyyFromfcNEQ*kyyFromfcNEQ + kzzFromfcNEQ*kzzFromfcNEQ);

	  //fac = fac * magS / (c1o3 * (one / om1 - c1o2));
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //doubflo *facAst = &QQ[dirZERO *sizeQ];

	  //fac = fac * alpha + facAst[k] * (one - alpha);
	  //facAst[k] = fac;
	  //(&QQ[dirZERO *sizeQ])[KQK] = fac;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ////doubflo uk = sqrtf(vx1*vx1 + vx2*vx2 + vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //doubflo phi = expf(magS/0.01f) - one;
	  //phi = (phi > one) ? one:phi;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //doubflo C = five;
	  //doubflo kappa = 0.41f;
	  //doubflo phi = (C * kappa * c1o2 * logf(magS / (c1o3 * (one / om1 - c1o2))) - one) / (C * kappa * c1o2 * logf(magS / (c1o3 * (one / om1 - c1o2))));
	  //phi = (phi < zero) ? zero:phi;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //doubflo sum = zero, count = zero;
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
	  //doubflo qMed = sum/count;
	  //doubflo phi = fac / (qMed + fac);
	  //phi = (phi > one) ? one:one;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  doubflo testQ = two;

      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = vx1 - (vx1 * nx_dirE[k] + vx2 * ny_dirE[k] + vx3 * nz_dirE[k]) * nx_dirE[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( nx_dirE[k]) + fac);
		 VeloX *= phi;
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q) - c2over27 * drho;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = vx1 - (vx1 * nx_dirW[k] + vx2 * ny_dirW[k] + vx3 * nz_dirW[k]) * nx_dirW[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-nx_dirW[k]) + fac);
		 VeloX *= phi;
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q) - c2over27 * drho;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
		 VeloY = vx2 - (vx1 * nx_dirN[k] + vx2 * ny_dirN[k] + vx3 * nz_dirN[k]) * ny_dirN[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( ny_dirN[k]) + fac);
		 VeloY *= phi;
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q) - c2over27 * drho;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
		 VeloY = vx2 - (vx1 * nx_dirS[k] + vx2 * ny_dirS[k] + vx3 * nz_dirS[k]) * ny_dirS[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-ny_dirS[k]) + fac);
		 VeloY *= phi;
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q) - c2over27 * drho;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
		 VeloZ = vx3 - (vx1 * nx_dirT[k] + vx2 * ny_dirT[k] + vx3 * nz_dirT[k]) * nz_dirT[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs( nz_dirT[k]) + fac);
		 VeloZ *= phi;
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
         (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q) - c2over27 * drho;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
		 VeloZ = vx3 - (vx1 * nx_dirB[k] + vx2 * ny_dirB[k] + vx3 * nz_dirB[k]) * nz_dirB[k];
		 //phi = fac * (one + magS / (Op0000002+uk) * (one-q));
		 //phi = phi > one ? one:phi;
		 //phi = fac; //Test
		 q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
		 phi = fac / (q * fabs(-nz_dirB[k]) + fac);
		 VeloZ *= phi;
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
         (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q) - c2over27 * drho;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
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
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q) - c1over54 * drho;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
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
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q) - c1over54 * drho;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
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
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q) - c1over54 * drho;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
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
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q) - c1over54 * drho;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
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
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q) - c1over54 * drho;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
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
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q) - c1over54 * drho;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
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
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q) - c1over54 * drho;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
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
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q) - c1over54 * drho;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
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
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q) - c1over54 * drho;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
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
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q) - c1over54 * drho;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
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
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q) - c1over54 * drho;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
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
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q) - c1over54 * drho;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
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
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
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
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
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
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
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
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
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
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
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
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
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
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
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
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QSlipNormDeviceComp27(doubflo* DD, 
												 int* k_Q, 
												 doubflo* QQ,
												 unsigned int sizeQ,
												 doubflo om1, 
												 doubflo* NormalX,
												 doubflo* NormalY,
												 doubflo* NormalZ,
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
      doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
      doubflo *nx_dirE,   *nx_dirW,   *nx_dirN,   *nx_dirS,   *nx_dirT,   *nx_dirB, 
              *nx_dirNE,  *nx_dirSW,  *nx_dirSE,  *nx_dirNW,  *nx_dirTE,  *nx_dirBW,
              *nx_dirBE,  *nx_dirTW,  *nx_dirTN,  *nx_dirBS,  *nx_dirBN,  *nx_dirTS,
              *nx_dirTNE, *nx_dirTSW, *nx_dirTSE, *nx_dirTNW, *nx_dirBNE, *nx_dirBSW,
              *nx_dirBSE, *nx_dirBNW; 
      nx_dirE   = &NormalX[dirE   *sizeQ];
      nx_dirW   = &NormalX[dirW   *sizeQ];
      nx_dirN   = &NormalX[dirN   *sizeQ];
      nx_dirS   = &NormalX[dirS   *sizeQ];
      nx_dirT   = &NormalX[dirT   *sizeQ];
      nx_dirB   = &NormalX[dirB   *sizeQ];
      nx_dirNE  = &NormalX[dirNE  *sizeQ];
      nx_dirSW  = &NormalX[dirSW  *sizeQ];
      nx_dirSE  = &NormalX[dirSE  *sizeQ];
      nx_dirNW  = &NormalX[dirNW  *sizeQ];
      nx_dirTE  = &NormalX[dirTE  *sizeQ];
      nx_dirBW  = &NormalX[dirBW  *sizeQ];
      nx_dirBE  = &NormalX[dirBE  *sizeQ];
      nx_dirTW  = &NormalX[dirTW  *sizeQ];
      nx_dirTN  = &NormalX[dirTN  *sizeQ];
      nx_dirBS  = &NormalX[dirBS  *sizeQ];
      nx_dirBN  = &NormalX[dirBN  *sizeQ];
      nx_dirTS  = &NormalX[dirTS  *sizeQ];
      nx_dirTNE = &NormalX[dirTNE *sizeQ];
      nx_dirTSW = &NormalX[dirTSW *sizeQ];
      nx_dirTSE = &NormalX[dirTSE *sizeQ];
      nx_dirTNW = &NormalX[dirTNW *sizeQ];
      nx_dirBNE = &NormalX[dirBNE *sizeQ];
      nx_dirBSW = &NormalX[dirBSW *sizeQ];
      nx_dirBSE = &NormalX[dirBSE *sizeQ];
      nx_dirBNW = &NormalX[dirBNW *sizeQ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo *ny_dirE,   *ny_dirW,   *ny_dirN,   *ny_dirS,   *ny_dirT,   *ny_dirB, 
              *ny_dirNE,  *ny_dirSW,  *ny_dirSE,  *ny_dirNW,  *ny_dirTE,  *ny_dirBW,
              *ny_dirBE,  *ny_dirTW,  *ny_dirTN,  *ny_dirBS,  *ny_dirBN,  *ny_dirTS,
              *ny_dirTNE, *ny_dirTSW, *ny_dirTSE, *ny_dirTNW, *ny_dirBNE, *ny_dirBSW,
              *ny_dirBSE, *ny_dirBNW; 
      ny_dirE   = &NormalY[dirE   *sizeQ];
      ny_dirW   = &NormalY[dirW   *sizeQ];
      ny_dirN   = &NormalY[dirN   *sizeQ];
      ny_dirS   = &NormalY[dirS   *sizeQ];
      ny_dirT   = &NormalY[dirT   *sizeQ];
      ny_dirB   = &NormalY[dirB   *sizeQ];
      ny_dirNE  = &NormalY[dirNE  *sizeQ];
      ny_dirSW  = &NormalY[dirSW  *sizeQ];
      ny_dirSE  = &NormalY[dirSE  *sizeQ];
      ny_dirNW  = &NormalY[dirNW  *sizeQ];
      ny_dirTE  = &NormalY[dirTE  *sizeQ];
      ny_dirBW  = &NormalY[dirBW  *sizeQ];
      ny_dirBE  = &NormalY[dirBE  *sizeQ];
      ny_dirTW  = &NormalY[dirTW  *sizeQ];
      ny_dirTN  = &NormalY[dirTN  *sizeQ];
      ny_dirBS  = &NormalY[dirBS  *sizeQ];
      ny_dirBN  = &NormalY[dirBN  *sizeQ];
      ny_dirTS  = &NormalY[dirTS  *sizeQ];
      ny_dirTNE = &NormalY[dirTNE *sizeQ];
      ny_dirTSW = &NormalY[dirTSW *sizeQ];
      ny_dirTSE = &NormalY[dirTSE *sizeQ];
      ny_dirTNW = &NormalY[dirTNW *sizeQ];
      ny_dirBNE = &NormalY[dirBNE *sizeQ];
      ny_dirBSW = &NormalY[dirBSW *sizeQ];
      ny_dirBSE = &NormalY[dirBSE *sizeQ];
      ny_dirBNW = &NormalY[dirBNW *sizeQ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo *nz_dirE,   *nz_dirW,   *nz_dirN,   *nz_dirS,   *nz_dirT,   *nz_dirB, 
              *nz_dirNE,  *nz_dirSW,  *nz_dirSE,  *nz_dirNW,  *nz_dirTE,  *nz_dirBW,
              *nz_dirBE,  *nz_dirTW,  *nz_dirTN,  *nz_dirBS,  *nz_dirBN,  *nz_dirTS,
              *nz_dirTNE, *nz_dirTSW, *nz_dirTSE, *nz_dirTNW, *nz_dirBNE, *nz_dirBSW,
              *nz_dirBSE, *nz_dirBNW; 
      nz_dirE   = &NormalZ[dirE   *sizeQ];
      nz_dirW   = &NormalZ[dirW   *sizeQ];
      nz_dirN   = &NormalZ[dirN   *sizeQ];
      nz_dirS   = &NormalZ[dirS   *sizeQ];
      nz_dirT   = &NormalZ[dirT   *sizeQ];
      nz_dirB   = &NormalZ[dirB   *sizeQ];
      nz_dirNE  = &NormalZ[dirNE  *sizeQ];
      nz_dirSW  = &NormalZ[dirSW  *sizeQ];
      nz_dirSE  = &NormalZ[dirSE  *sizeQ];
      nz_dirNW  = &NormalZ[dirNW  *sizeQ];
      nz_dirTE  = &NormalZ[dirTE  *sizeQ];
      nz_dirBW  = &NormalZ[dirBW  *sizeQ];
      nz_dirBE  = &NormalZ[dirBE  *sizeQ];
      nz_dirTW  = &NormalZ[dirTW  *sizeQ];
      nz_dirTN  = &NormalZ[dirTN  *sizeQ];
      nz_dirBS  = &NormalZ[dirBS  *sizeQ];
      nz_dirBN  = &NormalZ[dirBN  *sizeQ];
      nz_dirTS  = &NormalZ[dirTS  *sizeQ];
      nz_dirTNE = &NormalZ[dirTNE *sizeQ];
      nz_dirTSW = &NormalZ[dirTSW *sizeQ];
      nz_dirTSE = &NormalZ[dirTSE *sizeQ];
      nz_dirTNW = &NormalZ[dirTNW *sizeQ];
      nz_dirBNE = &NormalZ[dirBNE *sizeQ];
      nz_dirBSW = &NormalZ[dirBSW *sizeQ];
      nz_dirBSE = &NormalZ[dirBSE *sizeQ];
      nz_dirBNW = &NormalZ[dirBNW *sizeQ];
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
      doubflo f_W    = (D.f[dirE   ])[ke   ];
      doubflo f_E    = (D.f[dirW   ])[kw   ];
      doubflo f_S    = (D.f[dirN   ])[kn   ];
      doubflo f_N    = (D.f[dirS   ])[ks   ];
      doubflo f_B    = (D.f[dirT   ])[kt   ];
      doubflo f_T    = (D.f[dirB   ])[kb   ];
      doubflo f_SW   = (D.f[dirNE  ])[kne  ];
      doubflo f_NE   = (D.f[dirSW  ])[ksw  ];
      doubflo f_NW   = (D.f[dirSE  ])[kse  ];
      doubflo f_SE   = (D.f[dirNW  ])[knw  ];
      doubflo f_BW   = (D.f[dirTE  ])[kte  ];
      doubflo f_TE   = (D.f[dirBW  ])[kbw  ];
      doubflo f_TW   = (D.f[dirBE  ])[kbe  ];
      doubflo f_BE   = (D.f[dirTW  ])[ktw  ];
      doubflo f_BS   = (D.f[dirTN  ])[ktn  ];
      doubflo f_TN   = (D.f[dirBS  ])[kbs  ];
      doubflo f_TS   = (D.f[dirBN  ])[kbn  ];
      doubflo f_BN   = (D.f[dirTS  ])[kts  ];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (one + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (one + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (one + drho); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (one + drho);

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
	  doubflo VeloX = vx1;
	  doubflo VeloY = vx2;
	  doubflo VeloZ = vx3;
	  doubflo fac = c1o100;//0.5;
 	  doubflo phi = zero;
	  doubflo alpha = c1o100;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo kxyFromfcNEQ = -(three * om1 / (one-om1))*((f_SW+f_BSW+f_TSW-f_NW-f_BNW-f_TNW-f_SE-f_BSE-f_TSE+f_NE+f_BNE+f_TNE ) / (one + drho) - ((vx1*vx2)));
      doubflo kyzFromfcNEQ = -(three * om1 / (one-om1))*((f_BS+f_BSE+f_BSW-f_TS-f_TSE-f_TSW-f_BN-f_BNE-f_BNW+f_TN+f_TNE+f_TNW ) / (one + drho) - ((vx2*vx3)));
      doubflo kxzFromfcNEQ = -(three * om1 / (one-om1))*((f_BW+f_BSW+f_BNW-f_TW-f_TSW-f_TNW-f_BE-f_BSE-f_BNE+f_TE+f_TSE+f_TNE ) / (one + drho) - ((vx1*vx3)));

	  doubflo kxxFromfcNEQ = -(three * om1 / (one-om1))*((f_E+f_NE+f_SE+f_TE+f_BE+f_W+f_NW+f_SW+f_TW+f_BW+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (one + drho) - ((c1o3*drho + vx1*vx1)));
	  doubflo kyyFromfcNEQ = -(three * om1 / (one-om1))*((f_N+f_NE+f_NW+f_TN+f_BN+f_S+f_SE+f_SW+f_TS+f_BS+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (one + drho) - ((c1o3*drho + vx2*vx2)));
	  doubflo kzzFromfcNEQ = -(three * om1 / (one-om1))*((f_T+f_TE+f_TW+f_TN+f_BS+f_B+f_BE+f_BW+f_BN+f_BS+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (one + drho) - ((c1o3*drho + vx3*vx3)));

	  doubflo magS = sqrtf(kxyFromfcNEQ*kxyFromfcNEQ + kyzFromfcNEQ*kyzFromfcNEQ + kxzFromfcNEQ*kxzFromfcNEQ + kxxFromfcNEQ*kxxFromfcNEQ + kyyFromfcNEQ*kyyFromfcNEQ + kzzFromfcNEQ*kzzFromfcNEQ);

	  fac = fac * magS / (c1o3 * (one / om1 - c1o2));
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  doubflo *facAst = &QQ[dirZERO *sizeQ];

	  fac = fac * alpha + facAst[k] * (one - alpha);
	  facAst[k] = fac;
	  //(&QQ[dirZERO *sizeQ])[KQK] = fac;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ////doubflo uk = sqrtf(vx1*vx1 + vx2*vx2 + vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //doubflo phi = expf(magS/0.01f) - one;
	  //phi = (phi > one) ? one:phi;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //doubflo C = five;
	  //doubflo kappa = 0.41f;
	  //doubflo phi = (C * kappa * c1o2 * logf(magS / (c1o3 * (one / om1 - c1o2))) - one) / (C * kappa * c1o2 * logf(magS / (c1o3 * (one / om1 - c1o2))));
	  //phi = (phi < zero) ? zero:phi;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //doubflo sum = zero, count = zero;
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
	  //doubflo qMed = sum/count;
	  //doubflo phi = fac / (qMed + fac);
	  //phi = (phi > one) ? one:one;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  doubflo sliplength = 0.9f;//c1o2;
	  doubflo qSlip = zero;
	  doubflo un = zero;
	  doubflo ut = zero;
	  doubflo tangential = zero;
	  //doubflo smallSingle = Op0000002;

      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = vx1 - (vx1 * nx_dirE[k] + vx2 * ny_dirE[k] + vx3 * nz_dirE[k]) * nx_dirE[k];
		 un = fabs((vx1 * nx_dirE[k] + vx2 * ny_dirE[k] + vx3 * nz_dirE[k]) * nx_dirE[k]);
		 ut = fabs(VeloX);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W))/(one+q) - c2over27 * drho;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = vx1 - (vx1 * nx_dirW[k] + vx2 * ny_dirW[k] + vx3 * nz_dirW[k]) * nx_dirW[k];
		 un = fabs(-(vx1 * nx_dirW[k] + vx2 * ny_dirW[k] + vx3 * nz_dirW[k]) * nx_dirW[k]);
		 ut = fabs(-VeloX);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E))/(one+q) - c2over27 * drho;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
		 VeloY = vx2 - (vx1 * nx_dirN[k] + vx2 * ny_dirN[k] + vx3 * nz_dirN[k]) * ny_dirN[k];
		 un = fabs( (vx1 * nx_dirN[k] + vx2 * ny_dirN[k] + vx3 * nz_dirN[k]) * ny_dirN[k]);
		 ut = fabs( VeloY);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( ny_dirN[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S))/(one+q) - c2over27 * drho;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
		 VeloY = vx2 - (vx1 * nx_dirS[k] + vx2 * ny_dirS[k] + vx3 * nz_dirS[k]) * ny_dirS[k];
		 un = fabs(-(vx1 * nx_dirS[k] + vx2 * ny_dirS[k] + vx3 * nz_dirS[k]) * ny_dirS[k]);
		 ut = fabs(-VeloY);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-ny_dirS[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N))/(one+q) - c2over27 * drho;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
		 VeloZ = vx3 - (vx1 * nx_dirT[k] + vx2 * ny_dirT[k] + vx3 * nz_dirT[k]) * nz_dirT[k];
		 un = fabs( (vx1 * nx_dirT[k] + vx2 * ny_dirT[k] + vx3 * nz_dirT[k]) * nz_dirT[k]);
		 ut = fabs( VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nz_dirT[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
         (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B))/(one+q) - c2over27 * drho;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
		 VeloZ = vx3 - (vx1 * nx_dirB[k] + vx2 * ny_dirB[k] + vx3 * nz_dirB[k]) * nz_dirB[k];
		 un = fabs(-(vx1 * nx_dirB[k] + vx2 * ny_dirB[k] + vx3 * nz_dirB[k]) * nz_dirB[k]);
		 ut = fabs(-VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nz_dirB[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
         (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T))/(one+q) - c2over27 * drho;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = vx1 - (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * nx_dirNE[k];
		 VeloY = vx2 - (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * ny_dirNE[k];
		 un = fabs( (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * nx_dirNE[k] + (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * ny_dirNE[k]);
		 ut = fabs( VeloX + VeloY);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirNE[k]+ny_dirNE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW))/(one+q) - c1over54 * drho;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = vx1 - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * nx_dirSW[k];
		 VeloY = vx2 - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * ny_dirSW[k];
		 un = fabs(-(vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * nx_dirSW[k] - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * ny_dirSW[k]);
		 ut = fabs(-VeloX - VeloY);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirSW[k]-ny_dirSW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE))/(one+q) - c1over54 * drho;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = vx1 - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * nx_dirSE[k];
		 VeloY = vx2 - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * ny_dirSE[k];
		 un = fabs( (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * nx_dirSE[k] - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * ny_dirSE[k]);
		 ut = fabs( VeloX - VeloY);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirSE[k]-ny_dirSE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW))/(one+q) - c1over54 * drho;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = vx1 - (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * nx_dirNW[k];
		 VeloY = vx2 - (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * ny_dirNW[k];
		 un = fabs(-(vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * nx_dirNW[k] + (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * ny_dirNW[k]);
		 ut = fabs(-VeloX + VeloY);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirNW[k]+ny_dirNW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE))/(one+q) - c1over54 * drho;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = vx1 - (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nx_dirTE[k];
		 VeloZ = vx3 - (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nz_dirTE[k];
		 un = fabs( (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nx_dirTE[k] + (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nz_dirTE[k]);
		 ut = fabs( VeloX + VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirTE[k]+nz_dirTE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW))/(one+q) - c1over54 * drho;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = vx1 - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nx_dirBW[k];
		 VeloZ = vx3 - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nz_dirBW[k];
		 un = fabs(-(vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nx_dirBW[k] - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nz_dirBW[k]);
		 ut = fabs(-VeloX - VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirBW[k]-nz_dirBW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE))/(one+q) - c1over54 * drho;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
		 VeloX = vx1 - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nx_dirBE[k];
		 VeloZ = vx3 - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nz_dirBE[k];
		 un = fabs( (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nx_dirBE[k] - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nz_dirBE[k]);
		 ut = fabs( VeloX - VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( nx_dirBE[k]-nz_dirBE[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW))/(one+q) - c1over54 * drho;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
		 VeloX = vx1 - (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nx_dirTW[k];
		 VeloZ = vx3 - (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nz_dirTW[k];
		 un = fabs(-(vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nx_dirTW[k] + (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nz_dirTW[k]);
		 ut = fabs(-VeloX + VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-nx_dirTW[k]+nz_dirTW[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE))/(one+q) - c1over54 * drho;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
		 VeloY = vx2 - (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * ny_dirTN[k];
		 VeloZ = vx3 - (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * nz_dirTN[k];
		 un = fabs( (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * ny_dirTN[k] + (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * nz_dirTN[k]);
		 ut = fabs( VeloY + VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( ny_dirTN[k]+nz_dirTN[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS))/(one+q) - c1over54 * drho;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
		 VeloY = vx2 - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * ny_dirBS[k];
		 VeloZ = vx3 - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * nz_dirBS[k];
		 un = fabs(-(vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * ny_dirBS[k] - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * nz_dirBS[k]);
		 ut = fabs(-VeloY - VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-ny_dirBS[k]-nz_dirBS[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN))/(one+q) - c1over54 * drho;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
		 VeloY = vx2 - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * ny_dirBN[k];
		 VeloZ = vx3 - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * nz_dirBN[k];
		 un = fabs( (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * ny_dirBN[k] - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * nz_dirBN[k]);
		 ut = fabs( VeloY - VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs( ny_dirBN[k]-nz_dirBN[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS))/(one+q) - c1over54 * drho;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
		 VeloY = vx2 - (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * ny_dirTS[k];
		 VeloZ = vx3 - (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * nz_dirTS[k];
		 un = fabs(-(vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * ny_dirTS[k] + (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * nz_dirTS[k]);
		 ut = fabs(-VeloY + VeloZ);
		 tangential = ut / (ut + un + smallSingle);
		 qSlip = sliplength * fabs(-ny_dirTS[k]+nz_dirTS[k]);		//sliplength * e_i * n_i
		 //qSlip = (qSlip < zero) ? zero:qSlip;
		 //tangential = (tangential > one) ? one:tangential;
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN))/(one+q) - c1over54 * drho;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
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
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW))/(one+q) - c1over216 * drho;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
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
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE))/(one+q) - c1over216 * drho;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
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
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW))/(one+q) - c1over216 * drho;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
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
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE))/(one+q) - c1over216 * drho;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
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
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW))/(one+q) - c1over216 * drho;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
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
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE))/(one+q) - c1over216 * drho;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
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
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW))/(one+q) - c1over216 * drho;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
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
		 q = (q + qSlip)/(one + qSlip * (one - tangential) / (smallSingle + q));
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE))/(one+q) - c1over216 * drho;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


