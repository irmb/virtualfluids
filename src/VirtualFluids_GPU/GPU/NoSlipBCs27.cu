/* Device code */
#include "LBM/D3Q27.h"
#include "GPU/constant.h"

//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDeviceCompThinWallsPartOne27(int inx,
														 int iny,
														 doubflo* DD, 
														 int* k_Q, 
														 doubflo* QQ,
														 unsigned int sizeQ,
														 int kQ, 
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

   if(k<kQ)
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
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
            f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[dirE   ])[ke   ];
      f_E    = (D.f[dirW   ])[kw   ];
      f_S    = (D.f[dirN   ])[kn   ];
      f_N    = (D.f[dirS   ])[ks   ];
      f_B    = (D.f[dirT   ])[kt   ];
      f_T    = (D.f[dirB   ])[kb   ];
      f_SW   = (D.f[dirNE  ])[kne  ];
      f_NE   = (D.f[dirSW  ])[ksw  ];
      f_NW   = (D.f[dirSE  ])[kse  ];
      f_SE   = (D.f[dirNW  ])[knw  ];
      f_BW   = (D.f[dirTE  ])[kte  ];
      f_TE   = (D.f[dirBW  ])[kbw  ];
      f_TW   = (D.f[dirBE  ])[kbe  ];
      f_BE   = (D.f[dirTW  ])[ktw  ];
      f_BS   = (D.f[dirTN  ])[ktn  ];
      f_TN   = (D.f[dirBS  ])[kbs  ];
      f_TS   = (D.f[dirBN  ])[kbn  ];
      f_BN   = (D.f[dirTS  ])[kts  ];
      f_BSW  = (D.f[dirTNE ])[ktne ];
      f_BNE  = (D.f[dirTSW ])[ktsw ];
      f_BNW  = (D.f[dirTSE ])[ktse ];
      f_BSE  = (D.f[dirTNW ])[ktnw ];
      f_TSW  = (D.f[dirBNE ])[kbne ];
      f_TNE  = (D.f[dirBSW ])[kbsw ];
      f_TNW  = (D.f[dirBSE ])[kbse ];
      f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (one + drho); 


      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S)) / (one + drho); 

      vx3    =    (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                  (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                  (f_T - f_B)) / (one + drho); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (one + drho);

      //////////////////////////////////////////////////////////////////////////
      //if (evenOrOdd==false)
      //{
      //   D.f[dirE   ] = &DD[dirE   *size_Mat];
      //   D.f[dirW   ] = &DD[dirW   *size_Mat];
      //   D.f[dirN   ] = &DD[dirN   *size_Mat];
      //   D.f[dirS   ] = &DD[dirS   *size_Mat];
      //   D.f[dirT   ] = &DD[dirT   *size_Mat];
      //   D.f[dirB   ] = &DD[dirB   *size_Mat];
      //   D.f[dirNE  ] = &DD[dirNE  *size_Mat];
      //   D.f[dirSW  ] = &DD[dirSW  *size_Mat];
      //   D.f[dirSE  ] = &DD[dirSE  *size_Mat];
      //   D.f[dirNW  ] = &DD[dirNW  *size_Mat];
      //   D.f[dirTE  ] = &DD[dirTE  *size_Mat];
      //   D.f[dirBW  ] = &DD[dirBW  *size_Mat];
      //   D.f[dirBE  ] = &DD[dirBE  *size_Mat];
      //   D.f[dirTW  ] = &DD[dirTW  *size_Mat];
      //   D.f[dirTN  ] = &DD[dirTN  *size_Mat];
      //   D.f[dirBS  ] = &DD[dirBS  *size_Mat];
      //   D.f[dirBN  ] = &DD[dirBN  *size_Mat];
      //   D.f[dirTS  ] = &DD[dirTS  *size_Mat];
      //   D.f[dirZERO] = &DD[dirZERO*size_Mat];
      //   D.f[dirTNE ] = &DD[dirTNE *size_Mat];
      //   D.f[dirTSW ] = &DD[dirTSW *size_Mat];
      //   D.f[dirTSE ] = &DD[dirTSE *size_Mat];
      //   D.f[dirTNW ] = &DD[dirTNW *size_Mat];
      //   D.f[dirBNE ] = &DD[dirBNE *size_Mat];
      //   D.f[dirBSW ] = &DD[dirBSW *size_Mat];
      //   D.f[dirBSE ] = &DD[dirBSE *size_Mat];
      //   D.f[dirBNW ] = &DD[dirBNW *size_Mat];
      //} 
      //else
      //{
      //   D.f[dirW   ] = &DD[dirE   *size_Mat];
      //   D.f[dirE   ] = &DD[dirW   *size_Mat];
      //   D.f[dirS   ] = &DD[dirN   *size_Mat];
      //   D.f[dirN   ] = &DD[dirS   *size_Mat];
      //   D.f[dirB   ] = &DD[dirT   *size_Mat];
      //   D.f[dirT   ] = &DD[dirB   *size_Mat];
      //   D.f[dirSW  ] = &DD[dirNE  *size_Mat];
      //   D.f[dirNE  ] = &DD[dirSW  *size_Mat];
      //   D.f[dirNW  ] = &DD[dirSE  *size_Mat];
      //   D.f[dirSE  ] = &DD[dirNW  *size_Mat];
      //   D.f[dirBW  ] = &DD[dirTE  *size_Mat];
      //   D.f[dirTE  ] = &DD[dirBW  *size_Mat];
      //   D.f[dirTW  ] = &DD[dirBE  *size_Mat];
      //   D.f[dirBE  ] = &DD[dirTW  *size_Mat];
      //   D.f[dirBS  ] = &DD[dirTN  *size_Mat];
      //   D.f[dirTN  ] = &DD[dirBS  *size_Mat];
      //   D.f[dirTS  ] = &DD[dirBN  *size_Mat];
      //   D.f[dirBN  ] = &DD[dirTS  *size_Mat];
      //   D.f[dirZERO] = &DD[dirZERO*size_Mat];
      //   D.f[dirTNE ] = &DD[dirBSW *size_Mat];
      //   D.f[dirTSW ] = &DD[dirBNE *size_Mat];
      //   D.f[dirTSE ] = &DD[dirBNW *size_Mat];
      //   D.f[dirTNW ] = &DD[dirBSE *size_Mat];
      //   D.f[dirBNE ] = &DD[dirTSW *size_Mat];
      //   D.f[dirBSW ] = &DD[dirTNE *size_Mat];
      //   D.f[dirBSE ] = &DD[dirTNW *size_Mat];
      //   D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      //}
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
         //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W))/(one+q);
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E))/(one+q);
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S))/(one+q);
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N))/(one+q);
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
         (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B))/(one+q);
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
         (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T))/(one+q);
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW))/(one+q);
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE))/(one+q);
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW))/(one+q);
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE))/(one+q);
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW))/(one+q);
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE))/(one+q);
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW))/(one+q);
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE))/(one+q);
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS))/(one+q);
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN))/(one+q);
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS))/(one+q);
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN))/(one+q);
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW))/(one+q);
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE))/(one+q);
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW))/(one+q);
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE))/(one+q);
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW))/(one+q);
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE))/(one+q);
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW))/(one+q);
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE))/(one+q);
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDeviceCompThinWallsPartTwo27(int inx,
														 int iny,
														 doubflo* DD, 
														 int* k_Q, 
														 doubflo* QQ,
														 unsigned int sizeQ,
														 int kQ, 
														 doubflo om1, 
														 unsigned int* geom,
														 unsigned int* neighborX,
														 unsigned int* neighborY,
														 unsigned int* neighborZ,
														 unsigned int* neighborWSB,
														 unsigned int size_Mat, 
														 bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
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
	  //anti ET intermediate steps
	  unsigned int kmmm = neighborWSB[KQK]; // -1 -1 -1
	  unsigned int k0mm = neighborX[kmmm];  //  0 -1 -1
	  unsigned int km0m = neighborY[kmmm];  // -1  0 -1
	  unsigned int kmm0 = neighborZ[kmmm];  // -1 -1  0
	  unsigned int k0m0 = neighborX[kmm0];  //  0 -1  0
	  unsigned int km00 = neighborY[kmm0];  // -1  0  0
	  /////////////////////////////////////////////////
	  //final indices for anti ET
	  unsigned int kpmm = neighborX[k0mm];  //  1 -1 -1
	  unsigned int kmpm = neighborY[km0m];  // -1  1 -1
	  unsigned int kmmp = neighborZ[kmm0];  // -1 -1  1
	  unsigned int kmp0 = neighborY[km00];  // -1  1  0
	  unsigned int km0p = neighborZ[km00];  // -1  0  1
	  unsigned int k0mp = neighborZ[k0m0];  //  0 -1  1
	  ////////////////////////////////////////////////////////////////////////////////
	  Distributions27 D, DN;
	  if (evenOrOdd == true)
	  {
		  D.f[dirE] = &DD[dirE   *size_Mat];
		  D.f[dirW] = &DD[dirW   *size_Mat];
		  D.f[dirN] = &DD[dirN   *size_Mat];
		  D.f[dirS] = &DD[dirS   *size_Mat];
		  D.f[dirT] = &DD[dirT   *size_Mat];
		  D.f[dirB] = &DD[dirB   *size_Mat];
		  D.f[dirNE] = &DD[dirNE  *size_Mat];
		  D.f[dirSW] = &DD[dirSW  *size_Mat];
		  D.f[dirSE] = &DD[dirSE  *size_Mat];
		  D.f[dirNW] = &DD[dirNW  *size_Mat];
		  D.f[dirTE] = &DD[dirTE  *size_Mat];
		  D.f[dirBW] = &DD[dirBW  *size_Mat];
		  D.f[dirBE] = &DD[dirBE  *size_Mat];
		  D.f[dirTW] = &DD[dirTW  *size_Mat];
		  D.f[dirTN] = &DD[dirTN  *size_Mat];
		  D.f[dirBS] = &DD[dirBS  *size_Mat];
		  D.f[dirBN] = &DD[dirBN  *size_Mat];
		  D.f[dirTS] = &DD[dirTS  *size_Mat];
		  D.f[dirZERO] = &DD[dirZERO*size_Mat];
		  D.f[dirTNE] = &DD[dirTNE *size_Mat];
		  D.f[dirTSW] = &DD[dirTSW *size_Mat];
		  D.f[dirTSE] = &DD[dirTSE *size_Mat];
		  D.f[dirTNW] = &DD[dirTNW *size_Mat];
		  D.f[dirBNE] = &DD[dirBNE *size_Mat];
		  D.f[dirBSW] = &DD[dirBSW *size_Mat];
		  D.f[dirBSE] = &DD[dirBSE *size_Mat];
		  D.f[dirBNW] = &DD[dirBNW *size_Mat];
	  }
	  else
	  {
		  D.f[dirW] = &DD[dirE   *size_Mat];
		  D.f[dirE] = &DD[dirW   *size_Mat];
		  D.f[dirS] = &DD[dirN   *size_Mat];
		  D.f[dirN] = &DD[dirS   *size_Mat];
		  D.f[dirB] = &DD[dirT   *size_Mat];
		  D.f[dirT] = &DD[dirB   *size_Mat];
		  D.f[dirSW] = &DD[dirNE  *size_Mat];
		  D.f[dirNE] = &DD[dirSW  *size_Mat];
		  D.f[dirNW] = &DD[dirSE  *size_Mat];
		  D.f[dirSE] = &DD[dirNW  *size_Mat];
		  D.f[dirBW] = &DD[dirTE  *size_Mat];
		  D.f[dirTE] = &DD[dirBW  *size_Mat];
		  D.f[dirTW] = &DD[dirBE  *size_Mat];
		  D.f[dirBE] = &DD[dirTW  *size_Mat];
		  D.f[dirBS] = &DD[dirTN  *size_Mat];
		  D.f[dirTN] = &DD[dirBS  *size_Mat];
		  D.f[dirTS] = &DD[dirBN  *size_Mat];
		  D.f[dirBN] = &DD[dirTS  *size_Mat];
		  D.f[dirZERO] = &DD[dirZERO*size_Mat];
		  D.f[dirTNE] = &DD[dirBSW *size_Mat];
		  D.f[dirTSW] = &DD[dirBNE *size_Mat];
		  D.f[dirTSE] = &DD[dirBNW *size_Mat];
		  D.f[dirTNW] = &DD[dirBSE *size_Mat];
		  D.f[dirBNE] = &DD[dirTSW *size_Mat];
		  D.f[dirBSW] = &DD[dirTNE *size_Mat];
		  D.f[dirBSE] = &DD[dirTNW *size_Mat];
		  D.f[dirBNW] = &DD[dirTSE *size_Mat];
	  }
	  if (evenOrOdd==false)
      {
         DN.f[dirE   ] = &DD[dirE   *size_Mat];
         DN.f[dirW   ] = &DD[dirW   *size_Mat];
         DN.f[dirN   ] = &DD[dirN   *size_Mat];
         DN.f[dirS   ] = &DD[dirS   *size_Mat];
         DN.f[dirT   ] = &DD[dirT   *size_Mat];
         DN.f[dirB   ] = &DD[dirB   *size_Mat];
         DN.f[dirNE  ] = &DD[dirNE  *size_Mat];
         DN.f[dirSW  ] = &DD[dirSW  *size_Mat];
         DN.f[dirSE  ] = &DD[dirSE  *size_Mat];
         DN.f[dirNW  ] = &DD[dirNW  *size_Mat];
         DN.f[dirTE  ] = &DD[dirTE  *size_Mat];
         DN.f[dirBW  ] = &DD[dirBW  *size_Mat];
         DN.f[dirBE  ] = &DD[dirBE  *size_Mat];
         DN.f[dirTW  ] = &DD[dirTW  *size_Mat];
         DN.f[dirTN  ] = &DD[dirTN  *size_Mat];
         DN.f[dirBS  ] = &DD[dirBS  *size_Mat];
         DN.f[dirBN  ] = &DD[dirBN  *size_Mat];
         DN.f[dirTS  ] = &DD[dirTS  *size_Mat];
         DN.f[dirZERO] = &DD[dirZERO*size_Mat];
         DN.f[dirTNE ] = &DD[dirTNE *size_Mat];
         DN.f[dirTSW ] = &DD[dirTSW *size_Mat];
         DN.f[dirTSE ] = &DD[dirTSE *size_Mat];
         DN.f[dirTNW ] = &DD[dirTNW *size_Mat];
         DN.f[dirBNE ] = &DD[dirBNE *size_Mat];
         DN.f[dirBSW ] = &DD[dirBSW *size_Mat];
         DN.f[dirBSE ] = &DD[dirBSE *size_Mat];
         DN.f[dirBNW ] = &DD[dirBNW *size_Mat];
      } 
      else
      {
         DN.f[dirW   ] = &DD[dirE   *size_Mat];
         DN.f[dirE   ] = &DD[dirW   *size_Mat];
         DN.f[dirS   ] = &DD[dirN   *size_Mat];
         DN.f[dirN   ] = &DD[dirS   *size_Mat];
         DN.f[dirB   ] = &DD[dirT   *size_Mat];
         DN.f[dirT   ] = &DD[dirB   *size_Mat];
         DN.f[dirSW  ] = &DD[dirNE  *size_Mat];
         DN.f[dirNE  ] = &DD[dirSW  *size_Mat];
         DN.f[dirNW  ] = &DD[dirSE  *size_Mat];
         DN.f[dirSE  ] = &DD[dirNW  *size_Mat];
         DN.f[dirBW  ] = &DD[dirTE  *size_Mat];
         DN.f[dirTE  ] = &DD[dirBW  *size_Mat];
         DN.f[dirTW  ] = &DD[dirBE  *size_Mat];
         DN.f[dirBE  ] = &DD[dirTW  *size_Mat];
         DN.f[dirBS  ] = &DD[dirTN  *size_Mat];
         DN.f[dirTN  ] = &DD[dirBS  *size_Mat];
         DN.f[dirTS  ] = &DD[dirBN  *size_Mat];
         DN.f[dirBN  ] = &DD[dirTS  *size_Mat];
         DN.f[dirZERO] = &DD[dirZERO*size_Mat];
         DN.f[dirTNE ] = &DD[dirBSW *size_Mat];
         DN.f[dirTSW ] = &DD[dirBNE *size_Mat];
         DN.f[dirTSE ] = &DD[dirBNW *size_Mat];
         DN.f[dirTNW ] = &DD[dirBSE *size_Mat];
         DN.f[dirBNE ] = &DD[dirTSW *size_Mat];
         DN.f[dirBSW ] = &DD[dirTNE *size_Mat];
         DN.f[dirBSE ] = &DD[dirTNW *size_Mat];
         DN.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //directions allways exchange
	  //(-1 -1 -1) (-1  0  0) ( 0 -1  0) ( 0  0 -1) (-1 -1  0) (-1  0 -1) ( 0 -1 -1) ( 1  1 -1) ( 1 -1  1) (-1  1  1) ( 1 -1  0) ( 1  0 -1) ( 0  1 -1)
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //directions exchange if solid neighbor
	  //( 1  1  1) ( 1  0  0) ( 0  1  0) ( 0  0  1) ( 1  1  0) ( 1  0  1) ( 0  1  1) (-1 -1  1) (-1  1 -1) ( 1 -1 -1) (-1  1  0) (-1  0  1) ( 0 -1  1)
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  doubflo q, tmp;
      q = q_dirE[k];   if (q>=zero && q<=one){ if (geom[kw  ] < GEO_FLUID){tmp = (DN.f[dirW  ])[kw  ]; (DN.f[dirW  ])[kw  ]=(D.f[dirW  ])[kw  ]; (D.f[dirW  ])[kw  ]=tmp;}}
	  q = q_dirW[k];   if (q>=zero && q<=one){                            {tmp = (DN.f[dirE  ])[ke  ]; (DN.f[dirE  ])[ke  ]=(D.f[dirE  ])[ke  ]; (D.f[dirE  ])[ke  ]=tmp;}}
      q = q_dirN[k];   if (q>=zero && q<=one){ if (geom[ks  ] < GEO_FLUID){tmp = (DN.f[dirS  ])[ks  ]; (DN.f[dirS  ])[ks  ]=(D.f[dirS  ])[ks  ]; (D.f[dirS  ])[ks  ]=tmp;}}
      q = q_dirS[k];   if (q>=zero && q<=one){                            {tmp = (DN.f[dirN  ])[kn  ]; (DN.f[dirN  ])[kn  ]=(D.f[dirN  ])[kn  ]; (D.f[dirN  ])[kn  ]=tmp;}}
      q = q_dirT[k];   if (q>=zero && q<=one){ if (geom[kb  ] < GEO_FLUID){tmp = (DN.f[dirB  ])[kb  ]; (DN.f[dirB  ])[kb  ]=(D.f[dirB  ])[kb  ]; (D.f[dirB  ])[kb  ]=tmp;}}
      q = q_dirB[k];   if (q>=zero && q<=one){                            {tmp = (DN.f[dirT  ])[kt  ]; (DN.f[dirT  ])[kt  ]=(D.f[dirT  ])[kt  ]; (D.f[dirT  ])[kt  ]=tmp;}}
      q = q_dirNE[k];  if (q>=zero && q<=one){ if (geom[ksw ] < GEO_FLUID){tmp = (DN.f[dirSW ])[ksw ]; (DN.f[dirSW ])[ksw ]=(D.f[dirSW ])[ksw ]; (D.f[dirSW ])[ksw ]=tmp;}}
      q = q_dirSW[k];  if (q>=zero && q<=one){                            {tmp = (DN.f[dirNE ])[kne ]; (DN.f[dirNE ])[kne ]=(D.f[dirNE ])[kne ]; (D.f[dirNE ])[kne ]=tmp;}}
      q = q_dirSE[k];  if (q>=zero && q<=one){                            {tmp = (DN.f[dirNW ])[knw ]; (DN.f[dirNW ])[knw ]=(D.f[dirNW ])[knw ]; (D.f[dirNW ])[knw ]=tmp;}}
      q = q_dirNW[k];  if (q>=zero && q<=one){ if (geom[kmp0] < GEO_FLUID){tmp = (DN.f[dirSE ])[kse ]; (DN.f[dirSE ])[kse ]=(D.f[dirSE ])[kse ]; (D.f[dirSE ])[kse ]=tmp;}}
      q = q_dirTE[k];  if (q>=zero && q<=one){ if (geom[kbw ] < GEO_FLUID){tmp = (DN.f[dirBW ])[kbw ]; (DN.f[dirBW ])[kbw ]=(D.f[dirBW ])[kbw ]; (D.f[dirBW ])[kbw ]=tmp;}}
      q = q_dirBW[k];  if (q>=zero && q<=one){                            {tmp = (DN.f[dirTE ])[kte ]; (DN.f[dirTE ])[kte ]=(D.f[dirTE ])[kte ]; (D.f[dirTE ])[kte ]=tmp;}}
      q = q_dirBE[k];  if (q>=zero && q<=one){                            {tmp = (DN.f[dirTW ])[ktw ]; (DN.f[dirTW ])[ktw ]=(D.f[dirTW ])[ktw ]; (D.f[dirTW ])[ktw ]=tmp;}}
      q = q_dirTW[k];  if (q>=zero && q<=one){ if (geom[km0p] < GEO_FLUID){tmp = (DN.f[dirBE ])[kbe ]; (DN.f[dirBE ])[kbe ]=(D.f[dirBE ])[kbe ]; (D.f[dirBE ])[kbe ]=tmp;}}
      q = q_dirTN[k];  if (q>=zero && q<=one){ if (geom[kbs ] < GEO_FLUID){tmp = (DN.f[dirBS ])[kbs ]; (DN.f[dirBS ])[kbs ]=(D.f[dirBS ])[kbs ]; (D.f[dirBS ])[kbs ]=tmp;}}
      q = q_dirBS[k];  if (q>=zero && q<=one){                            {tmp = (DN.f[dirTN ])[ktn ]; (DN.f[dirTN ])[ktn ]=(D.f[dirTN ])[ktn ]; (D.f[dirTN ])[ktn ]=tmp;}}
      q = q_dirBN[k];  if (q>=zero && q<=one){                            {tmp = (DN.f[dirTS ])[kts ]; (DN.f[dirTS ])[kts ]=(D.f[dirTS ])[kts ]; (D.f[dirTS ])[kts ]=tmp;}}
      q = q_dirTS[k];  if (q>=zero && q<=one){ if (geom[k0mp] < GEO_FLUID){tmp = (DN.f[dirBN ])[kbn ]; (DN.f[dirBN ])[kbn ]=(D.f[dirBN ])[kbn ]; (D.f[dirBN ])[kbn ]=tmp;}}
      q = q_dirTNE[k]; if (q>=zero && q<=one){ if (geom[kbsw] < GEO_FLUID){tmp = (DN.f[dirBSW])[kbsw]; (DN.f[dirBSW])[kbsw]=(D.f[dirBSW])[kbsw]; (D.f[dirBSW])[kbsw]=tmp;}}
      q = q_dirBSW[k]; if (q>=zero && q<=one){                            {tmp = (DN.f[dirTNE])[ktne]; (DN.f[dirTNE])[ktne]=(D.f[dirTNE])[ktne]; (D.f[dirTNE])[ktne]=tmp;}}
      q = q_dirBNE[k]; if (q>=zero && q<=one){                            {tmp = (DN.f[dirTSW])[ktsw]; (DN.f[dirTSW])[ktsw]=(D.f[dirTSW])[ktsw]; (D.f[dirTSW])[ktsw]=tmp;}}
      q = q_dirTSW[k]; if (q>=zero && q<=one){ if (geom[kmmp] < GEO_FLUID){tmp = (DN.f[dirBNE])[kbne]; (DN.f[dirBNE])[kbne]=(D.f[dirBNE])[kbne]; (D.f[dirBNE])[kbne]=tmp;}}
      q = q_dirTSE[k]; if (q>=zero && q<=one){                            {tmp = (DN.f[dirBNW])[kbnw]; (DN.f[dirBNW])[kbnw]=(D.f[dirBNW])[kbnw]; (D.f[dirBNW])[kbnw]=tmp;}}
      q = q_dirBNW[k]; if (q>=zero && q<=one){ if (geom[kmpm] < GEO_FLUID){tmp = (DN.f[dirTSE])[ktse]; (DN.f[dirTSE])[ktse]=(D.f[dirTSE])[ktse]; (D.f[dirTSE])[ktse]=tmp;}}
      q = q_dirBSE[k]; if (q>=zero && q<=one){ if (geom[kpmm] < GEO_FLUID){tmp = (DN.f[dirTNW])[ktnw]; (DN.f[dirTNW])[ktnw]=(D.f[dirTNW])[ktnw]; (D.f[dirTNW])[ktnw]=tmp;}}
      q = q_dirTNW[k]; if (q>=zero && q<=one){                            {tmp = (DN.f[dirBSE])[kbse]; (DN.f[dirBSE])[kbse]=(D.f[dirBSE])[kbse]; (D.f[dirBSE])[kbse]=tmp;}}
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDevice3rdMomentsComp27(  int inx,
													 int iny,
													 doubflo* DD, 
													 int* k_Q, 
													 doubflo* QQ,
													 unsigned int sizeQ,
													 int kQ, 
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

   if(k<kQ)
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
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
            f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[dirE   ])[ke   ];
      f_E    = (D.f[dirW   ])[kw   ];
      f_S    = (D.f[dirN   ])[kn   ];
      f_N    = (D.f[dirS   ])[ks   ];
      f_B    = (D.f[dirT   ])[kt   ];
      f_T    = (D.f[dirB   ])[kb   ];
      f_SW   = (D.f[dirNE  ])[kne  ];
      f_NE   = (D.f[dirSW  ])[ksw  ];
      f_NW   = (D.f[dirSE  ])[kse  ];
      f_SE   = (D.f[dirNW  ])[knw  ];
      f_BW   = (D.f[dirTE  ])[kte  ];
      f_TE   = (D.f[dirBW  ])[kbw  ];
      f_TW   = (D.f[dirBE  ])[kbe  ];
      f_BE   = (D.f[dirTW  ])[ktw  ];
      f_BS   = (D.f[dirTN  ])[ktn  ];
      f_TN   = (D.f[dirBS  ])[kbs  ];
      f_TS   = (D.f[dirBN  ])[kbn  ];
      f_BN   = (D.f[dirTS  ])[kts  ];
      f_BSW  = (D.f[dirTNE ])[ktne ];
      f_BNE  = (D.f[dirTSW ])[ktsw ];
      f_BNW  = (D.f[dirTSE ])[ktse ];
      f_BSE  = (D.f[dirTNW ])[ktnw ];
      f_TSW  = (D.f[dirBNE ])[kbne ];
      f_TNE  = (D.f[dirBSW ])[kbsw ];
      f_TNW  = (D.f[dirBSE ])[kbse ];
      f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q, m3;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (one + drho); 


      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S)) / (one + drho); 

      vx3    =    (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
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
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_E - f_W - two * drho * c2over27 * (three*( vx1        ));
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W-m3+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W))/(one+q)+(m3*c1o2);
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_W - f_E - two * drho * c2over27 * (three*(-vx1        ));
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E-m3+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E))/(one+q)+(m3*c1o2);
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_N - f_S - two * drho * c2over27 * (three*( vx2        ));
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S-m3+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S))/(one+q)+(m3*c1o2);
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_S - f_N - two * drho * c2over27 * (three*(   -vx2     ));
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N-m3+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N))/(one+q)+(m3*c1o2);
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_T - f_B - two * drho * c2over27 * (three*(         vx3));
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
         (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B-m3+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B))/(one+q)+(m3*c1o2);
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_B - f_T - two * drho * c2over27 * (three*(        -vx3));
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
         (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T-m3+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T))/(one+q)+(m3*c1o2);
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_NE - f_SW - two * drho * c1over54 * (three*( vx1+vx2    ));
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW-m3+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW))/(one+q)+(m3*c1o2);
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_SW - f_NE - two * drho * c1over54 * (three*(-vx1-vx2    ));
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE-m3+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE))/(one+q)+(m3*c1o2);
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_SE - f_NW - two * drho * c1over54 * (three*( vx1-vx2    ));
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW-m3+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW))/(one+q)+(m3*c1o2);
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_NW - f_SE - two * drho * c1over54 * (three*(-vx1+vx2    ));
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE-m3+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE))/(one+q)+(m3*c1o2);
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_TE - f_BW - two * drho * c1over54 * (three*( vx1    +vx3));
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW-m3+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW))/(one+q)+(m3*c1o2);
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_BW - f_TE - two * drho * c1over54 * (three*(-vx1    -vx3));
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE-m3+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE))/(one+q)+(m3*c1o2);
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_BE - f_TW - two * drho * c1over54 * (three*( vx1    -vx3));
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW-m3+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW))/(one+q)+(m3*c1o2);
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_TW - f_BE - two * drho * c1over54 * (three*(-vx1    +vx3));
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE-m3+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE))/(one+q)+(m3*c1o2);
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_TN - f_BS - two * drho * c1over54 * (three*(     vx2+vx3));
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS-m3+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS))/(one+q)+(m3*c1o2);
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_BS - f_TN - two * drho * c1over54 * (three*(    -vx2-vx3));
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN-m3+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN))/(one+q)+(m3*c1o2);
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_BN - f_TS - two * drho * c1over54 * (three*(     vx2-vx3));
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS-m3+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS))/(one+q)+(m3*c1o2);
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_TS - f_BN - two * drho * c1over54 * (three*(    -vx2+vx3));
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN-m3+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN))/(one+q)+(m3*c1o2);
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_TNE - f_BSW - two * drho * c1over216 * (three*( vx1+vx2+vx3));
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW-m3+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW))/(one+q)+(m3*c1o2);
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_BSW - f_TNE - two * drho * c1over216 * (three*(-vx1-vx2-vx3));
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE-m3+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE))/(one+q)+(m3*c1o2);
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_BNE - f_TSW - two * drho * c1over216 * (three*( vx1+vx2-vx3));
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW-m3+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW))/(one+q)+(m3*c1o2);
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_TSW - f_BNE - two * drho * c1over216 * (three*(-vx1-vx2+vx3));
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE-m3+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE))/(one+q)+(m3*c1o2);
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_TSE - f_BNW - two * drho * c1over216 * (three*( vx1-vx2+vx3));
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW-m3+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW))/(one+q)+(m3*c1o2);
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_BNW - f_TSE - two * drho * c1over216 * (three*(-vx1+vx2-vx3));
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE-m3+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE))/(one+q)+(m3*c1o2);
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_BSE - f_TNW - two * drho * c1over216 * (three*( vx1-vx2-vx3));
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW-m3+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW))/(one+q)+(m3*c1o2);
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
		 m3 = f_TNW - f_BSE - two * drho * c1over216 * (three*(-vx1+vx2+vx3));
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE-m3+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE))/(one+q)+(m3*c1o2);
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDeviceIncompHighNu27(int inx,
												 int iny,
												 doubflo* DD, 
												 int* k_Q, 
												 doubflo* QQ,
												 unsigned int sizeQ,
												 int kQ, 
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

   if(k<kQ)
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
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
            f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_E   = (D.f[dirE   ])[ke   ];
      f_W   = (D.f[dirW   ])[kw   ];
      f_N   = (D.f[dirN   ])[kn   ];
      f_S   = (D.f[dirS   ])[ks   ];
      f_T   = (D.f[dirT   ])[kt   ];
      f_B   = (D.f[dirB   ])[kb   ];
      f_NE  = (D.f[dirNE  ])[kne  ];
      f_SW  = (D.f[dirSW  ])[ksw  ];
      f_SE  = (D.f[dirSE  ])[kse  ];
      f_NW  = (D.f[dirNW  ])[knw  ];
      f_TE  = (D.f[dirTE  ])[kte  ];
      f_BW  = (D.f[dirBW  ])[kbw  ];
      f_BE  = (D.f[dirBE  ])[kbe  ];
      f_TW  = (D.f[dirTW  ])[ktw  ];
      f_TN  = (D.f[dirTN  ])[ktn  ];
      f_BS  = (D.f[dirBS  ])[kbs  ];
      f_BN  = (D.f[dirBN  ])[kbn  ];
      f_TS  = (D.f[dirTS  ])[kts  ];
      f_TNE = (D.f[dirTNE ])[ktne ];
      f_TSW = (D.f[dirTSW ])[ktsw ];
      f_TSE = (D.f[dirTSE ])[ktse ];
      f_TNW = (D.f[dirTNW ])[ktnw ];
      f_BNE = (D.f[dirBNE ])[kbne ];
      f_BSW = (D.f[dirBSW ])[kbsw ];
      f_BSE = (D.f[dirBSE ])[kbse ];
      f_BNW = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W));// / (one + drho); 


      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S));// / (one + drho); 

      vx3    =    (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                  (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                  (f_T - f_B));// / (one + drho); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);// * (one + drho);

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
      doubflo VeloX = zero;
      doubflo VeloY = zero;
      doubflo VeloZ = zero;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) /** (one + drho)*/-cu_sq); 
         (D.f[dirW])[kw]=((one - q) * f_E + q * ((f_E + f_W) * (one - om1) + om1 * two * feq) - six*c2over27*( VeloX     )) / (q + one) ;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) /** (one + drho)*/-cu_sq); 
         (D.f[dirE])[ke]=((one - q) * f_W + q * ((f_W + f_E) * (one - om1) + om1 * two * feq) - six*c2over27*(-VeloX     )) / (q + one) ;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirS])[ks]=((one - q) * f_N + q * ((f_N + f_S) * (one - om1) + om1 * two * feq) - six*c2over27*( VeloY     )) / (q + one) ;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirN])[kn]=((one - q) * f_S + q * ((f_S + f_N) * (one - om1) + om1 * two * feq) - six*c2over27*(-VeloY     )) / (q + one) ;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirB])[kb]=((one - q) * f_T + q * ((f_T + f_B) * (one - om1) + om1 * two * feq) - six*c2over27*( VeloZ     )) / (q + one) ;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirT])[kt]=((one - q) * f_B + q * ((f_B + f_T) * (one - om1) + om1 * two * feq) - six*c2over27*(-VeloZ     )) / (q + one) ;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirSW])[ksw]=((one - q) * f_NE + q * ((f_NE + f_SW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX+VeloY)) / (q + one) ;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirNE])[kne]=((one - q) * f_SW + q * ((f_SW + f_NE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX-VeloY)) / (q + one) ;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirNW])[knw]=((one - q) * f_SE + q * ((f_SE + f_NW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX-VeloY)) / (q + one) ;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirSE])[kse]=((one - q) * f_NW + q * ((f_NW + f_SE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX+VeloY)) / (q + one) ;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBW])[kbw]=((one - q) * f_TE + q * ((f_TE + f_BW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX+VeloZ)) / (q + one) ;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTE])[kte]=((one - q) * f_BW + q * ((f_BW + f_TE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX-VeloZ)) / (q + one) ;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTW])[ktw]=((one - q) * f_BE + q * ((f_BE + f_TW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX-VeloZ)) / (q + one) ;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBE])[kbe]=((one - q) * f_TW + q * ((f_TW + f_BE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX+VeloZ)) / (q + one) ;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBS])[kbs]=((one - q) * f_TN + q * ((f_TN + f_BS) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloY+VeloZ)) / (q + one) ;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTN])[ktn]=((one - q) * f_BS + q * ((f_BS + f_TN) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloY-VeloZ)) / (q + one) ;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTS])[kts]=((one - q) * f_BN + q * ((f_BN + f_TS) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloY-VeloZ)) / (q + one) ;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBN])[kbn]=((one - q) * f_TS + q * ((f_TS + f_BN) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloY+VeloZ)) / (q + one) ;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBSW])[kbsw]=((one - q) * f_TNE + q * ((f_TNE + f_BSW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX+VeloY+VeloZ)) / (q + one) ;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTNE])[ktne]=((one - q) * f_BSW + q * ((f_BSW + f_TNE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX-VeloY-VeloZ)) / (q + one) ;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTSW])[ktsw]=((one - q) * f_BNE + q * ((f_BNE + f_TSW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX+VeloY-VeloZ)) / (q + one) ;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBNE])[kbne]=((one - q) * f_TSW + q * ((f_TSW + f_BNE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX-VeloY+VeloZ)) / (q + one) ;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBNW])[kbnw]=((one - q) * f_TSE + q * ((f_TSE + f_BNW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX-VeloY+VeloZ)) / (q + one) ;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTSE])[ktse]=((one - q) * f_BNW + q * ((f_BNW + f_TSE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX+VeloY-VeloZ)) / (q + one) ;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTNW])[ktnw]=((one - q) * f_BSE + q * ((f_BSE + f_TNW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX-VeloY-VeloZ)) / (q + one) ;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBSE])[kbse]=((one - q) * f_TNW + q * ((f_TNW + f_BSE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX+VeloY+VeloZ)) / (q + one) ;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDeviceCompHighNu27(  int inx,
												 int iny,
												 doubflo* DD, 
												 int* k_Q, 
												 doubflo* QQ,
												 unsigned int sizeQ,
												 int kQ, 
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

   if(k<kQ)
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
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
            f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_E   = (D.f[dirE   ])[ke   ];
      f_W   = (D.f[dirW   ])[kw   ];
      f_N   = (D.f[dirN   ])[kn   ];
      f_S   = (D.f[dirS   ])[ks   ];
      f_T   = (D.f[dirT   ])[kt   ];
      f_B   = (D.f[dirB   ])[kb   ];
      f_NE  = (D.f[dirNE  ])[kne  ];
      f_SW  = (D.f[dirSW  ])[ksw  ];
      f_SE  = (D.f[dirSE  ])[kse  ];
      f_NW  = (D.f[dirNW  ])[knw  ];
      f_TE  = (D.f[dirTE  ])[kte  ];
      f_BW  = (D.f[dirBW  ])[kbw  ];
      f_BE  = (D.f[dirBE  ])[kbe  ];
      f_TW  = (D.f[dirTW  ])[ktw  ];
      f_TN  = (D.f[dirTN  ])[ktn  ];
      f_BS  = (D.f[dirBS  ])[kbs  ];
      f_BN  = (D.f[dirBN  ])[kbn  ];
      f_TS  = (D.f[dirTS  ])[kts  ];
      f_TNE = (D.f[dirTNE ])[ktne ];
      f_TSW = (D.f[dirTSW ])[ktsw ];
      f_TSE = (D.f[dirTSE ])[ktse ];
      f_TNW = (D.f[dirTNW ])[ktnw ];
      f_BNE = (D.f[dirBNE ])[kbne ];
      f_BSW = (D.f[dirBSW ])[kbsw ];
      f_BSE = (D.f[dirBSE ])[kbse ];
      f_BNW = (D.f[dirBNW ])[kbnw ];
      //f_W    = (D.f[dirE   ])[ke   ];
      //f_E    = (D.f[dirW   ])[kw   ];
      //f_S    = (D.f[dirN   ])[kn   ];
      //f_N    = (D.f[dirS   ])[ks   ];
      //f_B    = (D.f[dirT   ])[kt   ];
      //f_T    = (D.f[dirB   ])[kb   ];
      //f_SW   = (D.f[dirNE  ])[kne  ];
      //f_NE   = (D.f[dirSW  ])[ksw  ];
      //f_NW   = (D.f[dirSE  ])[kse  ];
      //f_SE   = (D.f[dirNW  ])[knw  ];
      //f_BW   = (D.f[dirTE  ])[kte  ];
      //f_TE   = (D.f[dirBW  ])[kbw  ];
      //f_TW   = (D.f[dirBE  ])[kbe  ];
      //f_BE   = (D.f[dirTW  ])[ktw  ];
      //f_BS   = (D.f[dirTN  ])[ktn  ];
      //f_TN   = (D.f[dirBS  ])[kbs  ];
      //f_TS   = (D.f[dirBN  ])[kbn  ];
      //f_BN   = (D.f[dirTS  ])[kts  ];
      //f_BSW  = (D.f[dirTNE ])[ktne ];
      //f_BNE  = (D.f[dirTSW ])[ktsw ];
      //f_BNW  = (D.f[dirTSE ])[ktse ];
      //f_BSE  = (D.f[dirTNW ])[ktnw ];
      //f_TSW  = (D.f[dirBNE ])[kbne ];
      //f_TNE  = (D.f[dirBSW ])[kbsw ];
      //f_TNW  = (D.f[dirBSE ])[kbse ];
      //f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (one + drho); 


      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S)) / (one + drho); 

      vx3    =    (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
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
      doubflo VeloX = zero;
      doubflo VeloY = zero;
      doubflo VeloZ = zero;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirW])[kw]=((one - q) * f_E + q * ((f_E + f_W) * (one - om1) + om1 * two * feq) - six*c2over27*( VeloX     )) / (q + one) ;
         //(D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirE])[ke]=((one - q) * f_W + q * ((f_W + f_E) * (one - om1) + om1 * two * feq) - six*c2over27*(-VeloX     )) / (q + one) ;
         //(D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirS])[ks]=((one - q) * f_N + q * ((f_N + f_S) * (one - om1) + om1 * two * feq) - six*c2over27*( VeloY     )) / (q + one) ;
         //(D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirN])[kn]=((one - q) * f_S + q * ((f_S + f_N) * (one - om1) + om1 * two * feq) - six*c2over27*(-VeloY     )) / (q + one) ;
         //(D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
         (D.f[dirB])[kb]=((one - q) * f_T + q * ((f_T + f_B) * (one - om1) + om1 * two * feq) - six*c2over27*( VeloZ     )) / (q + one) ;
         //(D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
         (D.f[dirT])[kt]=((one - q) * f_B + q * ((f_B + f_T) * (one - om1) + om1 * two * feq) - six*c2over27*(-VeloZ     )) / (q + one) ;
         //(D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=((one - q) * f_NE + q * ((f_NE + f_SW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX+VeloY)) / (q + one) ;
         //(D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNE])[kne]=((one - q) * f_SW + q * ((f_SW + f_NE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX-VeloY)) / (q + one) ;
         //(D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNW])[knw]=((one - q) * f_SE + q * ((f_SE + f_NW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX-VeloY)) / (q + one) ;
         //(D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSE])[kse]=((one - q) * f_NW + q * ((f_NW + f_SE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX+VeloY)) / (q + one) ;
         //(D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=((one - q) * f_TE + q * ((f_TE + f_BW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX+VeloZ)) / (q + one) ;
         //(D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTE])[kte]=((one - q) * f_BW + q * ((f_BW + f_TE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX-VeloZ)) / (q + one) ;
         //(D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=((one - q) * f_BE + q * ((f_BE + f_TW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX-VeloZ)) / (q + one) ;
         //(D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=((one - q) * f_TW + q * ((f_TW + f_BE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX+VeloZ)) / (q + one) ;
         //(D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=((one - q) * f_TN + q * ((f_TN + f_BS) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloY+VeloZ)) / (q + one) ;
         //(D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=((one - q) * f_BS + q * ((f_BS + f_TN) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloY-VeloZ)) / (q + one) ;
         //(D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTS])[kts]=((one - q) * f_BN + q * ((f_BN + f_TS) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloY-VeloZ)) / (q + one) ;
         //(D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=((one - q) * f_TS + q * ((f_TS + f_BN) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloY+VeloZ)) / (q + one) ;
         //(D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=((one - q) * f_TNE + q * ((f_TNE + f_BSW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX+VeloY+VeloZ)) / (q + one) ;
         //(D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=((one - q) * f_BSW + q * ((f_BSW + f_TNE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX-VeloY-VeloZ)) / (q + one) ;
         //(D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=((one - q) * f_BNE + q * ((f_BNE + f_TSW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX+VeloY-VeloZ)) / (q + one) ;
         //(D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=((one - q) * f_TSW + q * ((f_TSW + f_BNE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX-VeloY+VeloZ)) / (q + one) ;
         //(D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=((one - q) * f_TSE + q * ((f_TSE + f_BNW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX-VeloY+VeloZ)) / (q + one) ;
         //(D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=((one - q) * f_BNW + q * ((f_BNW + f_TSE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX+VeloY-VeloZ)) / (q + one) ;
         //(D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=((one - q) * f_BSE + q * ((f_BSE + f_TNW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX-VeloY-VeloZ)) / (q + one) ;
         //(D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=((one - q) * f_TNW + q * ((f_TNW + f_BSE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX+VeloY+VeloZ)) / (q + one) ;
         //(D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDeviceComp27(int inx,
										 int iny,
										 doubflo* DD, 
										 int* k_Q, 
										 doubflo* QQ,
										 unsigned int sizeQ,
										 int kQ, 
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

   if(k<kQ)
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
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
            f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[dirE   ])[ke   ];
      f_E    = (D.f[dirW   ])[kw   ];
      f_S    = (D.f[dirN   ])[kn   ];
      f_N    = (D.f[dirS   ])[ks   ];
      f_B    = (D.f[dirT   ])[kt   ];
      f_T    = (D.f[dirB   ])[kb   ];
      f_SW   = (D.f[dirNE  ])[kne  ];
      f_NE   = (D.f[dirSW  ])[ksw  ];
      f_NW   = (D.f[dirSE  ])[kse  ];
      f_SE   = (D.f[dirNW  ])[knw  ];
      f_BW   = (D.f[dirTE  ])[kte  ];
      f_TE   = (D.f[dirBW  ])[kbw  ];
      f_TW   = (D.f[dirBE  ])[kbe  ];
      f_BE   = (D.f[dirTW  ])[ktw  ];
      f_BS   = (D.f[dirTN  ])[ktn  ];
      f_TN   = (D.f[dirBS  ])[kbs  ];
      f_TS   = (D.f[dirBN  ])[kbn  ];
      f_BN   = (D.f[dirTS  ])[kts  ];
      f_BSW  = (D.f[dirTNE ])[ktne ];
      f_BNE  = (D.f[dirTSW ])[ktsw ];
      f_BNW  = (D.f[dirTSE ])[ktse ];
      f_BSE  = (D.f[dirTNW ])[ktnw ];
      f_TSW  = (D.f[dirBNE ])[kbne ];
      f_TNE  = (D.f[dirBSW ])[kbsw ];
      f_TNW  = (D.f[dirBSE ])[kbse ];
      f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (one + drho); 


      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S)) / (one + drho); 

      vx3    =    (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
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
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W))/(one+q);
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E))/(one+q);
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S))/(one+q);
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N))/(one+q);
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
         (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B))/(one+q);
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
         (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T))/(one+q);
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW))/(one+q);
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE))/(one+q);
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW))/(one+q);
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE))/(one+q);
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW))/(one+q);
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE))/(one+q);
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW))/(one+q);
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE))/(one+q);
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS))/(one+q);
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN))/(one+q);
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS))/(one+q);
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN))/(one+q);
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW))/(one+q);
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE))/(one+q);
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW))/(one+q);
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE))/(one+q);
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW))/(one+q);
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE))/(one+q);
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW))/(one+q);
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE))/(one+q);
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDevice27(int inx,
                                     int iny,
                                     doubflo* DD, 
                                     int* k_Q, 
                                     doubflo* QQ,
                                     unsigned int sizeQ,
                                     int kQ, 
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

   if(k<kQ)
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
      //unsigned int nxny = nx*ny;
      //unsigned int kzero= KQK;
      //unsigned int ke   = KQK;
      //unsigned int kw   = KQK + 1;
      //unsigned int kn   = KQK;
      //unsigned int ks   = KQK + nx;
      //unsigned int kt   = KQK;
      //unsigned int kb   = KQK + nxny;
      //unsigned int ksw  = KQK + nx + 1;
      //unsigned int kne  = KQK;
      //unsigned int kse  = KQK + nx;
      //unsigned int knw  = KQK + 1;
      //unsigned int kbw  = KQK + nxny + 1;
      //unsigned int kte  = KQK;
      //unsigned int kbe  = KQK + nxny;
      //unsigned int ktw  = KQK + 1;
      //unsigned int kbs  = KQK + nxny + nx;
      //unsigned int ktn  = KQK;
      //unsigned int kbn  = KQK + nxny;
      //unsigned int kts  = KQK + nx;
      //unsigned int ktse = KQK + nx;
      //unsigned int kbnw = KQK + nxny + 1;
      //unsigned int ktnw = KQK + 1;
      //unsigned int kbse = KQK + nxny + nx;
      //unsigned int ktsw = KQK + nx + 1;
      //unsigned int kbne = KQK + nxny;
      //unsigned int ktne = KQK;
      //unsigned int kbsw = KQK + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
      //doubflo vx1, vx2, vx3, drho, feq, q;
      //drho    =   (D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
      //            (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
      //            (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
      //            (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
      //            (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
      //            (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
      //            (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
      //            (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
      //            (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
      //            (D.f[dirZERO])[kzero]+ 
      //            (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
      //            (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
      //            (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
      //            (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw];

      //vx1    =    (D.f[dirE   ])[ke  ]- (D.f[dirW   ])[kw  ]+ 
      //            (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]+
      //            (D.f[dirSE  ])[kse ]- (D.f[dirNW  ])[knw ]+
      //            (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]+
      //            (D.f[dirBE  ])[kbe ]- (D.f[dirTW  ])[ktw ]+
      //            (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]+ 
      //            (D.f[dirTSE ])[ktse]- (D.f[dirTNW ])[ktnw]+ 
      //            (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]+ 
      //            (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw];

      //vx2    =    (D.f[dirN   ])[kn  ]- (D.f[dirS   ])[ks  ]+
      //            (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]-
      //            (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
      //            (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]+
      //            (D.f[dirBN  ])[kbn ]- (D.f[dirTS  ])[kts ]+
      //            (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]- 
      //            (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
      //            (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
      //            (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw];

      //vx3    =    (D.f[dirT   ])[kt  ]- (D.f[dirB   ])[kb  ]+
      //            (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]-
      //            (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
      //            (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]-
      //            (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
      //            (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
      //            (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]- 
      //            (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
      //            (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw];

      //doubflo cu_sq=1.5f*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
            f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[dirE   ])[ke   ];
      f_E    = (D.f[dirW   ])[kw   ];
      f_S    = (D.f[dirN   ])[kn   ];
      f_N    = (D.f[dirS   ])[ks   ];
      f_B    = (D.f[dirT   ])[kt   ];
      f_T    = (D.f[dirB   ])[kb   ];
      f_SW   = (D.f[dirNE  ])[kne  ];
      f_NE   = (D.f[dirSW  ])[ksw  ];
      f_NW   = (D.f[dirSE  ])[kse  ];
      f_SE   = (D.f[dirNW  ])[knw  ];
      f_BW   = (D.f[dirTE  ])[kte  ];
      f_TE   = (D.f[dirBW  ])[kbw  ];
      f_TW   = (D.f[dirBE  ])[kbe  ];
      f_BE   = (D.f[dirTW  ])[ktw  ];
      f_BS   = (D.f[dirTN  ])[ktn  ];
      f_TN   = (D.f[dirBS  ])[kbs  ];
      f_TS   = (D.f[dirBN  ])[kbn  ];
      f_BN   = (D.f[dirTS  ])[kts  ];
      f_BSW  = (D.f[dirTNE ])[ktne ];
      f_BNE  = (D.f[dirTSW ])[ktsw ];
      f_BNW  = (D.f[dirTSE ])[ktse ];
      f_BSE  = (D.f[dirTNW ])[ktnw ];
      f_TSW  = (D.f[dirBNE ])[kbne ];
      f_TNE  = (D.f[dirBSW ])[kbsw ];
      f_TNW  = (D.f[dirBSE ])[kbse ];
      f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W); 


      vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S); 

      vx3    =    ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                  (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                  (f_T - f_B); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

	  //b�ser lecktest
	  //q = q_dirE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirW])[kw]=999.f;
   //   }

   //   q = q_dirW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirE])[ke]=999.f;
   //   }

   //   q = q_dirN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirS])[ks]=999.f;
   //   }

   //   q = q_dirS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirN])[kn]=999.f;
   //   }

   //   q = q_dirT[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirB])[kb]=999.f;
   //   }

   //   q = q_dirB[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirT])[kt]=999.f;
   //   }

   //   q = q_dirNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirSW])[ksw]=999.f;
   //   }

   //   q = q_dirSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirNE])[kne]=999.f;
   //   }

   //   q = q_dirSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirNW])[knw]=999.f;
   //   }

   //   q = q_dirNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirSE])[kse]=999.f;
   //   }

   //   q = q_dirTE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBW])[kbw]=999.f;
   //   }

   //   q = q_dirBW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTE])[kte]=999.f;
   //   }

   //   q = q_dirBE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTW])[ktw]=999.f;
   //   }

   //   q = q_dirTW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBE])[kbe]=999.f;
   //   }

   //   q = q_dirTN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBS])[kbs]=999.f;
   //   }

   //   q = q_dirBS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTN])[ktn]=999.f;
   //   }

   //   q = q_dirBN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTS])[kts]=999.f;
   //   }

   //   q = q_dirTS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBN])[kbn]=999.f;
   //   }

   //   q = q_dirTNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBSW])[kbsw]=999.f;
   //   }

   //   q = q_dirBSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTNE])[ktne]=999.f;
   //   }

   //   q = q_dirBNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTSW])[ktsw]=999.f;
   //   }

   //   q = q_dirTSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBNE])[kbne]=999.f;
   //   }

   //   q = q_dirTSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBNW])[kbnw]=999.f;
   //   }

   //   q = q_dirBNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTSE])[ktse]=999.f;
   //   }

   //   q = q_dirBSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTNW])[ktnw]=999.f;
   //   }

   //   q = q_dirTNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBSE])[kbse]=999.f;
   //   }

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
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        )-cu_sq); 
         (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W))/(one+q);
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
         (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E))/(one+q);
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
         (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S))/(one+q);
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
         (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N))/(one+q);
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3)-cu_sq); 
         (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B))/(one+q);
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
         (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T))/(one+q);
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW))/(one+q);
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE))/(one+q);
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW))/(one+q);
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE))/(one+q);
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW))/(one+q);
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE))/(one+q);
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW))/(one+q);
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE))/(one+q);
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS))/(one+q);
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN))/(one+q);
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS))/(one+q);
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN))/(one+q);
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW))/(one+q);
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE))/(one+q);
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW))/(one+q);
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE))/(one+q);
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW))/(one+q);
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE))/(one+q);
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW))/(one+q);
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE))/(one+q);
         //(D.f[dirBSE])[kbse]=zero;
      }

	 // q = q_dirE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq); 
  //       (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-feq*om1)/(one-om1)+q/(one+q)*(f_E+f_W);
		//// (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);
  //    }

  //    q = q_dirW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
  //       (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-feq*om1)/(one-om1)+q/(one+q)*(f_W+f_E);
  //    }

  //    q = q_dirN[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
  //       (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-feq*om1)/(one-om1)+q/(one+q)*(f_N+f_S);
  //    }

  //    q = q_dirS[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
  //       (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-feq*om1)/(one-om1)+q/(one+q)*(f_S+f_N);
  //    }

  //    q = q_dirT[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq); 
  //       (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-feq*om1)/(one-om1)+q/(one+q)*(f_T+f_B);
  //    }

  //    q = q_dirB[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
  //       (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-feq*om1)/(one-om1)+q/(one+q)*(f_B+f_T);
  //    }

  //    q = q_dirNE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
  //       (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-feq*om1)/(one-om1)+q/(one+q)*(f_NE+f_SW);
  //    }

  //    q = q_dirSW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
  //       (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-feq*om1)/(one-om1)+q/(one+q)*(f_SW+f_NE);
  //    }

  //    q = q_dirSE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
  //       (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-feq*om1)/(one-om1)+q/(one+q)*(f_SE+f_NW);
  //    }

  //    q = q_dirNW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
  //       (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-feq*om1)/(one-om1)+q/(one+q)*(f_NW+f_SE);
  //    }

  //    q = q_dirTE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
  //       (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-feq*om1)/(one-om1)+q/(one+q)*(f_TE+f_BW);
  //    }

  //    q = q_dirBW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
  //       (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-feq*om1)/(one-om1)+q/(one+q)*(f_BW+f_TE);
  //    }

  //    q = q_dirBE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
  //       (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-feq*om1)/(one-om1)+q/(one+q)*(f_BE+f_TW);
  //    }

  //    q = q_dirTW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
  //       (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-feq*om1)/(one-om1)+q/(one+q)*(f_TW+f_BE);
  //    }

  //    q = q_dirTN[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
  //       (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-feq*om1)/(one-om1)+q/(one+q)*(f_TN+f_BS);
  //    }

  //    q = q_dirBS[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
  //       (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-feq*om1)/(one-om1)+q/(one+q)*(f_BS+f_TN);
  //    }

  //    q = q_dirBN[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
  //       (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-feq*om1)/(one-om1)+q/(one+q)*(f_BN+f_TS);
  //    }

  //    q = q_dirTS[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
  //       (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-feq*om1)/(one-om1)+q/(one+q)*(f_TS+f_BN);
  //    }

  //    q = q_dirTNE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
  //       (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-feq*om1)/(one-om1)+q/(one+q)*(f_TNE+f_BSW);
  //    }

  //    q = q_dirBSW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
  //       (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-feq*om1)/(one-om1)+q/(one+q)*(f_BSW+f_TNE);
  //    }

  //    q = q_dirBNE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
  //       (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-feq*om1)/(one-om1)+q/(one+q)*(f_BNE+f_TSW);
  //    }

  //    q = q_dirTSW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
  //       (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-feq*om1)/(one-om1)+q/(one+q)*(f_TSW+f_BNE);
  //    }

  //    q = q_dirTSE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
  //       (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-feq*om1)/(one-om1)+q/(one+q)*(f_TSE+f_BNW);
  //    }

  //    q = q_dirBNW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
  //       (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-feq*om1)/(one-om1)+q/(one+q)*(f_BNW+f_TSE);
  //    }

  //    q = q_dirBSE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
  //       (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-feq*om1)/(one-om1)+q/(one+q)*(f_BSE+f_TNW);
  //    }

  //    q = q_dirTNW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
  //       (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-feq*om1)/(one-om1)+q/(one+q)*(f_TNW+f_BSE);
  //    }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void BBDevice27(int inx,
                                     int iny,
                                     doubflo* DD, 
                                     int* k_Q, 
                                     doubflo* QQ,
                                     unsigned int sizeQ,
                                     int kQ, 
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

   if(k<kQ)
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
      //unsigned int kzero= KQK;
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
      //unsigned int nxny = nx*ny;
      //unsigned int kzero= KQK;
      //unsigned int ke   = KQK;
      //unsigned int kw   = KQK + 1;
      //unsigned int kn   = KQK;
      //unsigned int ks   = KQK + nx;
      //unsigned int kt   = KQK;
      //unsigned int kb   = KQK + nxny;
      //unsigned int ksw  = KQK + nx + 1;
      //unsigned int kne  = KQK;
      //unsigned int kse  = KQK + nx;
      //unsigned int knw  = KQK + 1;
      //unsigned int kbw  = KQK + nxny + 1;
      //unsigned int kte  = KQK;
      //unsigned int kbe  = KQK + nxny;
      //unsigned int ktw  = KQK + 1;
      //unsigned int kbs  = KQK + nxny + nx;
      //unsigned int ktn  = KQK;
      //unsigned int kbn  = KQK + nxny;
      //unsigned int kts  = KQK + nx;
      //unsigned int ktse = KQK + nx;
      //unsigned int kbnw = KQK + nxny + 1;
      //unsigned int ktnw = KQK + 1;
      //unsigned int kbse = KQK + nxny + nx;
      //unsigned int ktsw = KQK + nx + 1;
      //unsigned int kbne = KQK + nxny;
      //unsigned int ktne = KQK;
      //unsigned int kbsw = KQK + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
     
      ////////////////////////////////////////////////////////////////////////////////
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[dirE   ])[ke   ];
      f_E    = (D.f[dirW   ])[kw   ];
      f_S    = (D.f[dirN   ])[kn   ];
      f_N    = (D.f[dirS   ])[ks   ];
      f_B    = (D.f[dirT   ])[kt   ];
      f_T    = (D.f[dirB   ])[kb   ];
      f_SW   = (D.f[dirNE  ])[kne  ];
      f_NE   = (D.f[dirSW  ])[ksw  ];
      f_NW   = (D.f[dirSE  ])[kse  ];
      f_SE   = (D.f[dirNW  ])[knw  ];
      f_BW   = (D.f[dirTE  ])[kte  ];
      f_TE   = (D.f[dirBW  ])[kbw  ];
      f_TW   = (D.f[dirBE  ])[kbe  ];
      f_BE   = (D.f[dirTW  ])[ktw  ];
      f_BS   = (D.f[dirTN  ])[ktn  ];
      f_TN   = (D.f[dirBS  ])[kbs  ];
      f_TS   = (D.f[dirBN  ])[kbn  ];
      f_BN   = (D.f[dirTS  ])[kts  ];
      f_BSW  = (D.f[dirTNE ])[ktne ];
      f_BNE  = (D.f[dirTSW ])[ktsw ];
      f_BNW  = (D.f[dirTSE ])[ktse ];
      f_BSE  = (D.f[dirTNW ])[ktnw ];
      f_TSW  = (D.f[dirBNE ])[kbne ];
      f_TNE  = (D.f[dirBSW ])[kbsw ];
      f_TNW  = (D.f[dirBSE ])[kbse ];
      f_TSE  = (D.f[dirBNW ])[kbnw ];

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
      doubflo q;
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirW])[kw]=f_E;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirE])[ke]=f_W;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirS])[ks]=f_N;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirN])[kn]=f_S;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirB])[kb]=f_T;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirT])[kt]=f_B;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirSW])[ksw]=f_NE;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirNE])[kne]=f_SW;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirNW])[knw]=f_SE;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirSE])[kse]=f_NW;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBW])[kbw]=f_TE;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTE])[kte]=f_BW;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTW])[ktw]=f_BE;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBE])[kbe]=f_TW;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBS])[kbs]=f_TN;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTN])[ktn]=f_BS;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTS])[kts]=f_BN;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBN])[kbn]=f_TS;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBSW])[kbsw]=f_TNE;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTNE])[ktne]=f_BSW;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTSW])[ktsw]=f_BNE;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBNE])[kbne]=f_TSW;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBNW])[kbnw]=f_TSE;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTSE])[ktse]=f_BNW;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTNW])[ktnw]=f_BSE;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBSE])[kbse]=f_TNW;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

