/* Device code */
#include "LBM/D3Q27.h"
#include "GPU/constant.h"

//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADPress7(  int inx,
                                       int iny,
                                       doubflo* DD, 
                                       doubflo* DD7, 
                                       doubflo* temp,
                                       doubflo* velo,
                                       doubflo diffusivity,
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

   Distributions7 D7;
   if (evenOrOdd==true)
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[1] = &DD7[1*size_Mat];
      D7.f[2] = &DD7[2*size_Mat];
      D7.f[3] = &DD7[3*size_Mat];
      D7.f[4] = &DD7[4*size_Mat];
      D7.f[5] = &DD7[5*size_Mat];
      D7.f[6] = &DD7[6*size_Mat];
   }
   else
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[2] = &DD7[1*size_Mat];
      D7.f[1] = &DD7[2*size_Mat];
      D7.f[4] = &DD7[3*size_Mat];
      D7.f[3] = &DD7[4*size_Mat];
      D7.f[6] = &DD7[5*size_Mat];
      D7.f[5] = &DD7[6*size_Mat];
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
      //////////////////////////////////////////////////////////////////////////////////
      doubflo  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB;//, 
      //         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
      //         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
      //         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
      //         *q_dirBSE, *q_dirBNW;

      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      //q_dirNE  = &QQ[dirNE  *sizeQ];
      //q_dirSW  = &QQ[dirSW  *sizeQ];
      //q_dirSE  = &QQ[dirSE  *sizeQ];
      //q_dirNW  = &QQ[dirNW  *sizeQ];
      //q_dirTE  = &QQ[dirTE  *sizeQ];
      //q_dirBW  = &QQ[dirBW  *sizeQ];
      //q_dirBE  = &QQ[dirBE  *sizeQ];
      //q_dirTW  = &QQ[dirTW  *sizeQ];
      //q_dirTN  = &QQ[dirTN  *sizeQ];
      //q_dirBS  = &QQ[dirBS  *sizeQ];
      //q_dirBN  = &QQ[dirBN  *sizeQ];
      //q_dirTS  = &QQ[dirTS  *sizeQ];
      //q_dirTNE = &QQ[dirTNE *sizeQ];
      //q_dirTSW = &QQ[dirTSW *sizeQ];
      //q_dirTSE = &QQ[dirTSE *sizeQ];
      //q_dirTNW = &QQ[dirTNW *sizeQ];
      //q_dirBNE = &QQ[dirBNE *sizeQ];
      //q_dirBSW = &QQ[dirBSW *sizeQ];
      //q_dirBSE = &QQ[dirBSE *sizeQ];
      //q_dirBNW = &QQ[dirBNW *sizeQ];
      //////////////////////////////////////////////////////////////////////////////////
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
      /*doubflo drho*/;
      //doubflo vx1_Inflow   = zero;
      //doubflo vx2_Inflow   = zero;
      //doubflo vx3_Inflow   = velo[k];
      //doubflo ux_sq_Inflow = vx1_Inflow * vx1_Inflow;
      //doubflo uy_sq_Inflow = vx2_Inflow * vx2_Inflow;
      //doubflo uz_sq_Inflow = vx3_Inflow * vx3_Inflow;


      //drho   =    f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
      //            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
      //            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      //doubflo vx1 =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //               ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //               (f_E - f_W); 

      //doubflo vx2 =  (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //               ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //               (f_N - f_S); 

      //doubflo vx3 =  ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //               (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //               (f_T - f_B); 

      doubflo rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ ((D.f[dirZERO])[kzero]);
      doubflo rho    =  rho0 + one;
      doubflo OORho  =  one/rho;
      doubflo vx1    =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      doubflo vx2    =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      doubflo vx3    =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));

	  //doubflo cu_sq       =1.5*(vx1*vx1+vx2*vx2+vx3*vx3);
      doubflo ux_sq       = vx1 * vx1;
      doubflo uy_sq       = vx2 * vx2;
      doubflo uz_sq       = vx3 * vx3;
      doubflo omegaD     = three - sqrt(three);
      doubflo Lam         = -(c1o2-one/omegaD);
      doubflo nue_d       = Lam/three;
      //doubflo ae          = zero;
      doubflo ae          = diffusivity/nue_d - one;

      doubflo f7_ZERO,f7_E,f7_W,f7_N,f7_S,f7_T,f7_B;
      doubflo /*feq7_ZERO,*/feq7_E,feq7_W,feq7_N,feq7_S,feq7_T,feq7_B;
      doubflo /*feqW7_ZERO,*/feqW7_E,feqW7_W,feqW7_N,feqW7_S,feqW7_T,feqW7_B;
      //doubflo TempD = temp[k];


      f7_ZERO =  (D7.f[0])[kzero];
      f7_W    =  (D7.f[1])[ke   ];
      f7_E    =  (D7.f[2])[kw   ];
      f7_S    =  (D7.f[3])[kn   ];
      f7_N    =  (D7.f[4])[ks   ];
      f7_B    =  (D7.f[5])[kt   ];
      f7_T    =  (D7.f[6])[kb   ];

      doubflo ConcD = f7_ZERO + f7_E + f7_W + f7_N + f7_S + f7_T + f7_B;

      //feq7_ZERO = ConcD*(c1o3*(ae*(-3.0))-(ux_sq+uy_sq+uz_sq));
      feq7_E    = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
      feq7_W    = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
      feq7_N    = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
      feq7_S    = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
      feq7_T    = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
      feq7_B    = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);

      //feq7_ZERO = TempD*(c1o3*(ae*(-3.0f))-(ux_sq+uy_sq+uz_sq));
      feqW7_E    = feq7_E;// TempD*(c1o6*(ae+one)+c1o2*(ux_sq_Inflow)+vx1_Inflow*c1o2);
      feqW7_W    = feq7_W;// TempD*(c1o6*(ae+one)+c1o2*(ux_sq_Inflow)-vx1_Inflow*c1o2);
      feqW7_N    = feq7_N;// TempD*(c1o6*(ae+one)+c1o2*(uy_sq_Inflow)+vx2_Inflow*c1o2);
      feqW7_S    = feq7_S;// TempD*(c1o6*(ae+one)+c1o2*(uy_sq_Inflow)-vx2_Inflow*c1o2);
      feqW7_T    = feq7_T;// TempD*(c1o6*(ae+one)+c1o2*(uz_sq_Inflow)+vx3_Inflow*c1o2);
      feqW7_B    = feq7_B;// TempD*(c1o6*(ae+one)+c1o2*(uz_sq_Inflow)-vx3_Inflow*c1o2);

      //////////////////////////////////////////////////////////////////////////
      //pointertausch
      if (evenOrOdd==false)
      {
         D7.f[0] = &DD7[0*size_Mat];
         D7.f[1] = &DD7[1*size_Mat];
         D7.f[2] = &DD7[2*size_Mat];
         D7.f[3] = &DD7[3*size_Mat];
         D7.f[4] = &DD7[4*size_Mat];
         D7.f[5] = &DD7[5*size_Mat];
         D7.f[6] = &DD7[6*size_Mat];
      }
      else
      {
         D7.f[0] = &DD7[0*size_Mat];
         D7.f[2] = &DD7[1*size_Mat];
         D7.f[1] = &DD7[2*size_Mat];
         D7.f[4] = &DD7[3*size_Mat];
         D7.f[3] = &DD7[4*size_Mat];
         D7.f[6] = &DD7[5*size_Mat];
         D7.f[5] = &DD7[6*size_Mat];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=0.1f;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //(D7.f[1])[ke   ] = f7_E - feq7_E + feqW7_W; //E
      //(D7.f[2])[kw   ] = f7_W - feq7_W + feqW7_E; //W
      //(D7.f[3])[kn   ] = f7_N - feq7_N + feqW7_S; //N
      //(D7.f[4])[ks   ] = f7_S - feq7_S + feqW7_N; //S
      //(D7.f[5])[kt   ] = f7_T - feq7_T + feqW7_B; //T
      //(D7.f[6])[kb   ] = f7_B - feq7_B + feqW7_T; //B

      //////////////////////////////////////////////////////////////////////////
      //mit Q's
      doubflo /*feq,*/ q;
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         //q=0.;
         (D7.f[2])[kw]=(two*feqW7_W-(f7_E*(q*omegaD-one)-omegaD*feq7_E*(q-one))/(omegaD-one)+f7_W*q)/(q+one);//f7_W - feq7_W + feqW7_E;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         //q=0.;
         (D7.f[1])[ke]=(two*feqW7_E-(f7_W*(q*omegaD-one)-omegaD*feq7_W*(q-one))/(omegaD-one)+f7_E*q)/(q+one);//f7_E - feq7_E + feqW7_W;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         //q=0.;
         (D7.f[4])[ks]=(two*feqW7_S-(f7_N*(q*omegaD-1.)-omegaD*feq7_N*(q-one))/(omegaD-one)+f7_S*q)/(q+one);//f7_S - feq7_S + feqW7_N;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         //q=0.;
         (D7.f[3])[kn]=(two*feqW7_N-(f7_S*(q*omegaD-one)-omegaD*feq7_S*(q-one))/(omegaD-one)+f7_N*q)/(q+one);//f7_N - feq7_N + feqW7_S;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         //q=0.;
         (D7.f[6])[kb]=(two*feqW7_B-(f7_T*(q*omegaD-1.)-omegaD*feq7_T*(q-one))/(omegaD-one)+f7_B*q)/(q+one);//f7_B - feq7_B + feqW7_T;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         //q=0.;
         (D7.f[5])[kt]=(two*feqW7_T-(f7_B*(q*omegaD-one)-omegaD*feq7_B*(q-one))/(omegaD-one)+f7_T*q)/(q+one);//f7_T - feq7_T + feqW7_B;
      }

      ////////////////////////////////////////////////////////////////////////////
      ////ohne Q's
      //doubflo /*feq,*/ q;
      //q = q_dirE[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[2])[kw]=f7_W - feq7_W + feqW7_E;
      //}

      //q = q_dirW[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[1])[ke]=f7_E - feq7_E + feqW7_W;
      //}

      //q = q_dirN[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[4])[ks]=f7_S - feq7_S + feqW7_N;
      //}

      //q = q_dirS[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[3])[kn]=f7_N - feq7_N + feqW7_S;
      //}

      //q = q_dirT[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[6])[kb]=f7_B - feq7_B + feqW7_T;
      //}

      //q = q_dirB[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[5])[kt]=f7_T - feq7_T + feqW7_B;
      //}
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADPress27( int inx,
                                       int iny,
                                       doubflo* DD, 
                                       doubflo* DD27, 
                                       doubflo* temp,
                                       doubflo* velo,
                                       doubflo diffusivity,
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

   Distributions27 D27;
   if (evenOrOdd==true)
   {
      D27.f[dirE   ] = &DD27[dirE   *size_Mat];
      D27.f[dirW   ] = &DD27[dirW   *size_Mat];
      D27.f[dirN   ] = &DD27[dirN   *size_Mat];
      D27.f[dirS   ] = &DD27[dirS   *size_Mat];
      D27.f[dirT   ] = &DD27[dirT   *size_Mat];
      D27.f[dirB   ] = &DD27[dirB   *size_Mat];
      D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
      D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
      D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
      D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
      D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
      D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
   } 
   else
   {
      D27.f[dirW   ] = &DD27[dirE   *size_Mat];
      D27.f[dirE   ] = &DD27[dirW   *size_Mat];
      D27.f[dirS   ] = &DD27[dirN   *size_Mat];
      D27.f[dirN   ] = &DD27[dirS   *size_Mat];
      D27.f[dirB   ] = &DD27[dirT   *size_Mat];
      D27.f[dirT   ] = &DD27[dirB   *size_Mat];
      D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
      D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
      D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
      D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
      D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
      D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
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
      doubflo  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
      doubflo f_ZERO = (D.f[dirZERO])[kzero];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, /*drho, feq,*/ q;
      //drho   = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
      //         f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
      //         f_T + f_B + f_N + f_S + f_E + f_W + f_ZERO; 

      //vx1    = ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //         (f_E - f_W); 


      //vx2    = (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //         (f_N - f_S); 

      //vx3    = ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //         (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //         (f_T - f_B); 
      doubflo rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ f_ZERO;
      doubflo rho    =  rho0 + one;
      doubflo OORho  =  one/rho;
      vx1            =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      vx2            =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      vx3            =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      doubflo f27_W    = (D27.f[dirE   ])[ke   ];
      doubflo f27_E    = (D27.f[dirW   ])[kw   ];
      doubflo f27_S    = (D27.f[dirN   ])[kn   ];
      doubflo f27_N    = (D27.f[dirS   ])[ks   ];
      doubflo f27_B    = (D27.f[dirT   ])[kt   ];
      doubflo f27_T    = (D27.f[dirB   ])[kb   ];
      doubflo f27_SW   = (D27.f[dirNE  ])[kne  ];
      doubflo f27_NE   = (D27.f[dirSW  ])[ksw  ];
      doubflo f27_NW   = (D27.f[dirSE  ])[kse  ];
      doubflo f27_SE   = (D27.f[dirNW  ])[knw  ];
      doubflo f27_BW   = (D27.f[dirTE  ])[kte  ];
      doubflo f27_TE   = (D27.f[dirBW  ])[kbw  ];
      doubflo f27_TW   = (D27.f[dirBE  ])[kbe  ];
      doubflo f27_BE   = (D27.f[dirTW  ])[ktw  ];
      doubflo f27_BS   = (D27.f[dirTN  ])[ktn  ];
      doubflo f27_TN   = (D27.f[dirBS  ])[kbs  ];
      doubflo f27_TS   = (D27.f[dirBN  ])[kbn  ];
      doubflo f27_BN   = (D27.f[dirTS  ])[kts  ];
      doubflo f27_ZERO = (D27.f[dirZERO])[kzero];
      doubflo f27_BSW  = (D27.f[dirTNE ])[ktne ];
      doubflo f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      doubflo f27_BNW  = (D27.f[dirTSE ])[ktse ];
      doubflo f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      doubflo f27_TSW  = (D27.f[dirBNE ])[kbne ];
      doubflo f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      doubflo f27_TNW  = (D27.f[dirBSE ])[kbse ];
      doubflo f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      doubflo ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
                        f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
                        f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      ////doubflo feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      //doubflo feq27_E    =   c2over27* (ConcD+(one+ConcD)*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq));
      //doubflo feq27_W    =   c2over27* (ConcD+(one+ConcD)*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq));
      //doubflo feq27_N    =   c2over27* (ConcD+(one+ConcD)*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq));
      //doubflo feq27_S    =   c2over27* (ConcD+(one+ConcD)*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq));
      //doubflo feq27_T    =   c2over27* (ConcD+(one+ConcD)*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq));
      //doubflo feq27_B    =   c2over27* (ConcD+(one+ConcD)*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq));
      //doubflo feq27_NE   =   c1over54* (ConcD+(one+ConcD)*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      //doubflo feq27_SW   =   c1over54* (ConcD+(one+ConcD)*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      //doubflo feq27_SE   =   c1over54* (ConcD+(one+ConcD)*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      //doubflo feq27_NW   =   c1over54* (ConcD+(one+ConcD)*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      //doubflo feq27_TE   =   c1over54* (ConcD+(one+ConcD)*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      //doubflo feq27_BW   =   c1over54* (ConcD+(one+ConcD)*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      //doubflo feq27_BE   =   c1over54* (ConcD+(one+ConcD)*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      //doubflo feq27_TW   =   c1over54* (ConcD+(one+ConcD)*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      //doubflo feq27_TN   =   c1over54* (ConcD+(one+ConcD)*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      //doubflo feq27_BS   =   c1over54* (ConcD+(one+ConcD)*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      //doubflo feq27_BN   =   c1over54* (ConcD+(one+ConcD)*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      //doubflo feq27_TS   =   c1over54* (ConcD+(one+ConcD)*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      //doubflo feq27_TNE  =   c1over216*(ConcD+(one+ConcD)*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      //doubflo feq27_BSW  =   c1over216*(ConcD+(one+ConcD)*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      //doubflo feq27_BNE  =   c1over216*(ConcD+(one+ConcD)*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      //doubflo feq27_TSW  =   c1over216*(ConcD+(one+ConcD)*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      //doubflo feq27_TSE  =   c1over216*(ConcD+(one+ConcD)*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      //doubflo feq27_BNW  =   c1over216*(ConcD+(one+ConcD)*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      //doubflo feq27_BSE  =   c1over216*(ConcD+(one+ConcD)*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      //doubflo feq27_TNW  =   c1over216*(ConcD+(one+ConcD)*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));
      //
      doubflo feq27_E    =   c2over27* ConcD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feq27_W    =   c2over27* ConcD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feq27_N    =   c2over27* ConcD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feq27_S    =   c2over27* ConcD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feq27_T    =   c2over27* ConcD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feq27_B    =   c2over27* ConcD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feq27_NE   =   c1over54* ConcD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feq27_SW   =   c1over54* ConcD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feq27_SE   =   c1over54* ConcD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feq27_NW   =   c1over54* ConcD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feq27_TE   =   c1over54* ConcD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feq27_BW   =   c1over54* ConcD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feq27_BE   =   c1over54* ConcD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feq27_TW   =   c1over54* ConcD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feq27_TN   =   c1over54* ConcD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feq27_BS   =   c1over54* ConcD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feq27_BN   =   c1over54* ConcD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feq27_TS   =   c1over54* ConcD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feq27_TNE  =   c1over216*ConcD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feq27_BSW  =   c1over216*ConcD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feq27_BNE  =   c1over216*ConcD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feq27_TSW  =   c1over216*ConcD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feq27_TSE  =   c1over216*ConcD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feq27_BNW  =   c1over216*ConcD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feq27_BSE  =   c1over216*ConcD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feq27_TNW  =   c1over216*ConcD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //doubflo TempD = temp[k];
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // velocity inflow
      //vx1   = zero;
      //vx2   = zero;
      //vx3   = velo[k];

      //doubflo feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      doubflo feqW27_E    =  feq27_E  ;// c2over27* TempD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feqW27_W    =  feq27_W  ;// c2over27* TempD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feqW27_N    =  feq27_N  ;// c2over27* TempD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feqW27_S    =  feq27_S  ;// c2over27* TempD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feqW27_T    =  feq27_T  ;// c2over27* TempD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feqW27_B    =  feq27_B  ;// c2over27* TempD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feqW27_NE   =  feq27_NE ;// c1over54* TempD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feqW27_SW   =  feq27_SW ;// c1over54* TempD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feqW27_SE   =  feq27_SE ;// c1over54* TempD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feqW27_NW   =  feq27_NW ;// c1over54* TempD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feqW27_TE   =  feq27_TE ;// c1over54* TempD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feqW27_BW   =  feq27_BW ;// c1over54* TempD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feqW27_BE   =  feq27_BE ;// c1over54* TempD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feqW27_TW   =  feq27_TW ;// c1over54* TempD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feqW27_TN   =  feq27_TN ;// c1over54* TempD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feqW27_BS   =  feq27_BS ;// c1over54* TempD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feqW27_BN   =  feq27_BN ;// c1over54* TempD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feqW27_TS   =  feq27_TS ;// c1over54* TempD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feqW27_TNE  =  feq27_TNE;// c1over216*TempD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feqW27_BSW  =  feq27_BSW;// c1over216*TempD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feqW27_BNE  =  feq27_BNE;// c1over216*TempD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feqW27_TSW  =  feq27_TSW;// c1over216*TempD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feqW27_TSE  =  feq27_TSE;// c1over216*TempD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feqW27_BNW  =  feq27_BNW;// c1over216*TempD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feqW27_BSE  =  feq27_BSE;// c1over216*TempD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feqW27_TNW  =  feq27_TNW;// c1over216*TempD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo omegaD     = three - sqrt(three);
      //doubflo Lam        = -(c1o2-one/omegaD);
      //doubflo nue_d      = Lam/three;
      //doubflo ae         = zero;
      //doubflo ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
      {
         D27.f[dirE   ] = &DD27[dirE   *size_Mat];
         D27.f[dirW   ] = &DD27[dirW   *size_Mat];
         D27.f[dirN   ] = &DD27[dirN   *size_Mat];
         D27.f[dirS   ] = &DD27[dirS   *size_Mat];
         D27.f[dirT   ] = &DD27[dirT   *size_Mat];
         D27.f[dirB   ] = &DD27[dirB   *size_Mat];
         D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
         D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
         D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
         D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
         D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
         D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
      } 
      else
      {
         D27.f[dirW   ] = &DD27[dirE   *size_Mat];
         D27.f[dirE   ] = &DD27[dirW   *size_Mat];
         D27.f[dirS   ] = &DD27[dirN   *size_Mat];
         D27.f[dirN   ] = &DD27[dirS   *size_Mat];
         D27.f[dirB   ] = &DD27[dirT   *size_Mat];
         D27.f[dirT   ] = &DD27[dirB   *size_Mat];
         D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
         D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
         D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
         D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
         D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
         D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[dirW  ])[kw  ]=(two*feqW27_W  -(f27_E  *(q*omegaD-one)-omegaD*feq27_E  *(q-one))/(omegaD-one)+f27_W  *q)/(q+one);
      q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[dirE  ])[ke  ]=(two*feqW27_E  -(f27_W  *(q*omegaD-one)-omegaD*feq27_W  *(q-one))/(omegaD-one)+f27_E  *q)/(q+one);
      q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[dirS  ])[ks  ]=(two*feqW27_S  -(f27_N  *(q*omegaD-one)-omegaD*feq27_N  *(q-one))/(omegaD-one)+f27_S  *q)/(q+one);
      q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[dirN  ])[kn  ]=(two*feqW27_N  -(f27_S  *(q*omegaD-one)-omegaD*feq27_S  *(q-one))/(omegaD-one)+f27_N  *q)/(q+one);
      q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[dirB  ])[kb  ]=(two*feqW27_B  -(f27_T  *(q*omegaD-one)-omegaD*feq27_T  *(q-one))/(omegaD-one)+f27_B  *q)/(q+one);
      q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[dirT  ])[kt  ]=(two*feqW27_T  -(f27_B  *(q*omegaD-one)-omegaD*feq27_B  *(q-one))/(omegaD-one)+f27_T  *q)/(q+one);
      q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[dirSW ])[ksw ]=(two*feqW27_SW -(f27_NE *(q*omegaD-one)-omegaD*feq27_NE *(q-one))/(omegaD-one)+f27_SW *q)/(q+one);
      q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[dirNE ])[kne ]=(two*feqW27_NE -(f27_SW *(q*omegaD-one)-omegaD*feq27_SW *(q-one))/(omegaD-one)+f27_NE *q)/(q+one);
      q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[dirNW ])[knw ]=(two*feqW27_NW -(f27_SE *(q*omegaD-one)-omegaD*feq27_SE *(q-one))/(omegaD-one)+f27_NW *q)/(q+one);
      q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[dirSE ])[kse ]=(two*feqW27_SE -(f27_NW *(q*omegaD-one)-omegaD*feq27_NW *(q-one))/(omegaD-one)+f27_SE *q)/(q+one);
      q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[dirBW ])[kbw ]=(two*feqW27_BW -(f27_TE *(q*omegaD-one)-omegaD*feq27_TE *(q-one))/(omegaD-one)+f27_BW *q)/(q+one);
      q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[dirTE ])[kte ]=(two*feqW27_TE -(f27_BW *(q*omegaD-one)-omegaD*feq27_BW *(q-one))/(omegaD-one)+f27_TE *q)/(q+one);
      q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[dirTW ])[ktw ]=(two*feqW27_TW -(f27_BE *(q*omegaD-one)-omegaD*feq27_BE *(q-one))/(omegaD-one)+f27_TW *q)/(q+one);
      q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[dirBE ])[kbe ]=(two*feqW27_BE -(f27_TW *(q*omegaD-one)-omegaD*feq27_TW *(q-one))/(omegaD-one)+f27_BE *q)/(q+one);
      q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[dirBS ])[kbs ]=(two*feqW27_BS -(f27_TN *(q*omegaD-one)-omegaD*feq27_TN *(q-one))/(omegaD-one)+f27_BS *q)/(q+one);
      q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[dirTN ])[ktn ]=(two*feqW27_TN -(f27_BS *(q*omegaD-one)-omegaD*feq27_BS *(q-one))/(omegaD-one)+f27_TN *q)/(q+one);
      q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[dirTS ])[kts ]=(two*feqW27_TS -(f27_BN *(q*omegaD-one)-omegaD*feq27_BN *(q-one))/(omegaD-one)+f27_TS *q)/(q+one);
      q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[dirBN ])[kbn ]=(two*feqW27_BN -(f27_TS *(q*omegaD-one)-omegaD*feq27_TS *(q-one))/(omegaD-one)+f27_BN *q)/(q+one);
      q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[dirBSW])[kbsw]=(two*feqW27_BSW-(f27_TNE*(q*omegaD-one)-omegaD*feq27_TNE*(q-one))/(omegaD-one)+f27_BSW*q)/(q+one);
      q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[dirTNE])[ktne]=(two*feqW27_TNE-(f27_BSW*(q*omegaD-one)-omegaD*feq27_BSW*(q-one))/(omegaD-one)+f27_TNE*q)/(q+one);
      q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[dirTSW])[ktsw]=(two*feqW27_TSW-(f27_BNE*(q*omegaD-one)-omegaD*feq27_BNE*(q-one))/(omegaD-one)+f27_TSW*q)/(q+one);
      q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[dirBNE])[kbne]=(two*feqW27_BNE-(f27_TSW*(q*omegaD-one)-omegaD*feq27_TSW*(q-one))/(omegaD-one)+f27_BNE*q)/(q+one);
      q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[dirBNW])[kbnw]=(two*feqW27_BNW-(f27_TSE*(q*omegaD-one)-omegaD*feq27_TSE*(q-one))/(omegaD-one)+f27_BNW*q)/(q+one);
      q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[dirTSE])[ktse]=(two*feqW27_TSE-(f27_BNW*(q*omegaD-one)-omegaD*feq27_BNW*(q-one))/(omegaD-one)+f27_TSE*q)/(q+one);
      q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[dirTNW])[ktnw]=(two*feqW27_TNW-(f27_BSE*(q*omegaD-one)-omegaD*feq27_BSE*(q-one))/(omegaD-one)+f27_TNW*q)/(q+one);
      q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[dirBSE])[kbse]=(two*feqW27_BSE-(f27_TNW*(q*omegaD-one)-omegaD*feq27_TNW*(q-one))/(omegaD-one)+f27_BSE*q)/(q+one);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADVel7( int inx,
                                    int iny,
                                    doubflo* DD, 
                                    doubflo* DD7, 
                                    doubflo* temp,
                                    doubflo* velo,
                                    doubflo diffusivity,
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

   Distributions7 D7;
   if (evenOrOdd==true)
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[1] = &DD7[1*size_Mat];
      D7.f[2] = &DD7[2*size_Mat];
      D7.f[3] = &DD7[3*size_Mat];
      D7.f[4] = &DD7[4*size_Mat];
      D7.f[5] = &DD7[5*size_Mat];
      D7.f[6] = &DD7[6*size_Mat];
   }
   else
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[2] = &DD7[1*size_Mat];
      D7.f[1] = &DD7[2*size_Mat];
      D7.f[4] = &DD7[3*size_Mat];
      D7.f[3] = &DD7[4*size_Mat];
      D7.f[6] = &DD7[5*size_Mat];
      D7.f[5] = &DD7[6*size_Mat];
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
      //////////////////////////////////////////////////////////////////////////////////
      doubflo  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB;//, 
      //         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
      //         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
      //         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
      //         *q_dirBSE, *q_dirBNW;

      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      //q_dirNE  = &QQ[dirNE  *sizeQ];
      //q_dirSW  = &QQ[dirSW  *sizeQ];
      //q_dirSE  = &QQ[dirSE  *sizeQ];
      //q_dirNW  = &QQ[dirNW  *sizeQ];
      //q_dirTE  = &QQ[dirTE  *sizeQ];
      //q_dirBW  = &QQ[dirBW  *sizeQ];
      //q_dirBE  = &QQ[dirBE  *sizeQ];
      //q_dirTW  = &QQ[dirTW  *sizeQ];
      //q_dirTN  = &QQ[dirTN  *sizeQ];
      //q_dirBS  = &QQ[dirBS  *sizeQ];
      //q_dirBN  = &QQ[dirBN  *sizeQ];
      //q_dirTS  = &QQ[dirTS  *sizeQ];
      //q_dirTNE = &QQ[dirTNE *sizeQ];
      //q_dirTSW = &QQ[dirTSW *sizeQ];
      //q_dirTSE = &QQ[dirTSE *sizeQ];
      //q_dirTNW = &QQ[dirTNW *sizeQ];
      //q_dirBNE = &QQ[dirBNE *sizeQ];
      //q_dirBSW = &QQ[dirBSW *sizeQ];
      //q_dirBSE = &QQ[dirBSE *sizeQ];
      //q_dirBNW = &QQ[dirBNW *sizeQ];
      //////////////////////////////////////////////////////////////////////////////////
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
      /*doubflo drho*/;
      doubflo vx1_Inflow   = zero;
      doubflo vx2_Inflow   = zero;
      doubflo vx3_Inflow   = velo[k];
      doubflo ux_sq_Inflow = vx1_Inflow * vx1_Inflow;
      doubflo uy_sq_Inflow = vx2_Inflow * vx2_Inflow;
      doubflo uz_sq_Inflow = vx3_Inflow * vx3_Inflow;


      ////drho   =    f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
      ////            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
      ////            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      //doubflo vx1 =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //               ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //               (f_E - f_W); 

      //doubflo vx2 =  (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //               ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //               (f_N - f_S); 

      //doubflo vx3 =  ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //               (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //               (f_T - f_B); 

      doubflo rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ ((D.f[dirZERO])[kzero]);
      doubflo rho    =  rho0 + one;
      doubflo OORho  =  one/rho;
      doubflo vx1    =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      doubflo vx2    =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      doubflo vx3    =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));

	  //doubflo cu_sq       =1.5f*(vx1*vx1+vx2*vx2+vx3*vx3);
      doubflo ux_sq       = vx1 * vx1;
      doubflo uy_sq       = vx2 * vx2;
      doubflo uz_sq       = vx3 * vx3;
      doubflo omegaD     = three - sqrt(three);
      doubflo Lam         = -(c1o2-one/omegaD);
      doubflo nue_d       = Lam/three;
      //doubflo ae          = zero;
      doubflo ae          = diffusivity/nue_d - one;

      doubflo f7_ZERO,f7_E,f7_W,f7_N,f7_S,f7_T,f7_B;
      doubflo /*feq7_ZERO,*/feq7_E,feq7_W,feq7_N,feq7_S,feq7_T,feq7_B;
      doubflo /*feqW7_ZERO,*/feqW7_E,feqW7_W,feqW7_N,feqW7_S,feqW7_T,feqW7_B;
      doubflo TempD = temp[k];


      f7_ZERO =  (D7.f[0])[kzero];
      f7_W    =  (D7.f[1])[ke   ];
      f7_E    =  (D7.f[2])[kw   ];
      f7_S    =  (D7.f[3])[kn   ];
      f7_N    =  (D7.f[4])[ks   ];
      f7_B    =  (D7.f[5])[kt   ];
      f7_T    =  (D7.f[6])[kb   ];

      doubflo ConcD = f7_ZERO + f7_E + f7_W + f7_N + f7_S + f7_T + f7_B;

      //feq7_ZERO = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
      feq7_E    = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
      feq7_W    = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
      feq7_N    = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
      feq7_S    = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
      feq7_T    = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
      feq7_B    = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);

      //feq7_ZERO = TempD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
      feqW7_E    = TempD*(c1o6*(ae+one)+c1o2*(ux_sq_Inflow)+vx1_Inflow*c1o2);
      feqW7_W    = TempD*(c1o6*(ae+one)+c1o2*(ux_sq_Inflow)-vx1_Inflow*c1o2);
      feqW7_N    = TempD*(c1o6*(ae+one)+c1o2*(uy_sq_Inflow)+vx2_Inflow*c1o2);
      feqW7_S    = TempD*(c1o6*(ae+one)+c1o2*(uy_sq_Inflow)-vx2_Inflow*c1o2);
      feqW7_T    = TempD*(c1o6*(ae+one)+c1o2*(uz_sq_Inflow)+vx3_Inflow*c1o2);
      feqW7_B    = TempD*(c1o6*(ae+one)+c1o2*(uz_sq_Inflow)-vx3_Inflow*c1o2);

      //////////////////////////////////////////////////////////////////////////
      //pointertausch
      if (evenOrOdd==false)
      {
         D7.f[0] = &DD7[0*size_Mat];
         D7.f[1] = &DD7[1*size_Mat];
         D7.f[2] = &DD7[2*size_Mat];
         D7.f[3] = &DD7[3*size_Mat];
         D7.f[4] = &DD7[4*size_Mat];
         D7.f[5] = &DD7[5*size_Mat];
         D7.f[6] = &DD7[6*size_Mat];
      }
      else
      {
         D7.f[0] = &DD7[0*size_Mat];
         D7.f[2] = &DD7[1*size_Mat];
         D7.f[1] = &DD7[2*size_Mat];
         D7.f[4] = &DD7[3*size_Mat];
         D7.f[3] = &DD7[4*size_Mat];
         D7.f[6] = &DD7[5*size_Mat];
         D7.f[5] = &DD7[6*size_Mat];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //(D7.f[1])[ke   ] = f7_E - feq7_E + feqW7_W; //E
      //(D7.f[2])[kw   ] = f7_W - feq7_W + feqW7_E; //W
      //(D7.f[3])[kn   ] = f7_N - feq7_N + feqW7_S; //N
      //(D7.f[4])[ks   ] = f7_S - feq7_S + feqW7_N; //S
      //(D7.f[5])[kt   ] = f7_T - feq7_T + feqW7_B; //T
      //(D7.f[6])[kb   ] = f7_B - feq7_B + feqW7_T; //B

      //////////////////////////////////////////////////////////////////////////
      //mit Q's
      doubflo /*feq,*/ q;
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         //q=0.;
         (D7.f[2])[kw]=(two*feqW7_W-(f7_E*(q*omegaD-one)-omegaD*feq7_E*(q-one))/(omegaD-one)+f7_W*q)/(q+one);//f7_W - feq7_W + feqW7_E;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         //q=0.;
         (D7.f[1])[ke]=(two*feqW7_E-(f7_W*(q*omegaD-one)-omegaD*feq7_W*(q-one))/(omegaD-one)+f7_E*q)/(q+one);//f7_E - feq7_E + feqW7_W;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         //q=0.;
         (D7.f[4])[ks]=(two*feqW7_S-(f7_N*(q*omegaD-one)-omegaD*feq7_N*(q-one))/(omegaD-one)+f7_S*q)/(q+one);//f7_S - feq7_S + feqW7_N;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         //q=0.;
         (D7.f[3])[kn]=(two*feqW7_N-(f7_S*(q*omegaD-one)-omegaD*feq7_S*(q-one))/(omegaD-one)+f7_N*q)/(q+one);//f7_N - feq7_N + feqW7_S;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         //q=0.;
         (D7.f[6])[kb]=(two*feqW7_B-(f7_T*(q*omegaD-one)-omegaD*feq7_T*(q-one))/(omegaD-one)+f7_B*q)/(q+one);//f7_B - feq7_B + feqW7_T;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         //q=0.;
         (D7.f[5])[kt]=(two*feqW7_T-(f7_B*(q*omegaD-one)-omegaD*feq7_B*(q-one))/(omegaD-one)+f7_T*q)/(q+one);//f7_T - feq7_T + feqW7_B;
      }

      ////////////////////////////////////////////////////////////////////////////
      ////ohne Q's
      //doubflo /*feq,*/ q;
      //q = q_dirE[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[2])[kw]=f7_W - feq7_W + feqW7_E;
      //}

      //q = q_dirW[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[1])[ke]=f7_E - feq7_E + feqW7_W;
      //}

      //q = q_dirN[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[4])[ks]=f7_S - feq7_S + feqW7_N;
      //}

      //q = q_dirS[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[3])[kn]=f7_N - feq7_N + feqW7_S;
      //}

      //q = q_dirT[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[6])[kb]=f7_B - feq7_B + feqW7_T;
      //}

      //q = q_dirB[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[5])[kt]=f7_T - feq7_T + feqW7_B;
      //}
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADVel27(int inx,
                                    int iny,
                                    doubflo* DD, 
                                    doubflo* DD27, 
                                    doubflo* temp,
                                    doubflo* velo,
                                    doubflo diffusivity,
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

   Distributions27 D27;
   if (evenOrOdd==true)
   {
      D27.f[dirE   ] = &DD27[dirE   *size_Mat];
      D27.f[dirW   ] = &DD27[dirW   *size_Mat];
      D27.f[dirN   ] = &DD27[dirN   *size_Mat];
      D27.f[dirS   ] = &DD27[dirS   *size_Mat];
      D27.f[dirT   ] = &DD27[dirT   *size_Mat];
      D27.f[dirB   ] = &DD27[dirB   *size_Mat];
      D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
      D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
      D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
      D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
      D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
      D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
   } 
   else
   {
      D27.f[dirW   ] = &DD27[dirE   *size_Mat];
      D27.f[dirE   ] = &DD27[dirW   *size_Mat];
      D27.f[dirS   ] = &DD27[dirN   *size_Mat];
      D27.f[dirN   ] = &DD27[dirS   *size_Mat];
      D27.f[dirB   ] = &DD27[dirT   *size_Mat];
      D27.f[dirT   ] = &DD27[dirB   *size_Mat];
      D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
      D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
      D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
      D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
      D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
      D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
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
      doubflo f_ZERO = (D.f[dirZERO])[kzero];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, /*drho, feq,*/ q;
      ////drho   = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
      ////         f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
      ////         f_T + f_B + f_N + f_S + f_E + f_W + f_ZERO; 

      //vx1    = ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //         (f_E - f_W); 


      //vx2    = (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //         (f_N - f_S); 

      //vx3    = ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //         (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //         (f_T - f_B); 
      doubflo rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ f_ZERO;
      doubflo rho    =  rho0 + one;
      doubflo OORho  =  one/rho;
      vx1     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      vx2     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      vx3     =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      doubflo f27_W    = (D27.f[dirE   ])[ke   ];
      doubflo f27_E    = (D27.f[dirW   ])[kw   ];
      doubflo f27_S    = (D27.f[dirN   ])[kn   ];
      doubflo f27_N    = (D27.f[dirS   ])[ks   ];
      doubflo f27_B    = (D27.f[dirT   ])[kt   ];
      doubflo f27_T    = (D27.f[dirB   ])[kb   ];
      doubflo f27_SW   = (D27.f[dirNE  ])[kne  ];
      doubflo f27_NE   = (D27.f[dirSW  ])[ksw  ];
      doubflo f27_NW   = (D27.f[dirSE  ])[kse  ];
      doubflo f27_SE   = (D27.f[dirNW  ])[knw  ];
      doubflo f27_BW   = (D27.f[dirTE  ])[kte  ];
      doubflo f27_TE   = (D27.f[dirBW  ])[kbw  ];
      doubflo f27_TW   = (D27.f[dirBE  ])[kbe  ];
      doubflo f27_BE   = (D27.f[dirTW  ])[ktw  ];
      doubflo f27_BS   = (D27.f[dirTN  ])[ktn  ];
      doubflo f27_TN   = (D27.f[dirBS  ])[kbs  ];
      doubflo f27_TS   = (D27.f[dirBN  ])[kbn  ];
      doubflo f27_BN   = (D27.f[dirTS  ])[kts  ];
      doubflo f27_ZERO = (D27.f[dirZERO])[kzero];
      doubflo f27_BSW  = (D27.f[dirTNE ])[ktne ];
      doubflo f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      doubflo f27_BNW  = (D27.f[dirTSE ])[ktse ];
      doubflo f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      doubflo f27_TSW  = (D27.f[dirBNE ])[kbne ];
      doubflo f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      doubflo f27_TNW  = (D27.f[dirBSE ])[kbse ];
      doubflo f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      doubflo ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
                        f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
                        f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      //doubflo feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      doubflo feq27_E    =   c2over27* ConcD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feq27_W    =   c2over27* ConcD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feq27_N    =   c2over27* ConcD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feq27_S    =   c2over27* ConcD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feq27_T    =   c2over27* ConcD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feq27_B    =   c2over27* ConcD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feq27_NE   =   c1over54* ConcD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feq27_SW   =   c1over54* ConcD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feq27_SE   =   c1over54* ConcD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feq27_NW   =   c1over54* ConcD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feq27_TE   =   c1over54* ConcD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feq27_BW   =   c1over54* ConcD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feq27_BE   =   c1over54* ConcD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feq27_TW   =   c1over54* ConcD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feq27_TN   =   c1over54* ConcD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feq27_BS   =   c1over54* ConcD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feq27_BN   =   c1over54* ConcD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feq27_TS   =   c1over54* ConcD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feq27_TNE  =   c1over216*ConcD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feq27_BSW  =   c1over216*ConcD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feq27_BNE  =   c1over216*ConcD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feq27_TSW  =   c1over216*ConcD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feq27_TSE  =   c1over216*ConcD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feq27_BNW  =   c1over216*ConcD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feq27_BSE  =   c1over216*ConcD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feq27_TNW  =   c1over216*ConcD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo TempD = temp[k];
      //doubflo TempD = four;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // velocity inflow
      //vx1   = velo[k]; //zero;
      //vx2   = zero; //velo[k];//zero;//velo[k];
      //vx3   = zero;

      ////doubflo feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      //doubflo feqW27_E    =   c2over27* (TempD+(one+TempD)*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq));
      //doubflo feqW27_W    =   c2over27* (TempD+(one+TempD)*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq));
      //doubflo feqW27_N    =   c2over27* (TempD+(one+TempD)*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq));
      //doubflo feqW27_S    =   c2over27* (TempD+(one+TempD)*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq));
      //doubflo feqW27_T    =   c2over27* (TempD+(one+TempD)*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq));
      //doubflo feqW27_B    =   c2over27* (TempD+(one+TempD)*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq));
      //doubflo feqW27_NE   =   c1over54* (TempD+(one+TempD)*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      //doubflo feqW27_SW   =   c1over54* (TempD+(one+TempD)*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      //doubflo feqW27_SE   =   c1over54* (TempD+(one+TempD)*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      //doubflo feqW27_NW   =   c1over54* (TempD+(one+TempD)*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      //doubflo feqW27_TE   =   c1over54* (TempD+(one+TempD)*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      //doubflo feqW27_BW   =   c1over54* (TempD+(one+TempD)*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      //doubflo feqW27_BE   =   c1over54* (TempD+(one+TempD)*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      //doubflo feqW27_TW   =   c1over54* (TempD+(one+TempD)*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      //doubflo feqW27_TN   =   c1over54* (TempD+(one+TempD)*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      //doubflo feqW27_BS   =   c1over54* (TempD+(one+TempD)*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      //doubflo feqW27_BN   =   c1over54* (TempD+(one+TempD)*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      //doubflo feqW27_TS   =   c1over54* (TempD+(one+TempD)*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      //doubflo feqW27_TNE  =   c1over216*(TempD+(one+TempD)*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      //doubflo feqW27_BSW  =   c1over216*(TempD+(one+TempD)*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      //doubflo feqW27_BNE  =   c1over216*(TempD+(one+TempD)*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      //doubflo feqW27_TSW  =   c1over216*(TempD+(one+TempD)*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      //doubflo feqW27_TSE  =   c1over216*(TempD+(one+TempD)*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      //doubflo feqW27_BNW  =   c1over216*(TempD+(one+TempD)*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      //doubflo feqW27_BSE  =   c1over216*(TempD+(one+TempD)*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      //doubflo feqW27_TNW  =   c1over216*(TempD+(one+TempD)*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));
      //
      doubflo feqW27_E    =   c2over27* TempD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feqW27_W    =   c2over27* TempD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feqW27_N    =   c2over27* TempD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feqW27_S    =   c2over27* TempD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feqW27_T    =   c2over27* TempD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feqW27_B    =   c2over27* TempD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feqW27_NE   =   c1over54* TempD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feqW27_SW   =   c1over54* TempD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feqW27_SE   =   c1over54* TempD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feqW27_NW   =   c1over54* TempD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feqW27_TE   =   c1over54* TempD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feqW27_BW   =   c1over54* TempD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feqW27_BE   =   c1over54* TempD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feqW27_TW   =   c1over54* TempD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feqW27_TN   =   c1over54* TempD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feqW27_BS   =   c1over54* TempD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feqW27_BN   =   c1over54* TempD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feqW27_TS   =   c1over54* TempD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feqW27_TNE  =   c1over216*TempD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feqW27_BSW  =   c1over216*TempD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feqW27_BNE  =   c1over216*TempD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feqW27_TSW  =   c1over216*TempD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feqW27_TSE  =   c1over216*TempD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feqW27_BNW  =   c1over216*TempD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feqW27_BSE  =   c1over216*TempD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feqW27_TNW  =   c1over216*TempD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo omegaD     = three - sqrt(three);
      //doubflo Lam        = -(c1o2-one/omegaD);
      //doubflo nue_d      = Lam/three;
      //doubflo ae         = zero;
      //doubflo ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
      {
         D27.f[dirE   ] = &DD27[dirE   *size_Mat];
         D27.f[dirW   ] = &DD27[dirW   *size_Mat];
         D27.f[dirN   ] = &DD27[dirN   *size_Mat];
         D27.f[dirS   ] = &DD27[dirS   *size_Mat];
         D27.f[dirT   ] = &DD27[dirT   *size_Mat];
         D27.f[dirB   ] = &DD27[dirB   *size_Mat];
         D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
         D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
         D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
         D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
         D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
         D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
      } 
      else
      {
         D27.f[dirW   ] = &DD27[dirE   *size_Mat];
         D27.f[dirE   ] = &DD27[dirW   *size_Mat];
         D27.f[dirS   ] = &DD27[dirN   *size_Mat];
         D27.f[dirN   ] = &DD27[dirS   *size_Mat];
         D27.f[dirB   ] = &DD27[dirT   *size_Mat];
         D27.f[dirT   ] = &DD27[dirB   *size_Mat];
         D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
         D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
         D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
         D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
         D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
         D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //(D27.f[dirW  ])[kw  ]= four;
      //(D27.f[dirE  ])[ke  ]= four;
      //(D27.f[dirS  ])[ks  ]= four;
      //(D27.f[dirN  ])[kn  ]= four;
      //(D27.f[dirB  ])[kb  ]= four;
      //(D27.f[dirT  ])[kt  ]= four;
      //(D27.f[dirSW ])[ksw ]= four;
      //(D27.f[dirNE ])[kne ]= four;
      //(D27.f[dirNW ])[knw ]= four;
      //(D27.f[dirSE ])[kse ]= four;
      //(D27.f[dirBW ])[kbw ]= four;
      //(D27.f[dirTE ])[kte ]= four;
      //(D27.f[dirTW ])[ktw ]= four;
      //(D27.f[dirBE ])[kbe ]= four;
      //(D27.f[dirBS ])[kbs ]= four;
      //(D27.f[dirTN ])[ktn ]= four;
      //(D27.f[dirTS ])[kts ]= four;
      //(D27.f[dirBN ])[kbn ]= four;
      //(D27.f[dirBSW])[kbsw]= four;
      //(D27.f[dirTNE])[ktne]= four;
      //(D27.f[dirTSW])[ktsw]= four;
      //(D27.f[dirBNE])[kbne]= four;
      //(D27.f[dirBNW])[kbnw]= four;
      //(D27.f[dirTSE])[ktse]= four;
      //(D27.f[dirTNW])[ktnw]= four;
      //(D27.f[dirBSE])[kbse]= four;
      q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[dirW  ])[kw  ]= -feqW27_W  + two * c2over27  * TempD;
      q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[dirE  ])[ke  ]= -feqW27_E  + two * c2over27  * TempD;
      q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[dirS  ])[ks  ]= -feqW27_S  + two * c2over27  * TempD;
      q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[dirN  ])[kn  ]= -feqW27_N  + two * c2over27  * TempD;
      q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[dirB  ])[kb  ]= -feqW27_B  + two * c2over27  * TempD;
      q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[dirT  ])[kt  ]= -feqW27_T  + two * c2over27  * TempD;
      q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[dirSW ])[ksw ]= -feqW27_SW + two * c1over54  * TempD;
      q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[dirNE ])[kne ]= -feqW27_NE + two * c1over54  * TempD;
      q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[dirNW ])[knw ]= -feqW27_NW + two * c1over54  * TempD;
      q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[dirSE ])[kse ]= -feqW27_SE + two * c1over54  * TempD;
      q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[dirBW ])[kbw ]= -feqW27_BW + two * c1over54  * TempD;
      q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[dirTE ])[kte ]= -feqW27_TE + two * c1over54  * TempD;
      q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[dirTW ])[ktw ]= -feqW27_TW + two * c1over54  * TempD;
      q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[dirBE ])[kbe ]= -feqW27_BE + two * c1over54  * TempD;
      q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[dirBS ])[kbs ]= -feqW27_BS + two * c1over54  * TempD;
      q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[dirTN ])[ktn ]= -feqW27_TN + two * c1over54  * TempD;
      q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[dirTS ])[kts ]= -feqW27_TS + two * c1over54  * TempD;
      q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[dirBN ])[kbn ]= -feqW27_BN + two * c1over54  * TempD;
      q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[dirBSW])[kbsw]= -feqW27_BSW+ two * c1over216 * TempD;
      q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[dirTNE])[ktne]= -feqW27_TNE+ two * c1over216 * TempD;
      q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[dirTSW])[ktsw]= -feqW27_TSW+ two * c1over216 * TempD;
      q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[dirBNE])[kbne]= -feqW27_BNE+ two * c1over216 * TempD;
      q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[dirBNW])[kbnw]= -feqW27_BNW+ two * c1over216 * TempD;
      q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[dirTSE])[ktse]= -feqW27_TSE+ two * c1over216 * TempD;
      q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[dirTNW])[ktnw]= -feqW27_TNW+ two * c1over216 * TempD;
      q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[dirBSE])[kbse]= -feqW27_BSE+ two * c1over216 * TempD;
      //q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[dirW  ])[kw  ]=(two*feqW27_W  -(f27_E  *(q*omegaD-one)-omegaD*feq27_E  *(q-one))/(omegaD-one)+f27_W  *q)/(q+one);
      //q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[dirE  ])[ke  ]=(two*feqW27_E  -(f27_W  *(q*omegaD-one)-omegaD*feq27_W  *(q-one))/(omegaD-one)+f27_E  *q)/(q+one);
      //q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[dirS  ])[ks  ]=(two*feqW27_S  -(f27_N  *(q*omegaD-one)-omegaD*feq27_N  *(q-one))/(omegaD-one)+f27_S  *q)/(q+one);
      //q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[dirN  ])[kn  ]=(two*feqW27_N  -(f27_S  *(q*omegaD-one)-omegaD*feq27_S  *(q-one))/(omegaD-one)+f27_N  *q)/(q+one);
      //q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[dirB  ])[kb  ]=(two*feqW27_B  -(f27_T  *(q*omegaD-one)-omegaD*feq27_T  *(q-one))/(omegaD-one)+f27_B  *q)/(q+one);
      //q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[dirT  ])[kt  ]=(two*feqW27_T  -(f27_B  *(q*omegaD-one)-omegaD*feq27_B  *(q-one))/(omegaD-one)+f27_T  *q)/(q+one);
      //q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[dirSW ])[ksw ]=(two*feqW27_SW -(f27_NE *(q*omegaD-one)-omegaD*feq27_NE *(q-one))/(omegaD-one)+f27_SW *q)/(q+one);
      //q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[dirNE ])[kne ]=(two*feqW27_NE -(f27_SW *(q*omegaD-one)-omegaD*feq27_SW *(q-one))/(omegaD-one)+f27_NE *q)/(q+one);
      //q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[dirNW ])[knw ]=(two*feqW27_NW -(f27_SE *(q*omegaD-one)-omegaD*feq27_SE *(q-one))/(omegaD-one)+f27_NW *q)/(q+one);
      //q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[dirSE ])[kse ]=(two*feqW27_SE -(f27_NW *(q*omegaD-one)-omegaD*feq27_NW *(q-one))/(omegaD-one)+f27_SE *q)/(q+one);
      //q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[dirBW ])[kbw ]=(two*feqW27_BW -(f27_TE *(q*omegaD-one)-omegaD*feq27_TE *(q-one))/(omegaD-one)+f27_BW *q)/(q+one);
      //q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[dirTE ])[kte ]=(two*feqW27_TE -(f27_BW *(q*omegaD-one)-omegaD*feq27_BW *(q-one))/(omegaD-one)+f27_TE *q)/(q+one);
      //q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[dirTW ])[ktw ]=(two*feqW27_TW -(f27_BE *(q*omegaD-one)-omegaD*feq27_BE *(q-one))/(omegaD-one)+f27_TW *q)/(q+one);
      //q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[dirBE ])[kbe ]=(two*feqW27_BE -(f27_TW *(q*omegaD-one)-omegaD*feq27_TW *(q-one))/(omegaD-one)+f27_BE *q)/(q+one);
      //q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[dirBS ])[kbs ]=(two*feqW27_BS -(f27_TN *(q*omegaD-one)-omegaD*feq27_TN *(q-one))/(omegaD-one)+f27_BS *q)/(q+one);
      //q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[dirTN ])[ktn ]=(two*feqW27_TN -(f27_BS *(q*omegaD-one)-omegaD*feq27_BS *(q-one))/(omegaD-one)+f27_TN *q)/(q+one);
      //q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[dirTS ])[kts ]=(two*feqW27_TS -(f27_BN *(q*omegaD-one)-omegaD*feq27_BN *(q-one))/(omegaD-one)+f27_TS *q)/(q+one);
      //q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[dirBN ])[kbn ]=(two*feqW27_BN -(f27_TS *(q*omegaD-one)-omegaD*feq27_TS *(q-one))/(omegaD-one)+f27_BN *q)/(q+one);
      //q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[dirBSW])[kbsw]=(two*feqW27_BSW-(f27_TNE*(q*omegaD-one)-omegaD*feq27_TNE*(q-one))/(omegaD-one)+f27_BSW*q)/(q+one);
      //q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[dirTNE])[ktne]=(two*feqW27_TNE-(f27_BSW*(q*omegaD-one)-omegaD*feq27_BSW*(q-one))/(omegaD-one)+f27_TNE*q)/(q+one);
      //q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[dirTSW])[ktsw]=(two*feqW27_TSW-(f27_BNE*(q*omegaD-one)-omegaD*feq27_BNE*(q-one))/(omegaD-one)+f27_TSW*q)/(q+one);
      //q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[dirBNE])[kbne]=(two*feqW27_BNE-(f27_TSW*(q*omegaD-one)-omegaD*feq27_TSW*(q-one))/(omegaD-one)+f27_BNE*q)/(q+one);
      //q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[dirBNW])[kbnw]=(two*feqW27_BNW-(f27_TSE*(q*omegaD-one)-omegaD*feq27_TSE*(q-one))/(omegaD-one)+f27_BNW*q)/(q+one);
      //q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[dirTSE])[ktse]=(two*feqW27_TSE-(f27_BNW*(q*omegaD-one)-omegaD*feq27_BNW*(q-one))/(omegaD-one)+f27_TSE*q)/(q+one);
      //q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[dirTNW])[ktnw]=(two*feqW27_TNW-(f27_BSE*(q*omegaD-one)-omegaD*feq27_BSE*(q-one))/(omegaD-one)+f27_TNW*q)/(q+one);
      //q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[dirBSE])[kbse]=(two*feqW27_BSE-(f27_TNW*(q*omegaD-one)-omegaD*feq27_TNW*(q-one))/(omegaD-one)+f27_BSE*q)/(q+one);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QAD7( int inx,
                                 int iny,
                                 doubflo* DD, 
                                 doubflo* DD7, 
                                 doubflo* temp,
                                 doubflo diffusivity,
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

   Distributions7 D7;
   if (evenOrOdd==true)
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[1] = &DD7[1*size_Mat];
      D7.f[2] = &DD7[2*size_Mat];
      D7.f[3] = &DD7[3*size_Mat];
      D7.f[4] = &DD7[4*size_Mat];
      D7.f[5] = &DD7[5*size_Mat];
      D7.f[6] = &DD7[6*size_Mat];
   }
   else
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[2] = &DD7[1*size_Mat];
      D7.f[1] = &DD7[2*size_Mat];
      D7.f[4] = &DD7[3*size_Mat];
      D7.f[3] = &DD7[4*size_Mat];
      D7.f[6] = &DD7[5*size_Mat];
      D7.f[5] = &DD7[6*size_Mat];
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
      //////////////////////////////////////////////////////////////////////////////////
      doubflo  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB;//, 
      //         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
      //         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
      //         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
      //         *q_dirBSE, *q_dirBNW;

      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      //q_dirNE  = &QQ[dirNE  *sizeQ];
      //q_dirSW  = &QQ[dirSW  *sizeQ];
      //q_dirSE  = &QQ[dirSE  *sizeQ];
      //q_dirNW  = &QQ[dirNW  *sizeQ];
      //q_dirTE  = &QQ[dirTE  *sizeQ];
      //q_dirBW  = &QQ[dirBW  *sizeQ];
      //q_dirBE  = &QQ[dirBE  *sizeQ];
      //q_dirTW  = &QQ[dirTW  *sizeQ];
      //q_dirTN  = &QQ[dirTN  *sizeQ];
      //q_dirBS  = &QQ[dirBS  *sizeQ];
      //q_dirBN  = &QQ[dirBN  *sizeQ];
      //q_dirTS  = &QQ[dirTS  *sizeQ];
      //q_dirTNE = &QQ[dirTNE *sizeQ];
      //q_dirTSW = &QQ[dirTSW *sizeQ];
      //q_dirTSE = &QQ[dirTSE *sizeQ];
      //q_dirTNW = &QQ[dirTNW *sizeQ];
      //q_dirBNE = &QQ[dirBNE *sizeQ];
      //q_dirBSW = &QQ[dirBSW *sizeQ];
      //q_dirBSE = &QQ[dirBSE *sizeQ];
      //q_dirBNW = &QQ[dirBNW *sizeQ];
      //////////////////////////////////////////////////////////////////////////////////
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
      doubflo vx1, vx2, vx3/*, drho*/;
      //drho   =    f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
      //            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
      //            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      //vx1    = ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //         (f_E - f_W); 


      //vx2    = (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //         (f_N - f_S); 

      //vx3    = ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //         (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //         (f_T - f_B); 
      doubflo rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ ((D.f[dirZERO])[kzero]);
      doubflo rho    =  rho0 + one;
      doubflo OORho  =  one/rho;
      vx1     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      vx2     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      vx3     =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));

      //doubflo cu_sq       =c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      doubflo ux_sq       = vx1 * vx1;
      doubflo uy_sq       = vx2 * vx2;
      doubflo uz_sq       = vx3 * vx3;
      doubflo omegaD     = three - sqrt(three);
      doubflo Lam         = -(c1o2-one/omegaD);
      doubflo nue_d       = Lam/three;
      //doubflo ae          = zero;
      doubflo ae          = diffusivity/nue_d - one;

      doubflo f7_ZERO,f7_E,f7_W,f7_N,f7_S,f7_T,f7_B;
      doubflo /*feq7_ZERO,*/feq7_E,feq7_W,feq7_N,feq7_S,feq7_T,feq7_B;
      doubflo /*feqW7_ZERO,*/feqW7_E,feqW7_W,feqW7_N,feqW7_S,feqW7_T,feqW7_B;
      doubflo TempD = temp[k];


      f7_ZERO =  (D7.f[0])[kzero];
      f7_W    =  (D7.f[1])[ke   ];
      f7_E    =  (D7.f[2])[kw   ];
      f7_S    =  (D7.f[3])[kn   ];
      f7_N    =  (D7.f[4])[ks   ];
      f7_B    =  (D7.f[5])[kt   ];
      f7_T    =  (D7.f[6])[kb   ];

      doubflo ConcD = f7_ZERO + f7_E + f7_W + f7_N + f7_S + f7_T + f7_B;

      //feq7_ZERO = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
      feq7_E    = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
      feq7_W    = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
      feq7_N    = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
      feq7_S    = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
      feq7_T    = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
      feq7_B    = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);

      //feq7_ZERO = TempD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
      feqW7_E    = TempD*(c1o6*(ae+one));//+c1o2*(ux_sq)+vx1*c1o2);
      feqW7_W    = TempD*(c1o6*(ae+one));//+c1o2*(ux_sq)-vx1*c1o2);
      feqW7_N    = TempD*(c1o6*(ae+one));//+c1o2*(uy_sq)+vx2*c1o2);
      feqW7_S    = TempD*(c1o6*(ae+one));//+c1o2*(uy_sq)-vx2*c1o2);
      feqW7_T    = TempD*(c1o6*(ae+one));//+c1o2*(uz_sq)+vx3*c1o2);
      feqW7_B    = TempD*(c1o6*(ae+one));//+c1o2*(uz_sq)-vx3*c1o2);

      //////////////////////////////////////////////////////////////////////////
      //pointertausch
      if (evenOrOdd==false)
      {
         D7.f[0] = &DD7[0*size_Mat];
         D7.f[1] = &DD7[1*size_Mat];
         D7.f[2] = &DD7[2*size_Mat];
         D7.f[3] = &DD7[3*size_Mat];
         D7.f[4] = &DD7[4*size_Mat];
         D7.f[5] = &DD7[5*size_Mat];
         D7.f[6] = &DD7[6*size_Mat];
      }
      else
      {
         D7.f[0] = &DD7[0*size_Mat];
         D7.f[2] = &DD7[1*size_Mat];
         D7.f[1] = &DD7[2*size_Mat];
         D7.f[4] = &DD7[3*size_Mat];
         D7.f[3] = &DD7[4*size_Mat];
         D7.f[6] = &DD7[5*size_Mat];
         D7.f[5] = &DD7[6*size_Mat];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //(D7.f[1])[ke   ] = f7_E - feq7_E + feqW7_W; //E
      //(D7.f[2])[kw   ] = f7_W - feq7_W + feqW7_E; //W
      //(D7.f[3])[kn   ] = f7_N - feq7_N + feqW7_S; //N
      //(D7.f[4])[ks   ] = f7_S - feq7_S + feqW7_N; //S
      //(D7.f[5])[kt   ] = f7_T - feq7_T + feqW7_B; //T
      //(D7.f[6])[kb   ] = f7_B - feq7_B + feqW7_T; //B

      ////////////////////////////////////////////////////////////////////////////
      ////mit Q's
      //doubflo /*feq,*/ q;
      //q = q_dirE[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=zero;
      //   (D7.f[2])[kw]=(two*feqW7_W-(f7_E*(q*omegaD-one)-omegaD*feq7_E*(q-one))/(omegaD-one)+f7_W*q)/(q+one);//f7_W - feq7_W + feqW7_E;
      //}

      //q = q_dirW[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=zero;
      //   (D7.f[1])[ke]=(two*feqW7_E-(f7_W*(q*omegaD-one)-omegaD*feq7_W*(q-one))/(omegaD-one)+f7_E*q)/(q+one);//f7_E - feq7_E + feqW7_W;
      //}

      //q = q_dirN[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=zero;
      //   (D7.f[4])[ks]=(two*feqW7_S-(f7_N*(q*omegaD-one)-omegaD*feq7_N*(q-one))/(omegaD-one)+f7_S*q)/(q+one);//f7_S - feq7_S + feqW7_N;
      //}

      //q = q_dirS[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=zero;
      //   (D7.f[3])[kn]=(two*feqW7_N-(f7_S*(q*omegaD-one)-omegaD*feq7_S*(q-one))/(omegaD-one)+f7_N*q)/(q+one);//f7_N - feq7_N + feqW7_S;
      //}

      //q = q_dirT[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=zero;
      //   (D7.f[6])[kb]=(two*feqW7_B-(f7_T*(q*omegaD-one)-omegaD*feq7_T*(q-one))/(omegaD-one)+f7_B*q)/(q+one);//f7_B - feq7_B + feqW7_T;
      //}

      //q = q_dirB[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=zero;
      //   (D7.f[5])[kt]=(two*feqW7_T-(f7_B*(q*omegaD-one)-omegaD*feq7_B*(q-one))/(omegaD-one)+f7_T*q)/(q+one);//f7_T - feq7_T + feqW7_B;
      //}

      //////////////////////////////////////////////////////////////////////////
      //ohne Q's
      doubflo /*feq,*/ q;
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         (D7.f[2])[kw]=f7_W - feq7_W + feqW7_E;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         (D7.f[1])[ke]=f7_E - feq7_E + feqW7_W;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         (D7.f[4])[ks]=f7_S - feq7_S + feqW7_N;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         (D7.f[3])[kn]=f7_N - feq7_N + feqW7_S;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         (D7.f[6])[kb]=f7_B - feq7_B + feqW7_T;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         (D7.f[5])[kt]=f7_T - feq7_T + feqW7_B;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADDirichlet27(	 int inx,
											 int iny,
											 doubflo* DD, 
											 doubflo* DD27, 
											 doubflo* temp,
											 doubflo diffusivity,
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

   Distributions27 D27;
   if (evenOrOdd==true)
   {
      D27.f[dirE   ] = &DD27[dirE   *size_Mat];
      D27.f[dirW   ] = &DD27[dirW   *size_Mat];
      D27.f[dirN   ] = &DD27[dirN   *size_Mat];
      D27.f[dirS   ] = &DD27[dirS   *size_Mat];
      D27.f[dirT   ] = &DD27[dirT   *size_Mat];
      D27.f[dirB   ] = &DD27[dirB   *size_Mat];
      D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
      D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
      D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
      D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
      D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
      D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
   } 
   else
   {
      D27.f[dirW   ] = &DD27[dirE   *size_Mat];
      D27.f[dirE   ] = &DD27[dirW   *size_Mat];
      D27.f[dirS   ] = &DD27[dirN   *size_Mat];
      D27.f[dirN   ] = &DD27[dirS   *size_Mat];
      D27.f[dirB   ] = &DD27[dirT   *size_Mat];
      D27.f[dirT   ] = &DD27[dirB   *size_Mat];
      D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
      D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
      D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
      D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
      D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
      D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
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
      doubflo  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
      doubflo f_ZERO = (D.f[dirZERO])[kzero];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, /*drho, feq,*/ q;
      ////drho   = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
      ////         f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
      ////         f_T + f_B + f_N + f_S + f_E + f_W + f_ZERO; 

      //vx1    = ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //         (f_E - f_W); 


      //vx2    = (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //         (f_N - f_S); 

      //vx3    = ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //         (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //         (f_T - f_B); 
      doubflo rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+f_ZERO;
      doubflo rho    =  rho0 + one;
      doubflo OORho  =  one/rho;
      vx1     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      vx2     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      vx3     =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      doubflo f27_W    = (D27.f[dirE   ])[ke   ];
      doubflo f27_E    = (D27.f[dirW   ])[kw   ];
      doubflo f27_S    = (D27.f[dirN   ])[kn   ];
      doubflo f27_N    = (D27.f[dirS   ])[ks   ];
      doubflo f27_B    = (D27.f[dirT   ])[kt   ];
      doubflo f27_T    = (D27.f[dirB   ])[kb   ];
      doubflo f27_SW   = (D27.f[dirNE  ])[kne  ];
      doubflo f27_NE   = (D27.f[dirSW  ])[ksw  ];
      doubflo f27_NW   = (D27.f[dirSE  ])[kse  ];
      doubflo f27_SE   = (D27.f[dirNW  ])[knw  ];
      doubflo f27_BW   = (D27.f[dirTE  ])[kte  ];
      doubflo f27_TE   = (D27.f[dirBW  ])[kbw  ];
      doubflo f27_TW   = (D27.f[dirBE  ])[kbe  ];
      doubflo f27_BE   = (D27.f[dirTW  ])[ktw  ];
      doubflo f27_BS   = (D27.f[dirTN  ])[ktn  ];
      doubflo f27_TN   = (D27.f[dirBS  ])[kbs  ];
      doubflo f27_TS   = (D27.f[dirBN  ])[kbn  ];
      doubflo f27_BN   = (D27.f[dirTS  ])[kts  ];
      doubflo f27_ZERO = (D27.f[dirZERO])[kzero];
      doubflo f27_BSW  = (D27.f[dirTNE ])[ktne ];
      doubflo f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      doubflo f27_BNW  = (D27.f[dirTSE ])[ktse ];
      doubflo f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      doubflo f27_TSW  = (D27.f[dirBNE ])[kbne ];
      doubflo f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      doubflo f27_TNW  = (D27.f[dirBSE ])[kbse ];
      doubflo f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      doubflo ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
         f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
         f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      //doubflo feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      doubflo feq27_E    =   c2over27* ConcD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feq27_W    =   c2over27* ConcD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feq27_N    =   c2over27* ConcD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feq27_S    =   c2over27* ConcD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feq27_T    =   c2over27* ConcD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feq27_B    =   c2over27* ConcD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feq27_NE   =   c1over54* ConcD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feq27_SW   =   c1over54* ConcD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feq27_SE   =   c1over54* ConcD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feq27_NW   =   c1over54* ConcD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feq27_TE   =   c1over54* ConcD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feq27_BW   =   c1over54* ConcD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feq27_BE   =   c1over54* ConcD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feq27_TW   =   c1over54* ConcD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feq27_TN   =   c1over54* ConcD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feq27_BS   =   c1over54* ConcD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feq27_BN   =   c1over54* ConcD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feq27_TS   =   c1over54* ConcD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feq27_TNE  =   c1over216*ConcD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feq27_BSW  =   c1over216*ConcD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feq27_BNE  =   c1over216*ConcD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feq27_TSW  =   c1over216*ConcD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feq27_TSE  =   c1over216*ConcD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feq27_BNW  =   c1over216*ConcD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feq27_BSE  =   c1over216*ConcD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feq27_TNW  =   c1over216*ConcD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo TempD = one;//temp[k];

      //doubflo feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      doubflo feqW27_E    =   c2over27* TempD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feqW27_W    =   c2over27* TempD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feqW27_N    =   c2over27* TempD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feqW27_S    =   c2over27* TempD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feqW27_T    =   c2over27* TempD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feqW27_B    =   c2over27* TempD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feqW27_NE   =   c1over54* TempD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feqW27_SW   =   c1over54* TempD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feqW27_SE   =   c1over54* TempD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feqW27_NW   =   c1over54* TempD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feqW27_TE   =   c1over54* TempD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feqW27_BW   =   c1over54* TempD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feqW27_BE   =   c1over54* TempD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feqW27_TW   =   c1over54* TempD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feqW27_TN   =   c1over54* TempD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feqW27_BS   =   c1over54* TempD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feqW27_BN   =   c1over54* TempD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feqW27_TS   =   c1over54* TempD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feqW27_TNE  =   c1over216*TempD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feqW27_BSW  =   c1over216*TempD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feqW27_BNE  =   c1over216*TempD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feqW27_TSW  =   c1over216*TempD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feqW27_TSE  =   c1over216*TempD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feqW27_BNW  =   c1over216*TempD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feqW27_BSE  =   c1over216*TempD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feqW27_TNW  =   c1over216*TempD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo omegaD     = three - sqrt(three);
      //doubflo Lam        = -(c1o2-one/omegaD);
      //doubflo nue_d      = Lam/three;
      //doubflo ae         = zero;
      //doubflo ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
      {
         D27.f[dirE   ] = &DD27[dirE   *size_Mat];
         D27.f[dirW   ] = &DD27[dirW   *size_Mat];
         D27.f[dirN   ] = &DD27[dirN   *size_Mat];
         D27.f[dirS   ] = &DD27[dirS   *size_Mat];
         D27.f[dirT   ] = &DD27[dirT   *size_Mat];
         D27.f[dirB   ] = &DD27[dirB   *size_Mat];
         D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
         D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
         D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
         D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
         D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
         D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
      } 
      else
      {
         D27.f[dirW   ] = &DD27[dirE   *size_Mat];
         D27.f[dirE   ] = &DD27[dirW   *size_Mat];
         D27.f[dirS   ] = &DD27[dirN   *size_Mat];
         D27.f[dirN   ] = &DD27[dirS   *size_Mat];
         D27.f[dirB   ] = &DD27[dirT   *size_Mat];
         D27.f[dirT   ] = &DD27[dirB   *size_Mat];
         D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
         D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
         D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
         D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
         D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
         D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=0.1f;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[  ke   ]; if (q>=zero && q<=one) (D27.f[dirW  ])[kw  ]=(two*feqW27_W  -(f27_E  *(q*omegaD-one)-omegaD*feq27_E  *(q-one))/(omegaD-one)+f27_W  *q)/(q+one);
      q = q_dirW[  kw   ]; if (q>=zero && q<=one) (D27.f[dirE  ])[ke  ]=(two*feqW27_E  -(f27_W  *(q*omegaD-one)-omegaD*feq27_W  *(q-one))/(omegaD-one)+f27_E  *q)/(q+one);
      q = q_dirN[  kn   ]; if (q>=zero && q<=one) (D27.f[dirS  ])[ks  ]=(two*feqW27_S  -(f27_N  *(q*omegaD-one)-omegaD*feq27_N  *(q-one))/(omegaD-one)+f27_S  *q)/(q+one);
      q = q_dirS[  ks   ]; if (q>=zero && q<=one) (D27.f[dirN  ])[kn  ]=(two*feqW27_N  -(f27_S  *(q*omegaD-one)-omegaD*feq27_S  *(q-one))/(omegaD-one)+f27_N  *q)/(q+one);
      q = q_dirT[  kt   ]; if (q>=zero && q<=one) (D27.f[dirB  ])[kb  ]=(two*feqW27_B  -(f27_T  *(q*omegaD-one)-omegaD*feq27_T  *(q-one))/(omegaD-one)+f27_B  *q)/(q+one);
      q = q_dirB[  kb   ]; if (q>=zero && q<=one) (D27.f[dirT  ])[kt  ]=(two*feqW27_T  -(f27_B  *(q*omegaD-one)-omegaD*feq27_B  *(q-one))/(omegaD-one)+f27_T  *q)/(q+one);
      q = q_dirNE[ kne  ]; if (q>=zero && q<=one) (D27.f[dirSW ])[ksw ]=(two*feqW27_SW -(f27_NE *(q*omegaD-one)-omegaD*feq27_NE *(q-one))/(omegaD-one)+f27_SW *q)/(q+one);
      q = q_dirSW[ ksw  ]; if (q>=zero && q<=one) (D27.f[dirNE ])[kne ]=(two*feqW27_NE -(f27_SW *(q*omegaD-one)-omegaD*feq27_SW *(q-one))/(omegaD-one)+f27_NE *q)/(q+one);
      q = q_dirSE[ kse  ]; if (q>=zero && q<=one) (D27.f[dirNW ])[knw ]=(two*feqW27_NW -(f27_SE *(q*omegaD-one)-omegaD*feq27_SE *(q-one))/(omegaD-one)+f27_NW *q)/(q+one);
      q = q_dirNW[ knw  ]; if (q>=zero && q<=one) (D27.f[dirSE ])[kse ]=(two*feqW27_SE -(f27_NW *(q*omegaD-one)-omegaD*feq27_NW *(q-one))/(omegaD-one)+f27_SE *q)/(q+one);
      q = q_dirTE[ kte  ]; if (q>=zero && q<=one) (D27.f[dirBW ])[kbw ]=(two*feqW27_BW -(f27_TE *(q*omegaD-one)-omegaD*feq27_TE *(q-one))/(omegaD-one)+f27_BW *q)/(q+one);
      q = q_dirBW[ kbw  ]; if (q>=zero && q<=one) (D27.f[dirTE ])[kte ]=(two*feqW27_TE -(f27_BW *(q*omegaD-one)-omegaD*feq27_BW *(q-one))/(omegaD-one)+f27_TE *q)/(q+one);
      q = q_dirBE[ kbe  ]; if (q>=zero && q<=one) (D27.f[dirTW ])[ktw ]=(two*feqW27_TW -(f27_BE *(q*omegaD-one)-omegaD*feq27_BE *(q-one))/(omegaD-one)+f27_TW *q)/(q+one);
      q = q_dirTW[ ktw  ]; if (q>=zero && q<=one) (D27.f[dirBE ])[kbe ]=(two*feqW27_BE -(f27_TW *(q*omegaD-one)-omegaD*feq27_TW *(q-one))/(omegaD-one)+f27_BE *q)/(q+one);
      q = q_dirTN[ ktn  ]; if (q>=zero && q<=one) (D27.f[dirBS ])[kbs ]=(two*feqW27_BS -(f27_TN *(q*omegaD-one)-omegaD*feq27_TN *(q-one))/(omegaD-one)+f27_BS *q)/(q+one);
      q = q_dirBS[ kbs  ]; if (q>=zero && q<=one) (D27.f[dirTN ])[ktn ]=(two*feqW27_TN -(f27_BS *(q*omegaD-one)-omegaD*feq27_BS *(q-one))/(omegaD-one)+f27_TN *q)/(q+one);
      q = q_dirBN[ kbn  ]; if (q>=zero && q<=one) (D27.f[dirTS ])[kts ]=(two*feqW27_TS -(f27_BN *(q*omegaD-one)-omegaD*feq27_BN *(q-one))/(omegaD-one)+f27_TS *q)/(q+one);
      q = q_dirTS[ kts  ]; if (q>=zero && q<=one) (D27.f[dirBN ])[kbn ]=(two*feqW27_BN -(f27_TS *(q*omegaD-one)-omegaD*feq27_TS *(q-one))/(omegaD-one)+f27_BN *q)/(q+one);
      q = q_dirTNE[ktne ]; if (q>=zero && q<=one) (D27.f[dirBSW])[kbsw]=(two*feqW27_BSW-(f27_TNE*(q*omegaD-one)-omegaD*feq27_TNE*(q-one))/(omegaD-one)+f27_BSW*q)/(q+one);
      q = q_dirBSW[kbsw ]; if (q>=zero && q<=one) (D27.f[dirTNE])[ktne]=(two*feqW27_TNE-(f27_BSW*(q*omegaD-one)-omegaD*feq27_BSW*(q-one))/(omegaD-one)+f27_TNE*q)/(q+one);
      q = q_dirBNE[kbne ]; if (q>=zero && q<=one) (D27.f[dirTSW])[ktsw]=(two*feqW27_TSW-(f27_BNE*(q*omegaD-one)-omegaD*feq27_BNE*(q-one))/(omegaD-one)+f27_TSW*q)/(q+one);
      q = q_dirTSW[ktsw ]; if (q>=zero && q<=one) (D27.f[dirBNE])[kbne]=(two*feqW27_BNE-(f27_TSW*(q*omegaD-one)-omegaD*feq27_TSW*(q-one))/(omegaD-one)+f27_BNE*q)/(q+one);
      q = q_dirTSE[ktse ]; if (q>=zero && q<=one) (D27.f[dirBNW])[kbnw]=(two*feqW27_BNW-(f27_TSE*(q*omegaD-one)-omegaD*feq27_TSE*(q-one))/(omegaD-one)+f27_BNW*q)/(q+one);
      q = q_dirBNW[kbnw ]; if (q>=zero && q<=one) (D27.f[dirTSE])[ktse]=(two*feqW27_TSE-(f27_BNW*(q*omegaD-one)-omegaD*feq27_BNW*(q-one))/(omegaD-one)+f27_TSE*q)/(q+one);
      q = q_dirBSE[kbse ]; if (q>=zero && q<=one) (D27.f[dirTNW])[ktnw]=(two*feqW27_TNW-(f27_BSE*(q*omegaD-one)-omegaD*feq27_BSE*(q-one))/(omegaD-one)+f27_TNW*q)/(q+one);
      q = q_dirTNW[ktnw ]; if (q>=zero && q<=one) (D27.f[dirBSE])[kbse]=(two*feqW27_BSE-(f27_TNW*(q*omegaD-one)-omegaD*feq27_TNW*(q-one))/(omegaD-one)+f27_BSE*q)/(q+one);
      //q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[dirW  ])[kw  ]=(two*feqW27_W  -(f27_E  *(q*omegaD-one)-omegaD*feq27_E  *(q-one))/(omegaD-one)+f27_W  *q)/(q+one);
      //q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[dirE  ])[ke  ]=(two*feqW27_E  -(f27_W  *(q*omegaD-one)-omegaD*feq27_W  *(q-one))/(omegaD-one)+f27_E  *q)/(q+one);
      //q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[dirS  ])[ks  ]=(two*feqW27_S  -(f27_N  *(q*omegaD-one)-omegaD*feq27_N  *(q-one))/(omegaD-one)+f27_S  *q)/(q+one);
      //q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[dirN  ])[kn  ]=(two*feqW27_N  -(f27_S  *(q*omegaD-one)-omegaD*feq27_S  *(q-one))/(omegaD-one)+f27_N  *q)/(q+one);
      //q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[dirB  ])[kb  ]=(two*feqW27_B  -(f27_T  *(q*omegaD-one)-omegaD*feq27_T  *(q-one))/(omegaD-one)+f27_B  *q)/(q+one);
      //q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[dirT  ])[kt  ]=(two*feqW27_T  -(f27_B  *(q*omegaD-one)-omegaD*feq27_B  *(q-one))/(omegaD-one)+f27_T  *q)/(q+one);
      //q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[dirSW ])[ksw ]=(two*feqW27_SW -(f27_NE *(q*omegaD-one)-omegaD*feq27_NE *(q-one))/(omegaD-one)+f27_SW *q)/(q+one);
      //q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[dirNE ])[kne ]=(two*feqW27_NE -(f27_SW *(q*omegaD-one)-omegaD*feq27_SW *(q-one))/(omegaD-one)+f27_NE *q)/(q+one);
      //q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[dirNW ])[knw ]=(two*feqW27_NW -(f27_SE *(q*omegaD-one)-omegaD*feq27_SE *(q-one))/(omegaD-one)+f27_NW *q)/(q+one);
      //q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[dirSE ])[kse ]=(two*feqW27_SE -(f27_NW *(q*omegaD-one)-omegaD*feq27_NW *(q-one))/(omegaD-one)+f27_SE *q)/(q+one);
      //q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[dirBW ])[kbw ]=(two*feqW27_BW -(f27_TE *(q*omegaD-one)-omegaD*feq27_TE *(q-one))/(omegaD-one)+f27_BW *q)/(q+one);
      //q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[dirTE ])[kte ]=(two*feqW27_TE -(f27_BW *(q*omegaD-one)-omegaD*feq27_BW *(q-one))/(omegaD-one)+f27_TE *q)/(q+one);
      //q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[dirTW ])[ktw ]=(two*feqW27_TW -(f27_BE *(q*omegaD-one)-omegaD*feq27_BE *(q-one))/(omegaD-one)+f27_TW *q)/(q+one);
      //q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[dirBE ])[kbe ]=(two*feqW27_BE -(f27_TW *(q*omegaD-one)-omegaD*feq27_TW *(q-one))/(omegaD-one)+f27_BE *q)/(q+one);
      //q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[dirBS ])[kbs ]=(two*feqW27_BS -(f27_TN *(q*omegaD-one)-omegaD*feq27_TN *(q-one))/(omegaD-one)+f27_BS *q)/(q+one);
      //q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[dirTN ])[ktn ]=(two*feqW27_TN -(f27_BS *(q*omegaD-one)-omegaD*feq27_BS *(q-one))/(omegaD-one)+f27_TN *q)/(q+one);
      //q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[dirTS ])[kts ]=(two*feqW27_TS -(f27_BN *(q*omegaD-one)-omegaD*feq27_BN *(q-one))/(omegaD-one)+f27_TS *q)/(q+one);
      //q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[dirBN ])[kbn ]=(two*feqW27_BN -(f27_TS *(q*omegaD-one)-omegaD*feq27_TS *(q-one))/(omegaD-one)+f27_BN *q)/(q+one);
      //q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[dirBSW])[kbsw]=(two*feqW27_BSW-(f27_TNE*(q*omegaD-one)-omegaD*feq27_TNE*(q-one))/(omegaD-one)+f27_BSW*q)/(q+one);
      //q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[dirTNE])[ktne]=(two*feqW27_TNE-(f27_BSW*(q*omegaD-one)-omegaD*feq27_BSW*(q-one))/(omegaD-one)+f27_TNE*q)/(q+one);
      //q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[dirTSW])[ktsw]=(two*feqW27_TSW-(f27_BNE*(q*omegaD-one)-omegaD*feq27_BNE*(q-one))/(omegaD-one)+f27_TSW*q)/(q+one);
      //q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[dirBNE])[kbne]=(two*feqW27_BNE-(f27_TSW*(q*omegaD-one)-omegaD*feq27_TSW*(q-one))/(omegaD-one)+f27_BNE*q)/(q+one);
      //q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[dirBNW])[kbnw]=(two*feqW27_BNW-(f27_TSE*(q*omegaD-one)-omegaD*feq27_TSE*(q-one))/(omegaD-one)+f27_BNW*q)/(q+one);
      //q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[dirTSE])[ktse]=(two*feqW27_TSE-(f27_BNW*(q*omegaD-one)-omegaD*feq27_BNW*(q-one))/(omegaD-one)+f27_TSE*q)/(q+one);
      //q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[dirTNW])[ktnw]=(two*feqW27_TNW-(f27_BSE*(q*omegaD-one)-omegaD*feq27_BSE*(q-one))/(omegaD-one)+f27_TNW*q)/(q+one);
      //q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[dirBSE])[kbse]=(two*feqW27_BSE-(f27_TNW*(q*omegaD-one)-omegaD*feq27_TNW*(q-one))/(omegaD-one)+f27_BSE*q)/(q+one);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADBB27(int inx,
                                   int iny,
                                   doubflo* DD, 
                                   doubflo* DD27, 
                                   doubflo* temp,
                                   doubflo diffusivity,
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

   Distributions27 D27;
   if (evenOrOdd==true)
   {
      D27.f[dirE   ] = &DD27[dirE   *size_Mat];
      D27.f[dirW   ] = &DD27[dirW   *size_Mat];
      D27.f[dirN   ] = &DD27[dirN   *size_Mat];
      D27.f[dirS   ] = &DD27[dirS   *size_Mat];
      D27.f[dirT   ] = &DD27[dirT   *size_Mat];
      D27.f[dirB   ] = &DD27[dirB   *size_Mat];
      D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
      D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
      D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
      D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
      D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
      D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
   } 
   else
   {
      D27.f[dirW   ] = &DD27[dirE   *size_Mat];
      D27.f[dirE   ] = &DD27[dirW   *size_Mat];
      D27.f[dirS   ] = &DD27[dirN   *size_Mat];
      D27.f[dirN   ] = &DD27[dirS   *size_Mat];
      D27.f[dirB   ] = &DD27[dirT   *size_Mat];
      D27.f[dirT   ] = &DD27[dirB   *size_Mat];
      D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
      D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
      D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
      D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
      D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
      D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
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
      doubflo  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
      doubflo f_ZERO = (D.f[dirZERO])[kzero];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, /*drho, feq,*/ q;
      ////drho   = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
      ////         f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
      ////         f_T + f_B + f_N + f_S + f_E + f_W + f_ZERO; 

      //vx1    = ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //         (f_E - f_W); 


      //vx2    = (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //         (f_N - f_S); 

      //vx3    = ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //         (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //         (f_T - f_B); 
      doubflo rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+f_ZERO;
      doubflo rho    =  rho0 + one;
      doubflo OORho  =  one/rho;
      vx1     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      vx2     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      vx3     =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      doubflo f27_W    = (D27.f[dirE   ])[ke   ];
      doubflo f27_E    = (D27.f[dirW   ])[kw   ];
      doubflo f27_S    = (D27.f[dirN   ])[kn   ];
      doubflo f27_N    = (D27.f[dirS   ])[ks   ];
      doubflo f27_B    = (D27.f[dirT   ])[kt   ];
      doubflo f27_T    = (D27.f[dirB   ])[kb   ];
      doubflo f27_SW   = (D27.f[dirNE  ])[kne  ];
      doubflo f27_NE   = (D27.f[dirSW  ])[ksw  ];
      doubflo f27_NW   = (D27.f[dirSE  ])[kse  ];
      doubflo f27_SE   = (D27.f[dirNW  ])[knw  ];
      doubflo f27_BW   = (D27.f[dirTE  ])[kte  ];
      doubflo f27_TE   = (D27.f[dirBW  ])[kbw  ];
      doubflo f27_TW   = (D27.f[dirBE  ])[kbe  ];
      doubflo f27_BE   = (D27.f[dirTW  ])[ktw  ];
      doubflo f27_BS   = (D27.f[dirTN  ])[ktn  ];
      doubflo f27_TN   = (D27.f[dirBS  ])[kbs  ];
      doubflo f27_TS   = (D27.f[dirBN  ])[kbn  ];
      doubflo f27_BN   = (D27.f[dirTS  ])[kts  ];
      doubflo f27_ZERO = (D27.f[dirZERO])[kzero];
      doubflo f27_BSW  = (D27.f[dirTNE ])[ktne ];
      doubflo f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      doubflo f27_BNW  = (D27.f[dirTSE ])[ktse ];
      doubflo f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      doubflo f27_TSW  = (D27.f[dirBNE ])[kbne ];
      doubflo f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      doubflo f27_TNW  = (D27.f[dirBSE ])[kbse ];
      doubflo f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      doubflo ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
         f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
         f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      //doubflo feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      doubflo feq27_E    =   c2over27* ConcD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feq27_W    =   c2over27* ConcD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feq27_N    =   c2over27* ConcD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feq27_S    =   c2over27* ConcD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feq27_T    =   c2over27* ConcD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feq27_B    =   c2over27* ConcD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feq27_NE   =   c1over54* ConcD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feq27_SW   =   c1over54* ConcD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feq27_SE   =   c1over54* ConcD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feq27_NW   =   c1over54* ConcD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feq27_TE   =   c1over54* ConcD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feq27_BW   =   c1over54* ConcD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feq27_BE   =   c1over54* ConcD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feq27_TW   =   c1over54* ConcD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feq27_TN   =   c1over54* ConcD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feq27_BS   =   c1over54* ConcD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feq27_BN   =   c1over54* ConcD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feq27_TS   =   c1over54* ConcD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feq27_TNE  =   c1over216*ConcD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feq27_BSW  =   c1over216*ConcD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feq27_BNE  =   c1over216*ConcD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feq27_TSW  =   c1over216*ConcD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feq27_TSE  =   c1over216*ConcD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feq27_BNW  =   c1over216*ConcD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feq27_BSE  =   c1over216*ConcD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feq27_TNW  =   c1over216*ConcD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo TempD = temp[k];

      //doubflo feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      doubflo feqW27_E    =   c2over27* TempD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feqW27_W    =   c2over27* TempD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feqW27_N    =   c2over27* TempD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feqW27_S    =   c2over27* TempD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feqW27_T    =   c2over27* TempD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feqW27_B    =   c2over27* TempD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feqW27_NE   =   c1over54* TempD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feqW27_SW   =   c1over54* TempD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feqW27_SE   =   c1over54* TempD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feqW27_NW   =   c1over54* TempD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feqW27_TE   =   c1over54* TempD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feqW27_BW   =   c1over54* TempD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feqW27_BE   =   c1over54* TempD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feqW27_TW   =   c1over54* TempD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feqW27_TN   =   c1over54* TempD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feqW27_BS   =   c1over54* TempD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feqW27_BN   =   c1over54* TempD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feqW27_TS   =   c1over54* TempD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feqW27_TNE  =   c1over216*TempD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feqW27_BSW  =   c1over216*TempD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feqW27_BNE  =   c1over216*TempD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feqW27_TSW  =   c1over216*TempD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feqW27_TSE  =   c1over216*TempD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feqW27_BNW  =   c1over216*TempD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feqW27_BSE  =   c1over216*TempD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feqW27_TNW  =   c1over216*TempD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo omegaD     = three - sqrt(three);
      //doubflo Lam        = -(c1o2-one/omegaD);
      //doubflo nue_d      = Lam/three;
      //doubflo ae         = zero;
      //doubflo ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
      {
         D27.f[dirE   ] = &DD27[dirE   *size_Mat];
         D27.f[dirW   ] = &DD27[dirW   *size_Mat];
         D27.f[dirN   ] = &DD27[dirN   *size_Mat];
         D27.f[dirS   ] = &DD27[dirS   *size_Mat];
         D27.f[dirT   ] = &DD27[dirT   *size_Mat];
         D27.f[dirB   ] = &DD27[dirB   *size_Mat];
         D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
         D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
         D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
         D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
         D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
         D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
      } 
      else
      {
         D27.f[dirW   ] = &DD27[dirE   *size_Mat];
         D27.f[dirE   ] = &DD27[dirW   *size_Mat];
         D27.f[dirS   ] = &DD27[dirN   *size_Mat];
         D27.f[dirN   ] = &DD27[dirS   *size_Mat];
         D27.f[dirB   ] = &DD27[dirT   *size_Mat];
         D27.f[dirT   ] = &DD27[dirB   *size_Mat];
         D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
         D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
         D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
         D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
         D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
         D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=0.1f;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[dirW  ])[kw  ]=f27_E  ;
      q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[dirE  ])[ke  ]=f27_W  ;
      q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[dirS  ])[ks  ]=f27_N  ;
      q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[dirN  ])[kn  ]=f27_S  ;
      q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[dirB  ])[kb  ]=f27_T  ;
      q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[dirT  ])[kt  ]=f27_B  ;
      q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[dirSW ])[ksw ]=f27_NE ;
      q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[dirNE ])[kne ]=f27_SW ;
      q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[dirNW ])[knw ]=f27_SE ;
      q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[dirSE ])[kse ]=f27_NW ;
      q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[dirBW ])[kbw ]=f27_TE ;
      q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[dirTE ])[kte ]=f27_BW ;
      q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[dirTW ])[ktw ]=f27_BE ;
      q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[dirBE ])[kbe ]=f27_TW ;
      q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[dirBS ])[kbs ]=f27_TN ;
      q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[dirTN ])[ktn ]=f27_BS ;
      q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[dirTS ])[kts ]=f27_BN ;
      q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[dirBN ])[kbn ]=f27_TS ;
      q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[dirBSW])[kbsw]=f27_TNE;
      q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[dirTNE])[ktne]=f27_BSW;
      q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[dirTSW])[ktsw]=f27_BNE;
      q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[dirBNE])[kbne]=f27_TSW;
      q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[dirBNW])[kbnw]=f27_TSE;
      q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[dirTSE])[ktse]=f27_BNW;
      q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[dirTNW])[ktnw]=f27_BSE;
      q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[dirBSE])[kbse]=f27_TNW;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////																																		  //////
//////                 										incomp   																		  //////
//////																																		  //////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QNoSlipADincomp7( int inx,
											 int iny,
											 doubflo* DD, 
											 doubflo* DD7, 
											 doubflo* temp,
											 doubflo diffusivity,
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

   Distributions7 D7;
   if (evenOrOdd==true)
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[1] = &DD7[1*size_Mat];
      D7.f[2] = &DD7[2*size_Mat];
      D7.f[3] = &DD7[3*size_Mat];
      D7.f[4] = &DD7[4*size_Mat];
      D7.f[5] = &DD7[5*size_Mat];
      D7.f[6] = &DD7[6*size_Mat];
   }
   else
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[2] = &DD7[1*size_Mat];
      D7.f[1] = &DD7[2*size_Mat];
      D7.f[4] = &DD7[3*size_Mat];
      D7.f[3] = &DD7[4*size_Mat];
      D7.f[6] = &DD7[5*size_Mat];
      D7.f[5] = &DD7[6*size_Mat];
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
      //////////////////////////////////////////////////////////////////////////////////
      doubflo  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB;

      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      //////////////////////////////////////////////////////////////////////////////////
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
      doubflo vx1 =  ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      doubflo vx2 =  ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      doubflo vx3 =  ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
		 ////dörrrrrty !!!!!!!!!!!!!
   //      doubflo vx1     =  ten * ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
   //      doubflo vx2     =  ten * ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
   //      doubflo vx3     =  ten * ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));

      //doubflo cu_sq       =c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      doubflo ux_sq       = vx1 * vx1;
      doubflo uy_sq       = vx2 * vx2;
      doubflo uz_sq       = vx3 * vx3;
      ////////////////////////////////////////////////////////////////////////////////
	  //BGK
      //doubflo omegaD     = three - sqrt(three);
      //doubflo Lam         = -(c1o2-one/omegaD);
      //doubflo nue_d       = Lam/three;
      //doubflo ae          = diffusivity/nue_d - one; //zero;

      doubflo f7_ZERO,f7_E,f7_W,f7_N,f7_S,f7_T,f7_B;
      doubflo /*feq7_ZERO,*/feq7_E,feq7_W,feq7_N,feq7_S,feq7_T,feq7_B;
      doubflo /*feqW7_ZERO,*/feqW7_E,feqW7_W,feqW7_N,feqW7_S,feqW7_T,feqW7_B;
      doubflo TempD = temp[k];


      f7_ZERO =  (D7.f[0])[kzero];
      f7_W    =  (D7.f[1])[ke   ];
      f7_E    =  (D7.f[2])[kw   ];
      f7_S    =  (D7.f[3])[kn   ];
      f7_N    =  (D7.f[4])[ks   ];
      f7_B    =  (D7.f[5])[kt   ];
      f7_T    =  (D7.f[6])[kb   ];

      doubflo ConcD = f7_ZERO + f7_E + f7_W + f7_N + f7_S + f7_T + f7_B;

      ////feq7_ZERO = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
      //feq7_E    = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
      //feq7_W    = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
      //feq7_N    = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
      //feq7_S    = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
      //feq7_T    = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
      //feq7_B    = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);

      ////feq7_ZERO = TempD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
      //feqW7_E    = TempD*(c1o6*(ae+one));//+c1o2*(ux_sq)+vx1*c1o2);
      //feqW7_W    = TempD*(c1o6*(ae+one));//+c1o2*(ux_sq)-vx1*c1o2);
      //feqW7_N    = TempD*(c1o6*(ae+one));//+c1o2*(uy_sq)+vx2*c1o2);
      //feqW7_S    = TempD*(c1o6*(ae+one));//+c1o2*(uy_sq)-vx2*c1o2);
      //feqW7_T    = TempD*(c1o6*(ae+one));//+c1o2*(uz_sq)+vx3*c1o2);
      //feqW7_B    = TempD*(c1o6*(ae+one));//+c1o2*(uz_sq)-vx3*c1o2);

      ////////////////////////////////////////////////////////////////////////////////
	  //TRT
      doubflo cs2     = c1o4;
      //////////////////////////////////////////////////////////////////////////
      //pointertausch
      if (evenOrOdd==false)
      {
         D7.f[0] = &DD7[0*size_Mat];
         D7.f[1] = &DD7[1*size_Mat];
         D7.f[2] = &DD7[2*size_Mat];
         D7.f[3] = &DD7[3*size_Mat];
         D7.f[4] = &DD7[4*size_Mat];
         D7.f[5] = &DD7[5*size_Mat];
         D7.f[6] = &DD7[6*size_Mat];
      }
      else
      {
         D7.f[0] = &DD7[0*size_Mat];
         D7.f[2] = &DD7[1*size_Mat];
         D7.f[1] = &DD7[2*size_Mat];
         D7.f[4] = &DD7[3*size_Mat];
         D7.f[3] = &DD7[4*size_Mat];
         D7.f[6] = &DD7[5*size_Mat];
         D7.f[5] = &DD7[6*size_Mat];
      }

      ////////////////////////////////////////////////////////////////////////////
      ////mit Q's
      //doubflo /*feq,*/ q;
      //q = q_dirE[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=zero;
      //   (D7.f[2])[kw]=(two*feqW7_W-(f7_E*(q*omegaD-one)-omegaD*feq7_E*(q-one))/(omegaD-one)+f7_W*q)/(q+one);//f7_W - feq7_W + feqW7_E;
      //}

      //q = q_dirW[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=zero;
      //   (D7.f[1])[ke]=(two*feqW7_E-(f7_W*(q*omegaD-one)-omegaD*feq7_W*(q-one))/(omegaD-one)+f7_E*q)/(q+one);//f7_E - feq7_E + feqW7_W;
      //}

      //q = q_dirN[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=zero;
      //   (D7.f[4])[ks]=(two*feqW7_S-(f7_N*(q*omegaD-one)-omegaD*feq7_N*(q-one))/(omegaD-one)+f7_S*q)/(q+one);//f7_S - feq7_S + feqW7_N;
      //}

      //q = q_dirS[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=zero;
      //   (D7.f[3])[kn]=(two*feqW7_N-(f7_S*(q*omegaD-one)-omegaD*feq7_S*(q-one))/(omegaD-one)+f7_N*q)/(q+one);//f7_N - feq7_N + feqW7_S;
      //}

      //q = q_dirT[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=zero;
      //   (D7.f[6])[kb]=(two*feqW7_B-(f7_T*(q*omegaD-one)-omegaD*feq7_T*(q-one))/(omegaD-one)+f7_B*q)/(q+one);//f7_B - feq7_B + feqW7_T;
      //}

      //q = q_dirB[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=zero;
      //   (D7.f[5])[kt]=(two*feqW7_T-(f7_B*(q*omegaD-one)-omegaD*feq7_B*(q-one))/(omegaD-one)+f7_T*q)/(q+one);//f7_T - feq7_T + feqW7_B;
      //}

      ////////////////////////////////////////////////////////////////////////////
      ////ohne Q's
      //doubflo /*feq,*/ q;
      //q = q_dirE[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[2])[kw]=f7_W;// - feq7_W + feqW7_E;
      //}

      //q = q_dirW[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[1])[ke]=f7_E;// - feq7_E + feqW7_W;
      //}

      //q = q_dirN[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[4])[ks]=f7_S;// - feq7_S + feqW7_N;
      //}

      //q = q_dirS[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[3])[kn]=f7_N;// - feq7_N + feqW7_S;
      //}

      //q = q_dirT[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[6])[kb]=f7_B;// - feq7_B + feqW7_T;
      //}

      //q = q_dirB[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[5])[kt]=f7_T;// - feq7_T + feqW7_B;
      //}


      //////////////////////////////////////////////////////////////////////////
      //ohne Q's aber mit TRT
      doubflo /*feq,*/ q;
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         (D7.f[2])[kw]= -f7_W + cs2 * TempD;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         (D7.f[1])[ke]= -f7_E + cs2 * TempD;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         (D7.f[4])[ks]= -f7_S + cs2 * TempD;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         (D7.f[3])[kn]= -f7_N + cs2 * TempD;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         (D7.f[6])[kb]= -f7_B + cs2 * TempD;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         (D7.f[5])[kt]= -f7_T + cs2 * TempD;
      }

   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QNoSlipADincomp27(int inx,
											 int iny,
											 doubflo* DD, 
											 doubflo* DD27, 
											 doubflo* temp,
											 doubflo diffusivity,
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

   Distributions27 D27;
   if (evenOrOdd==true)
   {
      D27.f[dirE   ] = &DD27[dirE   *size_Mat];
      D27.f[dirW   ] = &DD27[dirW   *size_Mat];
      D27.f[dirN   ] = &DD27[dirN   *size_Mat];
      D27.f[dirS   ] = &DD27[dirS   *size_Mat];
      D27.f[dirT   ] = &DD27[dirT   *size_Mat];
      D27.f[dirB   ] = &DD27[dirB   *size_Mat];
      D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
      D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
      D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
      D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
      D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
      D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
   } 
   else
   {
      D27.f[dirW   ] = &DD27[dirE   *size_Mat];
      D27.f[dirE   ] = &DD27[dirW   *size_Mat];
      D27.f[dirS   ] = &DD27[dirN   *size_Mat];
      D27.f[dirN   ] = &DD27[dirS   *size_Mat];
      D27.f[dirB   ] = &DD27[dirT   *size_Mat];
      D27.f[dirT   ] = &DD27[dirB   *size_Mat];
      D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
      D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
      D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
      D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
      D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
      D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
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
      doubflo  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
      doubflo f_ZERO = (D.f[dirZERO])[kzero];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1 =  ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      doubflo vx2 =  ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      doubflo vx3 =  ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      doubflo f27_W    = (D27.f[dirE   ])[ke   ];
      doubflo f27_E    = (D27.f[dirW   ])[kw   ];
      doubflo f27_S    = (D27.f[dirN   ])[kn   ];
      doubflo f27_N    = (D27.f[dirS   ])[ks   ];
      doubflo f27_B    = (D27.f[dirT   ])[kt   ];
      doubflo f27_T    = (D27.f[dirB   ])[kb   ];
      doubflo f27_SW   = (D27.f[dirNE  ])[kne  ];
      doubflo f27_NE   = (D27.f[dirSW  ])[ksw  ];
      doubflo f27_NW   = (D27.f[dirSE  ])[kse  ];
      doubflo f27_SE   = (D27.f[dirNW  ])[knw  ];
      doubflo f27_BW   = (D27.f[dirTE  ])[kte  ];
      doubflo f27_TE   = (D27.f[dirBW  ])[kbw  ];
      doubflo f27_TW   = (D27.f[dirBE  ])[kbe  ];
      doubflo f27_BE   = (D27.f[dirTW  ])[ktw  ];
      doubflo f27_BS   = (D27.f[dirTN  ])[ktn  ];
      doubflo f27_TN   = (D27.f[dirBS  ])[kbs  ];
      doubflo f27_TS   = (D27.f[dirBN  ])[kbn  ];
      doubflo f27_BN   = (D27.f[dirTS  ])[kts  ];
      doubflo f27_ZERO = (D27.f[dirZERO])[kzero];
      doubflo f27_BSW  = (D27.f[dirTNE ])[ktne ];
      doubflo f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      doubflo f27_BNW  = (D27.f[dirTSE ])[ktse ];
      doubflo f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      doubflo f27_TSW  = (D27.f[dirBNE ])[kbne ];
      doubflo f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      doubflo f27_TNW  = (D27.f[dirBSE ])[kbse ];
      doubflo f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      doubflo ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
         f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
         f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      //doubflo feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      doubflo feq27_E    =   c2over27* ConcD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feq27_W    =   c2over27* ConcD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feq27_N    =   c2over27* ConcD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feq27_S    =   c2over27* ConcD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feq27_T    =   c2over27* ConcD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feq27_B    =   c2over27* ConcD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feq27_NE   =   c1over54* ConcD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feq27_SW   =   c1over54* ConcD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feq27_SE   =   c1over54* ConcD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feq27_NW   =   c1over54* ConcD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feq27_TE   =   c1over54* ConcD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feq27_BW   =   c1over54* ConcD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feq27_BE   =   c1over54* ConcD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feq27_TW   =   c1over54* ConcD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feq27_TN   =   c1over54* ConcD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feq27_BS   =   c1over54* ConcD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feq27_BN   =   c1over54* ConcD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feq27_TS   =   c1over54* ConcD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feq27_TNE  =   c1over216*ConcD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feq27_BSW  =   c1over216*ConcD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feq27_BNE  =   c1over216*ConcD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feq27_TSW  =   c1over216*ConcD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feq27_TSE  =   c1over216*ConcD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feq27_BNW  =   c1over216*ConcD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feq27_BSE  =   c1over216*ConcD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feq27_TNW  =   c1over216*ConcD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo TempD = temp[k];

      //doubflo feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      doubflo feqW27_E    =   c2over27* TempD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feqW27_W    =   c2over27* TempD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feqW27_N    =   c2over27* TempD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feqW27_S    =   c2over27* TempD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feqW27_T    =   c2over27* TempD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feqW27_B    =   c2over27* TempD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feqW27_NE   =   c1over54* TempD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feqW27_SW   =   c1over54* TempD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feqW27_SE   =   c1over54* TempD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feqW27_NW   =   c1over54* TempD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feqW27_TE   =   c1over54* TempD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feqW27_BW   =   c1over54* TempD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feqW27_BE   =   c1over54* TempD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feqW27_TW   =   c1over54* TempD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feqW27_TN   =   c1over54* TempD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feqW27_BS   =   c1over54* TempD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feqW27_BN   =   c1over54* TempD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feqW27_TS   =   c1over54* TempD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feqW27_TNE  =   c1over216*TempD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feqW27_BSW  =   c1over216*TempD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feqW27_BNE  =   c1over216*TempD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feqW27_TSW  =   c1over216*TempD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feqW27_TSE  =   c1over216*TempD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feqW27_BNW  =   c1over216*TempD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feqW27_BSE  =   c1over216*TempD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feqW27_TNW  =   c1over216*TempD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo omegaD     = three - sqrt(three);
      //doubflo Lam        = -(c1o2-one/omegaD);
      //doubflo nue_d      = Lam/three;
      //doubflo ae         = zero;
      //doubflo ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
      {
         D27.f[dirE   ] = &DD27[dirE   *size_Mat];
         D27.f[dirW   ] = &DD27[dirW   *size_Mat];
         D27.f[dirN   ] = &DD27[dirN   *size_Mat];
         D27.f[dirS   ] = &DD27[dirS   *size_Mat];
         D27.f[dirT   ] = &DD27[dirT   *size_Mat];
         D27.f[dirB   ] = &DD27[dirB   *size_Mat];
         D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
         D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
         D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
         D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
         D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
         D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
      } 
      else
      {
         D27.f[dirW   ] = &DD27[dirE   *size_Mat];
         D27.f[dirE   ] = &DD27[dirW   *size_Mat];
         D27.f[dirS   ] = &DD27[dirN   *size_Mat];
         D27.f[dirN   ] = &DD27[dirS   *size_Mat];
         D27.f[dirB   ] = &DD27[dirT   *size_Mat];
         D27.f[dirT   ] = &DD27[dirB   *size_Mat];
         D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
         D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
         D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
         D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
         D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
         D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=0.1f;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  doubflo q;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[dirW  ])[kw  ]=(two*feqW27_W  -(f27_E  *(q*omegaD-one)-omegaD*feq27_E  *(q-one))/(omegaD-one)+f27_W  *q)/(q+one);
      q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[dirE  ])[ke  ]=(two*feqW27_E  -(f27_W  *(q*omegaD-one)-omegaD*feq27_W  *(q-one))/(omegaD-one)+f27_E  *q)/(q+one);
      q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[dirS  ])[ks  ]=(two*feqW27_S  -(f27_N  *(q*omegaD-one)-omegaD*feq27_N  *(q-one))/(omegaD-one)+f27_S  *q)/(q+one);
      q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[dirN  ])[kn  ]=(two*feqW27_N  -(f27_S  *(q*omegaD-one)-omegaD*feq27_S  *(q-one))/(omegaD-one)+f27_N  *q)/(q+one);
      q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[dirB  ])[kb  ]=(two*feqW27_B  -(f27_T  *(q*omegaD-one)-omegaD*feq27_T  *(q-one))/(omegaD-one)+f27_B  *q)/(q+one);
      q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[dirT  ])[kt  ]=(two*feqW27_T  -(f27_B  *(q*omegaD-one)-omegaD*feq27_B  *(q-one))/(omegaD-one)+f27_T  *q)/(q+one);
      q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[dirSW ])[ksw ]=(two*feqW27_SW -(f27_NE *(q*omegaD-one)-omegaD*feq27_NE *(q-one))/(omegaD-one)+f27_SW *q)/(q+one);
      q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[dirNE ])[kne ]=(two*feqW27_NE -(f27_SW *(q*omegaD-one)-omegaD*feq27_SW *(q-one))/(omegaD-one)+f27_NE *q)/(q+one);
      q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[dirNW ])[knw ]=(two*feqW27_NW -(f27_SE *(q*omegaD-one)-omegaD*feq27_SE *(q-one))/(omegaD-one)+f27_NW *q)/(q+one);
      q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[dirSE ])[kse ]=(two*feqW27_SE -(f27_NW *(q*omegaD-one)-omegaD*feq27_NW *(q-one))/(omegaD-one)+f27_SE *q)/(q+one);
      q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[dirBW ])[kbw ]=(two*feqW27_BW -(f27_TE *(q*omegaD-one)-omegaD*feq27_TE *(q-one))/(omegaD-one)+f27_BW *q)/(q+one);
      q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[dirTE ])[kte ]=(two*feqW27_TE -(f27_BW *(q*omegaD-one)-omegaD*feq27_BW *(q-one))/(omegaD-one)+f27_TE *q)/(q+one);
      q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[dirTW ])[ktw ]=(two*feqW27_TW -(f27_BE *(q*omegaD-one)-omegaD*feq27_BE *(q-one))/(omegaD-one)+f27_TW *q)/(q+one);
      q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[dirBE ])[kbe ]=(two*feqW27_BE -(f27_TW *(q*omegaD-one)-omegaD*feq27_TW *(q-one))/(omegaD-one)+f27_BE *q)/(q+one);
      q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[dirBS ])[kbs ]=(two*feqW27_BS -(f27_TN *(q*omegaD-one)-omegaD*feq27_TN *(q-one))/(omegaD-one)+f27_BS *q)/(q+one);
      q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[dirTN ])[ktn ]=(two*feqW27_TN -(f27_BS *(q*omegaD-one)-omegaD*feq27_BS *(q-one))/(omegaD-one)+f27_TN *q)/(q+one);
      q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[dirTS ])[kts ]=(two*feqW27_TS -(f27_BN *(q*omegaD-one)-omegaD*feq27_BN *(q-one))/(omegaD-one)+f27_TS *q)/(q+one);
      q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[dirBN ])[kbn ]=(two*feqW27_BN -(f27_TS *(q*omegaD-one)-omegaD*feq27_TS *(q-one))/(omegaD-one)+f27_BN *q)/(q+one);
      q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[dirBSW])[kbsw]=(two*feqW27_BSW-(f27_TNE*(q*omegaD-one)-omegaD*feq27_TNE*(q-one))/(omegaD-one)+f27_BSW*q)/(q+one);
      q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[dirTNE])[ktne]=(two*feqW27_TNE-(f27_BSW*(q*omegaD-one)-omegaD*feq27_BSW*(q-one))/(omegaD-one)+f27_TNE*q)/(q+one);
      q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[dirTSW])[ktsw]=(two*feqW27_TSW-(f27_BNE*(q*omegaD-one)-omegaD*feq27_BNE*(q-one))/(omegaD-one)+f27_TSW*q)/(q+one);
      q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[dirBNE])[kbne]=(two*feqW27_BNE-(f27_TSW*(q*omegaD-one)-omegaD*feq27_TSW*(q-one))/(omegaD-one)+f27_BNE*q)/(q+one);
      q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[dirBNW])[kbnw]=(two*feqW27_BNW-(f27_TSE*(q*omegaD-one)-omegaD*feq27_TSE*(q-one))/(omegaD-one)+f27_BNW*q)/(q+one);
      q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[dirTSE])[ktse]=(two*feqW27_TSE-(f27_BNW*(q*omegaD-one)-omegaD*feq27_BNW*(q-one))/(omegaD-one)+f27_TSE*q)/(q+one);
      q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[dirTNW])[ktnw]=(two*feqW27_TNW-(f27_BSE*(q*omegaD-one)-omegaD*feq27_BSE*(q-one))/(omegaD-one)+f27_TNW*q)/(q+one);
      q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[dirBSE])[kbse]=(two*feqW27_BSE-(f27_TNW*(q*omegaD-one)-omegaD*feq27_TNW*(q-one))/(omegaD-one)+f27_BSE*q)/(q+one);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADVeloIncomp7(  int inx,
											int iny,
											doubflo* DD, 
											doubflo* DD7, 
											doubflo* temp,
											doubflo* velo,
											doubflo diffusivity,
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

   Distributions7 D7;
   if (evenOrOdd==true)
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[1] = &DD7[1*size_Mat];
      D7.f[2] = &DD7[2*size_Mat];
      D7.f[3] = &DD7[3*size_Mat];
      D7.f[4] = &DD7[4*size_Mat];
      D7.f[5] = &DD7[5*size_Mat];
      D7.f[6] = &DD7[6*size_Mat];
   }
   else
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[2] = &DD7[1*size_Mat];
      D7.f[1] = &DD7[2*size_Mat];
      D7.f[4] = &DD7[3*size_Mat];
      D7.f[3] = &DD7[4*size_Mat];
      D7.f[6] = &DD7[5*size_Mat];
      D7.f[5] = &DD7[6*size_Mat];
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
      //////////////////////////////////////////////////////////////////////////////////
      doubflo  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB; 

      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      //////////////////////////////////////////////////////////////////////////////////
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
      doubflo vx1_Inflow   = zero;
      doubflo vx2_Inflow   = velo[k];
      doubflo vx3_Inflow   = zero;
      doubflo ux_sq_Inflow = vx1_Inflow * vx1_Inflow;
      doubflo uy_sq_Inflow = vx2_Inflow * vx2_Inflow;
      doubflo uz_sq_Inflow = vx3_Inflow * vx3_Inflow;
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      doubflo vx2 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      doubflo vx3 = ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
		 ////dörrrrrty !!!!!!!!!!!!!
   //      doubflo vx1     =  ten * ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
   //      doubflo vx2     =  ten * ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
   //      doubflo vx3     =  ten * ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
	  //doubflo cu_sq       =1.5f*(vx1*vx1+vx2*vx2+vx3*vx3);
      doubflo ux_sq       = vx1 * vx1;
      doubflo uy_sq       = vx2 * vx2;
      doubflo uz_sq       = vx3 * vx3;

      doubflo f7_ZERO,f7_E,f7_W,f7_N,f7_S,f7_T,f7_B;
      doubflo /*feq7_ZERO,*/feq7_E,feq7_W,feq7_N,feq7_S,feq7_T,feq7_B;
      doubflo /*feqW7_ZERO,*/feqW7_E,feqW7_W,feqW7_N,feqW7_S,feqW7_T,feqW7_B;
      doubflo TempD = temp[k];

      f7_ZERO =  (D7.f[0])[kzero];
      f7_W    =  (D7.f[1])[ke   ];
      f7_E    =  (D7.f[2])[kw   ];
      f7_S    =  (D7.f[3])[kn   ];
      f7_N    =  (D7.f[4])[ks   ];
      f7_B    =  (D7.f[5])[kt   ];
      f7_T    =  (D7.f[6])[kb   ];

      doubflo ConcD = f7_ZERO + f7_E + f7_W + f7_N + f7_S + f7_T + f7_B;

      ////////////////////////////////////////////////////////////////////////////////
	  //BGK
      //doubflo omegaD     = three - sqrt(three);
      //doubflo Lam         = -(c1o2-one/omegaD);
      //doubflo nue_d       = Lam/three;
      ////doubflo ae          = zero;
      //doubflo ae          = diffusivity/nue_d - one;

      ////feq7_ZERO = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
      //feq7_E    = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
      //feq7_W    = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
      //feq7_N    = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
      //feq7_S    = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
      //feq7_T    = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
      //feq7_B    = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);

      ////feq7_ZERO = TempD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
      //feqW7_E    = TempD*(c1o6*(ae+one)+c1o2*(ux_sq_Inflow)+vx1_Inflow*c1o2);
      //feqW7_W    = TempD*(c1o6*(ae+one)+c1o2*(ux_sq_Inflow)-vx1_Inflow*c1o2);
      //feqW7_N    = TempD*(c1o6*(ae+one)+c1o2*(uy_sq_Inflow)+vx2_Inflow*c1o2);
      //feqW7_S    = TempD*(c1o6*(ae+one)+c1o2*(uy_sq_Inflow)-vx2_Inflow*c1o2);
      //feqW7_T    = TempD*(c1o6*(ae+one)+c1o2*(uz_sq_Inflow)+vx3_Inflow*c1o2);
      //feqW7_B    = TempD*(c1o6*(ae+one)+c1o2*(uz_sq_Inflow)-vx3_Inflow*c1o2);

   //      ////////////////////////////////////////////////////////////////////////////////
		 ////TRT  Yoshida Kernel - based on Ying
         doubflo cs2         = c1o4;
   //      doubflo Lam         = diffusivity/(one)/cs2;
   //      doubflo omegaD      = - one / (Lam + c1o2);
   //      doubflo ae          = zero;
   //      ////////////////////////////////////////////////////////////////////////////////
		 //doubflo Mom000 = f7_ZERO + f7_W + f7_E + f7_N + f7_S + f7_T + f7_B; //1
   //      doubflo Mom100 = f7_E - f7_W;
   //      doubflo Mom010 = f7_N - f7_S;
   //      doubflo Mom001 = f7_T - f7_B;
   //      doubflo Mom222 = six*f7_ZERO - f7_W - f7_E - f7_N - f7_S - f7_T - f7_B;
   //      doubflo Mom200 = two*f7_W + two*f7_E - f7_N - f7_S - f7_T - f7_B;
   //      doubflo Mom022 = f7_N + f7_S - f7_T - f7_B;

   //      doubflo Meq000 = ConcD;
   //      doubflo Meq100 = ConcD*vx1;
   //      doubflo Meq010 = ConcD*vx2;
   //      doubflo Meq001 = ConcD*vx3;
   //      doubflo Meq222 = c3o4*ConcD;
   //      doubflo Meq200 = zero;
   //      doubflo Meq022 = zero;

   //      // relaxation TRT Yoshida

   //      // odd 
   //      Mom100 = omegaD * (Mom100-Meq100);
   //      Mom010 = omegaD * (Mom010-Meq010);
   //      Mom001 = omegaD * (Mom001-Meq001);
   //      
   //      // even
   //      Mom000 = -one*(Mom000-Meq000);
   //      Mom222 = -one*(Mom222-Meq222);
   //      Mom200 = -one*(Mom200-Meq200);
   //      Mom022 = -one*(Mom022-Meq022);
   //      
   //      //Back transformation to distributions
   //      f7_ZERO = f7_ZERO + c1o7*Mom000 + c1o7*Mom222;                                                  //1
   //      f7_E    = f7_E    + c1o7*Mom000 + c1o2*Mom100 - c1o6*c1o7*Mom222 + c1o6*Mom200;                 //2
   //      f7_W    = f7_W    + c1o7*Mom000 - c1o2*Mom100 - c1o6*c1o7*Mom222 + c1o6*Mom200;                 //3
   //      f7_N    = f7_N    + c1o7*Mom000 + c1o2*Mom010 - c1o6*c1o7*Mom222 - c1o12*Mom200 + c1o4 *Mom022; //4
   //      f7_S    = f7_S    + c1o7*Mom000 - c1o2*Mom010 - c1o6*c1o7*Mom222 - c1o12*Mom200 + c1o4 *Mom022; //5
   //      f7_T    = f7_T    + c1o7*Mom000 + c1o2*Mom001 - c1o6*c1o7*Mom222 - c1o12*Mom200 - c1o4 *Mom022; //6
   //      f7_B    = f7_B    + c1o7*Mom000 - c1o2*Mom001 - c1o6*c1o7*Mom222 - c1o12*Mom200 - c1o4 *Mom022; //7





      //////////////////////////////////////////////////////////////////////////
      //pointertausch
      if (evenOrOdd==false)
      {
         D7.f[0] = &DD7[0*size_Mat];
         D7.f[1] = &DD7[1*size_Mat];
         D7.f[2] = &DD7[2*size_Mat];
         D7.f[3] = &DD7[3*size_Mat];
         D7.f[4] = &DD7[4*size_Mat];
         D7.f[5] = &DD7[5*size_Mat];
         D7.f[6] = &DD7[6*size_Mat];
      }
      else
      {
         D7.f[0] = &DD7[0*size_Mat];
         D7.f[2] = &DD7[1*size_Mat];
         D7.f[1] = &DD7[2*size_Mat];
         D7.f[4] = &DD7[3*size_Mat];
         D7.f[3] = &DD7[4*size_Mat];
         D7.f[6] = &DD7[5*size_Mat];
         D7.f[5] = &DD7[6*size_Mat];
      }

      ////////////////////////////////////////////////////////////////////////////
      ////mit Q's
      //doubflo /*feq,*/ q;
      //q = q_dirE[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=0.;
      //   (D7.f[2])[kw]=(two*feqW7_W-(f7_E*(q*omegaD-one)-omegaD*feq7_E*(q-one))/(omegaD-one)+f7_W*q)/(q+one);//f7_W - feq7_W + feqW7_E;
      //}

      //q = q_dirW[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=0.;
      //   (D7.f[1])[ke]=(two*feqW7_E-(f7_W*(q*omegaD-one)-omegaD*feq7_W*(q-one))/(omegaD-one)+f7_E*q)/(q+one);//f7_E - feq7_E + feqW7_W;
      //}

      //q = q_dirN[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=0.;
      //   (D7.f[4])[ks]=(two*feqW7_S-(f7_N*(q*omegaD-one)-omegaD*feq7_N*(q-one))/(omegaD-one)+f7_S*q)/(q+one);//f7_S - feq7_S + feqW7_N;
      //}

      //q = q_dirS[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=0.;
      //   (D7.f[3])[kn]=(two*feqW7_N-(f7_S*(q*omegaD-one)-omegaD*feq7_S*(q-one))/(omegaD-one)+f7_N*q)/(q+one);//f7_N - feq7_N + feqW7_S;
      //}

      //q = q_dirT[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=0.;
      //   (D7.f[6])[kb]=(two*feqW7_B-(f7_T*(q*omegaD-one)-omegaD*feq7_T*(q-one))/(omegaD-one)+f7_B*q)/(q+one);//f7_B - feq7_B + feqW7_T;
      //}

      //q = q_dirB[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=0.;
      //   (D7.f[5])[kt]=(two*feqW7_T-(f7_B*(q*omegaD-one)-omegaD*feq7_B*(q-one))/(omegaD-one)+f7_T*q)/(q+one);//f7_T - feq7_T + feqW7_B;
      //}

      ////////////////////////////////////////////////////////////////////////////
      ////ohne Q's
      //doubflo /*feq,*/ q;
      //q = q_dirE[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[2])[kw]=f7_W - feq7_W + feqW7_E;
      //}

      //q = q_dirW[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[1])[ke]=f7_E - feq7_E + feqW7_W;
      //}

      //q = q_dirN[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[4])[ks]=f7_S - feq7_S + feqW7_N;
      //}

      //q = q_dirS[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[3])[kn]=f7_N - feq7_N + feqW7_S;
      //}

      //q = q_dirT[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[6])[kb]=f7_B - feq7_B + feqW7_T;
      //}

      //q = q_dirB[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[5])[kt]=f7_T - feq7_T + feqW7_B;
      //}

      //////////////////////////////////////////////////////////////////////////
      //ohne Q's aber mit TRT
      doubflo /*feq,*/ q;
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         (D7.f[2])[kw]= -f7_W + cs2 * TempD;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         (D7.f[1])[ke]= -f7_E + cs2 * TempD;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         (D7.f[4])[ks]= -f7_S + cs2 * TempD;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         (D7.f[3])[kn]= -f7_N + cs2 * TempD;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         (D7.f[6])[kb]= -f7_B + cs2 * TempD;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         (D7.f[5])[kt]= -f7_T + cs2 * TempD;
      }

   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADVeloIncomp27( int inx,
											int iny,
											doubflo* DD, 
											doubflo* DD27, 
											doubflo* temp,
											doubflo* velo,
											doubflo diffusivity,
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

   Distributions27 D27;
   if (evenOrOdd==true)
   {
      D27.f[dirE   ] = &DD27[dirE   *size_Mat];
      D27.f[dirW   ] = &DD27[dirW   *size_Mat];
      D27.f[dirN   ] = &DD27[dirN   *size_Mat];
      D27.f[dirS   ] = &DD27[dirS   *size_Mat];
      D27.f[dirT   ] = &DD27[dirT   *size_Mat];
      D27.f[dirB   ] = &DD27[dirB   *size_Mat];
      D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
      D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
      D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
      D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
      D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
      D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
   } 
   else
   {
      D27.f[dirW   ] = &DD27[dirE   *size_Mat];
      D27.f[dirE   ] = &DD27[dirW   *size_Mat];
      D27.f[dirS   ] = &DD27[dirN   *size_Mat];
      D27.f[dirN   ] = &DD27[dirS   *size_Mat];
      D27.f[dirB   ] = &DD27[dirT   *size_Mat];
      D27.f[dirT   ] = &DD27[dirB   *size_Mat];
      D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
      D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
      D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
      D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
      D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
      D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
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
      doubflo  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
      doubflo f_ZERO = (D.f[dirZERO])[kzero];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      doubflo vx2 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      doubflo vx3 = ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      doubflo f27_W    = (D27.f[dirE   ])[ke   ];
      doubflo f27_E    = (D27.f[dirW   ])[kw   ];
      doubflo f27_S    = (D27.f[dirN   ])[kn   ];
      doubflo f27_N    = (D27.f[dirS   ])[ks   ];
      doubflo f27_B    = (D27.f[dirT   ])[kt   ];
      doubflo f27_T    = (D27.f[dirB   ])[kb   ];
      doubflo f27_SW   = (D27.f[dirNE  ])[kne  ];
      doubflo f27_NE   = (D27.f[dirSW  ])[ksw  ];
      doubflo f27_NW   = (D27.f[dirSE  ])[kse  ];
      doubflo f27_SE   = (D27.f[dirNW  ])[knw  ];
      doubflo f27_BW   = (D27.f[dirTE  ])[kte  ];
      doubflo f27_TE   = (D27.f[dirBW  ])[kbw  ];
      doubflo f27_TW   = (D27.f[dirBE  ])[kbe  ];
      doubflo f27_BE   = (D27.f[dirTW  ])[ktw  ];
      doubflo f27_BS   = (D27.f[dirTN  ])[ktn  ];
      doubflo f27_TN   = (D27.f[dirBS  ])[kbs  ];
      doubflo f27_TS   = (D27.f[dirBN  ])[kbn  ];
      doubflo f27_BN   = (D27.f[dirTS  ])[kts  ];
      doubflo f27_ZERO = (D27.f[dirZERO])[kzero];
      doubflo f27_BSW  = (D27.f[dirTNE ])[ktne ];
      doubflo f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      doubflo f27_BNW  = (D27.f[dirTSE ])[ktse ];
      doubflo f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      doubflo f27_TSW  = (D27.f[dirBNE ])[kbne ];
      doubflo f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      doubflo f27_TNW  = (D27.f[dirBSE ])[kbse ];
      doubflo f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      doubflo ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
                        f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
                        f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      //doubflo feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      doubflo feq27_E    =   c2over27* ConcD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feq27_W    =   c2over27* ConcD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feq27_N    =   c2over27* ConcD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feq27_S    =   c2over27* ConcD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feq27_T    =   c2over27* ConcD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feq27_B    =   c2over27* ConcD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feq27_NE   =   c1over54* ConcD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feq27_SW   =   c1over54* ConcD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feq27_SE   =   c1over54* ConcD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feq27_NW   =   c1over54* ConcD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feq27_TE   =   c1over54* ConcD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feq27_BW   =   c1over54* ConcD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feq27_BE   =   c1over54* ConcD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feq27_TW   =   c1over54* ConcD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feq27_TN   =   c1over54* ConcD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feq27_BS   =   c1over54* ConcD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feq27_BN   =   c1over54* ConcD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feq27_TS   =   c1over54* ConcD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feq27_TNE  =   c1over216*ConcD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feq27_BSW  =   c1over216*ConcD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feq27_BNE  =   c1over216*ConcD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feq27_TSW  =   c1over216*ConcD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feq27_TSE  =   c1over216*ConcD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feq27_BNW  =   c1over216*ConcD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feq27_BSE  =   c1over216*ConcD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feq27_TNW  =   c1over216*ConcD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo TempD = temp[k];
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // velocity inflow
      vx1   = velo[k];//zero;
      vx2   = zero;//velo[k];
      vx3   = zero;

      //doubflo feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      doubflo feqW27_E    =   c2over27* TempD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feqW27_W    =   c2over27* TempD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feqW27_N    =   c2over27* TempD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feqW27_S    =   c2over27* TempD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feqW27_T    =   c2over27* TempD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feqW27_B    =   c2over27* TempD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feqW27_NE   =   c1over54* TempD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feqW27_SW   =   c1over54* TempD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feqW27_SE   =   c1over54* TempD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feqW27_NW   =   c1over54* TempD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feqW27_TE   =   c1over54* TempD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feqW27_BW   =   c1over54* TempD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feqW27_BE   =   c1over54* TempD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feqW27_TW   =   c1over54* TempD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feqW27_TN   =   c1over54* TempD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feqW27_BS   =   c1over54* TempD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feqW27_BN   =   c1over54* TempD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feqW27_TS   =   c1over54* TempD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feqW27_TNE  =   c1over216*TempD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feqW27_BSW  =   c1over216*TempD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feqW27_BNE  =   c1over216*TempD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feqW27_TSW  =   c1over216*TempD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feqW27_TSE  =   c1over216*TempD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feqW27_BNW  =   c1over216*TempD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feqW27_BSE  =   c1over216*TempD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feqW27_TNW  =   c1over216*TempD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo omegaD     = three - sqrt(three);
      //doubflo Lam        = -(c1o2-one/omegaD);
      //doubflo nue_d      = Lam/three;
      //doubflo ae         = zero;
      //doubflo ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
      {
         D27.f[dirE   ] = &DD27[dirE   *size_Mat];
         D27.f[dirW   ] = &DD27[dirW   *size_Mat];
         D27.f[dirN   ] = &DD27[dirN   *size_Mat];
         D27.f[dirS   ] = &DD27[dirS   *size_Mat];
         D27.f[dirT   ] = &DD27[dirT   *size_Mat];
         D27.f[dirB   ] = &DD27[dirB   *size_Mat];
         D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
         D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
         D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
         D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
         D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
         D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
      } 
      else
      {
         D27.f[dirW   ] = &DD27[dirE   *size_Mat];
         D27.f[dirE   ] = &DD27[dirW   *size_Mat];
         D27.f[dirS   ] = &DD27[dirN   *size_Mat];
         D27.f[dirN   ] = &DD27[dirS   *size_Mat];
         D27.f[dirB   ] = &DD27[dirT   *size_Mat];
         D27.f[dirT   ] = &DD27[dirB   *size_Mat];
         D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
         D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
         D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
         D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
         D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
         D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo q;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[dirW  ])[kw  ]= -feqW27_W  + two * c2over27  * TempD;
      q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[dirE  ])[ke  ]= -feqW27_E  + two * c2over27  * TempD;
      q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[dirS  ])[ks  ]= -feqW27_S  + two * c2over27  * TempD;
      q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[dirN  ])[kn  ]= -feqW27_N  + two * c2over27  * TempD;
      q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[dirB  ])[kb  ]= -feqW27_B  + two * c2over27  * TempD;
      q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[dirT  ])[kt  ]= -feqW27_T  + two * c2over27  * TempD;
      q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[dirSW ])[ksw ]= -feqW27_SW + two * c1over54  * TempD;
      q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[dirNE ])[kne ]= -feqW27_NE + two * c1over54  * TempD;
      q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[dirNW ])[knw ]= -feqW27_NW + two * c1over54  * TempD;
      q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[dirSE ])[kse ]= -feqW27_SE + two * c1over54  * TempD;
      q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[dirBW ])[kbw ]= -feqW27_BW + two * c1over54  * TempD;
      q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[dirTE ])[kte ]= -feqW27_TE + two * c1over54  * TempD;
      q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[dirTW ])[ktw ]= -feqW27_TW + two * c1over54  * TempD;
      q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[dirBE ])[kbe ]= -feqW27_BE + two * c1over54  * TempD;
      q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[dirBS ])[kbs ]= -feqW27_BS + two * c1over54  * TempD;
      q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[dirTN ])[ktn ]= -feqW27_TN + two * c1over54  * TempD;
      q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[dirTS ])[kts ]= -feqW27_TS + two * c1over54  * TempD;
      q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[dirBN ])[kbn ]= -feqW27_BN + two * c1over54  * TempD;
      q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[dirBSW])[kbsw]= -feqW27_BSW+ two * c1over216 * TempD;
      q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[dirTNE])[ktne]= -feqW27_TNE+ two * c1over216 * TempD;
      q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[dirTSW])[ktsw]= -feqW27_TSW+ two * c1over216 * TempD;
      q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[dirBNE])[kbne]= -feqW27_BNE+ two * c1over216 * TempD;
      q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[dirBNW])[kbnw]= -feqW27_BNW+ two * c1over216 * TempD;
      q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[dirTSE])[ktse]= -feqW27_TSE+ two * c1over216 * TempD;
      q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[dirTNW])[ktnw]= -feqW27_TNW+ two * c1over216 * TempD;
      q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[dirBSE])[kbse]= -feqW27_BSE+ two * c1over216 * TempD;
      //q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[dirW  ])[kw  ]=(two*feqW27_W  -(f27_E  *(q*omegaD-one)-omegaD*feq27_E  *(q-one))/(omegaD-one)+f27_W  *q)/(q+one);
      //q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[dirE  ])[ke  ]=(two*feqW27_E  -(f27_W  *(q*omegaD-one)-omegaD*feq27_W  *(q-one))/(omegaD-one)+f27_E  *q)/(q+one);
      //q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[dirS  ])[ks  ]=(two*feqW27_S  -(f27_N  *(q*omegaD-one)-omegaD*feq27_N  *(q-one))/(omegaD-one)+f27_S  *q)/(q+one);
      //q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[dirN  ])[kn  ]=(two*feqW27_N  -(f27_S  *(q*omegaD-one)-omegaD*feq27_S  *(q-one))/(omegaD-one)+f27_N  *q)/(q+one);
      //q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[dirB  ])[kb  ]=(two*feqW27_B  -(f27_T  *(q*omegaD-one)-omegaD*feq27_T  *(q-one))/(omegaD-one)+f27_B  *q)/(q+one);
      //q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[dirT  ])[kt  ]=(two*feqW27_T  -(f27_B  *(q*omegaD-one)-omegaD*feq27_B  *(q-one))/(omegaD-one)+f27_T  *q)/(q+one);
      //q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[dirSW ])[ksw ]=(two*feqW27_SW -(f27_NE *(q*omegaD-one)-omegaD*feq27_NE *(q-one))/(omegaD-one)+f27_SW *q)/(q+one);
      //q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[dirNE ])[kne ]=(two*feqW27_NE -(f27_SW *(q*omegaD-one)-omegaD*feq27_SW *(q-one))/(omegaD-one)+f27_NE *q)/(q+one);
      //q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[dirNW ])[knw ]=(two*feqW27_NW -(f27_SE *(q*omegaD-one)-omegaD*feq27_SE *(q-one))/(omegaD-one)+f27_NW *q)/(q+one);
      //q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[dirSE ])[kse ]=(two*feqW27_SE -(f27_NW *(q*omegaD-one)-omegaD*feq27_NW *(q-one))/(omegaD-one)+f27_SE *q)/(q+one);
      //q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[dirBW ])[kbw ]=(two*feqW27_BW -(f27_TE *(q*omegaD-one)-omegaD*feq27_TE *(q-one))/(omegaD-one)+f27_BW *q)/(q+one);
      //q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[dirTE ])[kte ]=(two*feqW27_TE -(f27_BW *(q*omegaD-one)-omegaD*feq27_BW *(q-one))/(omegaD-one)+f27_TE *q)/(q+one);
      //q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[dirTW ])[ktw ]=(two*feqW27_TW -(f27_BE *(q*omegaD-one)-omegaD*feq27_BE *(q-one))/(omegaD-one)+f27_TW *q)/(q+one);
      //q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[dirBE ])[kbe ]=(two*feqW27_BE -(f27_TW *(q*omegaD-one)-omegaD*feq27_TW *(q-one))/(omegaD-one)+f27_BE *q)/(q+one);
      //q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[dirBS ])[kbs ]=(two*feqW27_BS -(f27_TN *(q*omegaD-one)-omegaD*feq27_TN *(q-one))/(omegaD-one)+f27_BS *q)/(q+one);
      //q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[dirTN ])[ktn ]=(two*feqW27_TN -(f27_BS *(q*omegaD-one)-omegaD*feq27_BS *(q-one))/(omegaD-one)+f27_TN *q)/(q+one);
      //q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[dirTS ])[kts ]=(two*feqW27_TS -(f27_BN *(q*omegaD-one)-omegaD*feq27_BN *(q-one))/(omegaD-one)+f27_TS *q)/(q+one);
      //q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[dirBN ])[kbn ]=(two*feqW27_BN -(f27_TS *(q*omegaD-one)-omegaD*feq27_TS *(q-one))/(omegaD-one)+f27_BN *q)/(q+one);
      //q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[dirBSW])[kbsw]=(two*feqW27_BSW-(f27_TNE*(q*omegaD-one)-omegaD*feq27_TNE*(q-one))/(omegaD-one)+f27_BSW*q)/(q+one);
      //q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[dirTNE])[ktne]=(two*feqW27_TNE-(f27_BSW*(q*omegaD-one)-omegaD*feq27_BSW*(q-one))/(omegaD-one)+f27_TNE*q)/(q+one);
      //q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[dirTSW])[ktsw]=(two*feqW27_TSW-(f27_BNE*(q*omegaD-one)-omegaD*feq27_BNE*(q-one))/(omegaD-one)+f27_TSW*q)/(q+one);
      //q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[dirBNE])[kbne]=(two*feqW27_BNE-(f27_TSW*(q*omegaD-one)-omegaD*feq27_TSW*(q-one))/(omegaD-one)+f27_BNE*q)/(q+one);
      //q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[dirBNW])[kbnw]=(two*feqW27_BNW-(f27_TSE*(q*omegaD-one)-omegaD*feq27_TSE*(q-one))/(omegaD-one)+f27_BNW*q)/(q+one);
      //q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[dirTSE])[ktse]=(two*feqW27_TSE-(f27_BNW*(q*omegaD-one)-omegaD*feq27_BNW*(q-one))/(omegaD-one)+f27_TSE*q)/(q+one);
      //q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[dirTNW])[ktnw]=(two*feqW27_TNW-(f27_BSE*(q*omegaD-one)-omegaD*feq27_BSE*(q-one))/(omegaD-one)+f27_TNW*q)/(q+one);
      //q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[dirBSE])[kbse]=(two*feqW27_BSE-(f27_TNW*(q*omegaD-one)-omegaD*feq27_TNW*(q-one))/(omegaD-one)+f27_BSE*q)/(q+one);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADPressIncomp7(int inx,
										   int iny,
										   doubflo* DD, 
										   doubflo* DD7, 
										   doubflo* temp,
										   doubflo* velo,
										   doubflo diffusivity,
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

   Distributions7 D7;
   if (evenOrOdd==true)
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[1] = &DD7[1*size_Mat];
      D7.f[2] = &DD7[2*size_Mat];
      D7.f[3] = &DD7[3*size_Mat];
      D7.f[4] = &DD7[4*size_Mat];
      D7.f[5] = &DD7[5*size_Mat];
      D7.f[6] = &DD7[6*size_Mat];
   }
   else
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[2] = &DD7[1*size_Mat];
      D7.f[1] = &DD7[2*size_Mat];
      D7.f[4] = &DD7[3*size_Mat];
      D7.f[3] = &DD7[4*size_Mat];
      D7.f[6] = &DD7[5*size_Mat];
      D7.f[5] = &DD7[6*size_Mat];
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
      //////////////////////////////////////////////////////////////////////////////////
      doubflo  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB; 

      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      //////////////////////////////////////////////////////////////////////////////////
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
      doubflo vx1 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      doubflo vx2 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      doubflo vx3 = ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
		 ////dörrrrrty !!!!!!!!!!!!!
   //      doubflo vx1     =  ten * ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
   //      doubflo vx2     =  ten * ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
   //      doubflo vx3     =  ten * ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));

	  //doubflo cu_sq       =1.5*(vx1*vx1+vx2*vx2+vx3*vx3);
      doubflo ux_sq       = vx1 * vx1;
      doubflo uy_sq       = vx2 * vx2;
      doubflo uz_sq       = vx3 * vx3;
      //////////////////////////////////////////////////////////////////////////
	  //BGK
      //doubflo omegaD     = three - sqrt(three);
      //doubflo Lam         = -(c1o2-one/omegaD);
      //doubflo nue_d       = Lam/three;
      ////doubflo ae          = zero;
      //doubflo ae          = diffusivity/nue_d - one;

      doubflo f7_ZERO,f7_E,f7_W,f7_N,f7_S,f7_T,f7_B;
      doubflo /*feq7_ZERO,*/feq7_E,feq7_W,feq7_N,feq7_S,feq7_T,feq7_B;
      doubflo /*feqW7_ZERO,*/feqW7_E,feqW7_W,feqW7_N,feqW7_S,feqW7_T,feqW7_B;
      //doubflo TempD = temp[k];


      f7_ZERO =  (D7.f[0])[kzero];
      f7_W    =  (D7.f[1])[ke   ];
      f7_E    =  (D7.f[2])[kw   ];
      f7_S    =  (D7.f[3])[kn   ];
      f7_N    =  (D7.f[4])[ks   ];
      f7_B    =  (D7.f[5])[kt   ];
      f7_T    =  (D7.f[6])[kb   ];

      doubflo ConcD = f7_ZERO + f7_E + f7_W + f7_N + f7_S + f7_T + f7_B;

      ////feq7_ZERO = ConcD*(c1o3*(ae*(-3.0))-(ux_sq+uy_sq+uz_sq));
      //feq7_E    = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
      //feq7_W    = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
      //feq7_N    = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
      //feq7_S    = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
      //feq7_T    = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
      //feq7_B    = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);

      ////feq7_ZERO = TempD*(c1o3*(ae*(-3.0f))-(ux_sq+uy_sq+uz_sq));
      //feqW7_E    = feq7_E;// TempD*(c1o6*(ae+one)+c1o2*(ux_sq_Inflow)+vx1_Inflow*c1o2);
      //feqW7_W    = feq7_W;// TempD*(c1o6*(ae+one)+c1o2*(ux_sq_Inflow)-vx1_Inflow*c1o2);
      //feqW7_N    = feq7_N;// TempD*(c1o6*(ae+one)+c1o2*(uy_sq_Inflow)+vx2_Inflow*c1o2);
      //feqW7_S    = feq7_S;// TempD*(c1o6*(ae+one)+c1o2*(uy_sq_Inflow)-vx2_Inflow*c1o2);
      //feqW7_T    = feq7_T;// TempD*(c1o6*(ae+one)+c1o2*(uz_sq_Inflow)+vx3_Inflow*c1o2);
      //feqW7_B    = feq7_B;// TempD*(c1o6*(ae+one)+c1o2*(uz_sq_Inflow)-vx3_Inflow*c1o2);

      //////////////////////////////////////////////////////////////////////////
	  //TRT  Yoshida Kernel - based on Ying
      doubflo cs2         = c1o4;
      doubflo Lam         = diffusivity/(one)/cs2;
      doubflo omegaD      = - one / (Lam + c1o2);
      doubflo nue_d       = Lam/three;

      //////////////////////////////////////////////////////////////////////////
      //pointertausch
      if (evenOrOdd==false)
      {
         D7.f[0] = &DD7[0*size_Mat];
         D7.f[1] = &DD7[1*size_Mat];
         D7.f[2] = &DD7[2*size_Mat];
         D7.f[3] = &DD7[3*size_Mat];
         D7.f[4] = &DD7[4*size_Mat];
         D7.f[5] = &DD7[5*size_Mat];
         D7.f[6] = &DD7[6*size_Mat];
      }
      else
      {
         D7.f[0] = &DD7[0*size_Mat];
         D7.f[2] = &DD7[1*size_Mat];
         D7.f[1] = &DD7[2*size_Mat];
         D7.f[4] = &DD7[3*size_Mat];
         D7.f[3] = &DD7[4*size_Mat];
         D7.f[6] = &DD7[5*size_Mat];
         D7.f[5] = &DD7[6*size_Mat];
      }

      ////////////////////////////////////////////////////////////////////////////
      ////mit Q's
      //doubflo /*feq,*/ q;
      //q = q_dirE[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=0.;
      //   (D7.f[2])[kw]=(two*feqW7_W-(f7_E*(q*omegaD-one)-omegaD*feq7_E*(q-one))/(omegaD-one)+f7_W*q)/(q+one);//f7_W - feq7_W + feqW7_E;
      //}

      //q = q_dirW[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=0.;
      //   (D7.f[1])[ke]=(two*feqW7_E-(f7_W*(q*omegaD-one)-omegaD*feq7_W*(q-one))/(omegaD-one)+f7_E*q)/(q+one);//f7_E - feq7_E + feqW7_W;
      //}

      //q = q_dirN[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=0.;
      //   (D7.f[4])[ks]=(two*feqW7_S-(f7_N*(q*omegaD-1.)-omegaD*feq7_N*(q-one))/(omegaD-one)+f7_S*q)/(q+one);//f7_S - feq7_S + feqW7_N;
      //}

      //q = q_dirS[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=0.;
      //   (D7.f[3])[kn]=(two*feqW7_N-(f7_S*(q*omegaD-one)-omegaD*feq7_S*(q-one))/(omegaD-one)+f7_N*q)/(q+one);//f7_N - feq7_N + feqW7_S;
      //}

      //q = q_dirT[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=0.;
      //   (D7.f[6])[kb]=(two*feqW7_B-(f7_T*(q*omegaD-1.)-omegaD*feq7_T*(q-one))/(omegaD-one)+f7_B*q)/(q+one);//f7_B - feq7_B + feqW7_T;
      //}

      //q = q_dirB[k];
      //if (q>=zero && q<=one)
      //{
      //   //q=0.;
      //   (D7.f[5])[kt]=(two*feqW7_T-(f7_B*(q*omegaD-one)-omegaD*feq7_B*(q-one))/(omegaD-one)+f7_T*q)/(q+one);//f7_T - feq7_T + feqW7_B;
      //}

      ////////////////////////////////////////////////////////////////////////////
      ////ohne Q's
      //doubflo /*feq,*/ q;
      //q = q_dirE[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[2])[kw]= 0.1;//test
      //   //(D7.f[2])[kw]=f7_W - feq7_W + feqW7_E;
      //}

      //q = q_dirW[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[1])[ke]= 0.1;//test
      //   //(D7.f[1])[ke]=f7_E - feq7_E + feqW7_W;
      //}

      //q = q_dirN[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[4])[ks]= 0.1;//test
      //   //(D7.f[4])[ks]=f7_S - feq7_S + feqW7_N;
      //}

      //q = q_dirS[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[3])[kn]= 0.1;//test
      //   //(D7.f[3])[kn]=f7_N - feq7_N + feqW7_S;
      //}

      //q = q_dirT[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[6])[kb]= 0.1;//test
      //   //(D7.f[6])[kb]=f7_B - feq7_B + feqW7_T;
      //}

      //q = q_dirB[k];
      //if (q>=zero && q<=one)
      //{
      //   (D7.f[5])[kt]= 0.1;//test
      //   //(D7.f[5])[kt]=f7_T - feq7_T + feqW7_B;
      //}


      //////////////////////////////////////////////////////////////////////////
      //ohne Q's aber mit TRT
      doubflo /*feq,*/ q;
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         (D7.f[2])[kw]= f7_W + nue_d * ConcD;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         (D7.f[1])[ke]= f7_E + nue_d * ConcD;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         (D7.f[4])[ks]= f7_S + nue_d * ConcD;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         (D7.f[3])[kn]= f7_N + nue_d * ConcD;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         (D7.f[6])[kb]= f7_B + nue_d * ConcD;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         (D7.f[5])[kt]= f7_T + nue_d * ConcD;
      }

   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADPressIncomp27(   int inx,
											   int iny,
											   doubflo* DD, 
											   doubflo* DD27, 
											   doubflo* temp,
											   doubflo* velo,
											   doubflo diffusivity,
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

   Distributions27 D27;
   if (evenOrOdd==true)
   {
      D27.f[dirE   ] = &DD27[dirE   *size_Mat];
      D27.f[dirW   ] = &DD27[dirW   *size_Mat];
      D27.f[dirN   ] = &DD27[dirN   *size_Mat];
      D27.f[dirS   ] = &DD27[dirS   *size_Mat];
      D27.f[dirT   ] = &DD27[dirT   *size_Mat];
      D27.f[dirB   ] = &DD27[dirB   *size_Mat];
      D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
      D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
      D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
      D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
      D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
      D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
   } 
   else
   {
      D27.f[dirW   ] = &DD27[dirE   *size_Mat];
      D27.f[dirE   ] = &DD27[dirW   *size_Mat];
      D27.f[dirS   ] = &DD27[dirN   *size_Mat];
      D27.f[dirN   ] = &DD27[dirS   *size_Mat];
      D27.f[dirB   ] = &DD27[dirT   *size_Mat];
      D27.f[dirT   ] = &DD27[dirB   *size_Mat];
      D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
      D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
      D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
      D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
      D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
      D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
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
      doubflo  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
      doubflo f_ZERO = (D.f[dirZERO])[kzero];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1      = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      doubflo vx2      = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      doubflo vx3      = ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      doubflo f27_W    = (D27.f[dirE   ])[ke   ];
      doubflo f27_E    = (D27.f[dirW   ])[kw   ];
      doubflo f27_S    = (D27.f[dirN   ])[kn   ];
      doubflo f27_N    = (D27.f[dirS   ])[ks   ];
      doubflo f27_B    = (D27.f[dirT   ])[kt   ];
      doubflo f27_T    = (D27.f[dirB   ])[kb   ];
      doubflo f27_SW   = (D27.f[dirNE  ])[kne  ];
      doubflo f27_NE   = (D27.f[dirSW  ])[ksw  ];
      doubflo f27_NW   = (D27.f[dirSE  ])[kse  ];
      doubflo f27_SE   = (D27.f[dirNW  ])[knw  ];
      doubflo f27_BW   = (D27.f[dirTE  ])[kte  ];
      doubflo f27_TE   = (D27.f[dirBW  ])[kbw  ];
      doubflo f27_TW   = (D27.f[dirBE  ])[kbe  ];
      doubflo f27_BE   = (D27.f[dirTW  ])[ktw  ];
      doubflo f27_BS   = (D27.f[dirTN  ])[ktn  ];
      doubflo f27_TN   = (D27.f[dirBS  ])[kbs  ];
      doubflo f27_TS   = (D27.f[dirBN  ])[kbn  ];
      doubflo f27_BN   = (D27.f[dirTS  ])[kts  ];
      doubflo f27_ZERO = (D27.f[dirZERO])[kzero];
      doubflo f27_BSW  = (D27.f[dirTNE ])[ktne ];
      doubflo f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      doubflo f27_BNW  = (D27.f[dirTSE ])[ktse ];
      doubflo f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      doubflo f27_TSW  = (D27.f[dirBNE ])[kbne ];
      doubflo f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      doubflo f27_TNW  = (D27.f[dirBSE ])[kbse ];
      doubflo f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      doubflo ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
                        f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
                        f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      //doubflo feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      doubflo feq27_E    =   c2over27* ConcD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      doubflo feq27_W    =   c2over27* ConcD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      doubflo feq27_N    =   c2over27* ConcD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      doubflo feq27_S    =   c2over27* ConcD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      doubflo feq27_T    =   c2over27* ConcD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      doubflo feq27_B    =   c2over27* ConcD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      doubflo feq27_NE   =   c1over54* ConcD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      doubflo feq27_SW   =   c1over54* ConcD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      doubflo feq27_SE   =   c1over54* ConcD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      doubflo feq27_NW   =   c1over54* ConcD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      doubflo feq27_TE   =   c1over54* ConcD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      doubflo feq27_BW   =   c1over54* ConcD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      doubflo feq27_BE   =   c1over54* ConcD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      doubflo feq27_TW   =   c1over54* ConcD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      doubflo feq27_TN   =   c1over54* ConcD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      doubflo feq27_BS   =   c1over54* ConcD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      doubflo feq27_BN   =   c1over54* ConcD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      doubflo feq27_TS   =   c1over54* ConcD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      doubflo feq27_TNE  =   c1over216*ConcD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      doubflo feq27_BSW  =   c1over216*ConcD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      doubflo feq27_BNE  =   c1over216*ConcD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      doubflo feq27_TSW  =   c1over216*ConcD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      doubflo feq27_TSE  =   c1over216*ConcD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      doubflo feq27_BNW  =   c1over216*ConcD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      doubflo feq27_BSE  =   c1over216*ConcD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      doubflo feq27_TNW  =   c1over216*ConcD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo TempD = temp[k];
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // velocity inflow
      //vx1   = zero;
      //vx2   = zero;
      //vx3   = velo[k];

      //doubflo feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      doubflo feqW27_E    =  c2over27* TempD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq); //feq27_E  ;// 
      doubflo feqW27_W    =  c2over27* TempD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq); //feq27_W  ;// 
      doubflo feqW27_N    =  c2over27* TempD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq); //feq27_N  ;// 
      doubflo feqW27_S    =  c2over27* TempD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); //feq27_S  ;// 
      doubflo feqW27_T    =  c2over27* TempD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq); //feq27_T  ;// 
      doubflo feqW27_B    =  c2over27* TempD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq); //feq27_B  ;// 
      doubflo feqW27_NE   =  c1over54* TempD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); //feq27_NE ;// 
      doubflo feqW27_SW   =  c1over54* TempD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); //feq27_SW ;// 
      doubflo feqW27_SE   =  c1over54* TempD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); //feq27_SE ;// 
      doubflo feqW27_NW   =  c1over54* TempD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); //feq27_NW ;// 
      doubflo feqW27_TE   =  c1over54* TempD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); //feq27_TE ;// 
      doubflo feqW27_BW   =  c1over54* TempD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); //feq27_BW ;// 
      doubflo feqW27_BE   =  c1over54* TempD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); //feq27_BE ;// 
      doubflo feqW27_TW   =  c1over54* TempD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); //feq27_TW ;// 
      doubflo feqW27_TN   =  c1over54* TempD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); //feq27_TN ;// 
      doubflo feqW27_BS   =  c1over54* TempD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); //feq27_BS ;// 
      doubflo feqW27_BN   =  c1over54* TempD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); //feq27_BN ;// 
      doubflo feqW27_TS   =  c1over54* TempD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); //feq27_TS ;// 
      doubflo feqW27_TNE  =  c1over216*TempD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); //feq27_TNE;// 
      doubflo feqW27_BSW  =  c1over216*TempD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); //feq27_BSW;// 
      doubflo feqW27_BNE  =  c1over216*TempD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); //feq27_BNE;// 
      doubflo feqW27_TSW  =  c1over216*TempD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); //feq27_TSW;// 
      doubflo feqW27_TSE  =  c1over216*TempD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); //feq27_TSE;// 
      doubflo feqW27_BNW  =  c1over216*TempD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); //feq27_BNW;// 
      doubflo feqW27_BSE  =  c1over216*TempD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); //feq27_BSE;// 
      doubflo feqW27_TNW  =  c1over216*TempD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); //feq27_TNW;// 
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo omegaD     = three - sqrt(three);
      //doubflo Lam        = -(c1o2-one/omegaD);
      //doubflo nue_d      = Lam/three;
      //doubflo ae         = zero;
      //doubflo ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
      {
         D27.f[dirE   ] = &DD27[dirE   *size_Mat];
         D27.f[dirW   ] = &DD27[dirW   *size_Mat];
         D27.f[dirN   ] = &DD27[dirN   *size_Mat];
         D27.f[dirS   ] = &DD27[dirS   *size_Mat];
         D27.f[dirT   ] = &DD27[dirT   *size_Mat];
         D27.f[dirB   ] = &DD27[dirB   *size_Mat];
         D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
         D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
         D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
         D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
         D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
         D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
      } 
      else
      {
         D27.f[dirW   ] = &DD27[dirE   *size_Mat];
         D27.f[dirE   ] = &DD27[dirW   *size_Mat];
         D27.f[dirS   ] = &DD27[dirN   *size_Mat];
         D27.f[dirN   ] = &DD27[dirS   *size_Mat];
         D27.f[dirB   ] = &DD27[dirT   *size_Mat];
         D27.f[dirT   ] = &DD27[dirB   *size_Mat];
         D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
         D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
         D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
         D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
         D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
         D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
         D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
         D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
         D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
         D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
         D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
         D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
         D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
         D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
         D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
         D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
         D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
         D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
         D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
         D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
         D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo q;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[dirW  ])[kw  ]= -feqW27_W  + two * c2over27  * TempD;
      q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[dirE  ])[ke  ]= -feqW27_E  + two * c2over27  * TempD;
      q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[dirS  ])[ks  ]= -feqW27_S  + two * c2over27  * TempD;
      q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[dirN  ])[kn  ]= -feqW27_N  + two * c2over27  * TempD;
      q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[dirB  ])[kb  ]= -feqW27_B  + two * c2over27  * TempD;
      q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[dirT  ])[kt  ]= -feqW27_T  + two * c2over27  * TempD;
      q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[dirSW ])[ksw ]= -feqW27_SW + two * c1over54  * TempD;
      q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[dirNE ])[kne ]= -feqW27_NE + two * c1over54  * TempD;
      q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[dirNW ])[knw ]= -feqW27_NW + two * c1over54  * TempD;
      q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[dirSE ])[kse ]= -feqW27_SE + two * c1over54  * TempD;
      q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[dirBW ])[kbw ]= -feqW27_BW + two * c1over54  * TempD;
      q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[dirTE ])[kte ]= -feqW27_TE + two * c1over54  * TempD;
      q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[dirTW ])[ktw ]= -feqW27_TW + two * c1over54  * TempD;
      q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[dirBE ])[kbe ]= -feqW27_BE + two * c1over54  * TempD;
      q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[dirBS ])[kbs ]= -feqW27_BS + two * c1over54  * TempD;
      q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[dirTN ])[ktn ]= -feqW27_TN + two * c1over54  * TempD;
      q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[dirTS ])[kts ]= -feqW27_TS + two * c1over54  * TempD;
      q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[dirBN ])[kbn ]= -feqW27_BN + two * c1over54  * TempD;
      q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[dirBSW])[kbsw]= -feqW27_BSW+ two * c1over216 * TempD;
      q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[dirTNE])[ktne]= -feqW27_TNE+ two * c1over216 * TempD;
      q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[dirTSW])[ktsw]= -feqW27_TSW+ two * c1over216 * TempD;
      q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[dirBNE])[kbne]= -feqW27_BNE+ two * c1over216 * TempD;
      q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[dirBNW])[kbnw]= -feqW27_BNW+ two * c1over216 * TempD;
      q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[dirTSE])[ktse]= -feqW27_TSE+ two * c1over216 * TempD;
      q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[dirTNW])[ktnw]= -feqW27_TNW+ two * c1over216 * TempD;
      q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[dirBSE])[kbse]= -feqW27_BSE+ two * c1over216 * TempD;
      //q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[dirW  ])[kw  ]=(two*feqW27_W  -(f27_E  *(q*omegaD-one)-omegaD*feq27_E  *(q-one))/(omegaD-one)+f27_W  *q)/(q+one);
      //q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[dirE  ])[ke  ]=(two*feqW27_E  -(f27_W  *(q*omegaD-one)-omegaD*feq27_W  *(q-one))/(omegaD-one)+f27_E  *q)/(q+one);
      //q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[dirS  ])[ks  ]=(two*feqW27_S  -(f27_N  *(q*omegaD-one)-omegaD*feq27_N  *(q-one))/(omegaD-one)+f27_S  *q)/(q+one);
      //q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[dirN  ])[kn  ]=(two*feqW27_N  -(f27_S  *(q*omegaD-one)-omegaD*feq27_S  *(q-one))/(omegaD-one)+f27_N  *q)/(q+one);
      //q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[dirB  ])[kb  ]=(two*feqW27_B  -(f27_T  *(q*omegaD-one)-omegaD*feq27_T  *(q-one))/(omegaD-one)+f27_B  *q)/(q+one);
      //q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[dirT  ])[kt  ]=(two*feqW27_T  -(f27_B  *(q*omegaD-one)-omegaD*feq27_B  *(q-one))/(omegaD-one)+f27_T  *q)/(q+one);
      //q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[dirSW ])[ksw ]=(two*feqW27_SW -(f27_NE *(q*omegaD-one)-omegaD*feq27_NE *(q-one))/(omegaD-one)+f27_SW *q)/(q+one);
      //q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[dirNE ])[kne ]=(two*feqW27_NE -(f27_SW *(q*omegaD-one)-omegaD*feq27_SW *(q-one))/(omegaD-one)+f27_NE *q)/(q+one);
      //q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[dirNW ])[knw ]=(two*feqW27_NW -(f27_SE *(q*omegaD-one)-omegaD*feq27_SE *(q-one))/(omegaD-one)+f27_NW *q)/(q+one);
      //q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[dirSE ])[kse ]=(two*feqW27_SE -(f27_NW *(q*omegaD-one)-omegaD*feq27_NW *(q-one))/(omegaD-one)+f27_SE *q)/(q+one);
      //q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[dirBW ])[kbw ]=(two*feqW27_BW -(f27_TE *(q*omegaD-one)-omegaD*feq27_TE *(q-one))/(omegaD-one)+f27_BW *q)/(q+one);
      //q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[dirTE ])[kte ]=(two*feqW27_TE -(f27_BW *(q*omegaD-one)-omegaD*feq27_BW *(q-one))/(omegaD-one)+f27_TE *q)/(q+one);
      //q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[dirTW ])[ktw ]=(two*feqW27_TW -(f27_BE *(q*omegaD-one)-omegaD*feq27_BE *(q-one))/(omegaD-one)+f27_TW *q)/(q+one);
      //q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[dirBE ])[kbe ]=(two*feqW27_BE -(f27_TW *(q*omegaD-one)-omegaD*feq27_TW *(q-one))/(omegaD-one)+f27_BE *q)/(q+one);
      //q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[dirBS ])[kbs ]=(two*feqW27_BS -(f27_TN *(q*omegaD-one)-omegaD*feq27_TN *(q-one))/(omegaD-one)+f27_BS *q)/(q+one);
      //q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[dirTN ])[ktn ]=(two*feqW27_TN -(f27_BS *(q*omegaD-one)-omegaD*feq27_BS *(q-one))/(omegaD-one)+f27_TN *q)/(q+one);
      //q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[dirTS ])[kts ]=(two*feqW27_TS -(f27_BN *(q*omegaD-one)-omegaD*feq27_BN *(q-one))/(omegaD-one)+f27_TS *q)/(q+one);
      //q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[dirBN ])[kbn ]=(two*feqW27_BN -(f27_TS *(q*omegaD-one)-omegaD*feq27_TS *(q-one))/(omegaD-one)+f27_BN *q)/(q+one);
      //q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[dirBSW])[kbsw]=(two*feqW27_BSW-(f27_TNE*(q*omegaD-one)-omegaD*feq27_TNE*(q-one))/(omegaD-one)+f27_BSW*q)/(q+one);
      //q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[dirTNE])[ktne]=(two*feqW27_TNE-(f27_BSW*(q*omegaD-one)-omegaD*feq27_BSW*(q-one))/(omegaD-one)+f27_TNE*q)/(q+one);
      //q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[dirTSW])[ktsw]=(two*feqW27_TSW-(f27_BNE*(q*omegaD-one)-omegaD*feq27_BNE*(q-one))/(omegaD-one)+f27_TSW*q)/(q+one);
      //q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[dirBNE])[kbne]=(two*feqW27_BNE-(f27_TSW*(q*omegaD-one)-omegaD*feq27_TSW*(q-one))/(omegaD-one)+f27_BNE*q)/(q+one);
      //q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[dirBNW])[kbnw]=(two*feqW27_BNW-(f27_TSE*(q*omegaD-one)-omegaD*feq27_TSE*(q-one))/(omegaD-one)+f27_BNW*q)/(q+one);
      //q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[dirTSE])[ktse]=(two*feqW27_TSE-(f27_BNW*(q*omegaD-one)-omegaD*feq27_BNW*(q-one))/(omegaD-one)+f27_TSE*q)/(q+one);
      //q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[dirTNW])[ktnw]=(two*feqW27_TNW-(f27_BSE*(q*omegaD-one)-omegaD*feq27_BSE*(q-one))/(omegaD-one)+f27_TNW*q)/(q+one);
      //q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[dirBSE])[kbse]=(two*feqW27_BSE-(f27_TNW*(q*omegaD-one)-omegaD*feq27_TNW*(q-one))/(omegaD-one)+f27_BSE*q)/(q+one);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

