/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"

#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;

//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADPress7(  real* DD, 
                                       real* DD7, 
                                       real* temp,
                                       real* velo,
                                       real diffusivity,
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
   if (isEvenTimestep==true)
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

   if(k<numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB;//, 
      //         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
      //         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
      //         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
      //         *q_dirBSE, *q_dirBNW;

      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      //q_dirNE  = &QQ[dirNE  * numberOfBCnodes];
      //q_dirSW  = &QQ[dirSW  * numberOfBCnodes];
      //q_dirSE  = &QQ[dirSE  * numberOfBCnodes];
      //q_dirNW  = &QQ[dirNW  * numberOfBCnodes];
      //q_dirTE  = &QQ[dirTE  * numberOfBCnodes];
      //q_dirBW  = &QQ[dirBW  * numberOfBCnodes];
      //q_dirBE  = &QQ[dirBE  * numberOfBCnodes];
      //q_dirTW  = &QQ[dirTW  * numberOfBCnodes];
      //q_dirTN  = &QQ[dirTN  * numberOfBCnodes];
      //q_dirBS  = &QQ[dirBS  * numberOfBCnodes];
      //q_dirBN  = &QQ[dirBN  * numberOfBCnodes];
      //q_dirTS  = &QQ[dirTS  * numberOfBCnodes];
      //q_dirTNE = &QQ[dirTNE * numberOfBCnodes];
      //q_dirTSW = &QQ[dirTSW * numberOfBCnodes];
      //q_dirTSE = &QQ[dirTSE * numberOfBCnodes];
      //q_dirTNW = &QQ[dirTNW * numberOfBCnodes];
      //q_dirBNE = &QQ[dirBNE * numberOfBCnodes];
      //q_dirBSW = &QQ[dirBSW * numberOfBCnodes];
      //q_dirBSE = &QQ[dirBSE * numberOfBCnodes];
      //q_dirBNW = &QQ[dirBNW * numberOfBCnodes];
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
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
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
      /*real drho*/;
      //real vx1_Inflow   = zero;
      //real vx2_Inflow   = zero;
      //real vx3_Inflow   = velo[k];
      //real ux_sq_Inflow = vx1_Inflow * vx1_Inflow;
      //real uy_sq_Inflow = vx2_Inflow * vx2_Inflow;
      //real uz_sq_Inflow = vx3_Inflow * vx3_Inflow;


      //drho   =    f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
      //            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
      //            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]); 

      //real vx1 =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //               ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //               (f_E - f_W); 

      //real vx2 =  (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //               ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //               (f_N - f_S); 

      //real vx3 =  ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //               (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //               (f_T - f_B); 

      real rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ ((D.f[dirREST])[kzero]);
      real rho    =  rho0 + c1o1;
      real OORho  =  c1o1/rho;
      real vx1    =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      real vx2    =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      real vx3    =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));

	  //real cu_sq       =1.5*(vx1*vx1+vx2*vx2+vx3*vx3);
      real ux_sq       = vx1 * vx1;
      real uy_sq       = vx2 * vx2;
      real uz_sq       = vx3 * vx3;
      real omegaD     = c3o1 - sqrt(c3o1);
      real Lam         = -(c1o2-c1o1/omegaD);
      real nue_d       = Lam/c3o1;
      //real ae          = zero;
      real ae          = diffusivity/nue_d - c1o1;

      real f7_ZERO,f7_E,f7_W,f7_N,f7_S,f7_T,f7_B;
      real /*feq7_ZERO,*/feq7_E,feq7_W,feq7_N,feq7_S,feq7_T,feq7_B;
      real /*feqW7_ZERO,*/feqW7_E,feqW7_W,feqW7_N,feqW7_S,feqW7_T,feqW7_B;
      //real TempD = temp[k];


      f7_ZERO =  (D7.f[0])[kzero];
      f7_W    =  (D7.f[1])[ke   ];
      f7_E    =  (D7.f[2])[kw   ];
      f7_S    =  (D7.f[3])[kn   ];
      f7_N    =  (D7.f[4])[ks   ];
      f7_B    =  (D7.f[5])[kt   ];
      f7_T    =  (D7.f[6])[kb   ];

      real ConcD = f7_ZERO + f7_E + f7_W + f7_N + f7_S + f7_T + f7_B;

      //feq7_ZERO = ConcD*(c1o3*(ae*(-3.0))-(ux_sq+uy_sq+uz_sq));
      feq7_E    = ConcD*(c1o6*(ae+c1o1)+c1o2*(ux_sq)+vx1*c1o2);
      feq7_W    = ConcD*(c1o6*(ae+c1o1)+c1o2*(ux_sq)-vx1*c1o2);
      feq7_N    = ConcD*(c1o6*(ae+c1o1)+c1o2*(uy_sq)+vx2*c1o2);
      feq7_S    = ConcD*(c1o6*(ae+c1o1)+c1o2*(uy_sq)-vx2*c1o2);
      feq7_T    = ConcD*(c1o6*(ae+c1o1)+c1o2*(uz_sq)+vx3*c1o2);
      feq7_B    = ConcD*(c1o6*(ae+c1o1)+c1o2*(uz_sq)-vx3*c1o2);

      //feq7_ZERO = TempD*(c1o3*(ae*(-3.0f))-(ux_sq+uy_sq+uz_sq));
      feqW7_E    = feq7_E;// TempD*(c1o6*(ae+one)+c1o2*(ux_sq_Inflow)+vx1_Inflow*c1o2);
      feqW7_W    = feq7_W;// TempD*(c1o6*(ae+one)+c1o2*(ux_sq_Inflow)-vx1_Inflow*c1o2);
      feqW7_N    = feq7_N;// TempD*(c1o6*(ae+one)+c1o2*(uy_sq_Inflow)+vx2_Inflow*c1o2);
      feqW7_S    = feq7_S;// TempD*(c1o6*(ae+one)+c1o2*(uy_sq_Inflow)-vx2_Inflow*c1o2);
      feqW7_T    = feq7_T;// TempD*(c1o6*(ae+one)+c1o2*(uz_sq_Inflow)+vx3_Inflow*c1o2);
      feqW7_B    = feq7_B;// TempD*(c1o6*(ae+one)+c1o2*(uz_sq_Inflow)-vx3_Inflow*c1o2);

      //////////////////////////////////////////////////////////////////////////
      //pointertausch
      if (isEvenTimestep==false)
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
      //(D.f[dirREST])[k]=0.1f;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //(D7.f[1])[ke   ] = f7_E - feq7_E + feqW7_W; //E
      //(D7.f[2])[kw   ] = f7_W - feq7_W + feqW7_E; //W
      //(D7.f[3])[kn   ] = f7_N - feq7_N + feqW7_S; //N
      //(D7.f[4])[ks   ] = f7_S - feq7_S + feqW7_N; //S
      //(D7.f[5])[kt   ] = f7_T - feq7_T + feqW7_B; //T
      //(D7.f[6])[kb   ] = f7_B - feq7_B + feqW7_T; //B

      //////////////////////////////////////////////////////////////////////////
      //mit Q's
      real /*feq,*/ q;
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         //q=0.;
         (D7.f[2])[kw]=(c2o1*feqW7_W-(f7_E*(q*omegaD-c1o1)-omegaD*feq7_E*(q-c1o1))/(omegaD-c1o1)+f7_W*q)/(q+c1o1);//f7_W - feq7_W + feqW7_E;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         //q=0.;
         (D7.f[1])[ke]=(c2o1*feqW7_E-(f7_W*(q*omegaD-c1o1)-omegaD*feq7_W*(q-c1o1))/(omegaD-c1o1)+f7_E*q)/(q+c1o1);//f7_E - feq7_E + feqW7_W;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         //q=0.;
         (D7.f[4])[ks]=(c2o1*feqW7_S-(f7_N*(q*omegaD-1.)-omegaD*feq7_N*(q-c1o1))/(omegaD-c1o1)+f7_S*q)/(q+c1o1);//f7_S - feq7_S + feqW7_N;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         //q=0.;
         (D7.f[3])[kn]=(c2o1*feqW7_N-(f7_S*(q*omegaD-c1o1)-omegaD*feq7_S*(q-c1o1))/(omegaD-c1o1)+f7_N*q)/(q+c1o1);//f7_N - feq7_N + feqW7_S;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         //q=0.;
         (D7.f[6])[kb]=(c2o1*feqW7_B-(f7_T*(q*omegaD-1.)-omegaD*feq7_T*(q-c1o1))/(omegaD-c1o1)+f7_B*q)/(q+c1o1);//f7_B - feq7_B + feqW7_T;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         //q=0.;
         (D7.f[5])[kt]=(c2o1*feqW7_T-(f7_B*(q*omegaD-c1o1)-omegaD*feq7_B*(q-c1o1))/(omegaD-c1o1)+f7_T*q)/(q+c1o1);//f7_T - feq7_T + feqW7_B;
      }

      ////////////////////////////////////////////////////////////////////////////
      ////ohne Q's
      //real /*feq,*/ q;
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
extern "C" __global__ void QADPress27( real* DD, 
                                       real* DD27, 
                                       real* temp,
                                       real* velo,
                                       real diffusivity,
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
   if (isEvenTimestep==true)
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      q_dirNE  = &QQ[dirNE  * numberOfBCnodes];
      q_dirSW  = &QQ[dirSW  * numberOfBCnodes];
      q_dirSE  = &QQ[dirSE  * numberOfBCnodes];
      q_dirNW  = &QQ[dirNW  * numberOfBCnodes];
      q_dirTE  = &QQ[dirTE  * numberOfBCnodes];
      q_dirBW  = &QQ[dirBW  * numberOfBCnodes];
      q_dirBE  = &QQ[dirBE  * numberOfBCnodes];
      q_dirTW  = &QQ[dirTW  * numberOfBCnodes];
      q_dirTN  = &QQ[dirTN  * numberOfBCnodes];
      q_dirBS  = &QQ[dirBS  * numberOfBCnodes];
      q_dirBN  = &QQ[dirBN  * numberOfBCnodes];
      q_dirTS  = &QQ[dirTS  * numberOfBCnodes];
      q_dirTNE = &QQ[dirTNE * numberOfBCnodes];
      q_dirTSW = &QQ[dirTSW * numberOfBCnodes];
      q_dirTSE = &QQ[dirTSE * numberOfBCnodes];
      q_dirTNW = &QQ[dirTNW * numberOfBCnodes];
      q_dirBNE = &QQ[dirBNE * numberOfBCnodes];
      q_dirBSW = &QQ[dirBSW * numberOfBCnodes];
      q_dirBSE = &QQ[dirBSE * numberOfBCnodes];
      q_dirBNW = &QQ[dirBNW * numberOfBCnodes];
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
      real f_ZERO = (D.f[dirREST])[kzero];
      real f_BSW  = (D.f[dirTNE ])[ktne ];
      real f_BNE  = (D.f[dirTSW ])[ktsw ];
      real f_BNW  = (D.f[dirTSE ])[ktse ];
      real f_BSE  = (D.f[dirTNW ])[ktnw ];
      real f_TSW  = (D.f[dirBNE ])[kbne ];
      real f_TNE  = (D.f[dirBSW ])[kbsw ];
      real f_TNW  = (D.f[dirBSE ])[kbse ];
      real f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, /*drho, feq,*/ q;
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
      real rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ f_ZERO;
      real rho    =  rho0 + c1o1;
      real OORho  =  c1o1/rho;
      vx1            =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      vx2            =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      vx3            =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      real f27_W    = (D27.f[dirE   ])[ke   ];
      real f27_E    = (D27.f[dirW   ])[kw   ];
      real f27_S    = (D27.f[dirN   ])[kn   ];
      real f27_N    = (D27.f[dirS   ])[ks   ];
      real f27_B    = (D27.f[dirT   ])[kt   ];
      real f27_T    = (D27.f[dirB   ])[kb   ];
      real f27_SW   = (D27.f[dirNE  ])[kne  ];
      real f27_NE   = (D27.f[dirSW  ])[ksw  ];
      real f27_NW   = (D27.f[dirSE  ])[kse  ];
      real f27_SE   = (D27.f[dirNW  ])[knw  ];
      real f27_BW   = (D27.f[dirTE  ])[kte  ];
      real f27_TE   = (D27.f[dirBW  ])[kbw  ];
      real f27_TW   = (D27.f[dirBE  ])[kbe  ];
      real f27_BE   = (D27.f[dirTW  ])[ktw  ];
      real f27_BS   = (D27.f[dirTN  ])[ktn  ];
      real f27_TN   = (D27.f[dirBS  ])[kbs  ];
      real f27_TS   = (D27.f[dirBN  ])[kbn  ];
      real f27_BN   = (D27.f[dirTS  ])[kts  ];
      real f27_ZERO = (D27.f[dirREST])[kzero];
      real f27_BSW  = (D27.f[dirTNE ])[ktne ];
      real f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      real f27_BNW  = (D27.f[dirTSE ])[ktse ];
      real f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      real f27_TSW  = (D27.f[dirBNE ])[kbne ];
      real f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      real f27_TNW  = (D27.f[dirBSE ])[kbse ];
      real f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      real ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
                        f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
                        f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      ////real feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      //real feq27_E    =   c2over27* (ConcD+(one+ConcD)*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq));
      //real feq27_W    =   c2over27* (ConcD+(one+ConcD)*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq));
      //real feq27_N    =   c2over27* (ConcD+(one+ConcD)*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq));
      //real feq27_S    =   c2over27* (ConcD+(one+ConcD)*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq));
      //real feq27_T    =   c2over27* (ConcD+(one+ConcD)*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq));
      //real feq27_B    =   c2over27* (ConcD+(one+ConcD)*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq));
      //real feq27_NE   =   c1over54* (ConcD+(one+ConcD)*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      //real feq27_SW   =   c1over54* (ConcD+(one+ConcD)*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      //real feq27_SE   =   c1over54* (ConcD+(one+ConcD)*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      //real feq27_NW   =   c1over54* (ConcD+(one+ConcD)*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      //real feq27_TE   =   c1over54* (ConcD+(one+ConcD)*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      //real feq27_BW   =   c1over54* (ConcD+(one+ConcD)*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      //real feq27_BE   =   c1over54* (ConcD+(one+ConcD)*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      //real feq27_TW   =   c1over54* (ConcD+(one+ConcD)*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      //real feq27_TN   =   c1over54* (ConcD+(one+ConcD)*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      //real feq27_BS   =   c1over54* (ConcD+(one+ConcD)*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      //real feq27_BN   =   c1over54* (ConcD+(one+ConcD)*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      //real feq27_TS   =   c1over54* (ConcD+(one+ConcD)*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      //real feq27_TNE  =   c1over216*(ConcD+(one+ConcD)*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      //real feq27_BSW  =   c1over216*(ConcD+(one+ConcD)*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      //real feq27_BNE  =   c1over216*(ConcD+(one+ConcD)*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      //real feq27_TSW  =   c1over216*(ConcD+(one+ConcD)*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      //real feq27_TSE  =   c1over216*(ConcD+(one+ConcD)*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      //real feq27_BNW  =   c1over216*(ConcD+(one+ConcD)*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      //real feq27_BSE  =   c1over216*(ConcD+(one+ConcD)*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      //real feq27_TNW  =   c1over216*(ConcD+(one+ConcD)*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));
      //
      real feq27_E    =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      real feq27_W    =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      real feq27_N    =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      real feq27_S    =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      real feq27_T    =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      real feq27_B    =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      real feq27_NE   =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      real feq27_SW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      real feq27_SE   =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      real feq27_NW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      real feq27_TE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      real feq27_BW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      real feq27_BE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      real feq27_TW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      real feq27_TN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      real feq27_BS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      real feq27_BN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      real feq27_TS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      real feq27_TNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      real feq27_BSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      real feq27_BNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      real feq27_TSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      real feq27_TSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      real feq27_BNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      real feq27_BSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      real feq27_TNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //real TempD = temp[k];
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // velocity inflow
      //vx1   = zero;
      //vx2   = zero;
      //vx3   = velo[k];

      //real feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      real feqW27_E    =  feq27_E  ;// c2over27* TempD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      real feqW27_W    =  feq27_W  ;// c2over27* TempD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      real feqW27_N    =  feq27_N  ;// c2over27* TempD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      real feqW27_S    =  feq27_S  ;// c2over27* TempD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      real feqW27_T    =  feq27_T  ;// c2over27* TempD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      real feqW27_B    =  feq27_B  ;// c2over27* TempD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      real feqW27_NE   =  feq27_NE ;// c1over54* TempD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      real feqW27_SW   =  feq27_SW ;// c1over54* TempD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      real feqW27_SE   =  feq27_SE ;// c1over54* TempD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      real feqW27_NW   =  feq27_NW ;// c1over54* TempD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      real feqW27_TE   =  feq27_TE ;// c1over54* TempD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      real feqW27_BW   =  feq27_BW ;// c1over54* TempD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      real feqW27_BE   =  feq27_BE ;// c1over54* TempD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      real feqW27_TW   =  feq27_TW ;// c1over54* TempD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      real feqW27_TN   =  feq27_TN ;// c1over54* TempD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      real feqW27_BS   =  feq27_BS ;// c1over54* TempD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      real feqW27_BN   =  feq27_BN ;// c1over54* TempD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      real feqW27_TS   =  feq27_TS ;// c1over54* TempD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      real feqW27_TNE  =  feq27_TNE;// c1over216*TempD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      real feqW27_BSW  =  feq27_BSW;// c1over216*TempD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      real feqW27_BNE  =  feq27_BNE;// c1over216*TempD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      real feqW27_TSW  =  feq27_TSW;// c1over216*TempD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      real feqW27_TSE  =  feq27_TSE;// c1over216*TempD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      real feqW27_BNW  =  feq27_BNW;// c1over216*TempD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      real feqW27_BSE  =  feq27_BSE;// c1over216*TempD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      real feqW27_TNW  =  feq27_TNW;// c1over216*TempD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real omegaD     = c3o1 - sqrt(c3o1);
      //real Lam        = -(c1o2-one/omegaD);
      //real nue_d      = Lam/three;
      //real ae         = zero;
      //real ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      //(D.f[dirREST])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirW  ])[kw  ]=(c2o1*feqW27_W  -(f27_E  *(q*omegaD-c1o1)-omegaD*feq27_E  *(q-c1o1))/(omegaD-c1o1)+f27_W  *q)/(q+c1o1);
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirE  ])[ke  ]=(c2o1*feqW27_E  -(f27_W  *(q*omegaD-c1o1)-omegaD*feq27_W  *(q-c1o1))/(omegaD-c1o1)+f27_E  *q)/(q+c1o1);
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirS  ])[ks  ]=(c2o1*feqW27_S  -(f27_N  *(q*omegaD-c1o1)-omegaD*feq27_N  *(q-c1o1))/(omegaD-c1o1)+f27_S  *q)/(q+c1o1);
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirN  ])[kn  ]=(c2o1*feqW27_N  -(f27_S  *(q*omegaD-c1o1)-omegaD*feq27_S  *(q-c1o1))/(omegaD-c1o1)+f27_N  *q)/(q+c1o1);
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirB  ])[kb  ]=(c2o1*feqW27_B  -(f27_T  *(q*omegaD-c1o1)-omegaD*feq27_T  *(q-c1o1))/(omegaD-c1o1)+f27_B  *q)/(q+c1o1);
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirT  ])[kt  ]=(c2o1*feqW27_T  -(f27_B  *(q*omegaD-c1o1)-omegaD*feq27_B  *(q-c1o1))/(omegaD-c1o1)+f27_T  *q)/(q+c1o1);
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirSW ])[ksw ]=(c2o1*feqW27_SW -(f27_NE *(q*omegaD-c1o1)-omegaD*feq27_NE *(q-c1o1))/(omegaD-c1o1)+f27_SW *q)/(q+c1o1);
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirNE ])[kne ]=(c2o1*feqW27_NE -(f27_SW *(q*omegaD-c1o1)-omegaD*feq27_SW *(q-c1o1))/(omegaD-c1o1)+f27_NE *q)/(q+c1o1);
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirNW ])[knw ]=(c2o1*feqW27_NW -(f27_SE *(q*omegaD-c1o1)-omegaD*feq27_SE *(q-c1o1))/(omegaD-c1o1)+f27_NW *q)/(q+c1o1);
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirSE ])[kse ]=(c2o1*feqW27_SE -(f27_NW *(q*omegaD-c1o1)-omegaD*feq27_NW *(q-c1o1))/(omegaD-c1o1)+f27_SE *q)/(q+c1o1);
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBW ])[kbw ]=(c2o1*feqW27_BW -(f27_TE *(q*omegaD-c1o1)-omegaD*feq27_TE *(q-c1o1))/(omegaD-c1o1)+f27_BW *q)/(q+c1o1);
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTE ])[kte ]=(c2o1*feqW27_TE -(f27_BW *(q*omegaD-c1o1)-omegaD*feq27_BW *(q-c1o1))/(omegaD-c1o1)+f27_TE *q)/(q+c1o1);
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTW ])[ktw ]=(c2o1*feqW27_TW -(f27_BE *(q*omegaD-c1o1)-omegaD*feq27_BE *(q-c1o1))/(omegaD-c1o1)+f27_TW *q)/(q+c1o1);
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBE ])[kbe ]=(c2o1*feqW27_BE -(f27_TW *(q*omegaD-c1o1)-omegaD*feq27_TW *(q-c1o1))/(omegaD-c1o1)+f27_BE *q)/(q+c1o1);
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBS ])[kbs ]=(c2o1*feqW27_BS -(f27_TN *(q*omegaD-c1o1)-omegaD*feq27_TN *(q-c1o1))/(omegaD-c1o1)+f27_BS *q)/(q+c1o1);
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTN ])[ktn ]=(c2o1*feqW27_TN -(f27_BS *(q*omegaD-c1o1)-omegaD*feq27_BS *(q-c1o1))/(omegaD-c1o1)+f27_TN *q)/(q+c1o1);
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTS ])[kts ]=(c2o1*feqW27_TS -(f27_BN *(q*omegaD-c1o1)-omegaD*feq27_BN *(q-c1o1))/(omegaD-c1o1)+f27_TS *q)/(q+c1o1);
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBN ])[kbn ]=(c2o1*feqW27_BN -(f27_TS *(q*omegaD-c1o1)-omegaD*feq27_TS *(q-c1o1))/(omegaD-c1o1)+f27_BN *q)/(q+c1o1);
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSW])[kbsw]=(c2o1*feqW27_BSW-(f27_TNE*(q*omegaD-c1o1)-omegaD*feq27_TNE*(q-c1o1))/(omegaD-c1o1)+f27_BSW*q)/(q+c1o1);
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNE])[ktne]=(c2o1*feqW27_TNE-(f27_BSW*(q*omegaD-c1o1)-omegaD*feq27_BSW*(q-c1o1))/(omegaD-c1o1)+f27_TNE*q)/(q+c1o1);
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSW])[ktsw]=(c2o1*feqW27_TSW-(f27_BNE*(q*omegaD-c1o1)-omegaD*feq27_BNE*(q-c1o1))/(omegaD-c1o1)+f27_TSW*q)/(q+c1o1);
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNE])[kbne]=(c2o1*feqW27_BNE-(f27_TSW*(q*omegaD-c1o1)-omegaD*feq27_TSW*(q-c1o1))/(omegaD-c1o1)+f27_BNE*q)/(q+c1o1);
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNW])[kbnw]=(c2o1*feqW27_BNW-(f27_TSE*(q*omegaD-c1o1)-omegaD*feq27_TSE*(q-c1o1))/(omegaD-c1o1)+f27_BNW*q)/(q+c1o1);
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSE])[ktse]=(c2o1*feqW27_TSE-(f27_BNW*(q*omegaD-c1o1)-omegaD*feq27_BNW*(q-c1o1))/(omegaD-c1o1)+f27_TSE*q)/(q+c1o1);
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNW])[ktnw]=(c2o1*feqW27_TNW-(f27_BSE*(q*omegaD-c1o1)-omegaD*feq27_BSE*(q-c1o1))/(omegaD-c1o1)+f27_TNW*q)/(q+c1o1);
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSE])[kbse]=(c2o1*feqW27_BSE-(f27_TNW*(q*omegaD-c1o1)-omegaD*feq27_TNW*(q-c1o1))/(omegaD-c1o1)+f27_BSE*q)/(q+c1o1);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADPressNEQNeighbor27(
													real* DD,
													real* DD27,
													int* k_Q,
													int* k_N,
													int numberOfBCnodes,
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													unsigned int size_Mat,
													bool isEvenTimestep
												)
{
	Distributions27 D;
	if (isEvenTimestep == true)
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
		D.f[dirREST] = &DD[dirREST*size_Mat];
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
		D.f[dirREST] = &DD[dirREST*size_Mat];
		D.f[dirTNE] = &DD[dirBSW *size_Mat];
		D.f[dirTSW] = &DD[dirBNE *size_Mat];
		D.f[dirTSE] = &DD[dirBNW *size_Mat];
		D.f[dirTNW] = &DD[dirBSE *size_Mat];
		D.f[dirBNE] = &DD[dirTSW *size_Mat];
		D.f[dirBSW] = &DD[dirTNE *size_Mat];
		D.f[dirBSE] = &DD[dirTNW *size_Mat];
		D.f[dirBNW] = &DD[dirTSE *size_Mat];
	}

	Distributions27 D27;
	if (isEvenTimestep == true)
	{
		D27.f[dirE] = &DD27[dirE   *size_Mat];
		D27.f[dirW] = &DD27[dirW   *size_Mat];
		D27.f[dirN] = &DD27[dirN   *size_Mat];
		D27.f[dirS] = &DD27[dirS   *size_Mat];
		D27.f[dirT] = &DD27[dirT   *size_Mat];
		D27.f[dirB] = &DD27[dirB   *size_Mat];
		D27.f[dirNE] = &DD27[dirNE  *size_Mat];
		D27.f[dirSW] = &DD27[dirSW  *size_Mat];
		D27.f[dirSE] = &DD27[dirSE  *size_Mat];
		D27.f[dirNW] = &DD27[dirNW  *size_Mat];
		D27.f[dirTE] = &DD27[dirTE  *size_Mat];
		D27.f[dirBW] = &DD27[dirBW  *size_Mat];
		D27.f[dirBE] = &DD27[dirBE  *size_Mat];
		D27.f[dirTW] = &DD27[dirTW  *size_Mat];
		D27.f[dirTN] = &DD27[dirTN  *size_Mat];
		D27.f[dirBS] = &DD27[dirBS  *size_Mat];
		D27.f[dirBN] = &DD27[dirBN  *size_Mat];
		D27.f[dirTS] = &DD27[dirTS  *size_Mat];
		D27.f[dirREST] = &DD27[dirREST*size_Mat];
		D27.f[dirTNE] = &DD27[dirTNE *size_Mat];
		D27.f[dirTSW] = &DD27[dirTSW *size_Mat];
		D27.f[dirTSE] = &DD27[dirTSE *size_Mat];
		D27.f[dirTNW] = &DD27[dirTNW *size_Mat];
		D27.f[dirBNE] = &DD27[dirBNE *size_Mat];
		D27.f[dirBSW] = &DD27[dirBSW *size_Mat];
		D27.f[dirBSE] = &DD27[dirBSE *size_Mat];
		D27.f[dirBNW] = &DD27[dirBNW *size_Mat];
	}
	else
	{
		D27.f[dirW] = &DD27[dirE   *size_Mat];
		D27.f[dirE] = &DD27[dirW   *size_Mat];
		D27.f[dirS] = &DD27[dirN   *size_Mat];
		D27.f[dirN] = &DD27[dirS   *size_Mat];
		D27.f[dirB] = &DD27[dirT   *size_Mat];
		D27.f[dirT] = &DD27[dirB   *size_Mat];
		D27.f[dirSW] = &DD27[dirNE  *size_Mat];
		D27.f[dirNE] = &DD27[dirSW  *size_Mat];
		D27.f[dirNW] = &DD27[dirSE  *size_Mat];
		D27.f[dirSE] = &DD27[dirNW  *size_Mat];
		D27.f[dirBW] = &DD27[dirTE  *size_Mat];
		D27.f[dirTE] = &DD27[dirBW  *size_Mat];
		D27.f[dirTW] = &DD27[dirBE  *size_Mat];
		D27.f[dirBE] = &DD27[dirTW  *size_Mat];
		D27.f[dirBS] = &DD27[dirTN  *size_Mat];
		D27.f[dirTN] = &DD27[dirBS  *size_Mat];
		D27.f[dirTS] = &DD27[dirBN  *size_Mat];
		D27.f[dirBN] = &DD27[dirTS  *size_Mat];
		D27.f[dirREST] = &DD27[dirREST*size_Mat];
		D27.f[dirTNE] = &DD27[dirBSW *size_Mat];
		D27.f[dirTSW] = &DD27[dirBNE *size_Mat];
		D27.f[dirTSE] = &DD27[dirBNW *size_Mat];
		D27.f[dirTNW] = &DD27[dirBSE *size_Mat];
		D27.f[dirBNE] = &DD27[dirTSW *size_Mat];
		D27.f[dirBSW] = &DD27[dirTNE *size_Mat];
		D27.f[dirBSE] = &DD27[dirTNW *size_Mat];
		D27.f[dirBNW] = &DD27[dirTSE *size_Mat];
	}
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if (k < numberOfBCnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		//Fluid - BC Nodes
		////////////////////////////////////////////////////////////////////////////////
		//index
		unsigned int KQK = k_Q[k];
		unsigned int kzero = KQK;
		unsigned int ke = KQK;
		unsigned int kw = neighborX[KQK];
		unsigned int kn = KQK;
		unsigned int ks = neighborY[KQK];
		unsigned int kt = KQK;
		unsigned int kb = neighborZ[KQK];
		unsigned int ksw = neighborY[kw];
		unsigned int kne = KQK;
		unsigned int kse = ks;
		unsigned int knw = kw;
		unsigned int kbw = neighborZ[kw];
		unsigned int kte = KQK;
		unsigned int kbe = kb;
		unsigned int ktw = kw;
		unsigned int kbs = neighborZ[ks];
		unsigned int ktn = KQK;
		unsigned int kbn = kb;
		unsigned int kts = ks;
		unsigned int ktse = ks;
		unsigned int kbnw = kbw;
		unsigned int ktnw = kw;
		unsigned int kbse = kbs;
		unsigned int ktsw = ksw;
		unsigned int kbne = kb;
		unsigned int ktne = KQK;
		unsigned int kbsw = neighborZ[ksw];
		////////////////////////////////////////////////////////////////////////////////
		//distributions
		real f_W =    (D.f[dirE])[ke];
		real f_E =    (D.f[dirW])[kw];
		real f_S =    (D.f[dirN])[kn];
		real f_N =    (D.f[dirS])[ks];
		real f_B =    (D.f[dirT])[kt];
		real f_T =    (D.f[dirB])[kb];
		real f_SW =   (D.f[dirNE])[kne];
		real f_NE =   (D.f[dirSW])[ksw];
		real f_NW =   (D.f[dirSE])[kse];
		real f_SE =   (D.f[dirNW])[knw];
		real f_BW =   (D.f[dirTE])[kte];
		real f_TE =   (D.f[dirBW])[kbw];
		real f_TW =   (D.f[dirBE])[kbe];
		real f_BE =   (D.f[dirTW])[ktw];
		real f_BS =   (D.f[dirTN])[ktn];
		real f_TN =   (D.f[dirBS])[kbs];
		real f_TS =   (D.f[dirBN])[kbn];
		real f_BN =   (D.f[dirTS])[kts];
		real f_ZERO = (D.f[dirREST])[kzero];
		real f_BSW =  (D.f[dirTNE])[ktne];
		real f_BNE =  (D.f[dirTSW])[ktsw];
		real f_BNW =  (D.f[dirTSE])[ktse];
		real f_BSE =  (D.f[dirTNW])[ktnw];
		real f_TSW =  (D.f[dirBNE])[kbne];
		real f_TNE =  (D.f[dirBSW])[kbsw];
		real f_TNW =  (D.f[dirBSE])[kbse];
		real f_TSE =  (D.f[dirBNW])[kbnw];
		////////////////////////////////////////////////////////////////////////////////
		//macroscopic values
		real rho0 = 
			(f_TNE + f_BSW) + (f_TSW + f_BNE) + (f_TSE + f_BNW) + (f_TNW + f_BSE) + 
			(f_NE  + f_SW ) + (f_NW  + f_SE ) + (f_TE  + f_BW ) + (f_BE  + f_TW ) + 
			(f_TN  + f_BS ) + (f_BN  + f_TS ) + 
			(f_E   + f_W  ) + (f_N   + f_S  ) + (f_T   + f_B  ) +  f_ZERO;
		real rho = rho0 + c1o1;
		real OORho = c1o1 / rho;

		real vx1 = 
			OORho*((f_TNE - f_BSW) + (f_BNE - f_TSW) + (f_TSE - f_BNW) + (f_BSE - f_TNW) + 
			(f_NE - f_SW) + (f_SE - f_NW) + (f_TE - f_BW) + (f_BE - f_TW) + (f_E - f_W));
		real vx2 = 
			OORho*((f_TNE - f_BSW) + (f_BNE - f_TSW) + (f_BNW - f_TSE) + (f_TNW - f_BSE) + 
			(f_NE - f_SW) + (f_NW - f_SE) + (f_TN - f_BS) + (f_BN - f_TS) + (f_N - f_S));
		real vx3 = 
			OORho*((f_TNE - f_BSW) + (f_TSW - f_BNE) + (f_TSE - f_BNW) + (f_TNW - f_BSE) + 
			(f_TE - f_BW) + (f_TW - f_BE) + (f_TN - f_BS) + (f_TS - f_BN) + (f_T - f_B));
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//AD - BC Nodes
		////////////////////////////////////////////////////////////////////////////////
		//distributions
		real f27_W =    (D27.f[dirE])[ke];
		real f27_E =    (D27.f[dirW])[kw];
		real f27_S =    (D27.f[dirN])[kn];
		real f27_N =    (D27.f[dirS])[ks];
		real f27_B =    (D27.f[dirT])[kt];
		real f27_T =    (D27.f[dirB])[kb];
		real f27_SW =   (D27.f[dirNE])[kne];
		real f27_NE =   (D27.f[dirSW])[ksw];
		real f27_NW =   (D27.f[dirSE])[kse];
		real f27_SE =   (D27.f[dirNW])[knw];
		real f27_BW =   (D27.f[dirTE])[kte];
		real f27_TE =   (D27.f[dirBW])[kbw];
		real f27_TW =   (D27.f[dirBE])[kbe];
		real f27_BE =   (D27.f[dirTW])[ktw];
		real f27_BS =   (D27.f[dirTN])[ktn];
		real f27_TN =   (D27.f[dirBS])[kbs];
		real f27_TS =   (D27.f[dirBN])[kbn];
		real f27_BN =   (D27.f[dirTS])[kts];
		real f27_ZERO = (D27.f[dirREST])[kzero];
		real f27_BSW =  (D27.f[dirTNE])[ktne];
		real f27_BNE =  (D27.f[dirTSW])[ktsw];
		real f27_BNW =  (D27.f[dirTSE])[ktse];
		real f27_BSE =  (D27.f[dirTNW])[ktnw];
		real f27_TSW =  (D27.f[dirBNE])[kbne];
		real f27_TNE =  (D27.f[dirBSW])[kbsw];
		real f27_TNW =  (D27.f[dirBSE])[kbse];
		real f27_TSE =  (D27.f[dirBNW])[kbnw];
		////////////////////////////////////////////////////////////////////////////////
		real cusq = c3o2*(vx1*vx1 + vx2*vx2 + vx3*vx3);
		////////////////////////////////////////////////////////////////////////////////
		//concentration
		real ConcD = 
			(f27_TNE + f27_BSW) + (f27_TSW + f27_BNE) + (f27_TSE + f27_BNW) + (f27_TNW + f27_BSE) + 
			(f27_NE  + f27_SW ) + (f27_NW  + f27_SE ) + (f27_TE  + f27_BW ) + (f27_BE  + f27_TW ) + 
			(f27_TN  + f27_BS ) + (f27_BN  + f27_TS ) +   						
			(f27_E   + f27_W  ) + (f27_N   + f27_S  ) + (f27_T   + f27_B  ) +  f27_ZERO;
		////////////////////////////////////////////////////////////////////////////////
		//calculate non-equilibrium
		f27_ZERO  -=  c8o27* (ConcD-(ConcD+c1o1)*cusq);
		f27_E     -=  c2o27* (ConcD+(ConcD+c1o1)*(c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq));
		f27_W     -=  c2o27* (ConcD+(ConcD+c1o1)*(c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq));
		f27_N     -=  c2o27* (ConcD+(ConcD+c1o1)*(c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq));
		f27_S     -=  c2o27* (ConcD+(ConcD+c1o1)*(c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq));
		f27_T     -=  c2o27* (ConcD+(ConcD+c1o1)*(c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq));
		f27_B     -=  c2o27* (ConcD+(ConcD+c1o1)*(c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq));
		f27_NE    -=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq));
		f27_SW    -=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
		f27_SE    -=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq));
		f27_NW    -=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
		f27_TE    -=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq));
		f27_BW    -=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
		f27_BE    -=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq));
		f27_TW    -=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
		f27_TN    -=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq));
		f27_BS    -=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
		f27_BN    -=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq));
		f27_TS    -=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
		f27_TNE   -=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
		f27_BSW   -=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
		f27_BNE   -=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
		f27_TSW   -=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
		f27_TSE   -=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
		f27_BNW   -=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
		f27_BSE   -=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
		f27_TNW   -=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));
		////////////////////////////////////////////////////////////////////////////////
		ConcD = c0o1;
		////////////////////////////////////////////////////////////////////////////////
		//add BC equilibrium
		f27_ZERO  +=  c8o27* (ConcD-(ConcD+c1o1)*cusq);
		f27_E     +=  c2o27* (ConcD+(ConcD+c1o1)*(c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq));
		f27_W     +=  c2o27* (ConcD+(ConcD+c1o1)*(c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq));
		f27_N     +=  c2o27* (ConcD+(ConcD+c1o1)*(c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq));
		f27_S     +=  c2o27* (ConcD+(ConcD+c1o1)*(c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq));
		f27_T     +=  c2o27* (ConcD+(ConcD+c1o1)*(c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq));
		f27_B     +=  c2o27* (ConcD+(ConcD+c1o1)*(c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq));
		f27_NE    +=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq));
		f27_SW    +=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
		f27_SE    +=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq));
		f27_NW    +=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
		f27_TE    +=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq));
		f27_BW    +=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
		f27_BE    +=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq));
		f27_TW    +=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
		f27_TN    +=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq));
		f27_BS    +=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
		f27_BN    +=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq));
		f27_TS    +=  c1o54* (ConcD+(ConcD+c1o1)*(c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
		f27_TNE   +=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
		f27_BSW   +=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
		f27_BNE   +=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
		f27_TSW   +=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
		f27_TSE   +=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
		f27_BNW   +=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
		f27_BSE   +=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
		f27_TNW   +=  c1o216*(ConcD+(ConcD+c1o1)*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        __syncthreads();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Neighbors of BC Nodes
		////////////////////////////////////////////////////////////////////////////////
		//index neighbor
		unsigned int KNQK = k_N[k];
		unsigned int kNzero = KNQK;
		unsigned int kNe = KNQK;
		unsigned int kNw = neighborX[KNQK];
		unsigned int kNn = KNQK;
		unsigned int kNs = neighborY[KNQK];
		unsigned int kNt = KNQK;
		unsigned int kNb = neighborZ[KNQK];
		unsigned int kNsw = neighborY[kNw];
		unsigned int kNne = KNQK;
		unsigned int kNse = kNs;
		unsigned int kNnw = kNw;
		unsigned int kNbw = neighborZ[kNw];
		unsigned int kNte = KNQK;
		unsigned int kNbe = kNb;
		unsigned int kNtw = kNw;
		unsigned int kNbs = neighborZ[kNs];
		unsigned int kNtn = KNQK;
		unsigned int kNbn = kNb;
		unsigned int kNts = kNs;
		unsigned int kNtse = kNs;
		unsigned int kNbnw = kNbw;
		unsigned int kNtnw = kNw;
		unsigned int kNbse = kNbs;
		unsigned int kNtsw = kNsw;
		unsigned int kNbne = kNb;
		unsigned int kNtne = KNQK;
		unsigned int kNbsw = neighborZ[kNsw];
		////////////////////////////////////////////////////////////////////////////////
		//update distributions at neighbor nodes
        (D27.f[dirE   ])[kNe   ] = f27_W   ;  
        (D27.f[dirW   ])[kNw   ] = f27_E   ;	
        (D27.f[dirN   ])[kNn   ] = f27_S   ;	
        (D27.f[dirS   ])[kNs   ] = f27_N   ;	
        (D27.f[dirT   ])[kNt   ] = f27_B   ;	
        (D27.f[dirB   ])[kNb   ] = f27_T   ;	
        (D27.f[dirNE  ])[kNne  ] = f27_SW  ;	
        (D27.f[dirSW  ])[kNsw  ] = f27_NE  ;	
        (D27.f[dirSE  ])[kNse  ] = f27_NW  ;	
        (D27.f[dirNW  ])[kNnw  ] = f27_SE  ;	
        (D27.f[dirTE  ])[kNte  ] = f27_BW  ;	
        (D27.f[dirBW  ])[kNbw  ] = f27_TE  ;	
        (D27.f[dirBE  ])[kNbe  ] = f27_TW  ;	
        (D27.f[dirTW  ])[kNtw  ] = f27_BE  ;	
        (D27.f[dirTN  ])[kNtn  ] = f27_BS  ;	
        (D27.f[dirBS  ])[kNbs  ] = f27_TN  ;	
        (D27.f[dirBN  ])[kNbn  ] = f27_TS  ;	
        (D27.f[dirTS  ])[kNts  ] = f27_BN  ;	
        (D27.f[dirREST])[kNzero] = f27_ZERO;	
        (D27.f[dirTNE ])[kNtne ] = f27_BSW ;	
        (D27.f[dirTSW ])[kNtsw ] = f27_BNE ;	
        (D27.f[dirTSE ])[kNtse ] = f27_BNW ;	
        (D27.f[dirTNW ])[kNtnw ] = f27_BSE ;	
        (D27.f[dirBNE ])[kNbne ] = f27_TSW ;	
        (D27.f[dirBSW ])[kNbsw ] = f27_TNE ;	
        (D27.f[dirBSE ])[kNbse ] = f27_TNW ;	
        (D27.f[dirBNW ])[kNbnw ] = f27_TSE ;       
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADVel7( real* DD, 
                                    real* DD7, 
                                    real* temp,
                                    real* velo,
                                    real diffusivity,
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
   if (isEvenTimestep==true)
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

   if(k<numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB;//, 

      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
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
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
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
      /*real drho*/;
      real vx1_Inflow   = c0o1;
      real vx2_Inflow   = c0o1;
      real vx3_Inflow   = velo[k];
      real ux_sq_Inflow = vx1_Inflow * vx1_Inflow;
      real uy_sq_Inflow = vx2_Inflow * vx2_Inflow;
      real uz_sq_Inflow = vx3_Inflow * vx3_Inflow;


      ////drho   =    f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
      ////            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
      ////            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]); 

      //real vx1 =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //               ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //               (f_E - f_W); 

      //real vx2 =  (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //               ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //               (f_N - f_S); 

      //real vx3 =  ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //               (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //               (f_T - f_B); 

      real rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ ((D.f[dirREST])[kzero]);
      real rho    =  rho0 + c1o1;
      real OORho  =  c1o1/rho;
      real vx1    =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      real vx2    =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      real vx3    =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));

	  //real cu_sq       =1.5f*(vx1*vx1+vx2*vx2+vx3*vx3);
      real ux_sq       = vx1 * vx1;
      real uy_sq       = vx2 * vx2;
      real uz_sq       = vx3 * vx3;
      real omegaD     = c3o1 - sqrt(c3o1);
      real Lam         = -(c1o2-c1o1/omegaD);
      real nue_d       = Lam/c3o1;
      //real ae          = zero;
      real ae          = diffusivity/nue_d - c1o1;

      real f7_ZERO,f7_E,f7_W,f7_N,f7_S,f7_T,f7_B;
      real /*feq7_ZERO,*/feq7_E,feq7_W,feq7_N,feq7_S,feq7_T,feq7_B;
      real /*feqW7_ZERO,*/feqW7_E,feqW7_W,feqW7_N,feqW7_S,feqW7_T,feqW7_B;
      real TempD = temp[k];


      f7_ZERO =  (D7.f[0])[kzero];
      f7_W    =  (D7.f[1])[ke   ];
      f7_E    =  (D7.f[2])[kw   ];
      f7_S    =  (D7.f[3])[kn   ];
      f7_N    =  (D7.f[4])[ks   ];
      f7_B    =  (D7.f[5])[kt   ];
      f7_T    =  (D7.f[6])[kb   ];

      real ConcD = f7_ZERO + f7_E + f7_W + f7_N + f7_S + f7_T + f7_B;

      //feq7_ZERO = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
      feq7_E    = ConcD*(c1o6*(ae+c1o1)+c1o2*(ux_sq)+vx1*c1o2);
      feq7_W    = ConcD*(c1o6*(ae+c1o1)+c1o2*(ux_sq)-vx1*c1o2);
      feq7_N    = ConcD*(c1o6*(ae+c1o1)+c1o2*(uy_sq)+vx2*c1o2);
      feq7_S    = ConcD*(c1o6*(ae+c1o1)+c1o2*(uy_sq)-vx2*c1o2);
      feq7_T    = ConcD*(c1o6*(ae+c1o1)+c1o2*(uz_sq)+vx3*c1o2);
      feq7_B    = ConcD*(c1o6*(ae+c1o1)+c1o2*(uz_sq)-vx3*c1o2);

      //feq7_ZERO = TempD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
      feqW7_E    = TempD*(c1o6*(ae+c1o1)+c1o2*(ux_sq_Inflow)+vx1_Inflow*c1o2);
      feqW7_W    = TempD*(c1o6*(ae+c1o1)+c1o2*(ux_sq_Inflow)-vx1_Inflow*c1o2);
      feqW7_N    = TempD*(c1o6*(ae+c1o1)+c1o2*(uy_sq_Inflow)+vx2_Inflow*c1o2);
      feqW7_S    = TempD*(c1o6*(ae+c1o1)+c1o2*(uy_sq_Inflow)-vx2_Inflow*c1o2);
      feqW7_T    = TempD*(c1o6*(ae+c1o1)+c1o2*(uz_sq_Inflow)+vx3_Inflow*c1o2);
      feqW7_B    = TempD*(c1o6*(ae+c1o1)+c1o2*(uz_sq_Inflow)-vx3_Inflow*c1o2);

      //////////////////////////////////////////////////////////////////////////
      //pointertausch
      if (isEvenTimestep==false)
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
      //(D.f[dirREST])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //(D7.f[1])[ke   ] = f7_E - feq7_E + feqW7_W; //E
      //(D7.f[2])[kw   ] = f7_W - feq7_W + feqW7_E; //W
      //(D7.f[3])[kn   ] = f7_N - feq7_N + feqW7_S; //N
      //(D7.f[4])[ks   ] = f7_S - feq7_S + feqW7_N; //S
      //(D7.f[5])[kt   ] = f7_T - feq7_T + feqW7_B; //T
      //(D7.f[6])[kb   ] = f7_B - feq7_B + feqW7_T; //B

      //////////////////////////////////////////////////////////////////////////
      //mit Q's
      real /*feq,*/ q;
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         //q=0.;
         (D7.f[2])[kw]=(c2o1*feqW7_W-(f7_E*(q*omegaD-c1o1)-omegaD*feq7_E*(q-c1o1))/(omegaD-c1o1)+f7_W*q)/(q+c1o1);//f7_W - feq7_W + feqW7_E;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         //q=0.;
         (D7.f[1])[ke]=(c2o1*feqW7_E-(f7_W*(q*omegaD-c1o1)-omegaD*feq7_W*(q-c1o1))/(omegaD-c1o1)+f7_E*q)/(q+c1o1);//f7_E - feq7_E + feqW7_W;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         //q=0.;
         (D7.f[4])[ks]=(c2o1*feqW7_S-(f7_N*(q*omegaD-c1o1)-omegaD*feq7_N*(q-c1o1))/(omegaD-c1o1)+f7_S*q)/(q+c1o1);//f7_S - feq7_S + feqW7_N;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         //q=0.;
         (D7.f[3])[kn]=(c2o1*feqW7_N-(f7_S*(q*omegaD-c1o1)-omegaD*feq7_S*(q-c1o1))/(omegaD-c1o1)+f7_N*q)/(q+c1o1);//f7_N - feq7_N + feqW7_S;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         //q=0.;
         (D7.f[6])[kb]=(c2o1*feqW7_B-(f7_T*(q*omegaD-c1o1)-omegaD*feq7_T*(q-c1o1))/(omegaD-c1o1)+f7_B*q)/(q+c1o1);//f7_B - feq7_B + feqW7_T;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         //q=0.;
         (D7.f[5])[kt]=(c2o1*feqW7_T-(f7_B*(q*omegaD-c1o1)-omegaD*feq7_B*(q-c1o1))/(omegaD-c1o1)+f7_T*q)/(q+c1o1);//f7_T - feq7_T + feqW7_B;
      }

      ////////////////////////////////////////////////////////////////////////////
      ////ohne Q's
      //real /*feq,*/ q;
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
extern "C" __global__ void QADVel27(real* DD, 
                                    real* DD27, 
                                    real* temp,
                                    real* velo,
                                    real diffusivity,
                                    int* k_Q, 
                                    real* QQ,
                                    int numberOfBCnodes, 
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
   if (isEvenTimestep==true)
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      q_dirNE  = &QQ[dirNE  * numberOfBCnodes];
      q_dirSW  = &QQ[dirSW  * numberOfBCnodes];
      q_dirSE  = &QQ[dirSE  * numberOfBCnodes];
      q_dirNW  = &QQ[dirNW  * numberOfBCnodes];
      q_dirTE  = &QQ[dirTE  * numberOfBCnodes];
      q_dirBW  = &QQ[dirBW  * numberOfBCnodes];
      q_dirBE  = &QQ[dirBE  * numberOfBCnodes];
      q_dirTW  = &QQ[dirTW  * numberOfBCnodes];
      q_dirTN  = &QQ[dirTN  * numberOfBCnodes];
      q_dirBS  = &QQ[dirBS  * numberOfBCnodes];
      q_dirBN  = &QQ[dirBN  * numberOfBCnodes];
      q_dirTS  = &QQ[dirTS  * numberOfBCnodes];
      q_dirTNE = &QQ[dirTNE * numberOfBCnodes];
      q_dirTSW = &QQ[dirTSW * numberOfBCnodes];
      q_dirTSE = &QQ[dirTSE * numberOfBCnodes];
      q_dirTNW = &QQ[dirTNW * numberOfBCnodes];
      q_dirBNE = &QQ[dirBNE * numberOfBCnodes];
      q_dirBSW = &QQ[dirBSW * numberOfBCnodes];
      q_dirBSE = &QQ[dirBSE * numberOfBCnodes];
      q_dirBNW = &QQ[dirBNW * numberOfBCnodes];
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
      real f_ZERO = (D.f[dirREST])[kzero];
      real f_BSW  = (D.f[dirTNE ])[ktne ];
      real f_BNE  = (D.f[dirTSW ])[ktsw ];
      real f_BNW  = (D.f[dirTSE ])[ktse ];
      real f_BSE  = (D.f[dirTNW ])[ktnw ];
      real f_TSW  = (D.f[dirBNE ])[kbne ];
      real f_TNE  = (D.f[dirBSW ])[kbsw ];
      real f_TNW  = (D.f[dirBSE ])[kbse ];
      real f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, /*drho, feq,*/ q;
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
      real rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ f_ZERO;
      real rho    =  rho0 + c1o1;
      real OORho  =  c1o1/rho;
      vx1     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      vx2     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      vx3     =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      //real f27_W    = (D27.f[dirE   ])[ke   ];
      //real f27_E    = (D27.f[dirW   ])[kw   ];
      //real f27_S    = (D27.f[dirN   ])[kn   ];
      //real f27_N    = (D27.f[dirS   ])[ks   ];
      //real f27_B    = (D27.f[dirT   ])[kt   ];
      //real f27_T    = (D27.f[dirB   ])[kb   ];
      //real f27_SW   = (D27.f[dirNE  ])[kne  ];
      //real f27_NE   = (D27.f[dirSW  ])[ksw  ];
      //real f27_NW   = (D27.f[dirSE  ])[kse  ];
      //real f27_SE   = (D27.f[dirNW  ])[knw  ];
      //real f27_BW   = (D27.f[dirTE  ])[kte  ];
      //real f27_TE   = (D27.f[dirBW  ])[kbw  ];
      //real f27_TW   = (D27.f[dirBE  ])[kbe  ];
      //real f27_BE   = (D27.f[dirTW  ])[ktw  ];
      //real f27_BS   = (D27.f[dirTN  ])[ktn  ];
      //real f27_TN   = (D27.f[dirBS  ])[kbs  ];
      //real f27_TS   = (D27.f[dirBN  ])[kbn  ];
      //real f27_BN   = (D27.f[dirTS  ])[kts  ];
      //real f27_ZERO = (D27.f[dirREST])[kzero];
      //real f27_BSW  = (D27.f[dirTNE ])[ktne ];
      //real f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      //real f27_BNW  = (D27.f[dirTSE ])[ktse ];
      //real f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      //real f27_TSW  = (D27.f[dirBNE ])[kbne ];
      //real f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      //real f27_TNW  = (D27.f[dirBSE ])[kbse ];
      //real f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      //real ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
      //                  f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
      //                  f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      //real feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      //real feq27_E    =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      //real feq27_W    =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      //real feq27_N    =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      //real feq27_S    =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      //real feq27_T    =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      //real feq27_B    =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      //real feq27_NE   =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      //real feq27_SW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      //real feq27_SE   =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      //real feq27_NW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      //real feq27_TE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      //real feq27_BW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      //real feq27_BE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      //real feq27_TW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      //real feq27_TN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      //real feq27_BS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      //real feq27_BN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      //real feq27_TS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      //real feq27_TNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      //real feq27_BSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      //real feq27_BNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      //real feq27_TSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      //real feq27_TSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      //real feq27_BNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      //real feq27_BSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      //real feq27_TNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real TempD = temp[k];
      //real TempD = four;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // velocity inflow
      //vx1   = velo[k]; //zero;
      //vx2   = zero; //velo[k];//zero;//velo[k];
      //vx3   = zero;

      ////real feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      //real feqW27_E    =   c2over27* (TempD+(one+TempD)*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq));
      //real feqW27_W    =   c2over27* (TempD+(one+TempD)*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq));
      //real feqW27_N    =   c2over27* (TempD+(one+TempD)*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq));
      //real feqW27_S    =   c2over27* (TempD+(one+TempD)*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq));
      //real feqW27_T    =   c2over27* (TempD+(one+TempD)*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq));
      //real feqW27_B    =   c2over27* (TempD+(one+TempD)*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq));
      //real feqW27_NE   =   c1over54* (TempD+(one+TempD)*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      //real feqW27_SW   =   c1over54* (TempD+(one+TempD)*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      //real feqW27_SE   =   c1over54* (TempD+(one+TempD)*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      //real feqW27_NW   =   c1over54* (TempD+(one+TempD)*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      //real feqW27_TE   =   c1over54* (TempD+(one+TempD)*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      //real feqW27_BW   =   c1over54* (TempD+(one+TempD)*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      //real feqW27_BE   =   c1over54* (TempD+(one+TempD)*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      //real feqW27_TW   =   c1over54* (TempD+(one+TempD)*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      //real feqW27_TN   =   c1over54* (TempD+(one+TempD)*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      //real feqW27_BS   =   c1over54* (TempD+(one+TempD)*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      //real feqW27_BN   =   c1over54* (TempD+(one+TempD)*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      //real feqW27_TS   =   c1over54* (TempD+(one+TempD)*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      //real feqW27_TNE  =   c1over216*(TempD+(one+TempD)*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      //real feqW27_BSW  =   c1over216*(TempD+(one+TempD)*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      //real feqW27_BNE  =   c1over216*(TempD+(one+TempD)*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      //real feqW27_TSW  =   c1over216*(TempD+(one+TempD)*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      //real feqW27_TSE  =   c1over216*(TempD+(one+TempD)*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      //real feqW27_BNW  =   c1over216*(TempD+(one+TempD)*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      //real feqW27_BSE  =   c1over216*(TempD+(one+TempD)*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      //real feqW27_TNW  =   c1over216*(TempD+(one+TempD)*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));
      //
      real feqW27_E    =   c2o27* TempD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      real feqW27_W    =   c2o27* TempD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      real feqW27_N    =   c2o27* TempD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      real feqW27_S    =   c2o27* TempD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      real feqW27_T    =   c2o27* TempD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      real feqW27_B    =   c2o27* TempD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      real feqW27_NE   =   c1o54* TempD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      real feqW27_SW   =   c1o54* TempD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      real feqW27_SE   =   c1o54* TempD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      real feqW27_NW   =   c1o54* TempD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      real feqW27_TE   =   c1o54* TempD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      real feqW27_BW   =   c1o54* TempD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      real feqW27_BE   =   c1o54* TempD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      real feqW27_TW   =   c1o54* TempD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      real feqW27_TN   =   c1o54* TempD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      real feqW27_BS   =   c1o54* TempD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      real feqW27_BN   =   c1o54* TempD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      real feqW27_TS   =   c1o54* TempD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      real feqW27_TNE  =   c1o216*TempD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      real feqW27_BSW  =   c1o216*TempD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      real feqW27_BNE  =   c1o216*TempD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      real feqW27_TSW  =   c1o216*TempD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      real feqW27_TSE  =   c1o216*TempD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      real feqW27_BNW  =   c1o216*TempD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      real feqW27_BSE  =   c1o216*TempD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      real feqW27_TNW  =   c1o216*TempD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real omegaD     = c3o1 - sqrt(c3o1);
      //real Lam        = -(c1o2-one/omegaD);
      //real nue_d      = Lam/three;
      //real ae         = zero;
      //real ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      //(D.f[dirREST])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirREST])[k]=c1o10;
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
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirW  ])[kw  ]= -feqW27_W  + c2o1 * c2o27  * TempD;
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirE  ])[ke  ]= -feqW27_E  + c2o1 * c2o27  * TempD;
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirS  ])[ks  ]= -feqW27_S  + c2o1 * c2o27  * TempD;
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirN  ])[kn  ]= -feqW27_N  + c2o1 * c2o27  * TempD;
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirB  ])[kb  ]= -feqW27_B  + c2o1 * c2o27  * TempD;
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirT  ])[kt  ]= -feqW27_T  + c2o1 * c2o27  * TempD;
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirSW ])[ksw ]= -feqW27_SW + c2o1 * c1o54  * TempD;
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirNE ])[kne ]= -feqW27_NE + c2o1 * c1o54  * TempD;
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirNW ])[knw ]= -feqW27_NW + c2o1 * c1o54  * TempD;
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirSE ])[kse ]= -feqW27_SE + c2o1 * c1o54  * TempD;
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBW ])[kbw ]= -feqW27_BW + c2o1 * c1o54  * TempD;
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTE ])[kte ]= -feqW27_TE + c2o1 * c1o54  * TempD;
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTW ])[ktw ]= -feqW27_TW + c2o1 * c1o54  * TempD;
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBE ])[kbe ]= -feqW27_BE + c2o1 * c1o54  * TempD;
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBS ])[kbs ]= -feqW27_BS + c2o1 * c1o54  * TempD;
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTN ])[ktn ]= -feqW27_TN + c2o1 * c1o54  * TempD;
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTS ])[kts ]= -feqW27_TS + c2o1 * c1o54  * TempD;
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBN ])[kbn ]= -feqW27_BN + c2o1 * c1o54  * TempD;
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSW])[kbsw]= -feqW27_BSW+ c2o1 * c1o216 * TempD;
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNE])[ktne]= -feqW27_TNE+ c2o1 * c1o216 * TempD;
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSW])[ktsw]= -feqW27_TSW+ c2o1 * c1o216 * TempD;
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNE])[kbne]= -feqW27_BNE+ c2o1 * c1o216 * TempD;
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNW])[kbnw]= -feqW27_BNW+ c2o1 * c1o216 * TempD;
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSE])[ktse]= -feqW27_TSE+ c2o1 * c1o216 * TempD;
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNW])[ktnw]= -feqW27_TNW+ c2o1 * c1o216 * TempD;
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSE])[kbse]= -feqW27_BSE+ c2o1 * c1o216 * TempD;
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
extern "C" __global__ void QAD7( real* DD, 
                                 real* DD7, 
                                 real* temp,
                                 real diffusivity,
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
   if (isEvenTimestep==true)
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

   if(k<numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB;//, 
      //         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
      //         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
      //         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
      //         *q_dirBSE, *q_dirBNW;

      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      //q_dirNE  = &QQ[dirNE  * numberOfBCnodes];
      //q_dirSW  = &QQ[dirSW  * numberOfBCnodes];
      //q_dirSE  = &QQ[dirSE  * numberOfBCnodes];
      //q_dirNW  = &QQ[dirNW  * numberOfBCnodes];
      //q_dirTE  = &QQ[dirTE  * numberOfBCnodes];
      //q_dirBW  = &QQ[dirBW  * numberOfBCnodes];
      //q_dirBE  = &QQ[dirBE  * numberOfBCnodes];
      //q_dirTW  = &QQ[dirTW  * numberOfBCnodes];
      //q_dirTN  = &QQ[dirTN  * numberOfBCnodes];
      //q_dirBS  = &QQ[dirBS  * numberOfBCnodes];
      //q_dirBN  = &QQ[dirBN  * numberOfBCnodes];
      //q_dirTS  = &QQ[dirTS  * numberOfBCnodes];
      //q_dirTNE = &QQ[dirTNE * numberOfBCnodes];
      //q_dirTSW = &QQ[dirTSW * numberOfBCnodes];
      //q_dirTSE = &QQ[dirTSE * numberOfBCnodes];
      //q_dirTNW = &QQ[dirTNW * numberOfBCnodes];
      //q_dirBNE = &QQ[dirBNE * numberOfBCnodes];
      //q_dirBSW = &QQ[dirBSW * numberOfBCnodes];
      //q_dirBSE = &QQ[dirBSE * numberOfBCnodes];
      //q_dirBNW = &QQ[dirBNW * numberOfBCnodes];
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
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
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
      real vx1, vx2, vx3/*, drho*/;
      //drho   =    f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
      //            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
      //            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]); 

      //vx1    = ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //         (f_E - f_W); 


      //vx2    = (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //         (f_N - f_S); 

      //vx3    = ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //         (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //         (f_T - f_B); 
      real rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ ((D.f[dirREST])[kzero]);
      real rho    =  rho0 + c1o1;
      real OORho  =  c1o1/rho;
      vx1     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      vx2     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      vx3     =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));

      //real cu_sq       =c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      real ux_sq       = vx1 * vx1;
      real uy_sq       = vx2 * vx2;
      real uz_sq       = vx3 * vx3;
      real omegaD     = c3o1 - sqrt(c3o1);
      real Lam         = -(c1o2-c1o1/omegaD);
      real nue_d       = Lam/c3o1;
      //real ae          = zero;
      real ae          = diffusivity/nue_d - c1o1;

      real f7_ZERO,f7_E,f7_W,f7_N,f7_S,f7_T,f7_B;
      real /*feq7_ZERO,*/feq7_E,feq7_W,feq7_N,feq7_S,feq7_T,feq7_B;
      real /*feqW7_ZERO,*/feqW7_E,feqW7_W,feqW7_N,feqW7_S,feqW7_T,feqW7_B;
      real TempD = temp[k];


      f7_ZERO =  (D7.f[0])[kzero];
      f7_W    =  (D7.f[1])[ke   ];
      f7_E    =  (D7.f[2])[kw   ];
      f7_S    =  (D7.f[3])[kn   ];
      f7_N    =  (D7.f[4])[ks   ];
      f7_B    =  (D7.f[5])[kt   ];
      f7_T    =  (D7.f[6])[kb   ];

      real ConcD = f7_ZERO + f7_E + f7_W + f7_N + f7_S + f7_T + f7_B;

      //feq7_ZERO = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
      feq7_E    = ConcD*(c1o6*(ae+c1o1)+c1o2*(ux_sq)+vx1*c1o2);
      feq7_W    = ConcD*(c1o6*(ae+c1o1)+c1o2*(ux_sq)-vx1*c1o2);
      feq7_N    = ConcD*(c1o6*(ae+c1o1)+c1o2*(uy_sq)+vx2*c1o2);
      feq7_S    = ConcD*(c1o6*(ae+c1o1)+c1o2*(uy_sq)-vx2*c1o2);
      feq7_T    = ConcD*(c1o6*(ae+c1o1)+c1o2*(uz_sq)+vx3*c1o2);
      feq7_B    = ConcD*(c1o6*(ae+c1o1)+c1o2*(uz_sq)-vx3*c1o2);

      //feq7_ZERO = TempD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
      feqW7_E    = TempD*(c1o6*(ae+c1o1));//+c1o2*(ux_sq)+vx1*c1o2);
      feqW7_W    = TempD*(c1o6*(ae+c1o1));//+c1o2*(ux_sq)-vx1*c1o2);
      feqW7_N    = TempD*(c1o6*(ae+c1o1));//+c1o2*(uy_sq)+vx2*c1o2);
      feqW7_S    = TempD*(c1o6*(ae+c1o1));//+c1o2*(uy_sq)-vx2*c1o2);
      feqW7_T    = TempD*(c1o6*(ae+c1o1));//+c1o2*(uz_sq)+vx3*c1o2);
      feqW7_B    = TempD*(c1o6*(ae+c1o1));//+c1o2*(uz_sq)-vx3*c1o2);

      //////////////////////////////////////////////////////////////////////////
      //pointertausch
      if (isEvenTimestep==false)
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
      //(D.f[dirREST])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //(D7.f[1])[ke   ] = f7_E - feq7_E + feqW7_W; //E
      //(D7.f[2])[kw   ] = f7_W - feq7_W + feqW7_E; //W
      //(D7.f[3])[kn   ] = f7_N - feq7_N + feqW7_S; //N
      //(D7.f[4])[ks   ] = f7_S - feq7_S + feqW7_N; //S
      //(D7.f[5])[kt   ] = f7_T - feq7_T + feqW7_B; //T
      //(D7.f[6])[kb   ] = f7_B - feq7_B + feqW7_T; //B

      ////////////////////////////////////////////////////////////////////////////
      ////mit Q's
      //real /*feq,*/ q;
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
      real /*feq,*/ q;
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[2])[kw]=f7_W - feq7_W + feqW7_E;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[1])[ke]=f7_E - feq7_E + feqW7_W;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[4])[ks]=f7_S - feq7_S + feqW7_N;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[3])[kn]=f7_N - feq7_N + feqW7_S;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[6])[kb]=f7_B - feq7_B + feqW7_T;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[5])[kt]=f7_T - feq7_T + feqW7_B;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADDirichlet27(
											 real* DD, 
											 real* DD27, 
											 real* temp,
											 real diffusivity,
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
   if (isEvenTimestep==true)
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      q_dirNE  = &QQ[dirNE  * numberOfBCnodes];
      q_dirSW  = &QQ[dirSW  * numberOfBCnodes];
      q_dirSE  = &QQ[dirSE  * numberOfBCnodes];
      q_dirNW  = &QQ[dirNW  * numberOfBCnodes];
      q_dirTE  = &QQ[dirTE  * numberOfBCnodes];
      q_dirBW  = &QQ[dirBW  * numberOfBCnodes];
      q_dirBE  = &QQ[dirBE  * numberOfBCnodes];
      q_dirTW  = &QQ[dirTW  * numberOfBCnodes];
      q_dirTN  = &QQ[dirTN  * numberOfBCnodes];
      q_dirBS  = &QQ[dirBS  * numberOfBCnodes];
      q_dirBN  = &QQ[dirBN  * numberOfBCnodes];
      q_dirTS  = &QQ[dirTS  * numberOfBCnodes];
      q_dirTNE = &QQ[dirTNE * numberOfBCnodes];
      q_dirTSW = &QQ[dirTSW * numberOfBCnodes];
      q_dirTSE = &QQ[dirTSE * numberOfBCnodes];
      q_dirTNW = &QQ[dirTNW * numberOfBCnodes];
      q_dirBNE = &QQ[dirBNE * numberOfBCnodes];
      q_dirBSW = &QQ[dirBSW * numberOfBCnodes];
      q_dirBSE = &QQ[dirBSE * numberOfBCnodes];
      q_dirBNW = &QQ[dirBNW * numberOfBCnodes];
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
      real f_ZERO = (D.f[dirREST])[kzero];
      real f_BSW  = (D.f[dirTNE ])[ktne ];
      real f_BNE  = (D.f[dirTSW ])[ktsw ];
      real f_BNW  = (D.f[dirTSE ])[ktse ];
      real f_BSE  = (D.f[dirTNW ])[ktnw ];
      real f_TSW  = (D.f[dirBNE ])[kbne ];
      real f_TNE  = (D.f[dirBSW ])[kbsw ];
      real f_TNW  = (D.f[dirBSE ])[kbse ];
      real f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, /*drho, feq,*/ q;
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
      real rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+f_ZERO;
      real rho    =  rho0 + c1o1;
      real OORho  =  c1o1/rho;
      vx1     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      vx2     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      vx3     =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      real f27_W    = (D27.f[dirE   ])[ke   ];
      real f27_E    = (D27.f[dirW   ])[kw   ];
      real f27_S    = (D27.f[dirN   ])[kn   ];
      real f27_N    = (D27.f[dirS   ])[ks   ];
      real f27_B    = (D27.f[dirT   ])[kt   ];
      real f27_T    = (D27.f[dirB   ])[kb   ];
      real f27_SW   = (D27.f[dirNE  ])[kne  ];
      real f27_NE   = (D27.f[dirSW  ])[ksw  ];
      real f27_NW   = (D27.f[dirSE  ])[kse  ];
      real f27_SE   = (D27.f[dirNW  ])[knw  ];
      real f27_BW   = (D27.f[dirTE  ])[kte  ];
      real f27_TE   = (D27.f[dirBW  ])[kbw  ];
      real f27_TW   = (D27.f[dirBE  ])[kbe  ];
      real f27_BE   = (D27.f[dirTW  ])[ktw  ];
      real f27_BS   = (D27.f[dirTN  ])[ktn  ];
      real f27_TN   = (D27.f[dirBS  ])[kbs  ];
      real f27_TS   = (D27.f[dirBN  ])[kbn  ];
      real f27_BN   = (D27.f[dirTS  ])[kts  ];
      real f27_ZERO = (D27.f[dirREST])[kzero];
      real f27_BSW  = (D27.f[dirTNE ])[ktne ];
      real f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      real f27_BNW  = (D27.f[dirTSE ])[ktse ];
      real f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      real f27_TSW  = (D27.f[dirBNE ])[kbne ];
      real f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      real f27_TNW  = (D27.f[dirBSE ])[kbse ];
      real f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      real ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
         f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
         f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      //real feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      real feq27_E    =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      real feq27_W    =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      real feq27_N    =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      real feq27_S    =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      real feq27_T    =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      real feq27_B    =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      real feq27_NE   =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      real feq27_SW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      real feq27_SE   =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      real feq27_NW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      real feq27_TE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      real feq27_BW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      real feq27_BE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      real feq27_TW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      real feq27_TN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      real feq27_BS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      real feq27_BN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      real feq27_TS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      real feq27_TNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      real feq27_BSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      real feq27_BNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      real feq27_TSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      real feq27_TSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      real feq27_BNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      real feq27_BSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      real feq27_TNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real TempD = temp[k];//one;//temp[k];

      //real feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      real feqW27_E    =   c2o27* TempD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      real feqW27_W    =   c2o27* TempD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      real feqW27_N    =   c2o27* TempD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      real feqW27_S    =   c2o27* TempD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      real feqW27_T    =   c2o27* TempD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      real feqW27_B    =   c2o27* TempD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      real feqW27_NE   =   c1o54* TempD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      real feqW27_SW   =   c1o54* TempD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      real feqW27_SE   =   c1o54* TempD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      real feqW27_NW   =   c1o54* TempD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      real feqW27_TE   =   c1o54* TempD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      real feqW27_BW   =   c1o54* TempD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      real feqW27_BE   =   c1o54* TempD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      real feqW27_TW   =   c1o54* TempD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      real feqW27_TN   =   c1o54* TempD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      real feqW27_BS   =   c1o54* TempD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      real feqW27_BN   =   c1o54* TempD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      real feqW27_TS   =   c1o54* TempD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      real feqW27_TNE  =   c1o216*TempD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      real feqW27_BSW  =   c1o216*TempD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      real feqW27_BNE  =   c1o216*TempD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      real feqW27_TSW  =   c1o216*TempD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      real feqW27_TSE  =   c1o216*TempD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      real feqW27_BNW  =   c1o216*TempD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      real feqW27_BSE  =   c1o216*TempD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      real feqW27_TNW  =   c1o216*TempD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real omegaD     = c3o1 - sqrt(c3o1);
      //real Lam        = -(c1o2-one/omegaD);
      //real nue_d      = Lam/three;
      //real ae         = zero;
      //real ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      //(D.f[dirREST])[k]=0.1f;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[  ke   ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirW  ])[kw  ]=(c2o1*feqW27_W  -(f27_E  *(q*omegaD-c1o1)-omegaD*feq27_E  *(q-c1o1))/(omegaD-c1o1)+f27_W  *q)/(q+c1o1);
      q = q_dirW[  kw   ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirE  ])[ke  ]=(c2o1*feqW27_E  -(f27_W  *(q*omegaD-c1o1)-omegaD*feq27_W  *(q-c1o1))/(omegaD-c1o1)+f27_E  *q)/(q+c1o1);
      q = q_dirN[  kn   ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirS  ])[ks  ]=(c2o1*feqW27_S  -(f27_N  *(q*omegaD-c1o1)-omegaD*feq27_N  *(q-c1o1))/(omegaD-c1o1)+f27_S  *q)/(q+c1o1);
      q = q_dirS[  ks   ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirN  ])[kn  ]=(c2o1*feqW27_N  -(f27_S  *(q*omegaD-c1o1)-omegaD*feq27_S  *(q-c1o1))/(omegaD-c1o1)+f27_N  *q)/(q+c1o1);
      q = q_dirT[  kt   ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirB  ])[kb  ]=(c2o1*feqW27_B  -(f27_T  *(q*omegaD-c1o1)-omegaD*feq27_T  *(q-c1o1))/(omegaD-c1o1)+f27_B  *q)/(q+c1o1);
      q = q_dirB[  kb   ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirT  ])[kt  ]=(c2o1*feqW27_T  -(f27_B  *(q*omegaD-c1o1)-omegaD*feq27_B  *(q-c1o1))/(omegaD-c1o1)+f27_T  *q)/(q+c1o1);
      q = q_dirNE[ kne  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirSW ])[ksw ]=(c2o1*feqW27_SW -(f27_NE *(q*omegaD-c1o1)-omegaD*feq27_NE *(q-c1o1))/(omegaD-c1o1)+f27_SW *q)/(q+c1o1);
      q = q_dirSW[ ksw  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirNE ])[kne ]=(c2o1*feqW27_NE -(f27_SW *(q*omegaD-c1o1)-omegaD*feq27_SW *(q-c1o1))/(omegaD-c1o1)+f27_NE *q)/(q+c1o1);
      q = q_dirSE[ kse  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirNW ])[knw ]=(c2o1*feqW27_NW -(f27_SE *(q*omegaD-c1o1)-omegaD*feq27_SE *(q-c1o1))/(omegaD-c1o1)+f27_NW *q)/(q+c1o1);
      q = q_dirNW[ knw  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirSE ])[kse ]=(c2o1*feqW27_SE -(f27_NW *(q*omegaD-c1o1)-omegaD*feq27_NW *(q-c1o1))/(omegaD-c1o1)+f27_SE *q)/(q+c1o1);
      q = q_dirTE[ kte  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBW ])[kbw ]=(c2o1*feqW27_BW -(f27_TE *(q*omegaD-c1o1)-omegaD*feq27_TE *(q-c1o1))/(omegaD-c1o1)+f27_BW *q)/(q+c1o1);
      q = q_dirBW[ kbw  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTE ])[kte ]=(c2o1*feqW27_TE -(f27_BW *(q*omegaD-c1o1)-omegaD*feq27_BW *(q-c1o1))/(omegaD-c1o1)+f27_TE *q)/(q+c1o1);
      q = q_dirBE[ kbe  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTW ])[ktw ]=(c2o1*feqW27_TW -(f27_BE *(q*omegaD-c1o1)-omegaD*feq27_BE *(q-c1o1))/(omegaD-c1o1)+f27_TW *q)/(q+c1o1);
      q = q_dirTW[ ktw  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBE ])[kbe ]=(c2o1*feqW27_BE -(f27_TW *(q*omegaD-c1o1)-omegaD*feq27_TW *(q-c1o1))/(omegaD-c1o1)+f27_BE *q)/(q+c1o1);
      q = q_dirTN[ ktn  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBS ])[kbs ]=(c2o1*feqW27_BS -(f27_TN *(q*omegaD-c1o1)-omegaD*feq27_TN *(q-c1o1))/(omegaD-c1o1)+f27_BS *q)/(q+c1o1);
      q = q_dirBS[ kbs  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTN ])[ktn ]=(c2o1*feqW27_TN -(f27_BS *(q*omegaD-c1o1)-omegaD*feq27_BS *(q-c1o1))/(omegaD-c1o1)+f27_TN *q)/(q+c1o1);
      q = q_dirBN[ kbn  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTS ])[kts ]=(c2o1*feqW27_TS -(f27_BN *(q*omegaD-c1o1)-omegaD*feq27_BN *(q-c1o1))/(omegaD-c1o1)+f27_TS *q)/(q+c1o1);
      q = q_dirTS[ kts  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBN ])[kbn ]=(c2o1*feqW27_BN -(f27_TS *(q*omegaD-c1o1)-omegaD*feq27_TS *(q-c1o1))/(omegaD-c1o1)+f27_BN *q)/(q+c1o1);
      q = q_dirTNE[ktne ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSW])[kbsw]=(c2o1*feqW27_BSW-(f27_TNE*(q*omegaD-c1o1)-omegaD*feq27_TNE*(q-c1o1))/(omegaD-c1o1)+f27_BSW*q)/(q+c1o1);
      q = q_dirBSW[kbsw ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNE])[ktne]=(c2o1*feqW27_TNE-(f27_BSW*(q*omegaD-c1o1)-omegaD*feq27_BSW*(q-c1o1))/(omegaD-c1o1)+f27_TNE*q)/(q+c1o1);
      q = q_dirBNE[kbne ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSW])[ktsw]=(c2o1*feqW27_TSW-(f27_BNE*(q*omegaD-c1o1)-omegaD*feq27_BNE*(q-c1o1))/(omegaD-c1o1)+f27_TSW*q)/(q+c1o1);
      q = q_dirTSW[ktsw ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNE])[kbne]=(c2o1*feqW27_BNE-(f27_TSW*(q*omegaD-c1o1)-omegaD*feq27_TSW*(q-c1o1))/(omegaD-c1o1)+f27_BNE*q)/(q+c1o1);
      q = q_dirTSE[ktse ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNW])[kbnw]=(c2o1*feqW27_BNW-(f27_TSE*(q*omegaD-c1o1)-omegaD*feq27_TSE*(q-c1o1))/(omegaD-c1o1)+f27_BNW*q)/(q+c1o1);
      q = q_dirBNW[kbnw ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSE])[ktse]=(c2o1*feqW27_TSE-(f27_BNW*(q*omegaD-c1o1)-omegaD*feq27_BNW*(q-c1o1))/(omegaD-c1o1)+f27_TSE*q)/(q+c1o1);
      q = q_dirBSE[kbse ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNW])[ktnw]=(c2o1*feqW27_TNW-(f27_BSE*(q*omegaD-c1o1)-omegaD*feq27_BSE*(q-c1o1))/(omegaD-c1o1)+f27_TNW*q)/(q+c1o1);
      q = q_dirTNW[ktnw ]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSE])[kbse]=(c2o1*feqW27_BSE-(f27_TNW*(q*omegaD-c1o1)-omegaD*feq27_TNW*(q-c1o1))/(omegaD-c1o1)+f27_BSE*q)/(q+c1o1);
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
extern "C" __global__ void QADBB27( real* DD, 
                                   real* DD27, 
                                   real* temp,
                                   real diffusivity,
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
   //Distributions27 D;
   //if (isEvenTimestep==true)
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
   //   D.f[dirREST] = &DD[dirREST*size_Mat];
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
   //   D.f[dirREST] = &DD[dirREST*size_Mat];
   //   D.f[dirTNE ] = &DD[dirBSW *size_Mat];
   //   D.f[dirTSW ] = &DD[dirBNE *size_Mat];
   //   D.f[dirTSE ] = &DD[dirBNW *size_Mat];
   //   D.f[dirTNW ] = &DD[dirBSE *size_Mat];
   //   D.f[dirBNE ] = &DD[dirTSW *size_Mat];
   //   D.f[dirBSW ] = &DD[dirTNE *size_Mat];
   //   D.f[dirBSE ] = &DD[dirTNW *size_Mat];
   //   D.f[dirBNW ] = &DD[dirTSE *size_Mat];
   //}

   Distributions27 D27;
   if (isEvenTimestep==true)
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      q_dirNE  = &QQ[dirNE  * numberOfBCnodes];
      q_dirSW  = &QQ[dirSW  * numberOfBCnodes];
      q_dirSE  = &QQ[dirSE  * numberOfBCnodes];
      q_dirNW  = &QQ[dirNW  * numberOfBCnodes];
      q_dirTE  = &QQ[dirTE  * numberOfBCnodes];
      q_dirBW  = &QQ[dirBW  * numberOfBCnodes];
      q_dirBE  = &QQ[dirBE  * numberOfBCnodes];
      q_dirTW  = &QQ[dirTW  * numberOfBCnodes];
      q_dirTN  = &QQ[dirTN  * numberOfBCnodes];
      q_dirBS  = &QQ[dirBS  * numberOfBCnodes];
      q_dirBN  = &QQ[dirBN  * numberOfBCnodes];
      q_dirTS  = &QQ[dirTS  * numberOfBCnodes];
      q_dirTNE = &QQ[dirTNE * numberOfBCnodes];
      q_dirTSW = &QQ[dirTSW * numberOfBCnodes];
      q_dirTSE = &QQ[dirTSE * numberOfBCnodes];
      q_dirTNW = &QQ[dirTNW * numberOfBCnodes];
      q_dirBNE = &QQ[dirBNE * numberOfBCnodes];
      q_dirBSW = &QQ[dirBSW * numberOfBCnodes];
      q_dirBSE = &QQ[dirBSE * numberOfBCnodes];
      q_dirBNW = &QQ[dirBNW * numberOfBCnodes];
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
      ////////////////////////////////////////////////////////////////////////////////
      //real f_W    = (D.f[dirE   ])[ke   ];
      //real f_E    = (D.f[dirW   ])[kw   ];
      //real f_S    = (D.f[dirN   ])[kn   ];
      //real f_N    = (D.f[dirS   ])[ks   ];
      //real f_B    = (D.f[dirT   ])[kt   ];
      //real f_T    = (D.f[dirB   ])[kb   ];
      //real f_SW   = (D.f[dirNE  ])[kne  ];
      //real f_NE   = (D.f[dirSW  ])[ksw  ];
      //real f_NW   = (D.f[dirSE  ])[kse  ];
      //real f_SE   = (D.f[dirNW  ])[knw  ];
      //real f_BW   = (D.f[dirTE  ])[kte  ];
      //real f_TE   = (D.f[dirBW  ])[kbw  ];
      //real f_TW   = (D.f[dirBE  ])[kbe  ];
      //real f_BE   = (D.f[dirTW  ])[ktw  ];
      //real f_BS   = (D.f[dirTN  ])[ktn  ];
      //real f_TN   = (D.f[dirBS  ])[kbs  ];
      //real f_TS   = (D.f[dirBN  ])[kbn  ];
      //real f_BN   = (D.f[dirTS  ])[kts  ];
      //real f_ZERO = (D.f[dirREST])[kzero];
      //real f_BSW  = (D.f[dirTNE ])[ktne ];
      //real f_BNE  = (D.f[dirTSW ])[ktsw ];
      //real f_BNW  = (D.f[dirTSE ])[ktse ];
      //real f_BSE  = (D.f[dirTNW ])[ktnw ];
      //real f_TSW  = (D.f[dirBNE ])[kbne ];
      //real f_TNE  = (D.f[dirBSW ])[kbsw ];
      //real f_TNW  = (D.f[dirBSE ])[kbse ];
      //real f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      //real vx1, vx2, vx3, /*drho, feq,*/ q;
      real q;
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
      //real rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+f_ZERO;
      //real rho    =  rho0 + c1o1;
      //real OORho  =  c1o1/rho;
      //vx1     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      //vx2     =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      //vx3     =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      real f27_W    = (D27.f[dirE   ])[ke   ];
      real f27_E    = (D27.f[dirW   ])[kw   ];
      real f27_S    = (D27.f[dirN   ])[kn   ];
      real f27_N    = (D27.f[dirS   ])[ks   ];
      real f27_B    = (D27.f[dirT   ])[kt   ];
      real f27_T    = (D27.f[dirB   ])[kb   ];
      real f27_SW   = (D27.f[dirNE  ])[kne  ];
      real f27_NE   = (D27.f[dirSW  ])[ksw  ];
      real f27_NW   = (D27.f[dirSE  ])[kse  ];
      real f27_SE   = (D27.f[dirNW  ])[knw  ];
      real f27_BW   = (D27.f[dirTE  ])[kte  ];
      real f27_TE   = (D27.f[dirBW  ])[kbw  ];
      real f27_TW   = (D27.f[dirBE  ])[kbe  ];
      real f27_BE   = (D27.f[dirTW  ])[ktw  ];
      real f27_BS   = (D27.f[dirTN  ])[ktn  ];
      real f27_TN   = (D27.f[dirBS  ])[kbs  ];
      real f27_TS   = (D27.f[dirBN  ])[kbn  ];
      real f27_BN   = (D27.f[dirTS  ])[kts  ];
      //real f27_ZERO = (D27.f[dirREST])[kzero];
      real f27_BSW  = (D27.f[dirTNE ])[ktne ];
      real f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      real f27_BNW  = (D27.f[dirTSE ])[ktse ];
      real f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      real f27_TSW  = (D27.f[dirBNE ])[kbne ];
      real f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      real f27_TNW  = (D27.f[dirBSE ])[kbse ];
      real f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      //real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      //real ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
      //   f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
      //   f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      //real feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      //real feq27_E    =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      //real feq27_W    =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      //real feq27_N    =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      //real feq27_S    =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      //real feq27_T    =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      //real feq27_B    =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      //real feq27_NE   =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      //real feq27_SW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      //real feq27_SE   =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      //real feq27_NW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      //real feq27_TE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      //real feq27_BW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      //real feq27_BE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      //real feq27_TW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      //real feq27_TN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      //real feq27_BS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      //real feq27_BN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      //real feq27_TS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      //real feq27_TNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      //real feq27_BSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      //real feq27_BNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      //real feq27_TSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      //real feq27_TSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      //real feq27_BNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      //real feq27_BSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      //real feq27_TNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //real TempD = temp[k];

      //real feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      //real feqW27_E    =   c2o27* TempD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      //real feqW27_W    =   c2o27* TempD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      //real feqW27_N    =   c2o27* TempD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      //real feqW27_S    =   c2o27* TempD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      //real feqW27_T    =   c2o27* TempD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      //real feqW27_B    =   c2o27* TempD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      //real feqW27_NE   =   c1o54* TempD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      //real feqW27_SW   =   c1o54* TempD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      //real feqW27_SE   =   c1o54* TempD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      //real feqW27_NW   =   c1o54* TempD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      //real feqW27_TE   =   c1o54* TempD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      //real feqW27_BW   =   c1o54* TempD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      //real feqW27_BE   =   c1o54* TempD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      //real feqW27_TW   =   c1o54* TempD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      //real feqW27_TN   =   c1o54* TempD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      //real feqW27_BS   =   c1o54* TempD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      //real feqW27_BN   =   c1o54* TempD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      //real feqW27_TS   =   c1o54* TempD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      //real feqW27_TNE  =   c1o216*TempD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      //real feqW27_BSW  =   c1o216*TempD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      //real feqW27_BNE  =   c1o216*TempD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      //real feqW27_TSW  =   c1o216*TempD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      //real feqW27_TSE  =   c1o216*TempD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      //real feqW27_BNW  =   c1o216*TempD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      //real feqW27_BSE  =   c1o216*TempD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      //real feqW27_TNW  =   c1o216*TempD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real omegaD     = c3o1 - sqrt(c3o1);
      //real Lam        = -(c1o2-one/omegaD);
      //real nue_d      = Lam/three;
      //real ae         = zero;
      //real ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      //(D.f[dirREST])[k]=0.1f;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirW  ])[kw  ]=f27_E  ;
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirE  ])[ke  ]=f27_W  ;
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirS  ])[ks  ]=f27_N  ;
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirN  ])[kn  ]=f27_S  ;
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirB  ])[kb  ]=f27_T  ;
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirT  ])[kt  ]=f27_B  ;
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirSW ])[ksw ]=f27_NE ;
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirNE ])[kne ]=f27_SW ;
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirNW ])[knw ]=f27_SE ;
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirSE ])[kse ]=f27_NW ;
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBW ])[kbw ]=f27_TE ;
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTE ])[kte ]=f27_BW ;
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTW ])[ktw ]=f27_BE ;
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBE ])[kbe ]=f27_TW ;
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBS ])[kbs ]=f27_TN ;
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTN ])[ktn ]=f27_BS ;
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTS ])[kts ]=f27_BN ;
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBN ])[kbn ]=f27_TS ;
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSW])[kbsw]=f27_TNE;
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNE])[ktne]=f27_BSW;
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSW])[ktsw]=f27_BNE;
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNE])[kbne]=f27_TSW;
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNW])[kbnw]=f27_TSE;
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSE])[ktse]=f27_BNW;
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNW])[ktnw]=f27_BSE;
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSE])[kbse]=f27_TNW;
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
extern "C" __global__ void QNoSlipADincomp7(
											 real* DD, 
											 real* DD7, 
											 real* temp,
											 real diffusivity,
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
   //Distributions27 D;
   //if (isEvenTimestep==true)
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
   //   D.f[dirREST] = &DD[dirREST*size_Mat];
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
   //   D.f[dirREST] = &DD[dirREST*size_Mat];
   //   D.f[dirTNE ] = &DD[dirBSW *size_Mat];
   //   D.f[dirTSW ] = &DD[dirBNE *size_Mat];
   //   D.f[dirTSE ] = &DD[dirBNW *size_Mat];
   //   D.f[dirTNW ] = &DD[dirBSE *size_Mat];
   //   D.f[dirBNE ] = &DD[dirTSW *size_Mat];
   //   D.f[dirBSW ] = &DD[dirTNE *size_Mat];
   //   D.f[dirBSE ] = &DD[dirTNW *size_Mat];
   //   D.f[dirBNW ] = &DD[dirTSE *size_Mat];
   //}

   Distributions7 D7;
   if (isEvenTimestep==true)
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

   if(k<numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB;

      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      //////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      //unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      //unsigned int ksw  = neighborY[kw];
      //unsigned int kne  = KQK;
      //unsigned int kse  = ks;
      //unsigned int knw  = kw;
      //unsigned int kbw  = neighborZ[kw];
      //unsigned int kte  = KQK;
      //unsigned int kbe  = kb;
      //unsigned int ktw  = kw;
      //unsigned int kbs  = neighborZ[ks];
      //unsigned int ktn  = KQK;
      //unsigned int kbn  = kb;
      //unsigned int kts  = ks;
      //unsigned int ktse = ks;
      //unsigned int kbnw = kbw;
      //unsigned int ktnw = kw;
      //unsigned int kbse = kbs;
      //unsigned int ktsw = ksw;
      //unsigned int kbne = kb;
      //unsigned int ktne = KQK;
      //unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      //real f_W    = (D.f[dirE   ])[ke   ];
      //real f_E    = (D.f[dirW   ])[kw   ];
      //real f_S    = (D.f[dirN   ])[kn   ];
      //real f_N    = (D.f[dirS   ])[ks   ];
      //real f_B    = (D.f[dirT   ])[kt   ];
      //real f_T    = (D.f[dirB   ])[kb   ];
      //real f_SW   = (D.f[dirNE  ])[kne  ];
      //real f_NE   = (D.f[dirSW  ])[ksw  ];
      //real f_NW   = (D.f[dirSE  ])[kse  ];
      //real f_SE   = (D.f[dirNW  ])[knw  ];
      //real f_BW   = (D.f[dirTE  ])[kte  ];
      //real f_TE   = (D.f[dirBW  ])[kbw  ];
      //real f_TW   = (D.f[dirBE  ])[kbe  ];
      //real f_BE   = (D.f[dirTW  ])[ktw  ];
      //real f_BS   = (D.f[dirTN  ])[ktn  ];
      //real f_TN   = (D.f[dirBS  ])[kbs  ];
      //real f_TS   = (D.f[dirBN  ])[kbn  ];
      //real f_BN   = (D.f[dirTS  ])[kts  ];
      //real f_BSW  = (D.f[dirTNE ])[ktne ];
      //real f_BNE  = (D.f[dirTSW ])[ktsw ];
      //real f_BNW  = (D.f[dirTSE ])[ktse ];
      //real f_BSE  = (D.f[dirTNW ])[ktnw ];
      //real f_TSW  = (D.f[dirBNE ])[kbne ];
      //real f_TNE  = (D.f[dirBSW ])[kbsw ];
      //real f_TNW  = (D.f[dirBSE ])[kbse ];
      //real f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      //real vx1 =  ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      //real vx2 =  ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      //real vx3 =  ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
		 ////d�rrrrrty !!!!!!!!!!!!!
   //      real vx1     =  ten * ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
   //      real vx2     =  ten * ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
   //      real vx3     =  ten * ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));

      //real cu_sq       =c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      //real ux_sq       = vx1 * vx1;
      //real uy_sq       = vx2 * vx2;
      //real uz_sq       = vx3 * vx3;
      ////////////////////////////////////////////////////////////////////////////////
	  //BGK
      //real omegaD     = three - sqrt(three);
      //real Lam         = -(c1o2-one/omegaD);
      //real nue_d       = Lam/three;
      //real ae          = diffusivity/nue_d - one; //zero;

      real f7_E,f7_W,f7_N,f7_S,f7_T,f7_B;
      //real /*feq7_ZERO,*/feq7_E,feq7_W,feq7_N,feq7_S,feq7_T,feq7_B;
      //real /*feqW7_ZERO,*/feqW7_E,feqW7_W,feqW7_N,feqW7_S,feqW7_T,feqW7_B;
      real TempD = temp[k];


      //f7_ZERO =  (D7.f[0])[kzero];
      f7_W    =  (D7.f[1])[ke   ];
      f7_E    =  (D7.f[2])[kw   ];
      f7_S    =  (D7.f[3])[kn   ];
      f7_N    =  (D7.f[4])[ks   ];
      f7_B    =  (D7.f[5])[kt   ];
      f7_T    =  (D7.f[6])[kb   ];

      //real ConcD = f7_ZERO + f7_E + f7_W + f7_N + f7_S + f7_T + f7_B;

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
      real cs2     = c1o4;
      //////////////////////////////////////////////////////////////////////////
      //pointertausch
      if (isEvenTimestep==false)
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
      //real /*feq,*/ q;
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
      //real /*feq,*/ q;
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
      real /*feq,*/ q;
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[2])[kw]= -f7_W + cs2 * TempD;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[1])[ke]= -f7_E + cs2 * TempD;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[4])[ks]= -f7_S + cs2 * TempD;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[3])[kn]= -f7_N + cs2 * TempD;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[6])[kb]= -f7_B + cs2 * TempD;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[5])[kt]= -f7_T + cs2 * TempD;
      }

   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QNoSlipADincomp27(
											 real* DD, 
											 real* DD27, 
											 real* temp,
											 real diffusivity,
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
   if (isEvenTimestep==true)
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      q_dirNE  = &QQ[dirNE  * numberOfBCnodes];
      q_dirSW  = &QQ[dirSW  * numberOfBCnodes];
      q_dirSE  = &QQ[dirSE  * numberOfBCnodes];
      q_dirNW  = &QQ[dirNW  * numberOfBCnodes];
      q_dirTE  = &QQ[dirTE  * numberOfBCnodes];
      q_dirBW  = &QQ[dirBW  * numberOfBCnodes];
      q_dirBE  = &QQ[dirBE  * numberOfBCnodes];
      q_dirTW  = &QQ[dirTW  * numberOfBCnodes];
      q_dirTN  = &QQ[dirTN  * numberOfBCnodes];
      q_dirBS  = &QQ[dirBS  * numberOfBCnodes];
      q_dirBN  = &QQ[dirBN  * numberOfBCnodes];
      q_dirTS  = &QQ[dirTS  * numberOfBCnodes];
      q_dirTNE = &QQ[dirTNE * numberOfBCnodes];
      q_dirTSW = &QQ[dirTSW * numberOfBCnodes];
      q_dirTSE = &QQ[dirTSE * numberOfBCnodes];
      q_dirTNW = &QQ[dirTNW * numberOfBCnodes];
      q_dirBNE = &QQ[dirBNE * numberOfBCnodes];
      q_dirBSW = &QQ[dirBSW * numberOfBCnodes];
      q_dirBSE = &QQ[dirBSE * numberOfBCnodes];
      q_dirBNW = &QQ[dirBNW * numberOfBCnodes];
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
      //real f_ZERO = (D.f[dirREST])[kzero];
      real f_BSW  = (D.f[dirTNE ])[ktne ];
      real f_BNE  = (D.f[dirTSW ])[ktsw ];
      real f_BNW  = (D.f[dirTSE ])[ktse ];
      real f_BSE  = (D.f[dirTNW ])[ktnw ];
      real f_TSW  = (D.f[dirBNE ])[kbne ];
      real f_TNE  = (D.f[dirBSW ])[kbsw ];
      real f_TNW  = (D.f[dirBSE ])[kbse ];
      real f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1 =  ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      real vx2 =  ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      real vx3 =  ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      real f27_W    = (D27.f[dirE   ])[ke   ];
      real f27_E    = (D27.f[dirW   ])[kw   ];
      real f27_S    = (D27.f[dirN   ])[kn   ];
      real f27_N    = (D27.f[dirS   ])[ks   ];
      real f27_B    = (D27.f[dirT   ])[kt   ];
      real f27_T    = (D27.f[dirB   ])[kb   ];
      real f27_SW   = (D27.f[dirNE  ])[kne  ];
      real f27_NE   = (D27.f[dirSW  ])[ksw  ];
      real f27_NW   = (D27.f[dirSE  ])[kse  ];
      real f27_SE   = (D27.f[dirNW  ])[knw  ];
      real f27_BW   = (D27.f[dirTE  ])[kte  ];
      real f27_TE   = (D27.f[dirBW  ])[kbw  ];
      real f27_TW   = (D27.f[dirBE  ])[kbe  ];
      real f27_BE   = (D27.f[dirTW  ])[ktw  ];
      real f27_BS   = (D27.f[dirTN  ])[ktn  ];
      real f27_TN   = (D27.f[dirBS  ])[kbs  ];
      real f27_TS   = (D27.f[dirBN  ])[kbn  ];
      real f27_BN   = (D27.f[dirTS  ])[kts  ];
      real f27_ZERO = (D27.f[dirREST])[kzero];
      real f27_BSW  = (D27.f[dirTNE ])[ktne ];
      real f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      real f27_BNW  = (D27.f[dirTSE ])[ktse ];
      real f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      real f27_TSW  = (D27.f[dirBNE ])[kbne ];
      real f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      real f27_TNW  = (D27.f[dirBSE ])[kbse ];
      real f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      real ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
         f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
         f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      //real feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      real feq27_E    =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      real feq27_W    =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      real feq27_N    =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      real feq27_S    =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      real feq27_T    =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      real feq27_B    =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      real feq27_NE   =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      real feq27_SW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      real feq27_SE   =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      real feq27_NW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      real feq27_TE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      real feq27_BW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      real feq27_BE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      real feq27_TW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      real feq27_TN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      real feq27_BS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      real feq27_BN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      real feq27_TS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      real feq27_TNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      real feq27_BSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      real feq27_BNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      real feq27_TSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      real feq27_TSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      real feq27_BNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      real feq27_BSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      real feq27_TNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real TempD = temp[k];

      //real feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      real feqW27_E    =   c2o27* TempD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      real feqW27_W    =   c2o27* TempD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      real feqW27_N    =   c2o27* TempD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      real feqW27_S    =   c2o27* TempD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      real feqW27_T    =   c2o27* TempD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      real feqW27_B    =   c2o27* TempD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      real feqW27_NE   =   c1o54* TempD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      real feqW27_SW   =   c1o54* TempD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      real feqW27_SE   =   c1o54* TempD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      real feqW27_NW   =   c1o54* TempD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      real feqW27_TE   =   c1o54* TempD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      real feqW27_BW   =   c1o54* TempD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      real feqW27_BE   =   c1o54* TempD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      real feqW27_TW   =   c1o54* TempD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      real feqW27_TN   =   c1o54* TempD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      real feqW27_BS   =   c1o54* TempD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      real feqW27_BN   =   c1o54* TempD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      real feqW27_TS   =   c1o54* TempD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      real feqW27_TNE  =   c1o216*TempD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      real feqW27_BSW  =   c1o216*TempD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      real feqW27_BNE  =   c1o216*TempD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      real feqW27_TSW  =   c1o216*TempD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      real feqW27_TSE  =   c1o216*TempD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      real feqW27_BNW  =   c1o216*TempD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      real feqW27_BSE  =   c1o216*TempD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      real feqW27_TNW  =   c1o216*TempD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real omegaD     = c3o1 - sqrt(c3o1);
      //real Lam        = -(c1o2-one/omegaD);
      //real nue_d      = Lam/three;
      //real ae         = zero;
      //real ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      //(D.f[dirREST])[k]=0.1f;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real q;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirW  ])[kw  ]=(c2o1*feqW27_W  -(f27_E  *(q*omegaD-c1o1)-omegaD*feq27_E  *(q-c1o1))/(omegaD-c1o1)+f27_W  *q)/(q+c1o1);
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirE  ])[ke  ]=(c2o1*feqW27_E  -(f27_W  *(q*omegaD-c1o1)-omegaD*feq27_W  *(q-c1o1))/(omegaD-c1o1)+f27_E  *q)/(q+c1o1);
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirS  ])[ks  ]=(c2o1*feqW27_S  -(f27_N  *(q*omegaD-c1o1)-omegaD*feq27_N  *(q-c1o1))/(omegaD-c1o1)+f27_S  *q)/(q+c1o1);
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirN  ])[kn  ]=(c2o1*feqW27_N  -(f27_S  *(q*omegaD-c1o1)-omegaD*feq27_S  *(q-c1o1))/(omegaD-c1o1)+f27_N  *q)/(q+c1o1);
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirB  ])[kb  ]=(c2o1*feqW27_B  -(f27_T  *(q*omegaD-c1o1)-omegaD*feq27_T  *(q-c1o1))/(omegaD-c1o1)+f27_B  *q)/(q+c1o1);
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirT  ])[kt  ]=(c2o1*feqW27_T  -(f27_B  *(q*omegaD-c1o1)-omegaD*feq27_B  *(q-c1o1))/(omegaD-c1o1)+f27_T  *q)/(q+c1o1);
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirSW ])[ksw ]=(c2o1*feqW27_SW -(f27_NE *(q*omegaD-c1o1)-omegaD*feq27_NE *(q-c1o1))/(omegaD-c1o1)+f27_SW *q)/(q+c1o1);
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirNE ])[kne ]=(c2o1*feqW27_NE -(f27_SW *(q*omegaD-c1o1)-omegaD*feq27_SW *(q-c1o1))/(omegaD-c1o1)+f27_NE *q)/(q+c1o1);
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirNW ])[knw ]=(c2o1*feqW27_NW -(f27_SE *(q*omegaD-c1o1)-omegaD*feq27_SE *(q-c1o1))/(omegaD-c1o1)+f27_NW *q)/(q+c1o1);
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirSE ])[kse ]=(c2o1*feqW27_SE -(f27_NW *(q*omegaD-c1o1)-omegaD*feq27_NW *(q-c1o1))/(omegaD-c1o1)+f27_SE *q)/(q+c1o1);
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBW ])[kbw ]=(c2o1*feqW27_BW -(f27_TE *(q*omegaD-c1o1)-omegaD*feq27_TE *(q-c1o1))/(omegaD-c1o1)+f27_BW *q)/(q+c1o1);
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTE ])[kte ]=(c2o1*feqW27_TE -(f27_BW *(q*omegaD-c1o1)-omegaD*feq27_BW *(q-c1o1))/(omegaD-c1o1)+f27_TE *q)/(q+c1o1);
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTW ])[ktw ]=(c2o1*feqW27_TW -(f27_BE *(q*omegaD-c1o1)-omegaD*feq27_BE *(q-c1o1))/(omegaD-c1o1)+f27_TW *q)/(q+c1o1);
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBE ])[kbe ]=(c2o1*feqW27_BE -(f27_TW *(q*omegaD-c1o1)-omegaD*feq27_TW *(q-c1o1))/(omegaD-c1o1)+f27_BE *q)/(q+c1o1);
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBS ])[kbs ]=(c2o1*feqW27_BS -(f27_TN *(q*omegaD-c1o1)-omegaD*feq27_TN *(q-c1o1))/(omegaD-c1o1)+f27_BS *q)/(q+c1o1);
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTN ])[ktn ]=(c2o1*feqW27_TN -(f27_BS *(q*omegaD-c1o1)-omegaD*feq27_BS *(q-c1o1))/(omegaD-c1o1)+f27_TN *q)/(q+c1o1);
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTS ])[kts ]=(c2o1*feqW27_TS -(f27_BN *(q*omegaD-c1o1)-omegaD*feq27_BN *(q-c1o1))/(omegaD-c1o1)+f27_TS *q)/(q+c1o1);
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBN ])[kbn ]=(c2o1*feqW27_BN -(f27_TS *(q*omegaD-c1o1)-omegaD*feq27_TS *(q-c1o1))/(omegaD-c1o1)+f27_BN *q)/(q+c1o1);
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSW])[kbsw]=(c2o1*feqW27_BSW-(f27_TNE*(q*omegaD-c1o1)-omegaD*feq27_TNE*(q-c1o1))/(omegaD-c1o1)+f27_BSW*q)/(q+c1o1);
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNE])[ktne]=(c2o1*feqW27_TNE-(f27_BSW*(q*omegaD-c1o1)-omegaD*feq27_BSW*(q-c1o1))/(omegaD-c1o1)+f27_TNE*q)/(q+c1o1);
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSW])[ktsw]=(c2o1*feqW27_TSW-(f27_BNE*(q*omegaD-c1o1)-omegaD*feq27_BNE*(q-c1o1))/(omegaD-c1o1)+f27_TSW*q)/(q+c1o1);
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNE])[kbne]=(c2o1*feqW27_BNE-(f27_TSW*(q*omegaD-c1o1)-omegaD*feq27_TSW*(q-c1o1))/(omegaD-c1o1)+f27_BNE*q)/(q+c1o1);
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNW])[kbnw]=(c2o1*feqW27_BNW-(f27_TSE*(q*omegaD-c1o1)-omegaD*feq27_TSE*(q-c1o1))/(omegaD-c1o1)+f27_BNW*q)/(q+c1o1);
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSE])[ktse]=(c2o1*feqW27_TSE-(f27_BNW*(q*omegaD-c1o1)-omegaD*feq27_BNW*(q-c1o1))/(omegaD-c1o1)+f27_TSE*q)/(q+c1o1);
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNW])[ktnw]=(c2o1*feqW27_TNW-(f27_BSE*(q*omegaD-c1o1)-omegaD*feq27_BSE*(q-c1o1))/(omegaD-c1o1)+f27_TNW*q)/(q+c1o1);
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSE])[kbse]=(c2o1*feqW27_BSE-(f27_TNW*(q*omegaD-c1o1)-omegaD*feq27_TNW*(q-c1o1))/(omegaD-c1o1)+f27_BSE*q)/(q+c1o1);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADVeloIncomp7(
											real* DD, 
											real* DD7, 
											real* temp,
											real* velo,
											real diffusivity,
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
   //Distributions27 D;
   //if (isEvenTimestep==true)
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
   //   D.f[dirREST] = &DD[dirREST*size_Mat];
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
   //   D.f[dirREST] = &DD[dirREST*size_Mat];
   //   D.f[dirTNE ] = &DD[dirBSW *size_Mat];
   //   D.f[dirTSW ] = &DD[dirBNE *size_Mat];
   //   D.f[dirTSE ] = &DD[dirBNW *size_Mat];
   //   D.f[dirTNW ] = &DD[dirBSE *size_Mat];
   //   D.f[dirBNE ] = &DD[dirTSW *size_Mat];
   //   D.f[dirBSW ] = &DD[dirTNE *size_Mat];
   //   D.f[dirBSE ] = &DD[dirTNW *size_Mat];
   //   D.f[dirBNW ] = &DD[dirTSE *size_Mat];
   //}

   Distributions7 D7;
   if (isEvenTimestep==true)
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

   if(k<numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB; 

      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      //////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      //unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      //unsigned int ksw  = neighborY[kw];
      //unsigned int kne  = KQK;
      //unsigned int kse  = ks;
      //unsigned int knw  = kw;
      //unsigned int kbw  = neighborZ[kw];
      //unsigned int kte  = KQK;
      //unsigned int kbe  = kb;
      //unsigned int ktw  = kw;
      //unsigned int kbs  = neighborZ[ks];
      //unsigned int ktn  = KQK;
      //unsigned int kbn  = kb;
      //unsigned int kts  = ks;
      //unsigned int ktse = ks;
      //unsigned int kbnw = kbw;
      //unsigned int ktnw = kw;
      //unsigned int kbse = kbs;
      //unsigned int ktsw = ksw;
      //unsigned int kbne = kb;
      //unsigned int ktne = KQK;
      //unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      //real f_W    = (D.f[dirE   ])[ke   ];
      //real f_E    = (D.f[dirW   ])[kw   ];
      //real f_S    = (D.f[dirN   ])[kn   ];
      //real f_N    = (D.f[dirS   ])[ks   ];
      //real f_B    = (D.f[dirT   ])[kt   ];
      //real f_T    = (D.f[dirB   ])[kb   ];
      //real f_SW   = (D.f[dirNE  ])[kne  ];
      //real f_NE   = (D.f[dirSW  ])[ksw  ];
      //real f_NW   = (D.f[dirSE  ])[kse  ];
      //real f_SE   = (D.f[dirNW  ])[knw  ];
      //real f_BW   = (D.f[dirTE  ])[kte  ];
      //real f_TE   = (D.f[dirBW  ])[kbw  ];
      //real f_TW   = (D.f[dirBE  ])[kbe  ];
      //real f_BE   = (D.f[dirTW  ])[ktw  ];
      //real f_BS   = (D.f[dirTN  ])[ktn  ];
      //real f_TN   = (D.f[dirBS  ])[kbs  ];
      //real f_TS   = (D.f[dirBN  ])[kbn  ];
      //real f_BN   = (D.f[dirTS  ])[kts  ];
      //real f_BSW  = (D.f[dirTNE ])[ktne ];
      //real f_BNE  = (D.f[dirTSW ])[ktsw ];
      //real f_BNW  = (D.f[dirTSE ])[ktse ];
      //real f_BSE  = (D.f[dirTNW ])[ktnw ];
      //real f_TSW  = (D.f[dirBNE ])[kbne ];
      //real f_TNE  = (D.f[dirBSW ])[kbsw ];
      //real f_TNW  = (D.f[dirBSE ])[kbse ];
      //real f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      //real vx1_Inflow   = c0o1;
      //real vx2_Inflow   = velo[k];
      //real vx3_Inflow   = c0o1;
      //real ux_sq_Inflow = vx1_Inflow * vx1_Inflow;
      //real uy_sq_Inflow = vx2_Inflow * vx2_Inflow;
      //real uz_sq_Inflow = vx3_Inflow * vx3_Inflow;
      ////////////////////////////////////////////////////////////////////////////////
      //real vx1 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      //real vx2 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      //real vx3 = ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
		 ////d�rrrrrty !!!!!!!!!!!!!
   //      real vx1     =  ten * ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
   //      real vx2     =  ten * ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
   //      real vx3     =  ten * ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
	  //real cu_sq       =1.5f*(vx1*vx1+vx2*vx2+vx3*vx3);
      //real ux_sq       = vx1 * vx1;
      //real uy_sq       = vx2 * vx2;
      //real uz_sq       = vx3 * vx3;

      real f7_E,f7_W,f7_N,f7_S,f7_T,f7_B;
      //real /*feq7_ZERO,*/feq7_E,feq7_W,feq7_N,feq7_S,feq7_T,feq7_B;
      //real /*feqW7_ZERO,*/feqW7_E,feqW7_W,feqW7_N,feqW7_S,feqW7_T,feqW7_B;
      real TempD = temp[k];

      //f7_ZERO =  (D7.f[0])[kzero];
      f7_W    =  (D7.f[1])[ke   ];
      f7_E    =  (D7.f[2])[kw   ];
      f7_S    =  (D7.f[3])[kn   ];
      f7_N    =  (D7.f[4])[ks   ];
      f7_B    =  (D7.f[5])[kt   ];
      f7_T    =  (D7.f[6])[kb   ];

      //real ConcD = f7_ZERO + f7_E + f7_W + f7_N + f7_S + f7_T + f7_B;

      ////////////////////////////////////////////////////////////////////////////////
	  //BGK
      //real omegaD     = three - sqrt(three);
      //real Lam         = -(c1o2-one/omegaD);
      //real nue_d       = Lam/three;
      ////real ae          = zero;
      //real ae          = diffusivity/nue_d - one;

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
         real cs2         = c1o4;
   //      real Lam         = diffusivity/(one)/cs2;
   //      real omegaD      = - one / (Lam + c1o2);
   //      real ae          = zero;
   //      ////////////////////////////////////////////////////////////////////////////////
		 //real Mom000 = f7_ZERO + f7_W + f7_E + f7_N + f7_S + f7_T + f7_B; //1
   //      real Mom100 = f7_E - f7_W;
   //      real Mom010 = f7_N - f7_S;
   //      real Mom001 = f7_T - f7_B;
   //      real Mom222 = six*f7_ZERO - f7_W - f7_E - f7_N - f7_S - f7_T - f7_B;
   //      real Mom200 = two*f7_W + two*f7_E - f7_N - f7_S - f7_T - f7_B;
   //      real Mom022 = f7_N + f7_S - f7_T - f7_B;

   //      real Meq000 = ConcD;
   //      real Meq100 = ConcD*vx1;
   //      real Meq010 = ConcD*vx2;
   //      real Meq001 = ConcD*vx3;
   //      real Meq222 = c3o4*ConcD;
   //      real Meq200 = zero;
   //      real Meq022 = zero;

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
      if (isEvenTimestep==false)
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
      //real /*feq,*/ q;
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
      //real /*feq,*/ q;
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
      real /*feq,*/ q;
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[2])[kw]= -f7_W + cs2 * TempD;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[1])[ke]= -f7_E + cs2 * TempD;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[4])[ks]= -f7_S + cs2 * TempD;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[3])[kn]= -f7_N + cs2 * TempD;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[6])[kb]= -f7_B + cs2 * TempD;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[5])[kt]= -f7_T + cs2 * TempD;
      }

   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADVeloIncomp27(
											real* DD, 
											real* DD27, 
											real* temp,
											real* velo,
											real diffusivity,
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
   if (isEvenTimestep==true)
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      q_dirNE  = &QQ[dirNE  * numberOfBCnodes];
      q_dirSW  = &QQ[dirSW  * numberOfBCnodes];
      q_dirSE  = &QQ[dirSE  * numberOfBCnodes];
      q_dirNW  = &QQ[dirNW  * numberOfBCnodes];
      q_dirTE  = &QQ[dirTE  * numberOfBCnodes];
      q_dirBW  = &QQ[dirBW  * numberOfBCnodes];
      q_dirBE  = &QQ[dirBE  * numberOfBCnodes];
      q_dirTW  = &QQ[dirTW  * numberOfBCnodes];
      q_dirTN  = &QQ[dirTN  * numberOfBCnodes];
      q_dirBS  = &QQ[dirBS  * numberOfBCnodes];
      q_dirBN  = &QQ[dirBN  * numberOfBCnodes];
      q_dirTS  = &QQ[dirTS  * numberOfBCnodes];
      q_dirTNE = &QQ[dirTNE * numberOfBCnodes];
      q_dirTSW = &QQ[dirTSW * numberOfBCnodes];
      q_dirTSE = &QQ[dirTSE * numberOfBCnodes];
      q_dirTNW = &QQ[dirTNW * numberOfBCnodes];
      q_dirBNE = &QQ[dirBNE * numberOfBCnodes];
      q_dirBSW = &QQ[dirBSW * numberOfBCnodes];
      q_dirBSE = &QQ[dirBSE * numberOfBCnodes];
      q_dirBNW = &QQ[dirBNW * numberOfBCnodes];
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
      //real f_ZERO = (D.f[dirREST])[kzero];
      real f_BSW  = (D.f[dirTNE ])[ktne ];
      real f_BNE  = (D.f[dirTSW ])[ktsw ];
      real f_BNW  = (D.f[dirTSE ])[ktse ];
      real f_BSE  = (D.f[dirTNW ])[ktnw ];
      real f_TSW  = (D.f[dirBNE ])[kbne ];
      real f_TNE  = (D.f[dirBSW ])[kbsw ];
      real f_TNW  = (D.f[dirBSE ])[kbse ];
      real f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      real vx2 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      real vx3 = ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      //real f27_W    = (D27.f[dirE   ])[ke   ];
      //real f27_E    = (D27.f[dirW   ])[kw   ];
      //real f27_S    = (D27.f[dirN   ])[kn   ];
      //real f27_N    = (D27.f[dirS   ])[ks   ];
      //real f27_B    = (D27.f[dirT   ])[kt   ];
      //real f27_T    = (D27.f[dirB   ])[kb   ];
      //real f27_SW   = (D27.f[dirNE  ])[kne  ];
      //real f27_NE   = (D27.f[dirSW  ])[ksw  ];
      //real f27_NW   = (D27.f[dirSE  ])[kse  ];
      //real f27_SE   = (D27.f[dirNW  ])[knw  ];
      //real f27_BW   = (D27.f[dirTE  ])[kte  ];
      //real f27_TE   = (D27.f[dirBW  ])[kbw  ];
      //real f27_TW   = (D27.f[dirBE  ])[kbe  ];
      //real f27_BE   = (D27.f[dirTW  ])[ktw  ];
      //real f27_BS   = (D27.f[dirTN  ])[ktn  ];
      //real f27_TN   = (D27.f[dirBS  ])[kbs  ];
      //real f27_TS   = (D27.f[dirBN  ])[kbn  ];
      //real f27_BN   = (D27.f[dirTS  ])[kts  ];
      //real f27_ZERO = (D27.f[dirREST])[kzero];
      //real f27_BSW  = (D27.f[dirTNE ])[ktne ];
      //real f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      //real f27_BNW  = (D27.f[dirTSE ])[ktse ];
      //real f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      //real f27_TSW  = (D27.f[dirBNE ])[kbne ];
      //real f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      //real f27_TNW  = (D27.f[dirBSE ])[kbse ];
      //real f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      //real ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
      //                  f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
      //                  f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      //real feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      //real feq27_E    =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      //real feq27_W    =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      //real feq27_N    =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      //real feq27_S    =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      //real feq27_T    =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      //real feq27_B    =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      //real feq27_NE   =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      //real feq27_SW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      //real feq27_SE   =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      //real feq27_NW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      //real feq27_TE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      //real feq27_BW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      //real feq27_BE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      //real feq27_TW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      //real feq27_TN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      //real feq27_BS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      //real feq27_BN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      //real feq27_TS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      //real feq27_TNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      //real feq27_BSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      //real feq27_BNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      //real feq27_TSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      //real feq27_TSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      //real feq27_BNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      //real feq27_BSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      //real feq27_TNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real TempD = temp[k];
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // velocity inflow
      vx1   = velo[k];//zero;
      vx2   = c0o1;//velo[k];
      vx3   = c0o1;

      //real feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      real feqW27_E    =   c2o27* TempD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      real feqW27_W    =   c2o27* TempD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      real feqW27_N    =   c2o27* TempD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      real feqW27_S    =   c2o27* TempD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      real feqW27_T    =   c2o27* TempD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      real feqW27_B    =   c2o27* TempD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      real feqW27_NE   =   c1o54* TempD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      real feqW27_SW   =   c1o54* TempD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      real feqW27_SE   =   c1o54* TempD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      real feqW27_NW   =   c1o54* TempD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      real feqW27_TE   =   c1o54* TempD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      real feqW27_BW   =   c1o54* TempD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      real feqW27_BE   =   c1o54* TempD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      real feqW27_TW   =   c1o54* TempD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      real feqW27_TN   =   c1o54* TempD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      real feqW27_BS   =   c1o54* TempD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      real feqW27_BN   =   c1o54* TempD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      real feqW27_TS   =   c1o54* TempD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      real feqW27_TNE  =   c1o216*TempD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      real feqW27_BSW  =   c1o216*TempD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      real feqW27_BNE  =   c1o216*TempD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      real feqW27_TSW  =   c1o216*TempD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      real feqW27_TSE  =   c1o216*TempD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      real feqW27_BNW  =   c1o216*TempD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      real feqW27_BSE  =   c1o216*TempD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      real feqW27_TNW  =   c1o216*TempD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real omegaD     = c3o1 - sqrt(c3o1);
      //real Lam        = -(c1o2-one/omegaD);
      //real nue_d      = Lam/three;
      //real ae         = zero;
      //real ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      //(D.f[dirREST])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real q;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirW  ])[kw  ]= -feqW27_W  + c2o1 * c2o27  * TempD;
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirE  ])[ke  ]= -feqW27_E  + c2o1 * c2o27  * TempD;
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirS  ])[ks  ]= -feqW27_S  + c2o1 * c2o27  * TempD;
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirN  ])[kn  ]= -feqW27_N  + c2o1 * c2o27  * TempD;
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirB  ])[kb  ]= -feqW27_B  + c2o1 * c2o27  * TempD;
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirT  ])[kt  ]= -feqW27_T  + c2o1 * c2o27  * TempD;
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirSW ])[ksw ]= -feqW27_SW + c2o1 * c1o54  * TempD;
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirNE ])[kne ]= -feqW27_NE + c2o1 * c1o54  * TempD;
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirNW ])[knw ]= -feqW27_NW + c2o1 * c1o54  * TempD;
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirSE ])[kse ]= -feqW27_SE + c2o1 * c1o54  * TempD;
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBW ])[kbw ]= -feqW27_BW + c2o1 * c1o54  * TempD;
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTE ])[kte ]= -feqW27_TE + c2o1 * c1o54  * TempD;
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTW ])[ktw ]= -feqW27_TW + c2o1 * c1o54  * TempD;
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBE ])[kbe ]= -feqW27_BE + c2o1 * c1o54  * TempD;
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBS ])[kbs ]= -feqW27_BS + c2o1 * c1o54  * TempD;
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTN ])[ktn ]= -feqW27_TN + c2o1 * c1o54  * TempD;
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTS ])[kts ]= -feqW27_TS + c2o1 * c1o54  * TempD;
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBN ])[kbn ]= -feqW27_BN + c2o1 * c1o54  * TempD;
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSW])[kbsw]= -feqW27_BSW+ c2o1 * c1o216 * TempD;
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNE])[ktne]= -feqW27_TNE+ c2o1 * c1o216 * TempD;
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSW])[ktsw]= -feqW27_TSW+ c2o1 * c1o216 * TempD;
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNE])[kbne]= -feqW27_BNE+ c2o1 * c1o216 * TempD;
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNW])[kbnw]= -feqW27_BNW+ c2o1 * c1o216 * TempD;
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSE])[ktse]= -feqW27_TSE+ c2o1 * c1o216 * TempD;
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNW])[ktnw]= -feqW27_TNW+ c2o1 * c1o216 * TempD;
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSE])[kbse]= -feqW27_BSE+ c2o1 * c1o216 * TempD;
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
										   real* DD, 
										   real* DD7, 
										   real* temp,
										   real* velo,
										   real diffusivity,
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
  /* Distributions27 D;
   if (isEvenTimestep==true)
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
      D.f[dirTNE ] = &DD[dirBSW *size_Mat];
      D.f[dirTSW ] = &DD[dirBNE *size_Mat];
      D.f[dirTSE ] = &DD[dirBNW *size_Mat];
      D.f[dirTNW ] = &DD[dirBSE *size_Mat];
      D.f[dirBNE ] = &DD[dirTSW *size_Mat];
      D.f[dirBSW ] = &DD[dirTNE *size_Mat];
      D.f[dirBSE ] = &DD[dirTNW *size_Mat];
      D.f[dirBNW ] = &DD[dirTSE *size_Mat];
   }*/

   Distributions7 D7;
   if (isEvenTimestep==true)
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

   if(k<numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB; 

      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
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
      //unsigned int ksw  = neighborY[kw];
      //unsigned int kne  = KQK;
      //unsigned int kse  = ks;
      //unsigned int knw  = kw;
      //unsigned int kbw  = neighborZ[kw];
      //unsigned int kte  = KQK;
      //unsigned int kbe  = kb;
      //unsigned int ktw  = kw;
      //unsigned int kbs  = neighborZ[ks];
      //unsigned int ktn  = KQK;
      //unsigned int kbn  = kb;
      //unsigned int kts  = ks;
      //unsigned int ktse = ks;
      //unsigned int kbnw = kbw;
      //unsigned int ktnw = kw;
      //unsigned int kbse = kbs;
      //unsigned int ktsw = ksw;
      //unsigned int kbne = kb;
      //unsigned int ktne = KQK;
      //unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
    /*  real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
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
      f_TSE  = (D.f[dirBNW ])[kbnw ];*/
      ////////////////////////////////////////////////////////////////////////////////
      //real vx1 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      //real vx2 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      //real vx3 = ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
		 ////d�rrrrrty !!!!!!!!!!!!!
   //      real vx1     =  ten * ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
   //      real vx2     =  ten * ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
   //      real vx3     =  ten * ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));

	  //real cu_sq       =1.5*(vx1*vx1+vx2*vx2+vx3*vx3);
      //real ux_sq       = vx1 * vx1;
      //real uy_sq       = vx2 * vx2;
      //real uz_sq       = vx3 * vx3;
      //////////////////////////////////////////////////////////////////////////
	  //BGK
      //real omegaD     = three - sqrt(three);
      //real Lam         = -(c1o2-one/omegaD);
      //real nue_d       = Lam/three;
      ////real ae          = zero;
      //real ae          = diffusivity/nue_d - one;

      real f7_ZERO,f7_E,f7_W,f7_N,f7_S,f7_T,f7_B;
      //real /*feq7_ZERO,*/feq7_E,feq7_W,feq7_N,feq7_S,feq7_T,feq7_B;
      //real /*feqW7_ZERO,*/feqW7_E,feqW7_W,feqW7_N,feqW7_S,feqW7_T,feqW7_B;
      //real TempD = temp[k];


      f7_ZERO =  (D7.f[0])[kzero];
      f7_W    =  (D7.f[1])[ke   ];
      f7_E    =  (D7.f[2])[kw   ];
      f7_S    =  (D7.f[3])[kn   ];
      f7_N    =  (D7.f[4])[ks   ];
      f7_B    =  (D7.f[5])[kt   ];
      f7_T    =  (D7.f[6])[kb   ];

      real ConcD = f7_ZERO + f7_E + f7_W + f7_N + f7_S + f7_T + f7_B;

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
      real cs2         = c1o4;
      real Lam         = diffusivity/(c1o1)/cs2;
      //real omegaD      = - c1o1 / (Lam + c1o2);
      real nue_d       = Lam/c3o1;

      //////////////////////////////////////////////////////////////////////////
      //pointertausch
      if (isEvenTimestep==false)
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
      //real /*feq,*/ q;
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
      //real /*feq,*/ q;
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
      real /*feq,*/ q;
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[2])[kw]= f7_W + nue_d * ConcD;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[1])[ke]= f7_E + nue_d * ConcD;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[4])[ks]= f7_S + nue_d * ConcD;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[3])[kn]= f7_N + nue_d * ConcD;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[6])[kb]= f7_B + nue_d * ConcD;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D7.f[5])[kt]= f7_T + nue_d * ConcD;
      }

   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QADPressIncomp27(
											   real* DD,
											   real* DD27,
											   real* temp,
											   real* velo,
											   real diffusivity,
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
   if (isEvenTimestep==true)
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
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

   if(k < numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      q_dirNE  = &QQ[dirNE  * numberOfBCnodes];
      q_dirSW  = &QQ[dirSW  * numberOfBCnodes];
      q_dirSE  = &QQ[dirSE  * numberOfBCnodes];
      q_dirNW  = &QQ[dirNW  * numberOfBCnodes];
      q_dirTE  = &QQ[dirTE  * numberOfBCnodes];
      q_dirBW  = &QQ[dirBW  * numberOfBCnodes];
      q_dirBE  = &QQ[dirBE  * numberOfBCnodes];
      q_dirTW  = &QQ[dirTW  * numberOfBCnodes];
      q_dirTN  = &QQ[dirTN  * numberOfBCnodes];
      q_dirBS  = &QQ[dirBS  * numberOfBCnodes];
      q_dirBN  = &QQ[dirBN  * numberOfBCnodes];
      q_dirTS  = &QQ[dirTS  * numberOfBCnodes];
      q_dirTNE = &QQ[dirTNE * numberOfBCnodes];
      q_dirTSW = &QQ[dirTSW * numberOfBCnodes];
      q_dirTSE = &QQ[dirTSE * numberOfBCnodes];
      q_dirTNW = &QQ[dirTNW * numberOfBCnodes];
      q_dirBNE = &QQ[dirBNE * numberOfBCnodes];
      q_dirBSW = &QQ[dirBSW * numberOfBCnodes];
      q_dirBSE = &QQ[dirBSE * numberOfBCnodes];
      q_dirBNW = &QQ[dirBNW * numberOfBCnodes];
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
      //real f_ZERO = (D.f[dirREST])[kzero];
      real f_BSW  = (D.f[dirTNE ])[ktne ];
      real f_BNE  = (D.f[dirTSW ])[ktsw ];
      real f_BNW  = (D.f[dirTSE ])[ktse ];
      real f_BSE  = (D.f[dirTNW ])[ktnw ];
      real f_TSW  = (D.f[dirBNE ])[kbne ];
      real f_TNE  = (D.f[dirBSW ])[kbsw ];
      real f_TNW  = (D.f[dirBSE ])[kbse ];
      real f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1      = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      real vx2      = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      real vx3      = ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      //real f27_W    = (D27.f[dirE   ])[ke   ];
      //real f27_E    = (D27.f[dirW   ])[kw   ];
      //real f27_S    = (D27.f[dirN   ])[kn   ];
      //real f27_N    = (D27.f[dirS   ])[ks   ];
      //real f27_B    = (D27.f[dirT   ])[kt   ];
      //real f27_T    = (D27.f[dirB   ])[kb   ];
      //real f27_SW   = (D27.f[dirNE  ])[kne  ];
      //real f27_NE   = (D27.f[dirSW  ])[ksw  ];
      //real f27_NW   = (D27.f[dirSE  ])[kse  ];
      //real f27_SE   = (D27.f[dirNW  ])[knw  ];
      //real f27_BW   = (D27.f[dirTE  ])[kte  ];
      //real f27_TE   = (D27.f[dirBW  ])[kbw  ];
      //real f27_TW   = (D27.f[dirBE  ])[kbe  ];
      //real f27_BE   = (D27.f[dirTW  ])[ktw  ];
      //real f27_BS   = (D27.f[dirTN  ])[ktn  ];
      //real f27_TN   = (D27.f[dirBS  ])[kbs  ];
      //real f27_TS   = (D27.f[dirBN  ])[kbn  ];
      //real f27_BN   = (D27.f[dirTS  ])[kts  ];
      //real f27_ZERO = (D27.f[dirREST])[kzero];
      //real f27_BSW  = (D27.f[dirTNE ])[ktne ];
      //real f27_BNE  = (D27.f[dirTSW ])[ktsw ];
      //real f27_BNW  = (D27.f[dirTSE ])[ktse ];
      //real f27_BSE  = (D27.f[dirTNW ])[ktnw ];
      //real f27_TSW  = (D27.f[dirBNE ])[kbne ];
      //real f27_TNE  = (D27.f[dirBSW ])[kbsw ];
      //real f27_TNW  = (D27.f[dirBSE ])[kbse ];
      //real f27_TSE  = (D27.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      //real ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
                        //f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
                        //f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      //real feq27_ZERO =   c8over27* ConcD*(one-cu_sq);
      /*real feq27_E    =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      real feq27_W    =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      real feq27_N    =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      real feq27_S    =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      real feq27_T    =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      real feq27_B    =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      real feq27_NE   =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      real feq27_SW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      real feq27_SE   =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      real feq27_NW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      real feq27_TE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      real feq27_BW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      real feq27_BE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      real feq27_TW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      real feq27_TN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      real feq27_BS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      real feq27_BN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      real feq27_TS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      real feq27_TNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      real feq27_BSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      real feq27_BNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      real feq27_TSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      real feq27_TSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      real feq27_BNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      real feq27_BSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      real feq27_TNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);*/
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real TempD = temp[k];
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // velocity inflow
      //vx1   = zero;
      //vx2   = zero;
      //vx3   = velo[k];

      //real feqW27_ZERO =   c8over27* TempD*(one-cu_sq);
      real feqW27_E    =  c2o27* TempD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq); //feq27_E  ;// 
      real feqW27_W    =  c2o27* TempD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq); //feq27_W  ;// 
      real feqW27_N    =  c2o27* TempD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq); //feq27_N  ;// 
      real feqW27_S    =  c2o27* TempD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq); //feq27_S  ;// 
      real feqW27_T    =  c2o27* TempD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq); //feq27_T  ;// 
      real feqW27_B    =  c2o27* TempD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq); //feq27_B  ;// 
      real feqW27_NE   =  c1o54* TempD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); //feq27_NE ;// 
      real feqW27_SW   =  c1o54* TempD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); //feq27_SW ;// 
      real feqW27_SE   =  c1o54* TempD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); //feq27_SE ;// 
      real feqW27_NW   =  c1o54* TempD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); //feq27_NW ;// 
      real feqW27_TE   =  c1o54* TempD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); //feq27_TE ;// 
      real feqW27_BW   =  c1o54* TempD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); //feq27_BW ;// 
      real feqW27_BE   =  c1o54* TempD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); //feq27_BE ;// 
      real feqW27_TW   =  c1o54* TempD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); //feq27_TW ;// 
      real feqW27_TN   =  c1o54* TempD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); //feq27_TN ;// 
      real feqW27_BS   =  c1o54* TempD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); //feq27_BS ;// 
      real feqW27_BN   =  c1o54* TempD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); //feq27_BN ;// 
      real feqW27_TS   =  c1o54* TempD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); //feq27_TS ;// 
      real feqW27_TNE  =  c1o216*TempD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); //feq27_TNE;// 
      real feqW27_BSW  =  c1o216*TempD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); //feq27_BSW;// 
      real feqW27_BNE  =  c1o216*TempD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); //feq27_BNE;// 
      real feqW27_TSW  =  c1o216*TempD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); //feq27_TSW;// 
      real feqW27_TSE  =  c1o216*TempD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); //feq27_TSE;// 
      real feqW27_BNW  =  c1o216*TempD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); //feq27_BNW;// 
      real feqW27_BSE  =  c1o216*TempD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); //feq27_BSE;// 
      real feqW27_TNW  =  c1o216*TempD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); //feq27_TNW;// 
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real omegaD     = c3o1 - sqrt(c3o1);
      //real Lam        = -(c1o2-one/omegaD);
      //real nue_d      = Lam/three;
      //real ae         = zero;
      //real ae         = diffusivity/nue_d - one;


      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
         D27.f[dirREST] = &DD27[dirREST*size_Mat];
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
      //(D.f[dirREST])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real q;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirW  ])[kw  ]= -feqW27_W  + c2o1 * c2o27  * TempD;
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirE  ])[ke  ]= -feqW27_E  + c2o1 * c2o27  * TempD;
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirS  ])[ks  ]= -feqW27_S  + c2o1 * c2o27  * TempD;
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirN  ])[kn  ]= -feqW27_N  + c2o1 * c2o27  * TempD;
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirB  ])[kb  ]= -feqW27_B  + c2o1 * c2o27  * TempD;
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dirT  ])[kt  ]= -feqW27_T  + c2o1 * c2o27  * TempD;
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirSW ])[ksw ]= -feqW27_SW + c2o1 * c1o54  * TempD;
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirNE ])[kne ]= -feqW27_NE + c2o1 * c1o54  * TempD;
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirNW ])[knw ]= -feqW27_NW + c2o1 * c1o54  * TempD;
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirSE ])[kse ]= -feqW27_SE + c2o1 * c1o54  * TempD;
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBW ])[kbw ]= -feqW27_BW + c2o1 * c1o54  * TempD;
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTE ])[kte ]= -feqW27_TE + c2o1 * c1o54  * TempD;
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTW ])[ktw ]= -feqW27_TW + c2o1 * c1o54  * TempD;
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBE ])[kbe ]= -feqW27_BE + c2o1 * c1o54  * TempD;
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBS ])[kbs ]= -feqW27_BS + c2o1 * c1o54  * TempD;
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTN ])[ktn ]= -feqW27_TN + c2o1 * c1o54  * TempD;
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirTS ])[kts ]= -feqW27_TS + c2o1 * c1o54  * TempD;
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dirBN ])[kbn ]= -feqW27_BN + c2o1 * c1o54  * TempD;
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSW])[kbsw]= -feqW27_BSW+ c2o1 * c1o216 * TempD;
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNE])[ktne]= -feqW27_TNE+ c2o1 * c1o216 * TempD;
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSW])[ktsw]= -feqW27_TSW+ c2o1 * c1o216 * TempD;
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNE])[kbne]= -feqW27_BNE+ c2o1 * c1o216 * TempD;
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBNW])[kbnw]= -feqW27_BNW+ c2o1 * c1o216 * TempD;
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTSE])[ktse]= -feqW27_TSE+ c2o1 * c1o216 * TempD;
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirTNW])[ktnw]= -feqW27_TNW+ c2o1 * c1o216 * TempD;
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dirBSE])[kbse]= -feqW27_BSE+ c2o1 * c1o216 * TempD;
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













////////////////////////////////////////////////////////////////////////////////
inline __device__ real calcDistributionBC_AD_interpol(real q, real weight, real v, real v_sq, real f, real finf, real omegaDiffusivity, real jTangential, real concentration) {
    real feq = weight * concentration * (c1o1 + c3o1 * v + c9o2 * v * v * concentration - v_sq * concentration);
    return (c1o1 - q) / (c1o1 + q) * ((f - feq * omegaDiffusivity) / (c1o1 - omegaDiffusivity)) + (q * (f + finf) - c6o1 * weight * (jTangential)) / (c1o1 + q);
}
////////////////////////////////////////////////////////////////////////////////
inline __device__ real calcDistributionBC_AD(real q, real weight, real v, real v_sq, real f, real finf, real omegaDiffusivity, real jTangential, real concentration) {
    return f - c6o1 * weight * jTangential;
}


// has to be excecuted before Fluid BCs
//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void AD_SlipVelDeviceComp(
    real *normalX,
    real *normalY,
    real *normalZ,
    real *distributions,
    real *distributionsAD,
    int *QindexArray,
    real *Qarrays,
    uint numberOfBCnodes,
    real omegaDiffusivity,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    uint size_Mat,
    bool isEvenTimestep)
{
    Distributions27 D;
    if (isEvenTimestep)
    {
        D.f[dirE   ] = &distributions[dirE    * size_Mat];
        D.f[dirW   ] = &distributions[dirW    * size_Mat];
        D.f[dirN   ] = &distributions[dirN    * size_Mat];
        D.f[dirS   ] = &distributions[dirS    * size_Mat];
        D.f[dirT   ] = &distributions[dirT    * size_Mat];
        D.f[dirB   ] = &distributions[dirB    * size_Mat];
        D.f[dirNE  ] = &distributions[dirNE   * size_Mat];
        D.f[dirSW  ] = &distributions[dirSW   * size_Mat];
        D.f[dirSE  ] = &distributions[dirSE   * size_Mat];
        D.f[dirNW  ] = &distributions[dirNW   * size_Mat];
        D.f[dirTE  ] = &distributions[dirTE   * size_Mat];
        D.f[dirBW  ] = &distributions[dirBW   * size_Mat];
        D.f[dirBE  ] = &distributions[dirBE   * size_Mat];
        D.f[dirTW  ] = &distributions[dirTW   * size_Mat];
        D.f[dirTN  ] = &distributions[dirTN   * size_Mat];
        D.f[dirBS  ] = &distributions[dirBS   * size_Mat];
        D.f[dirBN  ] = &distributions[dirBN   * size_Mat];
        D.f[dirTS  ] = &distributions[dirTS   * size_Mat];
        D.f[dirREST] = &distributions[dirREST * size_Mat];
        D.f[dirTNE ] = &distributions[dirTNE  * size_Mat];
        D.f[dirTSW ] = &distributions[dirTSW  * size_Mat];
        D.f[dirTSE ] = &distributions[dirTSE  * size_Mat];
        D.f[dirTNW ] = &distributions[dirTNW  * size_Mat];
        D.f[dirBNE ] = &distributions[dirBNE  * size_Mat];
        D.f[dirBSW ] = &distributions[dirBSW  * size_Mat];
        D.f[dirBSE ] = &distributions[dirBSE  * size_Mat];
        D.f[dirBNW ] = &distributions[dirBNW  * size_Mat];
    }
    else
    {
        D.f[dirW   ] = &distributions[dirE    * size_Mat];
        D.f[dirE   ] = &distributions[dirW    * size_Mat];
        D.f[dirS   ] = &distributions[dirN    * size_Mat];
        D.f[dirN   ] = &distributions[dirS    * size_Mat];
        D.f[dirB   ] = &distributions[dirT    * size_Mat];
        D.f[dirT   ] = &distributions[dirB    * size_Mat];
        D.f[dirSW  ] = &distributions[dirNE   * size_Mat];
        D.f[dirNE  ] = &distributions[dirSW   * size_Mat];
        D.f[dirNW  ] = &distributions[dirSE   * size_Mat];
        D.f[dirSE  ] = &distributions[dirNW   * size_Mat];
        D.f[dirBW  ] = &distributions[dirTE   * size_Mat];
        D.f[dirTE  ] = &distributions[dirBW   * size_Mat];
        D.f[dirTW  ] = &distributions[dirBE   * size_Mat];
        D.f[dirBE  ] = &distributions[dirTW   * size_Mat];
        D.f[dirBS  ] = &distributions[dirTN   * size_Mat];
        D.f[dirTN  ] = &distributions[dirBS   * size_Mat];
        D.f[dirTS  ] = &distributions[dirBN   * size_Mat];
        D.f[dirBN  ] = &distributions[dirTS   * size_Mat];
        D.f[dirREST] = &distributions[dirREST * size_Mat];
        D.f[dirTNE ] = &distributions[dirBSW  * size_Mat];
        D.f[dirTSW ] = &distributions[dirBNE  * size_Mat];
        D.f[dirTSE ] = &distributions[dirBNW  * size_Mat];
        D.f[dirTNW ] = &distributions[dirBSE  * size_Mat];
        D.f[dirBNE ] = &distributions[dirTSW  * size_Mat];
        D.f[dirBSW ] = &distributions[dirTNE  * size_Mat];
        D.f[dirBSE ] = &distributions[dirTNW  * size_Mat];
        D.f[dirBNW ] = &distributions[dirTSE  * size_Mat];
    }
    ////////////////////////////////////////////////////////////////////////////////
    Distributions27 DAD;
    if (isEvenTimestep)
    {
        DAD.f[dirE   ] = &distributionsAD[dirE    * size_Mat];
        DAD.f[dirW   ] = &distributionsAD[dirW    * size_Mat];
        DAD.f[dirN   ] = &distributionsAD[dirN    * size_Mat];
        DAD.f[dirS   ] = &distributionsAD[dirS    * size_Mat];
        DAD.f[dirT   ] = &distributionsAD[dirT    * size_Mat];
        DAD.f[dirB   ] = &distributionsAD[dirB    * size_Mat];
        DAD.f[dirNE  ] = &distributionsAD[dirNE   * size_Mat];
        DAD.f[dirSW  ] = &distributionsAD[dirSW   * size_Mat];
        DAD.f[dirSE  ] = &distributionsAD[dirSE   * size_Mat];
        DAD.f[dirNW  ] = &distributionsAD[dirNW   * size_Mat];
        DAD.f[dirTE  ] = &distributionsAD[dirTE   * size_Mat];
        DAD.f[dirBW  ] = &distributionsAD[dirBW   * size_Mat];
        DAD.f[dirBE  ] = &distributionsAD[dirBE   * size_Mat];
        DAD.f[dirTW  ] = &distributionsAD[dirTW   * size_Mat];
        DAD.f[dirTN  ] = &distributionsAD[dirTN   * size_Mat];
        DAD.f[dirBS  ] = &distributionsAD[dirBS   * size_Mat];
        DAD.f[dirBN  ] = &distributionsAD[dirBN   * size_Mat];
        DAD.f[dirTS  ] = &distributionsAD[dirTS   * size_Mat];
        DAD.f[dirREST] = &distributionsAD[dirREST * size_Mat];
        DAD.f[dirTNE ] = &distributionsAD[dirTNE  * size_Mat];
        DAD.f[dirTSW ] = &distributionsAD[dirTSW  * size_Mat];
        DAD.f[dirTSE ] = &distributionsAD[dirTSE  * size_Mat];
        DAD.f[dirTNW ] = &distributionsAD[dirTNW  * size_Mat];
        DAD.f[dirBNE ] = &distributionsAD[dirBNE  * size_Mat];
        DAD.f[dirBSW ] = &distributionsAD[dirBSW  * size_Mat];
        DAD.f[dirBSE ] = &distributionsAD[dirBSE  * size_Mat];
        DAD.f[dirBNW ] = &distributionsAD[dirBNW  * size_Mat];
    }
    else
    {
        DAD.f[dirW   ] = &distributionsAD[dirE    * size_Mat];
        DAD.f[dirE   ] = &distributionsAD[dirW    * size_Mat];
        DAD.f[dirS   ] = &distributionsAD[dirN    * size_Mat];
        DAD.f[dirN   ] = &distributionsAD[dirS    * size_Mat];
        DAD.f[dirB   ] = &distributionsAD[dirT    * size_Mat];
        DAD.f[dirT   ] = &distributionsAD[dirB    * size_Mat];
        DAD.f[dirSW  ] = &distributionsAD[dirNE   * size_Mat];
        DAD.f[dirNE  ] = &distributionsAD[dirSW   * size_Mat];
        DAD.f[dirNW  ] = &distributionsAD[dirSE   * size_Mat];
        DAD.f[dirSE  ] = &distributionsAD[dirNW   * size_Mat];
        DAD.f[dirBW  ] = &distributionsAD[dirTE   * size_Mat];
        DAD.f[dirTE  ] = &distributionsAD[dirBW   * size_Mat];
        DAD.f[dirTW  ] = &distributionsAD[dirBE   * size_Mat];
        DAD.f[dirBE  ] = &distributionsAD[dirTW   * size_Mat];
        DAD.f[dirBS  ] = &distributionsAD[dirTN   * size_Mat];
        DAD.f[dirTN  ] = &distributionsAD[dirBS   * size_Mat];
        DAD.f[dirTS  ] = &distributionsAD[dirBN   * size_Mat];
        DAD.f[dirBN  ] = &distributionsAD[dirTS   * size_Mat];
        DAD.f[dirREST] = &distributionsAD[dirREST * size_Mat];
        DAD.f[dirTNE ] = &distributionsAD[dirBSW  * size_Mat];
        DAD.f[dirTSW ] = &distributionsAD[dirBNE  * size_Mat];
        DAD.f[dirTSE ] = &distributionsAD[dirBNW  * size_Mat];
        DAD.f[dirTNW ] = &distributionsAD[dirBSE  * size_Mat];
        DAD.f[dirBNE ] = &distributionsAD[dirTSW  * size_Mat];
        DAD.f[dirBSW ] = &distributionsAD[dirTNE  * size_Mat];
        DAD.f[dirBSE ] = &distributionsAD[dirTNW  * size_Mat];
        DAD.f[dirBNW ] = &distributionsAD[dirTSE  * size_Mat];
    }
    ////////////////////////////////////////////////////////////////////////////////
    const unsigned  x = threadIdx.x;  // Globaler x-Index
    const unsigned  y = blockIdx.x;   // Globaler y-Index
    const unsigned  z = blockIdx.y;   // Globaler z-Index

    const unsigned nx = blockDim.x;
    const unsigned ny = gridDim.x;

    const unsigned k = nx * (ny * z + y) + x;
    //////////////////////////////////////////////////////////////////////////

    if (k < numberOfBCnodes)
    {
        ////////////////////////////////////////////////////////////////////////////////
        real NormX = normalX[k];
        real NormY = normalY[k];
        real NormZ = normalZ[k];
        ////////////////////////////////////////////////////////////////////////////////
        real* q_dirE, * q_dirW, * q_dirN, * q_dirS, * q_dirT, * q_dirB,
            * q_dirNE, * q_dirSW, * q_dirSE, * q_dirNW, * q_dirTE, * q_dirBW,
            * q_dirBE, * q_dirTW, * q_dirTN, * q_dirBS, * q_dirBN, * q_dirTS,
            * q_dirTNE, * q_dirTSW, * q_dirTSE, * q_dirTNW, * q_dirBNE, * q_dirBSW,
            * q_dirBSE, * q_dirBNW;
        q_dirE   = &Qarrays[dirE   * numberOfBCnodes];
        q_dirW   = &Qarrays[dirW   * numberOfBCnodes];
        q_dirN   = &Qarrays[dirN   * numberOfBCnodes];
        q_dirS   = &Qarrays[dirS   * numberOfBCnodes];
        q_dirT   = &Qarrays[dirT   * numberOfBCnodes];
        q_dirB   = &Qarrays[dirB   * numberOfBCnodes];
        q_dirNE  = &Qarrays[dirNE  * numberOfBCnodes];
        q_dirSW  = &Qarrays[dirSW  * numberOfBCnodes];
        q_dirSE  = &Qarrays[dirSE  * numberOfBCnodes];
        q_dirNW  = &Qarrays[dirNW  * numberOfBCnodes];
        q_dirTE  = &Qarrays[dirTE  * numberOfBCnodes];
        q_dirBW  = &Qarrays[dirBW  * numberOfBCnodes];
        q_dirBE  = &Qarrays[dirBE  * numberOfBCnodes];
        q_dirTW  = &Qarrays[dirTW  * numberOfBCnodes];
        q_dirTN  = &Qarrays[dirTN  * numberOfBCnodes];
        q_dirBS  = &Qarrays[dirBS  * numberOfBCnodes];
        q_dirBN  = &Qarrays[dirBN  * numberOfBCnodes];
        q_dirTS  = &Qarrays[dirTS  * numberOfBCnodes];
        q_dirTNE = &Qarrays[dirTNE * numberOfBCnodes];
        q_dirTSW = &Qarrays[dirTSW * numberOfBCnodes];
        q_dirTSE = &Qarrays[dirTSE * numberOfBCnodes];
        q_dirTNW = &Qarrays[dirTNW * numberOfBCnodes];
        q_dirBNE = &Qarrays[dirBNE * numberOfBCnodes];
        q_dirBSW = &Qarrays[dirBSW * numberOfBCnodes];
        q_dirBSE = &Qarrays[dirBSE * numberOfBCnodes];
        q_dirBNW = &Qarrays[dirBNW * numberOfBCnodes];
        ////////////////////////////////////////////////////////////////////////////////
        //index
        unsigned int KQK   = QindexArray[k];
        unsigned int kzero = KQK;
        unsigned int ke    = KQK;
        unsigned int kw    = neighborX[KQK];
        unsigned int kn    = KQK;
        unsigned int ks    = neighborY[KQK];
        unsigned int kt    = KQK;
        unsigned int kb    = neighborZ[KQK];
        unsigned int ksw   = neighborY[kw];
        unsigned int kne   = KQK;
        unsigned int kse   = ks;
        unsigned int knw   = kw;
        unsigned int kbw   = neighborZ[kw];
        unsigned int kte   = KQK;
        unsigned int kbe   = kb;
        unsigned int ktw   = kw;
        unsigned int kbs   = neighborZ[ks];
        unsigned int ktn   = KQK;
        unsigned int kbn   = kb;
        unsigned int kts   = ks;
        unsigned int ktse  = ks;
        unsigned int kbnw  = kbw;
        unsigned int ktnw  = kw;
        unsigned int kbse  = kbs;
        unsigned int ktsw  = ksw;
        unsigned int kbne  = kb;
        unsigned int ktne  = KQK;
        unsigned int kbsw  = neighborZ[ksw];
        ////////////////////////////////////////////////////////////////////////////////
        real f_E, f_W, f_N, f_S, f_T, f_B, f_NE, f_SW, f_SE, f_NW, f_TE, f_BW, f_BE,
            f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

        f_W   = (D.f[dirE])[ke];
        f_E   = (D.f[dirW])[kw];
        f_S   = (D.f[dirN])[kn];
        f_N   = (D.f[dirS])[ks];
        f_B   = (D.f[dirT])[kt];
        f_T   = (D.f[dirB])[kb];
        f_SW  = (D.f[dirNE])[kne];
        f_NE  = (D.f[dirSW])[ksw];
        f_NW  = (D.f[dirSE])[kse];
        f_SE  = (D.f[dirNW])[knw];
        f_BW  = (D.f[dirTE])[kte];
        f_TE  = (D.f[dirBW])[kbw];
        f_TW  = (D.f[dirBE])[kbe];
        f_BE  = (D.f[dirTW])[ktw];
        f_BS  = (D.f[dirTN])[ktn];
        f_TN  = (D.f[dirBS])[kbs];
        f_TS  = (D.f[dirBN])[kbn];
        f_BN  = (D.f[dirTS])[kts];
        f_BSW = (D.f[dirTNE])[ktne];
        f_BNE = (D.f[dirTSW])[ktsw];
        f_BNW = (D.f[dirTSE])[ktse];
        f_BSE = (D.f[dirTNW])[ktnw];
        f_TSW = (D.f[dirBNE])[kbne];
        f_TNE = (D.f[dirBSW])[kbsw];
        f_TNW = (D.f[dirBSE])[kbse];
        f_TSE = (D.f[dirBNW])[kbnw];
        ////////////////////////////////////////////////////////////////////////////////
        real vx1, vx2, vx3, drho, q;
        drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]);

        vx1 = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
            ((f_BE - f_TW) + (f_TE - f_BW)) + ((f_SE - f_NW) + (f_NE - f_SW)) +
            (f_E - f_W)) / (c1o1 + drho);


        vx2 = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
            ((f_BN - f_TS) + (f_TN - f_BS)) + (-(f_SE - f_NW) + (f_NE - f_SW)) +
            (f_N - f_S)) / (c1o1 + drho);

        vx3 = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
            (-(f_BN - f_TS) + (f_TN - f_BS)) + ((f_TE - f_BW) - (f_BE - f_TW)) +
            (f_T - f_B)) / (c1o1 + drho);

        real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3) * (c1o1 + drho);

        ////////////////////////////////////////////////////////////////////////////////
        f_W   = (DAD.f[dirE])[ke];
        f_E   = (DAD.f[dirW])[kw];
        f_S   = (DAD.f[dirN])[kn];
        f_N   = (DAD.f[dirS])[ks];
        f_B   = (DAD.f[dirT])[kt];
        f_T   = (DAD.f[dirB])[kb];
        f_SW  = (DAD.f[dirNE])[kne];
        f_NE  = (DAD.f[dirSW])[ksw];
        f_NW  = (DAD.f[dirSE])[kse];
        f_SE  = (DAD.f[dirNW])[knw];
        f_BW  = (DAD.f[dirTE])[kte];
        f_TE  = (DAD.f[dirBW])[kbw];
        f_TW  = (DAD.f[dirBE])[kbe];
        f_BE  = (DAD.f[dirTW])[ktw];
        f_BS  = (DAD.f[dirTN])[ktn];
        f_TN  = (DAD.f[dirBS])[kbs];
        f_TS  = (DAD.f[dirBN])[kbn];
        f_BN  = (DAD.f[dirTS])[kts];
        f_BSW = (DAD.f[dirTNE])[ktne];
        f_BNE = (DAD.f[dirTSW])[ktsw];
        f_BNW = (DAD.f[dirTSE])[ktse];
        f_BSE = (DAD.f[dirTNW])[ktnw];
        f_TSW = (DAD.f[dirBNE])[kbne];
        f_TNE = (DAD.f[dirBSW])[kbsw];
        f_TNW = (DAD.f[dirBSE])[kbse];
        f_TSE = (DAD.f[dirBNW])[kbnw];
        //////////////////////////////////////////////////////////////////////////
        if (!isEvenTimestep)
        {
            DAD.f[dirE   ] = &distributionsAD[dirE    * size_Mat];
            DAD.f[dirW   ] = &distributionsAD[dirW    * size_Mat];
            DAD.f[dirN   ] = &distributionsAD[dirN    * size_Mat];
            DAD.f[dirS   ] = &distributionsAD[dirS    * size_Mat];
            DAD.f[dirT   ] = &distributionsAD[dirT    * size_Mat];
            DAD.f[dirB   ] = &distributionsAD[dirB    * size_Mat];
            DAD.f[dirNE  ] = &distributionsAD[dirNE   * size_Mat];
            DAD.f[dirSW  ] = &distributionsAD[dirSW   * size_Mat];
            DAD.f[dirSE  ] = &distributionsAD[dirSE   * size_Mat];
            DAD.f[dirNW  ] = &distributionsAD[dirNW   * size_Mat];
            DAD.f[dirTE  ] = &distributionsAD[dirTE   * size_Mat];
            DAD.f[dirBW  ] = &distributionsAD[dirBW   * size_Mat];
            DAD.f[dirBE  ] = &distributionsAD[dirBE   * size_Mat];
            DAD.f[dirTW  ] = &distributionsAD[dirTW   * size_Mat];
            DAD.f[dirTN  ] = &distributionsAD[dirTN   * size_Mat];
            DAD.f[dirBS  ] = &distributionsAD[dirBS   * size_Mat];
            DAD.f[dirBN  ] = &distributionsAD[dirBN   * size_Mat];
            DAD.f[dirTS  ] = &distributionsAD[dirTS   * size_Mat];
            DAD.f[dirREST] = &distributionsAD[dirREST * size_Mat];
            DAD.f[dirTNE ] = &distributionsAD[dirTNE  * size_Mat];
            DAD.f[dirTSW ] = &distributionsAD[dirTSW  * size_Mat];
            DAD.f[dirTSE ] = &distributionsAD[dirTSE  * size_Mat];
            DAD.f[dirTNW ] = &distributionsAD[dirTNW  * size_Mat];
            DAD.f[dirBNE ] = &distributionsAD[dirBNE  * size_Mat];
            DAD.f[dirBSW ] = &distributionsAD[dirBSW  * size_Mat];
            DAD.f[dirBSE ] = &distributionsAD[dirBSE  * size_Mat];
            DAD.f[dirBNW ] = &distributionsAD[dirBNW  * size_Mat];
        }
        else
        {
            DAD.f[dirW   ] = &distributionsAD[dirE    * size_Mat];
            DAD.f[dirE   ] = &distributionsAD[dirW    * size_Mat];
            DAD.f[dirS   ] = &distributionsAD[dirN    * size_Mat];
            DAD.f[dirN   ] = &distributionsAD[dirS    * size_Mat];
            DAD.f[dirB   ] = &distributionsAD[dirT    * size_Mat];
            DAD.f[dirT   ] = &distributionsAD[dirB    * size_Mat];
            DAD.f[dirSW  ] = &distributionsAD[dirNE   * size_Mat];
            DAD.f[dirNE  ] = &distributionsAD[dirSW   * size_Mat];
            DAD.f[dirNW  ] = &distributionsAD[dirSE   * size_Mat];
            DAD.f[dirSE  ] = &distributionsAD[dirNW   * size_Mat];
            DAD.f[dirBW  ] = &distributionsAD[dirTE   * size_Mat];
            DAD.f[dirTE  ] = &distributionsAD[dirBW   * size_Mat];
            DAD.f[dirTW  ] = &distributionsAD[dirBE   * size_Mat];
            DAD.f[dirBE  ] = &distributionsAD[dirTW   * size_Mat];
            DAD.f[dirBS  ] = &distributionsAD[dirTN   * size_Mat];
            DAD.f[dirTN  ] = &distributionsAD[dirBS   * size_Mat];
            DAD.f[dirTS  ] = &distributionsAD[dirBN   * size_Mat];
            DAD.f[dirBN  ] = &distributionsAD[dirTS   * size_Mat];
            DAD.f[dirREST] = &distributionsAD[dirREST * size_Mat];
            DAD.f[dirTNE ] = &distributionsAD[dirBSW  * size_Mat];
            DAD.f[dirTSW ] = &distributionsAD[dirBNE  * size_Mat];
            DAD.f[dirTSE ] = &distributionsAD[dirBNW  * size_Mat];
            DAD.f[dirTNW ] = &distributionsAD[dirBSE  * size_Mat];
            DAD.f[dirBNE ] = &distributionsAD[dirTSW  * size_Mat];
            DAD.f[dirBSW ] = &distributionsAD[dirTNE  * size_Mat];
            DAD.f[dirBSE ] = &distributionsAD[dirTNW  * size_Mat];
            DAD.f[dirBNW ] = &distributionsAD[dirTSE  * size_Mat];
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        real concentration =
            f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]);

        real jx1 =
            (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
            ((f_BE - f_TW) + (f_TE - f_BW)) + ((f_SE - f_NW) + (f_NE - f_SW)) +
            (f_E - f_W)) - (vx1 * concentration);

        real jx2 =
            ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
            ((f_BN - f_TS) + (f_TN - f_BS)) + (-(f_SE - f_NW) + (f_NE - f_SW)) +
            (f_N - f_S)) - (vx2 * concentration);

        real jx3 =
            (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
            (-(f_BN - f_TS) + (f_TN - f_BS)) + ((f_TE - f_BW) - (f_BE - f_TW)) +
            (f_T - f_B)) - (vx3 * concentration);

        //jx1 *= (c2o1 - omegaDiffusivity) / (c2o1 - c2o1 * omegaDiffusivity);
        //jx2 *= (c2o1 - omegaDiffusivity) / (c2o1 - c2o1 * omegaDiffusivity);
        //jx3 *= (c2o1 - omegaDiffusivity) / (c2o1 - c2o1 * omegaDiffusivity);

        real NormJ = jx1 * NormX + jx2 * NormY + jx3 * NormZ;

        real jTan1 = jx1 - NormJ * NormX;
        real jTan2 = jx2 - NormJ * NormY;
        real jTan3 = jx3 - NormJ * NormZ;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        q = q_dirE[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[dirW  ])[kw  ] = calcDistributionBC_AD(q, c2o27,   vx1,         cu_sq, f_E,   f_W,   omegaDiffusivity,        jTan1,       concentration); }
        q = q_dirW[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[dirE  ])[ke  ] = calcDistributionBC_AD(q, c2o27,  -vx1,         cu_sq, f_W,   f_E,   omegaDiffusivity,       -jTan1,       concentration); }
        q = q_dirN[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[dirS  ])[ks  ] = calcDistributionBC_AD(q, c2o27,   vx2,         cu_sq, f_N,   f_S,   omegaDiffusivity,        jTan2,       concentration); }
        q = q_dirS[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[dirN  ])[kn  ] = calcDistributionBC_AD(q, c2o27,  -vx2,         cu_sq, f_S,   f_N,   omegaDiffusivity,       -jTan2,       concentration); }
        q = q_dirT[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[dirB  ])[kb  ] = calcDistributionBC_AD(q, c2o27,   vx3,         cu_sq, f_T,   f_B,   omegaDiffusivity,        jTan3,       concentration); }
        q = q_dirB[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[dirT  ])[kt  ] = calcDistributionBC_AD(q, c2o27,  -vx3,         cu_sq, f_B,   f_T,   omegaDiffusivity,       -jTan3,       concentration); }
        q = q_dirNE[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dirSW ])[ksw ] = calcDistributionBC_AD(q, c1o54,   vx1+vx2,     cu_sq, f_NE,  f_SW,  omegaDiffusivity,  jTan1+jTan2,       concentration); }
        q = q_dirSW[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dirNE ])[kne ] = calcDistributionBC_AD(q, c1o54,  -vx1-vx2,     cu_sq, f_SW,  f_NE,  omegaDiffusivity, -jTan1-jTan2,       concentration); }
        q = q_dirSE[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dirNW ])[knw ] = calcDistributionBC_AD(q, c1o54,   vx1-vx2,     cu_sq, f_SE,  f_NW,  omegaDiffusivity,  jTan1-jTan2,       concentration); }
        q = q_dirNW[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dirSE ])[kse ] = calcDistributionBC_AD(q, c1o54,  -vx1+vx2,     cu_sq, f_NW,  f_SE,  omegaDiffusivity, -jTan1+jTan2,       concentration); }
        q = q_dirTE[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dirBW ])[kbw ] = calcDistributionBC_AD(q, c1o54,   vx1    +vx3, cu_sq, f_TE,  f_BW,  omegaDiffusivity,  jTan1      +jTan3, concentration); }
        q = q_dirBW[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dirTE ])[kte ] = calcDistributionBC_AD(q, c1o54,  -vx1    -vx3, cu_sq, f_BW,  f_TE,  omegaDiffusivity, -jTan1      -jTan3, concentration); }
        q = q_dirBE[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dirTW ])[ktw ] = calcDistributionBC_AD(q, c1o54,   vx1    -vx3, cu_sq, f_BE,  f_TW,  omegaDiffusivity,  jTan1      -jTan3, concentration); }
        q = q_dirTW[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dirBE ])[kbe ] = calcDistributionBC_AD(q, c1o54,  -vx1    +vx3, cu_sq, f_TW,  f_BE,  omegaDiffusivity, -jTan1      +jTan3, concentration); }
        q = q_dirTN[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dirBS ])[kbs ] = calcDistributionBC_AD(q, c1o54,       vx2+vx3, cu_sq, f_TN,  f_BS,  omegaDiffusivity,        jTan2+jTan3, concentration); }
        q = q_dirBS[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dirTN ])[ktn ] = calcDistributionBC_AD(q, c1o54,      -vx2-vx3, cu_sq, f_BS,  f_TN,  omegaDiffusivity,       -jTan2-jTan3, concentration); }
        q = q_dirBN[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dirTS ])[kts ] = calcDistributionBC_AD(q, c1o54,       vx2-vx3, cu_sq, f_BN,  f_TS,  omegaDiffusivity,        jTan2-jTan3, concentration); }
        q = q_dirTS[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dirBN ])[kbn ] = calcDistributionBC_AD(q, c1o54,      -vx2+vx3, cu_sq, f_TS,  f_BN,  omegaDiffusivity,       -jTan2+jTan3, concentration); }
        q = q_dirTNE[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dirBSW])[kbsw] = calcDistributionBC_AD(q, c1o216,  vx1+vx2+vx3, cu_sq, f_TNE, f_BSW, omegaDiffusivity,  jTan1+jTan2+jTan3, concentration); }
        q = q_dirBSW[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dirTNE])[ktne] = calcDistributionBC_AD(q, c1o216, -vx1-vx2-vx3, cu_sq, f_BSW, f_TNE, omegaDiffusivity, -jTan1-jTan2-jTan3, concentration); }
        q = q_dirBNE[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dirTSW])[ktsw] = calcDistributionBC_AD(q, c1o216,  vx1+vx2-vx3, cu_sq, f_BNE, f_TSW, omegaDiffusivity,  jTan1+jTan2-jTan3, concentration); }
        q = q_dirTSW[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dirBNE])[kbne] = calcDistributionBC_AD(q, c1o216, -vx1-vx2+vx3, cu_sq, f_TSW, f_BNE, omegaDiffusivity, -jTan1-jTan2+jTan3, concentration); }
        q = q_dirTSE[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dirBNW])[kbnw] = calcDistributionBC_AD(q, c1o216,  vx1-vx2+vx3, cu_sq, f_TSE, f_BNW, omegaDiffusivity,  jTan1-jTan2+jTan3, concentration); }
        q = q_dirBNW[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dirTSE])[ktse] = calcDistributionBC_AD(q, c1o216, -vx1+vx2-vx3, cu_sq, f_BNW, f_TSE, omegaDiffusivity, -jTan1+jTan2-jTan3, concentration); }
        q = q_dirBSE[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dirTNW])[ktnw] = calcDistributionBC_AD(q, c1o216,  vx1-vx2-vx3, cu_sq, f_BSE, f_TNW, omegaDiffusivity,  jTan1-jTan2-jTan3, concentration); }
        q = q_dirTNW[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dirBSE])[kbse] = calcDistributionBC_AD(q, c1o216, -vx1+vx2+vx3, cu_sq, f_TNW, f_BSE, omegaDiffusivity, -jTan1+jTan2+jTan3, concentration); }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
