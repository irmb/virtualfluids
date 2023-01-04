/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"

#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

//////////////////////////////////////////////////////////////////////////////
__global__ void QADPress7(  real* DD, 
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
                                       unsigned long long numberOfLBnodes, 
                                       bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
   }

   Distributions7 D7;
   if (isEvenTimestep==true)
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[1] = &DD7[1*numberOfLBnodes];
      D7.f[2] = &DD7[2*numberOfLBnodes];
      D7.f[3] = &DD7[3*numberOfLBnodes];
      D7.f[4] = &DD7[4*numberOfLBnodes];
      D7.f[5] = &DD7[5*numberOfLBnodes];
      D7.f[6] = &DD7[6*numberOfLBnodes];
   }
   else
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[2] = &DD7[1*numberOfLBnodes];
      D7.f[1] = &DD7[2*numberOfLBnodes];
      D7.f[4] = &DD7[3*numberOfLBnodes];
      D7.f[3] = &DD7[4*numberOfLBnodes];
      D7.f[6] = &DD7[5*numberOfLBnodes];
      D7.f[5] = &DD7[6*numberOfLBnodes];
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

      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      //q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      //q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      //q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      //q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      //q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      //q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      //q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      //q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      //q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      //q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      //q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      //q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      //q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      //q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      //q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      //q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      //q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      //q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      //q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      //q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
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

      f_W    = (D.f[DIR_P00])[ke   ];
      f_E    = (D.f[DIR_M00])[kw   ];
      f_S    = (D.f[DIR_0P0])[kn   ];
      f_N    = (D.f[DIR_0M0])[ks   ];
      f_B    = (D.f[DIR_00P])[kt   ];
      f_T    = (D.f[DIR_00M])[kb   ];
      f_SW   = (D.f[DIR_PP0])[kne  ];
      f_NE   = (D.f[DIR_MM0])[ksw  ];
      f_NW   = (D.f[DIR_PM0])[kse  ];
      f_SE   = (D.f[DIR_MP0])[knw  ];
      f_BW   = (D.f[DIR_P0P])[kte  ];
      f_TE   = (D.f[DIR_M0M])[kbw  ];
      f_TW   = (D.f[DIR_P0M])[kbe  ];
      f_BE   = (D.f[DIR_M0P])[ktw  ];
      f_BS   = (D.f[DIR_0PP])[ktn  ];
      f_TN   = (D.f[DIR_0MM])[kbs  ];
      f_TS   = (D.f[DIR_0PM])[kbn  ];
      f_BN   = (D.f[DIR_0MP])[kts  ];
      f_BSW  = (D.f[DIR_PPP])[ktne ];
      f_BNE  = (D.f[DIR_MMP])[ktsw ];
      f_BNW  = (D.f[DIR_PMP])[ktse ];
      f_BSE  = (D.f[DIR_MPP])[ktnw ];
      f_TSW  = (D.f[DIR_PPM])[kbne ];
      f_TNE  = (D.f[DIR_MMM])[kbsw ];
      f_TNW  = (D.f[DIR_PMM])[kbse ];
      f_TSE  = (D.f[DIR_MPM])[kbnw ];
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
      //            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]); 

      //real vx1 =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //               ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //               (f_E - f_W); 

      //real vx2 =  (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //               ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //               (f_N - f_S); 

      //real vx3 =  ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //               (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //               (f_T - f_B); 

      real rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ ((D.f[DIR_000])[kzero]);
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
         D7.f[0] = &DD7[0*numberOfLBnodes];
         D7.f[1] = &DD7[1*numberOfLBnodes];
         D7.f[2] = &DD7[2*numberOfLBnodes];
         D7.f[3] = &DD7[3*numberOfLBnodes];
         D7.f[4] = &DD7[4*numberOfLBnodes];
         D7.f[5] = &DD7[5*numberOfLBnodes];
         D7.f[6] = &DD7[6*numberOfLBnodes];
      }
      else
      {
         D7.f[0] = &DD7[0*numberOfLBnodes];
         D7.f[2] = &DD7[1*numberOfLBnodes];
         D7.f[1] = &DD7[2*numberOfLBnodes];
         D7.f[4] = &DD7[3*numberOfLBnodes];
         D7.f[3] = &DD7[4*numberOfLBnodes];
         D7.f[6] = &DD7[5*numberOfLBnodes];
         D7.f[5] = &DD7[6*numberOfLBnodes];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[DIR_000])[k]=0.1f;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //(D7.f[1])[ke   ] = f7_E - feq7_E + feqW7_W; //DIR_P00
      //(D7.f[2])[kw   ] = f7_W - feq7_W + feqW7_E; //DIR_M00
      //(D7.f[3])[kn   ] = f7_N - feq7_N + feqW7_S; //DIR_0P0
      //(D7.f[4])[ks   ] = f7_S - feq7_S + feqW7_N; //DIR_0M0
      //(D7.f[5])[kt   ] = f7_T - feq7_T + feqW7_B; //DIR_00P
      //(D7.f[6])[kb   ] = f7_B - feq7_B + feqW7_T; //DIR_00M

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
__global__ void QADPress27( real* DD, 
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
                                       unsigned long long numberOfLBnodes, 
                                       bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
   }

   Distributions27 D27;
   if (isEvenTimestep==true)
   {
      D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
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
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
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
      real f_W    = (D.f[DIR_P00])[ke   ];
      real f_E    = (D.f[DIR_M00])[kw   ];
      real f_S    = (D.f[DIR_0P0])[kn   ];
      real f_N    = (D.f[DIR_0M0])[ks   ];
      real f_B    = (D.f[DIR_00P])[kt   ];
      real f_T    = (D.f[DIR_00M])[kb   ];
      real f_SW   = (D.f[DIR_PP0])[kne  ];
      real f_NE   = (D.f[DIR_MM0])[ksw  ];
      real f_NW   = (D.f[DIR_PM0])[kse  ];
      real f_SE   = (D.f[DIR_MP0])[knw  ];
      real f_BW   = (D.f[DIR_P0P])[kte  ];
      real f_TE   = (D.f[DIR_M0M])[kbw  ];
      real f_TW   = (D.f[DIR_P0M])[kbe  ];
      real f_BE   = (D.f[DIR_M0P])[ktw  ];
      real f_BS   = (D.f[DIR_0PP])[ktn  ];
      real f_TN   = (D.f[DIR_0MM])[kbs  ];
      real f_TS   = (D.f[DIR_0PM])[kbn  ];
      real f_BN   = (D.f[DIR_0MP])[kts  ];
      real f_ZERO = (D.f[DIR_000])[kzero];
      real f_BSW  = (D.f[DIR_PPP])[ktne ];
      real f_BNE  = (D.f[DIR_MMP])[ktsw ];
      real f_BNW  = (D.f[DIR_PMP])[ktse ];
      real f_BSE  = (D.f[DIR_MPP])[ktnw ];
      real f_TSW  = (D.f[DIR_PPM])[kbne ];
      real f_TNE  = (D.f[DIR_MMM])[kbsw ];
      real f_TNW  = (D.f[DIR_PMM])[kbse ];
      real f_TSE  = (D.f[DIR_MPM])[kbnw ];
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
      real f27_W    = (D27.f[DIR_P00])[ke   ];
      real f27_E    = (D27.f[DIR_M00])[kw   ];
      real f27_S    = (D27.f[DIR_0P0])[kn   ];
      real f27_N    = (D27.f[DIR_0M0])[ks   ];
      real f27_B    = (D27.f[DIR_00P])[kt   ];
      real f27_T    = (D27.f[DIR_00M])[kb   ];
      real f27_SW   = (D27.f[DIR_PP0])[kne  ];
      real f27_NE   = (D27.f[DIR_MM0])[ksw  ];
      real f27_NW   = (D27.f[DIR_PM0])[kse  ];
      real f27_SE   = (D27.f[DIR_MP0])[knw  ];
      real f27_BW   = (D27.f[DIR_P0P])[kte  ];
      real f27_TE   = (D27.f[DIR_M0M])[kbw  ];
      real f27_TW   = (D27.f[DIR_P0M])[kbe  ];
      real f27_BE   = (D27.f[DIR_M0P])[ktw  ];
      real f27_BS   = (D27.f[DIR_0PP])[ktn  ];
      real f27_TN   = (D27.f[DIR_0MM])[kbs  ];
      real f27_TS   = (D27.f[DIR_0PM])[kbn  ];
      real f27_BN   = (D27.f[DIR_0MP])[kts  ];
      real f27_ZERO = (D27.f[DIR_000])[kzero];
      real f27_BSW  = (D27.f[DIR_PPP])[ktne ];
      real f27_BNE  = (D27.f[DIR_MMP])[ktsw ];
      real f27_BNW  = (D27.f[DIR_PMP])[ktse ];
      real f27_BSE  = (D27.f[DIR_MPP])[ktnw ];
      real f27_TSW  = (D27.f[DIR_PPM])[kbne ];
      real f27_TNE  = (D27.f[DIR_MMM])[kbsw ];
      real f27_TNW  = (D27.f[DIR_PMM])[kbse ];
      real f27_TSE  = (D27.f[DIR_MPM])[kbnw ];
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
         D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[DIR_000])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M00])[kw  ]=(c2o1*feqW27_W  -(f27_E  *(q*omegaD-c1o1)-omegaD*feq27_E  *(q-c1o1))/(omegaD-c1o1)+f27_W  *q)/(q+c1o1);
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P00])[ke  ]=(c2o1*feqW27_E  -(f27_W  *(q*omegaD-c1o1)-omegaD*feq27_W  *(q-c1o1))/(omegaD-c1o1)+f27_E  *q)/(q+c1o1);
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0M0])[ks  ]=(c2o1*feqW27_S  -(f27_N  *(q*omegaD-c1o1)-omegaD*feq27_N  *(q-c1o1))/(omegaD-c1o1)+f27_S  *q)/(q+c1o1);
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0P0])[kn  ]=(c2o1*feqW27_N  -(f27_S  *(q*omegaD-c1o1)-omegaD*feq27_S  *(q-c1o1))/(omegaD-c1o1)+f27_N  *q)/(q+c1o1);
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00M])[kb  ]=(c2o1*feqW27_B  -(f27_T  *(q*omegaD-c1o1)-omegaD*feq27_T  *(q-c1o1))/(omegaD-c1o1)+f27_B  *q)/(q+c1o1);
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00P])[kt  ]=(c2o1*feqW27_T  -(f27_B  *(q*omegaD-c1o1)-omegaD*feq27_B  *(q-c1o1))/(omegaD-c1o1)+f27_T  *q)/(q+c1o1);
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MM0])[ksw ]=(c2o1*feqW27_SW -(f27_NE *(q*omegaD-c1o1)-omegaD*feq27_NE *(q-c1o1))/(omegaD-c1o1)+f27_SW *q)/(q+c1o1);
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PP0])[kne ]=(c2o1*feqW27_NE -(f27_SW *(q*omegaD-c1o1)-omegaD*feq27_SW *(q-c1o1))/(omegaD-c1o1)+f27_NE *q)/(q+c1o1);
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MP0])[knw ]=(c2o1*feqW27_NW -(f27_SE *(q*omegaD-c1o1)-omegaD*feq27_SE *(q-c1o1))/(omegaD-c1o1)+f27_NW *q)/(q+c1o1);
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PM0])[kse ]=(c2o1*feqW27_SE -(f27_NW *(q*omegaD-c1o1)-omegaD*feq27_NW *(q-c1o1))/(omegaD-c1o1)+f27_SE *q)/(q+c1o1);
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0M])[kbw ]=(c2o1*feqW27_BW -(f27_TE *(q*omegaD-c1o1)-omegaD*feq27_TE *(q-c1o1))/(omegaD-c1o1)+f27_BW *q)/(q+c1o1);
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0P])[kte ]=(c2o1*feqW27_TE -(f27_BW *(q*omegaD-c1o1)-omegaD*feq27_BW *(q-c1o1))/(omegaD-c1o1)+f27_TE *q)/(q+c1o1);
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0P])[ktw ]=(c2o1*feqW27_TW -(f27_BE *(q*omegaD-c1o1)-omegaD*feq27_BE *(q-c1o1))/(omegaD-c1o1)+f27_TW *q)/(q+c1o1);
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0M])[kbe ]=(c2o1*feqW27_BE -(f27_TW *(q*omegaD-c1o1)-omegaD*feq27_TW *(q-c1o1))/(omegaD-c1o1)+f27_BE *q)/(q+c1o1);
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MM])[kbs ]=(c2o1*feqW27_BS -(f27_TN *(q*omegaD-c1o1)-omegaD*feq27_TN *(q-c1o1))/(omegaD-c1o1)+f27_BS *q)/(q+c1o1);
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PP])[ktn ]=(c2o1*feqW27_TN -(f27_BS *(q*omegaD-c1o1)-omegaD*feq27_BS *(q-c1o1))/(omegaD-c1o1)+f27_TN *q)/(q+c1o1);
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MP])[kts ]=(c2o1*feqW27_TS -(f27_BN *(q*omegaD-c1o1)-omegaD*feq27_BN *(q-c1o1))/(omegaD-c1o1)+f27_TS *q)/(q+c1o1);
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PM])[kbn ]=(c2o1*feqW27_BN -(f27_TS *(q*omegaD-c1o1)-omegaD*feq27_TS *(q-c1o1))/(omegaD-c1o1)+f27_BN *q)/(q+c1o1);
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMM])[kbsw]=(c2o1*feqW27_BSW-(f27_TNE*(q*omegaD-c1o1)-omegaD*feq27_TNE*(q-c1o1))/(omegaD-c1o1)+f27_BSW*q)/(q+c1o1);
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPP])[ktne]=(c2o1*feqW27_TNE-(f27_BSW*(q*omegaD-c1o1)-omegaD*feq27_BSW*(q-c1o1))/(omegaD-c1o1)+f27_TNE*q)/(q+c1o1);
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMP])[ktsw]=(c2o1*feqW27_TSW-(f27_BNE*(q*omegaD-c1o1)-omegaD*feq27_BNE*(q-c1o1))/(omegaD-c1o1)+f27_TSW*q)/(q+c1o1);
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPM])[kbne]=(c2o1*feqW27_BNE-(f27_TSW*(q*omegaD-c1o1)-omegaD*feq27_TSW*(q-c1o1))/(omegaD-c1o1)+f27_BNE*q)/(q+c1o1);
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPM])[kbnw]=(c2o1*feqW27_BNW-(f27_TSE*(q*omegaD-c1o1)-omegaD*feq27_TSE*(q-c1o1))/(omegaD-c1o1)+f27_BNW*q)/(q+c1o1);
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMP])[ktse]=(c2o1*feqW27_TSE-(f27_BNW*(q*omegaD-c1o1)-omegaD*feq27_BNW*(q-c1o1))/(omegaD-c1o1)+f27_TSE*q)/(q+c1o1);
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPP])[ktnw]=(c2o1*feqW27_TNW-(f27_BSE*(q*omegaD-c1o1)-omegaD*feq27_BSE*(q-c1o1))/(omegaD-c1o1)+f27_TNW*q)/(q+c1o1);
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMM])[kbse]=(c2o1*feqW27_BSE-(f27_TNW*(q*omegaD-c1o1)-omegaD*feq27_TNW*(q-c1o1))/(omegaD-c1o1)+f27_BSE*q)/(q+c1o1);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QADPressNEQNeighbor27(
													real* DD,
													real* DD27,
													int* k_Q,
													int* k_N,
													int numberOfBCnodes,
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													unsigned long long numberOfLBnodes,
													bool isEvenTimestep
												)
{
	Distributions27 D;
	if (isEvenTimestep == true)
	{
		D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
		D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
		D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
		D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
		D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
		D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
		D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
		D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
		D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
		D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
		D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
		D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
		D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
		D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
		D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
		D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
		D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
		D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
		D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
		D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
		D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
		D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
		D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
		D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
		D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
		D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
		D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
	}
	else
	{
		D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
		D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
		D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
		D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
		D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
		D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
		D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
		D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
		D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
		D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
		D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
		D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
		D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
		D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
		D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
		D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
		D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
		D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
		D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
		D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
		D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
		D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
		D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
		D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
		D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
		D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
		D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
	}

	Distributions27 D27;
	if (isEvenTimestep == true)
	{
		D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
		D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
		D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
		D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
		D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
		D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
		D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
		D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
		D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
		D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
		D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
		D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
		D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
		D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
		D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
		D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
		D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
		D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
		D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
		D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
		D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
		D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
		D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
		D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
		D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
		D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
		D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
	}
	else
	{
		D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
		D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
		D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
		D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
		D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
		D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
		D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
		D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
		D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
		D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
		D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
		D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
		D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
		D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
		D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
		D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
		D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
		D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
		D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
		D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
		D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
		D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
		D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
		D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
		D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
		D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
		D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
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
		real f_W =    (D.f[DIR_P00])[ke];
		real f_E =    (D.f[DIR_M00])[kw];
		real f_S =    (D.f[DIR_0P0])[kn];
		real f_N =    (D.f[DIR_0M0])[ks];
		real f_B =    (D.f[DIR_00P])[kt];
		real f_T =    (D.f[DIR_00M])[kb];
		real f_SW =   (D.f[DIR_PP0])[kne];
		real f_NE =   (D.f[DIR_MM0])[ksw];
		real f_NW =   (D.f[DIR_PM0])[kse];
		real f_SE =   (D.f[DIR_MP0])[knw];
		real f_BW =   (D.f[DIR_P0P])[kte];
		real f_TE =   (D.f[DIR_M0M])[kbw];
		real f_TW =   (D.f[DIR_P0M])[kbe];
		real f_BE =   (D.f[DIR_M0P])[ktw];
		real f_BS =   (D.f[DIR_0PP])[ktn];
		real f_TN =   (D.f[DIR_0MM])[kbs];
		real f_TS =   (D.f[DIR_0PM])[kbn];
		real f_BN =   (D.f[DIR_0MP])[kts];
		real f_ZERO = (D.f[DIR_000])[kzero];
		real f_BSW =  (D.f[DIR_PPP])[ktne];
		real f_BNE =  (D.f[DIR_MMP])[ktsw];
		real f_BNW =  (D.f[DIR_PMP])[ktse];
		real f_BSE =  (D.f[DIR_MPP])[ktnw];
		real f_TSW =  (D.f[DIR_PPM])[kbne];
		real f_TNE =  (D.f[DIR_MMM])[kbsw];
		real f_TNW =  (D.f[DIR_PMM])[kbse];
		real f_TSE =  (D.f[DIR_MPM])[kbnw];
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
		real f27_W =    (D27.f[DIR_P00])[ke];
		real f27_E =    (D27.f[DIR_M00])[kw];
		real f27_S =    (D27.f[DIR_0P0])[kn];
		real f27_N =    (D27.f[DIR_0M0])[ks];
		real f27_B =    (D27.f[DIR_00P])[kt];
		real f27_T =    (D27.f[DIR_00M])[kb];
		real f27_SW =   (D27.f[DIR_PP0])[kne];
		real f27_NE =   (D27.f[DIR_MM0])[ksw];
		real f27_NW =   (D27.f[DIR_PM0])[kse];
		real f27_SE =   (D27.f[DIR_MP0])[knw];
		real f27_BW =   (D27.f[DIR_P0P])[kte];
		real f27_TE =   (D27.f[DIR_M0M])[kbw];
		real f27_TW =   (D27.f[DIR_P0M])[kbe];
		real f27_BE =   (D27.f[DIR_M0P])[ktw];
		real f27_BS =   (D27.f[DIR_0PP])[ktn];
		real f27_TN =   (D27.f[DIR_0MM])[kbs];
		real f27_TS =   (D27.f[DIR_0PM])[kbn];
		real f27_BN =   (D27.f[DIR_0MP])[kts];
		real f27_ZERO = (D27.f[DIR_000])[kzero];
		real f27_BSW =  (D27.f[DIR_PPP])[ktne];
		real f27_BNE =  (D27.f[DIR_MMP])[ktsw];
		real f27_BNW =  (D27.f[DIR_PMP])[ktse];
		real f27_BSE =  (D27.f[DIR_MPP])[ktnw];
		real f27_TSW =  (D27.f[DIR_PPM])[kbne];
		real f27_TNE =  (D27.f[DIR_MMM])[kbsw];
		real f27_TNW =  (D27.f[DIR_PMM])[kbse];
		real f27_TSE =  (D27.f[DIR_MPM])[kbnw];
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
        (D27.f[DIR_P00])[kNe   ] = f27_W   ;  
        (D27.f[DIR_M00])[kNw   ] = f27_E   ;	
        (D27.f[DIR_0P0])[kNn   ] = f27_S   ;	
        (D27.f[DIR_0M0])[kNs   ] = f27_N   ;	
        (D27.f[DIR_00P])[kNt   ] = f27_B   ;	
        (D27.f[DIR_00M])[kNb   ] = f27_T   ;	
        (D27.f[DIR_PP0])[kNne  ] = f27_SW  ;	
        (D27.f[DIR_MM0])[kNsw  ] = f27_NE  ;	
        (D27.f[DIR_PM0])[kNse  ] = f27_NW  ;	
        (D27.f[DIR_MP0])[kNnw  ] = f27_SE  ;	
        (D27.f[DIR_P0P])[kNte  ] = f27_BW  ;	
        (D27.f[DIR_M0M])[kNbw  ] = f27_TE  ;	
        (D27.f[DIR_P0M])[kNbe  ] = f27_TW  ;	
        (D27.f[DIR_M0P])[kNtw  ] = f27_BE  ;	
        (D27.f[DIR_0PP])[kNtn  ] = f27_BS  ;	
        (D27.f[DIR_0MM])[kNbs  ] = f27_TN  ;	
        (D27.f[DIR_0PM])[kNbn  ] = f27_TS  ;	
        (D27.f[DIR_0MP])[kNts  ] = f27_BN  ;	
        (D27.f[DIR_000])[kNzero] = f27_ZERO;	
        (D27.f[DIR_PPP])[kNtne ] = f27_BSW ;	
        (D27.f[DIR_MMP])[kNtsw ] = f27_BNE ;	
        (D27.f[DIR_PMP])[kNtse ] = f27_BNW ;	
        (D27.f[DIR_MPP])[kNtnw ] = f27_BSE ;	
        (D27.f[DIR_PPM])[kNbne ] = f27_TSW ;	
        (D27.f[DIR_MMM])[kNbsw ] = f27_TNE ;	
        (D27.f[DIR_PMM])[kNbse ] = f27_TNW ;	
        (D27.f[DIR_MPM])[kNbnw ] = f27_TSE ;       
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QADVel7( real* DD, 
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
                                    unsigned long long numberOfLBnodes, 
                                    bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
   }

   Distributions7 D7;
   if (isEvenTimestep==true)
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[1] = &DD7[1*numberOfLBnodes];
      D7.f[2] = &DD7[2*numberOfLBnodes];
      D7.f[3] = &DD7[3*numberOfLBnodes];
      D7.f[4] = &DD7[4*numberOfLBnodes];
      D7.f[5] = &DD7[5*numberOfLBnodes];
      D7.f[6] = &DD7[6*numberOfLBnodes];
   }
   else
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[2] = &DD7[1*numberOfLBnodes];
      D7.f[1] = &DD7[2*numberOfLBnodes];
      D7.f[4] = &DD7[3*numberOfLBnodes];
      D7.f[3] = &DD7[4*numberOfLBnodes];
      D7.f[6] = &DD7[5*numberOfLBnodes];
      D7.f[5] = &DD7[6*numberOfLBnodes];
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

      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
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

      f_W    = (D.f[DIR_P00])[ke   ];
      f_E    = (D.f[DIR_M00])[kw   ];
      f_S    = (D.f[DIR_0P0])[kn   ];
      f_N    = (D.f[DIR_0M0])[ks   ];
      f_B    = (D.f[DIR_00P])[kt   ];
      f_T    = (D.f[DIR_00M])[kb   ];
      f_SW   = (D.f[DIR_PP0])[kne  ];
      f_NE   = (D.f[DIR_MM0])[ksw  ];
      f_NW   = (D.f[DIR_PM0])[kse  ];
      f_SE   = (D.f[DIR_MP0])[knw  ];
      f_BW   = (D.f[DIR_P0P])[kte  ];
      f_TE   = (D.f[DIR_M0M])[kbw  ];
      f_TW   = (D.f[DIR_P0M])[kbe  ];
      f_BE   = (D.f[DIR_M0P])[ktw  ];
      f_BS   = (D.f[DIR_0PP])[ktn  ];
      f_TN   = (D.f[DIR_0MM])[kbs  ];
      f_TS   = (D.f[DIR_0PM])[kbn  ];
      f_BN   = (D.f[DIR_0MP])[kts  ];
      f_BSW  = (D.f[DIR_PPP])[ktne ];
      f_BNE  = (D.f[DIR_MMP])[ktsw ];
      f_BNW  = (D.f[DIR_PMP])[ktse ];
      f_BSE  = (D.f[DIR_MPP])[ktnw ];
      f_TSW  = (D.f[DIR_PPM])[kbne ];
      f_TNE  = (D.f[DIR_MMM])[kbsw ];
      f_TNW  = (D.f[DIR_PMM])[kbse ];
      f_TSE  = (D.f[DIR_MPM])[kbnw ];
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
      ////            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]); 

      //real vx1 =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //               ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //               (f_E - f_W); 

      //real vx2 =  (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //               ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //               (f_N - f_S); 

      //real vx3 =  ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //               (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //               (f_T - f_B); 

      real rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ ((D.f[DIR_000])[kzero]);
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
         D7.f[0] = &DD7[0*numberOfLBnodes];
         D7.f[1] = &DD7[1*numberOfLBnodes];
         D7.f[2] = &DD7[2*numberOfLBnodes];
         D7.f[3] = &DD7[3*numberOfLBnodes];
         D7.f[4] = &DD7[4*numberOfLBnodes];
         D7.f[5] = &DD7[5*numberOfLBnodes];
         D7.f[6] = &DD7[6*numberOfLBnodes];
      }
      else
      {
         D7.f[0] = &DD7[0*numberOfLBnodes];
         D7.f[2] = &DD7[1*numberOfLBnodes];
         D7.f[1] = &DD7[2*numberOfLBnodes];
         D7.f[4] = &DD7[3*numberOfLBnodes];
         D7.f[3] = &DD7[4*numberOfLBnodes];
         D7.f[6] = &DD7[5*numberOfLBnodes];
         D7.f[5] = &DD7[6*numberOfLBnodes];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[DIR_000])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //(D7.f[1])[ke   ] = f7_E - feq7_E + feqW7_W; //DIR_P00
      //(D7.f[2])[kw   ] = f7_W - feq7_W + feqW7_E; //DIR_M00
      //(D7.f[3])[kn   ] = f7_N - feq7_N + feqW7_S; //DIR_0P0
      //(D7.f[4])[ks   ] = f7_S - feq7_S + feqW7_N; //DIR_0M0
      //(D7.f[5])[kt   ] = f7_T - feq7_T + feqW7_B; //DIR_00P
      //(D7.f[6])[kb   ] = f7_B - feq7_B + feqW7_T; //DIR_00M

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
__global__ void QADVel27(real* DD, 
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
                                    unsigned long long numberOfLBnodes, 
                                    bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
   }

   Distributions27 D27;
   if (isEvenTimestep==true)
   {
      D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
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
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
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
      real f_W    = (D.f[DIR_P00])[ke   ];
      real f_E    = (D.f[DIR_M00])[kw   ];
      real f_S    = (D.f[DIR_0P0])[kn   ];
      real f_N    = (D.f[DIR_0M0])[ks   ];
      real f_B    = (D.f[DIR_00P])[kt   ];
      real f_T    = (D.f[DIR_00M])[kb   ];
      real f_SW   = (D.f[DIR_PP0])[kne  ];
      real f_NE   = (D.f[DIR_MM0])[ksw  ];
      real f_NW   = (D.f[DIR_PM0])[kse  ];
      real f_SE   = (D.f[DIR_MP0])[knw  ];
      real f_BW   = (D.f[DIR_P0P])[kte  ];
      real f_TE   = (D.f[DIR_M0M])[kbw  ];
      real f_TW   = (D.f[DIR_P0M])[kbe  ];
      real f_BE   = (D.f[DIR_M0P])[ktw  ];
      real f_BS   = (D.f[DIR_0PP])[ktn  ];
      real f_TN   = (D.f[DIR_0MM])[kbs  ];
      real f_TS   = (D.f[DIR_0PM])[kbn  ];
      real f_BN   = (D.f[DIR_0MP])[kts  ];
      real f_ZERO = (D.f[DIR_000])[kzero];
      real f_BSW  = (D.f[DIR_PPP])[ktne ];
      real f_BNE  = (D.f[DIR_MMP])[ktsw ];
      real f_BNW  = (D.f[DIR_PMP])[ktse ];
      real f_BSE  = (D.f[DIR_MPP])[ktnw ];
      real f_TSW  = (D.f[DIR_PPM])[kbne ];
      real f_TNE  = (D.f[DIR_MMM])[kbsw ];
      real f_TNW  = (D.f[DIR_PMM])[kbse ];
      real f_TSE  = (D.f[DIR_MPM])[kbnw ];
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
      //real f27_W    = (D27.f[DIR_P00])[ke   ];
      //real f27_E    = (D27.f[DIR_M00])[kw   ];
      //real f27_S    = (D27.f[DIR_0P0])[kn   ];
      //real f27_N    = (D27.f[DIR_0M0])[ks   ];
      //real f27_B    = (D27.f[DIR_00P])[kt   ];
      //real f27_T    = (D27.f[DIR_00M])[kb   ];
      //real f27_SW   = (D27.f[DIR_PP0])[kne  ];
      //real f27_NE   = (D27.f[DIR_MM0])[ksw  ];
      //real f27_NW   = (D27.f[DIR_PM0])[kse  ];
      //real f27_SE   = (D27.f[DIR_MP0])[knw  ];
      //real f27_BW   = (D27.f[DIR_P0P])[kte  ];
      //real f27_TE   = (D27.f[DIR_M0M])[kbw  ];
      //real f27_TW   = (D27.f[DIR_P0M])[kbe  ];
      //real f27_BE   = (D27.f[DIR_M0P])[ktw  ];
      //real f27_BS   = (D27.f[DIR_0PP])[ktn  ];
      //real f27_TN   = (D27.f[DIR_0MM])[kbs  ];
      //real f27_TS   = (D27.f[DIR_0PM])[kbn  ];
      //real f27_BN   = (D27.f[DIR_0MP])[kts  ];
      //real f27_ZERO = (D27.f[DIR_000])[kzero];
      //real f27_BSW  = (D27.f[DIR_PPP])[ktne ];
      //real f27_BNE  = (D27.f[DIR_MMP])[ktsw ];
      //real f27_BNW  = (D27.f[DIR_PMP])[ktse ];
      //real f27_BSE  = (D27.f[DIR_MPP])[ktnw ];
      //real f27_TSW  = (D27.f[DIR_PPM])[kbne ];
      //real f27_TNE  = (D27.f[DIR_MMM])[kbsw ];
      //real f27_TNW  = (D27.f[DIR_PMM])[kbse ];
      //real f27_TSE  = (D27.f[DIR_MPM])[kbnw ];
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
         D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[DIR_000])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[DIR_000])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //(D27.f[DIR_M00])[kw  ]= four;
      //(D27.f[DIR_P00])[ke  ]= four;
      //(D27.f[DIR_0M0])[ks  ]= four;
      //(D27.f[DIR_0P0])[kn  ]= four;
      //(D27.f[DIR_00M])[kb  ]= four;
      //(D27.f[DIR_00P])[kt  ]= four;
      //(D27.f[DIR_MM0])[ksw ]= four;
      //(D27.f[DIR_PP0])[kne ]= four;
      //(D27.f[DIR_MP0])[knw ]= four;
      //(D27.f[DIR_PM0])[kse ]= four;
      //(D27.f[DIR_M0M])[kbw ]= four;
      //(D27.f[DIR_P0P])[kte ]= four;
      //(D27.f[DIR_M0P])[ktw ]= four;
      //(D27.f[DIR_P0M])[kbe ]= four;
      //(D27.f[DIR_0MM])[kbs ]= four;
      //(D27.f[DIR_0PP])[ktn ]= four;
      //(D27.f[DIR_0MP])[kts ]= four;
      //(D27.f[DIR_0PM])[kbn ]= four;
      //(D27.f[DIR_MMM])[kbsw]= four;
      //(D27.f[DIR_PPP])[ktne]= four;
      //(D27.f[DIR_MMP])[ktsw]= four;
      //(D27.f[DIR_PPM])[kbne]= four;
      //(D27.f[DIR_MPM])[kbnw]= four;
      //(D27.f[DIR_PMP])[ktse]= four;
      //(D27.f[DIR_MPP])[ktnw]= four;
      //(D27.f[DIR_PMM])[kbse]= four;
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M00])[kw  ]= -feqW27_W  + c2o1 * c2o27  * TempD;
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P00])[ke  ]= -feqW27_E  + c2o1 * c2o27  * TempD;
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0M0])[ks  ]= -feqW27_S  + c2o1 * c2o27  * TempD;
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0P0])[kn  ]= -feqW27_N  + c2o1 * c2o27  * TempD;
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00M])[kb  ]= -feqW27_B  + c2o1 * c2o27  * TempD;
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00P])[kt  ]= -feqW27_T  + c2o1 * c2o27  * TempD;
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MM0])[ksw ]= -feqW27_SW + c2o1 * c1o54  * TempD;
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PP0])[kne ]= -feqW27_NE + c2o1 * c1o54  * TempD;
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MP0])[knw ]= -feqW27_NW + c2o1 * c1o54  * TempD;
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PM0])[kse ]= -feqW27_SE + c2o1 * c1o54  * TempD;
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0M])[kbw ]= -feqW27_BW + c2o1 * c1o54  * TempD;
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0P])[kte ]= -feqW27_TE + c2o1 * c1o54  * TempD;
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0P])[ktw ]= -feqW27_TW + c2o1 * c1o54  * TempD;
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0M])[kbe ]= -feqW27_BE + c2o1 * c1o54  * TempD;
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MM])[kbs ]= -feqW27_BS + c2o1 * c1o54  * TempD;
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PP])[ktn ]= -feqW27_TN + c2o1 * c1o54  * TempD;
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MP])[kts ]= -feqW27_TS + c2o1 * c1o54  * TempD;
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PM])[kbn ]= -feqW27_BN + c2o1 * c1o54  * TempD;
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMM])[kbsw]= -feqW27_BSW+ c2o1 * c1o216 * TempD;
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPP])[ktne]= -feqW27_TNE+ c2o1 * c1o216 * TempD;
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMP])[ktsw]= -feqW27_TSW+ c2o1 * c1o216 * TempD;
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPM])[kbne]= -feqW27_BNE+ c2o1 * c1o216 * TempD;
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPM])[kbnw]= -feqW27_BNW+ c2o1 * c1o216 * TempD;
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMP])[ktse]= -feqW27_TSE+ c2o1 * c1o216 * TempD;
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPP])[ktnw]= -feqW27_TNW+ c2o1 * c1o216 * TempD;
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMM])[kbse]= -feqW27_BSE+ c2o1 * c1o216 * TempD;
      //q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[DIR_M00])[kw  ]=(two*feqW27_W  -(f27_E  *(q*omegaD-one)-omegaD*feq27_E  *(q-one))/(omegaD-one)+f27_W  *q)/(q+one);
      //q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[DIR_P00])[ke  ]=(two*feqW27_E  -(f27_W  *(q*omegaD-one)-omegaD*feq27_W  *(q-one))/(omegaD-one)+f27_E  *q)/(q+one);
      //q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[DIR_0M0])[ks  ]=(two*feqW27_S  -(f27_N  *(q*omegaD-one)-omegaD*feq27_N  *(q-one))/(omegaD-one)+f27_S  *q)/(q+one);
      //q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[DIR_0P0])[kn  ]=(two*feqW27_N  -(f27_S  *(q*omegaD-one)-omegaD*feq27_S  *(q-one))/(omegaD-one)+f27_N  *q)/(q+one);
      //q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[DIR_00M])[kb  ]=(two*feqW27_B  -(f27_T  *(q*omegaD-one)-omegaD*feq27_T  *(q-one))/(omegaD-one)+f27_B  *q)/(q+one);
      //q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[DIR_00P])[kt  ]=(two*feqW27_T  -(f27_B  *(q*omegaD-one)-omegaD*feq27_B  *(q-one))/(omegaD-one)+f27_T  *q)/(q+one);
      //q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[DIR_MM0])[ksw ]=(two*feqW27_SW -(f27_NE *(q*omegaD-one)-omegaD*feq27_NE *(q-one))/(omegaD-one)+f27_SW *q)/(q+one);
      //q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[DIR_PP0])[kne ]=(two*feqW27_NE -(f27_SW *(q*omegaD-one)-omegaD*feq27_SW *(q-one))/(omegaD-one)+f27_NE *q)/(q+one);
      //q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[DIR_MP0])[knw ]=(two*feqW27_NW -(f27_SE *(q*omegaD-one)-omegaD*feq27_SE *(q-one))/(omegaD-one)+f27_NW *q)/(q+one);
      //q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[DIR_PM0])[kse ]=(two*feqW27_SE -(f27_NW *(q*omegaD-one)-omegaD*feq27_NW *(q-one))/(omegaD-one)+f27_SE *q)/(q+one);
      //q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[DIR_M0M])[kbw ]=(two*feqW27_BW -(f27_TE *(q*omegaD-one)-omegaD*feq27_TE *(q-one))/(omegaD-one)+f27_BW *q)/(q+one);
      //q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[DIR_P0P])[kte ]=(two*feqW27_TE -(f27_BW *(q*omegaD-one)-omegaD*feq27_BW *(q-one))/(omegaD-one)+f27_TE *q)/(q+one);
      //q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[DIR_M0P])[ktw ]=(two*feqW27_TW -(f27_BE *(q*omegaD-one)-omegaD*feq27_BE *(q-one))/(omegaD-one)+f27_TW *q)/(q+one);
      //q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[DIR_P0M])[kbe ]=(two*feqW27_BE -(f27_TW *(q*omegaD-one)-omegaD*feq27_TW *(q-one))/(omegaD-one)+f27_BE *q)/(q+one);
      //q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[DIR_0MM])[kbs ]=(two*feqW27_BS -(f27_TN *(q*omegaD-one)-omegaD*feq27_TN *(q-one))/(omegaD-one)+f27_BS *q)/(q+one);
      //q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[DIR_0PP])[ktn ]=(two*feqW27_TN -(f27_BS *(q*omegaD-one)-omegaD*feq27_BS *(q-one))/(omegaD-one)+f27_TN *q)/(q+one);
      //q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[DIR_0MP])[kts ]=(two*feqW27_TS -(f27_BN *(q*omegaD-one)-omegaD*feq27_BN *(q-one))/(omegaD-one)+f27_TS *q)/(q+one);
      //q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[DIR_0PM])[kbn ]=(two*feqW27_BN -(f27_TS *(q*omegaD-one)-omegaD*feq27_TS *(q-one))/(omegaD-one)+f27_BN *q)/(q+one);
      //q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[DIR_MMM])[kbsw]=(two*feqW27_BSW-(f27_TNE*(q*omegaD-one)-omegaD*feq27_TNE*(q-one))/(omegaD-one)+f27_BSW*q)/(q+one);
      //q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[DIR_PPP])[ktne]=(two*feqW27_TNE-(f27_BSW*(q*omegaD-one)-omegaD*feq27_BSW*(q-one))/(omegaD-one)+f27_TNE*q)/(q+one);
      //q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[DIR_MMP])[ktsw]=(two*feqW27_TSW-(f27_BNE*(q*omegaD-one)-omegaD*feq27_BNE*(q-one))/(omegaD-one)+f27_TSW*q)/(q+one);
      //q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[DIR_PPM])[kbne]=(two*feqW27_BNE-(f27_TSW*(q*omegaD-one)-omegaD*feq27_TSW*(q-one))/(omegaD-one)+f27_BNE*q)/(q+one);
      //q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[DIR_MPM])[kbnw]=(two*feqW27_BNW-(f27_TSE*(q*omegaD-one)-omegaD*feq27_TSE*(q-one))/(omegaD-one)+f27_BNW*q)/(q+one);
      //q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[DIR_PMP])[ktse]=(two*feqW27_TSE-(f27_BNW*(q*omegaD-one)-omegaD*feq27_BNW*(q-one))/(omegaD-one)+f27_TSE*q)/(q+one);
      //q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[DIR_MPP])[ktnw]=(two*feqW27_TNW-(f27_BSE*(q*omegaD-one)-omegaD*feq27_BSE*(q-one))/(omegaD-one)+f27_TNW*q)/(q+one);
      //q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[DIR_PMM])[kbse]=(two*feqW27_BSE-(f27_TNW*(q*omegaD-one)-omegaD*feq27_TNW*(q-one))/(omegaD-one)+f27_BSE*q)/(q+one);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QAD7( real* DD, 
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
                                 unsigned long long numberOfLBnodes, 
                                 bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
   }

   Distributions7 D7;
   if (isEvenTimestep==true)
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[1] = &DD7[1*numberOfLBnodes];
      D7.f[2] = &DD7[2*numberOfLBnodes];
      D7.f[3] = &DD7[3*numberOfLBnodes];
      D7.f[4] = &DD7[4*numberOfLBnodes];
      D7.f[5] = &DD7[5*numberOfLBnodes];
      D7.f[6] = &DD7[6*numberOfLBnodes];
   }
   else
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[2] = &DD7[1*numberOfLBnodes];
      D7.f[1] = &DD7[2*numberOfLBnodes];
      D7.f[4] = &DD7[3*numberOfLBnodes];
      D7.f[3] = &DD7[4*numberOfLBnodes];
      D7.f[6] = &DD7[5*numberOfLBnodes];
      D7.f[5] = &DD7[6*numberOfLBnodes];
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

      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      //q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      //q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      //q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      //q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      //q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      //q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      //q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      //q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      //q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      //q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      //q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      //q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      //q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      //q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      //q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      //q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      //q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      //q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      //q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      //q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
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

      f_W    = (D.f[DIR_P00])[ke   ];
      f_E    = (D.f[DIR_M00])[kw   ];
      f_S    = (D.f[DIR_0P0])[kn   ];
      f_N    = (D.f[DIR_0M0])[ks   ];
      f_B    = (D.f[DIR_00P])[kt   ];
      f_T    = (D.f[DIR_00M])[kb   ];
      f_SW   = (D.f[DIR_PP0])[kne  ];
      f_NE   = (D.f[DIR_MM0])[ksw  ];
      f_NW   = (D.f[DIR_PM0])[kse  ];
      f_SE   = (D.f[DIR_MP0])[knw  ];
      f_BW   = (D.f[DIR_P0P])[kte  ];
      f_TE   = (D.f[DIR_M0M])[kbw  ];
      f_TW   = (D.f[DIR_P0M])[kbe  ];
      f_BE   = (D.f[DIR_M0P])[ktw  ];
      f_BS   = (D.f[DIR_0PP])[ktn  ];
      f_TN   = (D.f[DIR_0MM])[kbs  ];
      f_TS   = (D.f[DIR_0PM])[kbn  ];
      f_BN   = (D.f[DIR_0MP])[kts  ];
      f_BSW  = (D.f[DIR_PPP])[ktne ];
      f_BNE  = (D.f[DIR_MMP])[ktsw ];
      f_BNW  = (D.f[DIR_PMP])[ktse ];
      f_BSE  = (D.f[DIR_MPP])[ktnw ];
      f_TSW  = (D.f[DIR_PPM])[kbne ];
      f_TNE  = (D.f[DIR_MMM])[kbsw ];
      f_TNW  = (D.f[DIR_PMM])[kbse ];
      f_TSE  = (D.f[DIR_MPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3/*, drho*/;
      //drho   =    f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
      //            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
      //            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]); 

      //vx1    = ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //         (f_E - f_W); 


      //vx2    = (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //         ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //         (f_N - f_S); 

      //vx3    = ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //         (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //         (f_T - f_B); 
      real rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+ ((D.f[DIR_000])[kzero]);
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
         D7.f[0] = &DD7[0*numberOfLBnodes];
         D7.f[1] = &DD7[1*numberOfLBnodes];
         D7.f[2] = &DD7[2*numberOfLBnodes];
         D7.f[3] = &DD7[3*numberOfLBnodes];
         D7.f[4] = &DD7[4*numberOfLBnodes];
         D7.f[5] = &DD7[5*numberOfLBnodes];
         D7.f[6] = &DD7[6*numberOfLBnodes];
      }
      else
      {
         D7.f[0] = &DD7[0*numberOfLBnodes];
         D7.f[2] = &DD7[1*numberOfLBnodes];
         D7.f[1] = &DD7[2*numberOfLBnodes];
         D7.f[4] = &DD7[3*numberOfLBnodes];
         D7.f[3] = &DD7[4*numberOfLBnodes];
         D7.f[6] = &DD7[5*numberOfLBnodes];
         D7.f[5] = &DD7[6*numberOfLBnodes];
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[DIR_000])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //(D7.f[1])[ke   ] = f7_E - feq7_E + feqW7_W; //DIR_P00
      //(D7.f[2])[kw   ] = f7_W - feq7_W + feqW7_E; //DIR_M00
      //(D7.f[3])[kn   ] = f7_N - feq7_N + feqW7_S; //DIR_0P0
      //(D7.f[4])[ks   ] = f7_S - feq7_S + feqW7_N; //DIR_0M0
      //(D7.f[5])[kt   ] = f7_T - feq7_T + feqW7_B; //DIR_00P
      //(D7.f[6])[kb   ] = f7_B - feq7_B + feqW7_T; //DIR_00M

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
__global__ void QADDirichlet27(
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
											 unsigned long long numberOfLBnodes, 
											 bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
   }

   Distributions27 D27;
   if (isEvenTimestep==true)
   {
      D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
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
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
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
      real f_W    = (D.f[DIR_P00])[ke   ];
      real f_E    = (D.f[DIR_M00])[kw   ];
      real f_S    = (D.f[DIR_0P0])[kn   ];
      real f_N    = (D.f[DIR_0M0])[ks   ];
      real f_B    = (D.f[DIR_00P])[kt   ];
      real f_T    = (D.f[DIR_00M])[kb   ];
      real f_SW   = (D.f[DIR_PP0])[kne  ];
      real f_NE   = (D.f[DIR_MM0])[ksw  ];
      real f_NW   = (D.f[DIR_PM0])[kse  ];
      real f_SE   = (D.f[DIR_MP0])[knw  ];
      real f_BW   = (D.f[DIR_P0P])[kte  ];
      real f_TE   = (D.f[DIR_M0M])[kbw  ];
      real f_TW   = (D.f[DIR_P0M])[kbe  ];
      real f_BE   = (D.f[DIR_M0P])[ktw  ];
      real f_BS   = (D.f[DIR_0PP])[ktn  ];
      real f_TN   = (D.f[DIR_0MM])[kbs  ];
      real f_TS   = (D.f[DIR_0PM])[kbn  ];
      real f_BN   = (D.f[DIR_0MP])[kts  ];
      real f_ZERO = (D.f[DIR_000])[kzero];
      real f_BSW  = (D.f[DIR_PPP])[ktne ];
      real f_BNE  = (D.f[DIR_MMP])[ktsw ];
      real f_BNW  = (D.f[DIR_PMP])[ktse ];
      real f_BSE  = (D.f[DIR_MPP])[ktnw ];
      real f_TSW  = (D.f[DIR_PPM])[kbne ];
      real f_TNE  = (D.f[DIR_MMM])[kbsw ];
      real f_TNW  = (D.f[DIR_PMM])[kbse ];
      real f_TSE  = (D.f[DIR_MPM])[kbnw ];
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
      real f27_W    = (D27.f[DIR_P00])[ke   ];
      real f27_E    = (D27.f[DIR_M00])[kw   ];
      real f27_S    = (D27.f[DIR_0P0])[kn   ];
      real f27_N    = (D27.f[DIR_0M0])[ks   ];
      real f27_B    = (D27.f[DIR_00P])[kt   ];
      real f27_T    = (D27.f[DIR_00M])[kb   ];
      real f27_SW   = (D27.f[DIR_PP0])[kne  ];
      real f27_NE   = (D27.f[DIR_MM0])[ksw  ];
      real f27_NW   = (D27.f[DIR_PM0])[kse  ];
      real f27_SE   = (D27.f[DIR_MP0])[knw  ];
      real f27_BW   = (D27.f[DIR_P0P])[kte  ];
      real f27_TE   = (D27.f[DIR_M0M])[kbw  ];
      real f27_TW   = (D27.f[DIR_P0M])[kbe  ];
      real f27_BE   = (D27.f[DIR_M0P])[ktw  ];
      real f27_BS   = (D27.f[DIR_0PP])[ktn  ];
      real f27_TN   = (D27.f[DIR_0MM])[kbs  ];
      real f27_TS   = (D27.f[DIR_0PM])[kbn  ];
      real f27_BN   = (D27.f[DIR_0MP])[kts  ];
      real f27_ZERO = (D27.f[DIR_000])[kzero];
      real f27_BSW  = (D27.f[DIR_PPP])[ktne ];
      real f27_BNE  = (D27.f[DIR_MMP])[ktsw ];
      real f27_BNW  = (D27.f[DIR_PMP])[ktse ];
      real f27_BSE  = (D27.f[DIR_MPP])[ktnw ];
      real f27_TSW  = (D27.f[DIR_PPM])[kbne ];
      real f27_TNE  = (D27.f[DIR_MMM])[kbsw ];
      real f27_TNW  = (D27.f[DIR_PMM])[kbse ];
      real f27_TSE  = (D27.f[DIR_MPM])[kbnw ];
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
         D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[DIR_000])[k]=0.1f;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[  ke   ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M00])[kw  ]=(c2o1*feqW27_W  -(f27_E  *(q*omegaD-c1o1)-omegaD*feq27_E  *(q-c1o1))/(omegaD-c1o1)+f27_W  *q)/(q+c1o1);
      q = q_dirW[  kw   ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P00])[ke  ]=(c2o1*feqW27_E  -(f27_W  *(q*omegaD-c1o1)-omegaD*feq27_W  *(q-c1o1))/(omegaD-c1o1)+f27_E  *q)/(q+c1o1);
      q = q_dirN[  kn   ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0M0])[ks  ]=(c2o1*feqW27_S  -(f27_N  *(q*omegaD-c1o1)-omegaD*feq27_N  *(q-c1o1))/(omegaD-c1o1)+f27_S  *q)/(q+c1o1);
      q = q_dirS[  ks   ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0P0])[kn  ]=(c2o1*feqW27_N  -(f27_S  *(q*omegaD-c1o1)-omegaD*feq27_S  *(q-c1o1))/(omegaD-c1o1)+f27_N  *q)/(q+c1o1);
      q = q_dirT[  kt   ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00M])[kb  ]=(c2o1*feqW27_B  -(f27_T  *(q*omegaD-c1o1)-omegaD*feq27_T  *(q-c1o1))/(omegaD-c1o1)+f27_B  *q)/(q+c1o1);
      q = q_dirB[  kb   ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00P])[kt  ]=(c2o1*feqW27_T  -(f27_B  *(q*omegaD-c1o1)-omegaD*feq27_B  *(q-c1o1))/(omegaD-c1o1)+f27_T  *q)/(q+c1o1);
      q = q_dirNE[ kne  ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MM0])[ksw ]=(c2o1*feqW27_SW -(f27_NE *(q*omegaD-c1o1)-omegaD*feq27_NE *(q-c1o1))/(omegaD-c1o1)+f27_SW *q)/(q+c1o1);
      q = q_dirSW[ ksw  ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PP0])[kne ]=(c2o1*feqW27_NE -(f27_SW *(q*omegaD-c1o1)-omegaD*feq27_SW *(q-c1o1))/(omegaD-c1o1)+f27_NE *q)/(q+c1o1);
      q = q_dirSE[ kse  ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MP0])[knw ]=(c2o1*feqW27_NW -(f27_SE *(q*omegaD-c1o1)-omegaD*feq27_SE *(q-c1o1))/(omegaD-c1o1)+f27_NW *q)/(q+c1o1);
      q = q_dirNW[ knw  ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PM0])[kse ]=(c2o1*feqW27_SE -(f27_NW *(q*omegaD-c1o1)-omegaD*feq27_NW *(q-c1o1))/(omegaD-c1o1)+f27_SE *q)/(q+c1o1);
      q = q_dirTE[ kte  ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0M])[kbw ]=(c2o1*feqW27_BW -(f27_TE *(q*omegaD-c1o1)-omegaD*feq27_TE *(q-c1o1))/(omegaD-c1o1)+f27_BW *q)/(q+c1o1);
      q = q_dirBW[ kbw  ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0P])[kte ]=(c2o1*feqW27_TE -(f27_BW *(q*omegaD-c1o1)-omegaD*feq27_BW *(q-c1o1))/(omegaD-c1o1)+f27_TE *q)/(q+c1o1);
      q = q_dirBE[ kbe  ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0P])[ktw ]=(c2o1*feqW27_TW -(f27_BE *(q*omegaD-c1o1)-omegaD*feq27_BE *(q-c1o1))/(omegaD-c1o1)+f27_TW *q)/(q+c1o1);
      q = q_dirTW[ ktw  ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0M])[kbe ]=(c2o1*feqW27_BE -(f27_TW *(q*omegaD-c1o1)-omegaD*feq27_TW *(q-c1o1))/(omegaD-c1o1)+f27_BE *q)/(q+c1o1);
      q = q_dirTN[ ktn  ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MM])[kbs ]=(c2o1*feqW27_BS -(f27_TN *(q*omegaD-c1o1)-omegaD*feq27_TN *(q-c1o1))/(omegaD-c1o1)+f27_BS *q)/(q+c1o1);
      q = q_dirBS[ kbs  ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PP])[ktn ]=(c2o1*feqW27_TN -(f27_BS *(q*omegaD-c1o1)-omegaD*feq27_BS *(q-c1o1))/(omegaD-c1o1)+f27_TN *q)/(q+c1o1);
      q = q_dirBN[ kbn  ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MP])[kts ]=(c2o1*feqW27_TS -(f27_BN *(q*omegaD-c1o1)-omegaD*feq27_BN *(q-c1o1))/(omegaD-c1o1)+f27_TS *q)/(q+c1o1);
      q = q_dirTS[ kts  ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PM])[kbn ]=(c2o1*feqW27_BN -(f27_TS *(q*omegaD-c1o1)-omegaD*feq27_TS *(q-c1o1))/(omegaD-c1o1)+f27_BN *q)/(q+c1o1);
      q = q_dirTNE[ktne ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMM])[kbsw]=(c2o1*feqW27_BSW-(f27_TNE*(q*omegaD-c1o1)-omegaD*feq27_TNE*(q-c1o1))/(omegaD-c1o1)+f27_BSW*q)/(q+c1o1);
      q = q_dirBSW[kbsw ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPP])[ktne]=(c2o1*feqW27_TNE-(f27_BSW*(q*omegaD-c1o1)-omegaD*feq27_BSW*(q-c1o1))/(omegaD-c1o1)+f27_TNE*q)/(q+c1o1);
      q = q_dirBNE[kbne ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMP])[ktsw]=(c2o1*feqW27_TSW-(f27_BNE*(q*omegaD-c1o1)-omegaD*feq27_BNE*(q-c1o1))/(omegaD-c1o1)+f27_TSW*q)/(q+c1o1);
      q = q_dirTSW[ktsw ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPM])[kbne]=(c2o1*feqW27_BNE-(f27_TSW*(q*omegaD-c1o1)-omegaD*feq27_TSW*(q-c1o1))/(omegaD-c1o1)+f27_BNE*q)/(q+c1o1);
      q = q_dirTSE[ktse ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPM])[kbnw]=(c2o1*feqW27_BNW-(f27_TSE*(q*omegaD-c1o1)-omegaD*feq27_TSE*(q-c1o1))/(omegaD-c1o1)+f27_BNW*q)/(q+c1o1);
      q = q_dirBNW[kbnw ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMP])[ktse]=(c2o1*feqW27_TSE-(f27_BNW*(q*omegaD-c1o1)-omegaD*feq27_BNW*(q-c1o1))/(omegaD-c1o1)+f27_TSE*q)/(q+c1o1);
      q = q_dirBSE[kbse ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPP])[ktnw]=(c2o1*feqW27_TNW-(f27_BSE*(q*omegaD-c1o1)-omegaD*feq27_BSE*(q-c1o1))/(omegaD-c1o1)+f27_TNW*q)/(q+c1o1);
      q = q_dirTNW[ktnw ]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMM])[kbse]=(c2o1*feqW27_BSE-(f27_TNW*(q*omegaD-c1o1)-omegaD*feq27_TNW*(q-c1o1))/(omegaD-c1o1)+f27_BSE*q)/(q+c1o1);
      //q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[DIR_M00])[kw  ]=(two*feqW27_W  -(f27_E  *(q*omegaD-one)-omegaD*feq27_E  *(q-one))/(omegaD-one)+f27_W  *q)/(q+one);
      //q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[DIR_P00])[ke  ]=(two*feqW27_E  -(f27_W  *(q*omegaD-one)-omegaD*feq27_W  *(q-one))/(omegaD-one)+f27_E  *q)/(q+one);
      //q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[DIR_0M0])[ks  ]=(two*feqW27_S  -(f27_N  *(q*omegaD-one)-omegaD*feq27_N  *(q-one))/(omegaD-one)+f27_S  *q)/(q+one);
      //q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[DIR_0P0])[kn  ]=(two*feqW27_N  -(f27_S  *(q*omegaD-one)-omegaD*feq27_S  *(q-one))/(omegaD-one)+f27_N  *q)/(q+one);
      //q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[DIR_00M])[kb  ]=(two*feqW27_B  -(f27_T  *(q*omegaD-one)-omegaD*feq27_T  *(q-one))/(omegaD-one)+f27_B  *q)/(q+one);
      //q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[DIR_00P])[kt  ]=(two*feqW27_T  -(f27_B  *(q*omegaD-one)-omegaD*feq27_B  *(q-one))/(omegaD-one)+f27_T  *q)/(q+one);
      //q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[DIR_MM0])[ksw ]=(two*feqW27_SW -(f27_NE *(q*omegaD-one)-omegaD*feq27_NE *(q-one))/(omegaD-one)+f27_SW *q)/(q+one);
      //q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[DIR_PP0])[kne ]=(two*feqW27_NE -(f27_SW *(q*omegaD-one)-omegaD*feq27_SW *(q-one))/(omegaD-one)+f27_NE *q)/(q+one);
      //q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[DIR_MP0])[knw ]=(two*feqW27_NW -(f27_SE *(q*omegaD-one)-omegaD*feq27_SE *(q-one))/(omegaD-one)+f27_NW *q)/(q+one);
      //q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[DIR_PM0])[kse ]=(two*feqW27_SE -(f27_NW *(q*omegaD-one)-omegaD*feq27_NW *(q-one))/(omegaD-one)+f27_SE *q)/(q+one);
      //q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[DIR_M0M])[kbw ]=(two*feqW27_BW -(f27_TE *(q*omegaD-one)-omegaD*feq27_TE *(q-one))/(omegaD-one)+f27_BW *q)/(q+one);
      //q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[DIR_P0P])[kte ]=(two*feqW27_TE -(f27_BW *(q*omegaD-one)-omegaD*feq27_BW *(q-one))/(omegaD-one)+f27_TE *q)/(q+one);
      //q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[DIR_M0P])[ktw ]=(two*feqW27_TW -(f27_BE *(q*omegaD-one)-omegaD*feq27_BE *(q-one))/(omegaD-one)+f27_TW *q)/(q+one);
      //q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[DIR_P0M])[kbe ]=(two*feqW27_BE -(f27_TW *(q*omegaD-one)-omegaD*feq27_TW *(q-one))/(omegaD-one)+f27_BE *q)/(q+one);
      //q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[DIR_0MM])[kbs ]=(two*feqW27_BS -(f27_TN *(q*omegaD-one)-omegaD*feq27_TN *(q-one))/(omegaD-one)+f27_BS *q)/(q+one);
      //q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[DIR_0PP])[ktn ]=(two*feqW27_TN -(f27_BS *(q*omegaD-one)-omegaD*feq27_BS *(q-one))/(omegaD-one)+f27_TN *q)/(q+one);
      //q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[DIR_0MP])[kts ]=(two*feqW27_TS -(f27_BN *(q*omegaD-one)-omegaD*feq27_BN *(q-one))/(omegaD-one)+f27_TS *q)/(q+one);
      //q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[DIR_0PM])[kbn ]=(two*feqW27_BN -(f27_TS *(q*omegaD-one)-omegaD*feq27_TS *(q-one))/(omegaD-one)+f27_BN *q)/(q+one);
      //q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[DIR_MMM])[kbsw]=(two*feqW27_BSW-(f27_TNE*(q*omegaD-one)-omegaD*feq27_TNE*(q-one))/(omegaD-one)+f27_BSW*q)/(q+one);
      //q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[DIR_PPP])[ktne]=(two*feqW27_TNE-(f27_BSW*(q*omegaD-one)-omegaD*feq27_BSW*(q-one))/(omegaD-one)+f27_TNE*q)/(q+one);
      //q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[DIR_MMP])[ktsw]=(two*feqW27_TSW-(f27_BNE*(q*omegaD-one)-omegaD*feq27_BNE*(q-one))/(omegaD-one)+f27_TSW*q)/(q+one);
      //q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[DIR_PPM])[kbne]=(two*feqW27_BNE-(f27_TSW*(q*omegaD-one)-omegaD*feq27_TSW*(q-one))/(omegaD-one)+f27_BNE*q)/(q+one);
      //q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[DIR_MPM])[kbnw]=(two*feqW27_BNW-(f27_TSE*(q*omegaD-one)-omegaD*feq27_TSE*(q-one))/(omegaD-one)+f27_BNW*q)/(q+one);
      //q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[DIR_PMP])[ktse]=(two*feqW27_TSE-(f27_BNW*(q*omegaD-one)-omegaD*feq27_BNW*(q-one))/(omegaD-one)+f27_TSE*q)/(q+one);
      //q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[DIR_MPP])[ktnw]=(two*feqW27_TNW-(f27_BSE*(q*omegaD-one)-omegaD*feq27_BSE*(q-one))/(omegaD-one)+f27_TNW*q)/(q+one);
      //q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[DIR_PMM])[kbse]=(two*feqW27_BSE-(f27_TNW*(q*omegaD-one)-omegaD*feq27_TNW*(q-one))/(omegaD-one)+f27_BSE*q)/(q+one);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








































//////////////////////////////////////////////////////////////////////////////
__global__ void QADBB27( real* DD, 
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
                                   unsigned long long numberOfLBnodes, 
                                   bool isEvenTimestep)
{
   //Distributions27 D;
   //if (isEvenTimestep==true)
   //{
   //   D.f[DIR_P00] = &DD[DIR_P00 * size_Mat];
   //   D.f[DIR_M00] = &DD[DIR_M00 * size_Mat];
   //   D.f[DIR_0P0] = &DD[DIR_0P0 * size_Mat];
   //   D.f[DIR_0M0] = &DD[DIR_0M0 * size_Mat];
   //   D.f[DIR_00P] = &DD[DIR_00P * size_Mat];
   //   D.f[DIR_00M] = &DD[DIR_00M * size_Mat];
   //   D.f[DIR_PP0] = &DD[DIR_PP0 * size_Mat];
   //   D.f[DIR_MM0] = &DD[DIR_MM0 * size_Mat];
   //   D.f[DIR_PM0] = &DD[DIR_PM0 * size_Mat];
   //   D.f[DIR_MP0] = &DD[DIR_MP0 * size_Mat];
   //   D.f[DIR_P0P] = &DD[DIR_P0P * size_Mat];
   //   D.f[DIR_M0M] = &DD[DIR_M0M * size_Mat];
   //   D.f[DIR_P0M] = &DD[DIR_P0M * size_Mat];
   //   D.f[DIR_M0P] = &DD[DIR_M0P * size_Mat];
   //   D.f[DIR_0PP] = &DD[DIR_0PP * size_Mat];
   //   D.f[DIR_0MM] = &DD[DIR_0MM * size_Mat];
   //   D.f[DIR_0PM] = &DD[DIR_0PM * size_Mat];
   //   D.f[DIR_0MP] = &DD[DIR_0MP * size_Mat];
   //   D.f[DIR_000] = &DD[DIR_000 * size_Mat];
   //   D.f[DIR_PPP] = &DD[DIR_PPP * size_Mat];
   //   D.f[DIR_MMP] = &DD[DIR_MMP * size_Mat];
   //   D.f[DIR_PMP] = &DD[DIR_PMP * size_Mat];
   //   D.f[DIR_MPP] = &DD[DIR_MPP * size_Mat];
   //   D.f[DIR_PPM] = &DD[DIR_PPM * size_Mat];
   //   D.f[DIR_MMM] = &DD[DIR_MMM * size_Mat];
   //   D.f[DIR_PMM] = &DD[DIR_PMM * size_Mat];
   //   D.f[DIR_MPM] = &DD[DIR_MPM * size_Mat];
   //} 
   //else
   //{
   //   D.f[DIR_M00] = &DD[DIR_P00 * size_Mat];
   //   D.f[DIR_P00] = &DD[DIR_M00 * size_Mat];
   //   D.f[DIR_0M0] = &DD[DIR_0P0 * size_Mat];
   //   D.f[DIR_0P0] = &DD[DIR_0M0 * size_Mat];
   //   D.f[DIR_00M] = &DD[DIR_00P * size_Mat];
   //   D.f[DIR_00P] = &DD[DIR_00M * size_Mat];
   //   D.f[DIR_MM0] = &DD[DIR_PP0 * size_Mat];
   //   D.f[DIR_PP0] = &DD[DIR_MM0 * size_Mat];
   //   D.f[DIR_MP0] = &DD[DIR_PM0 * size_Mat];
   //   D.f[DIR_PM0] = &DD[DIR_MP0 * size_Mat];
   //   D.f[DIR_M0M] = &DD[DIR_P0P * size_Mat];
   //   D.f[DIR_P0P] = &DD[DIR_M0M * size_Mat];
   //   D.f[DIR_M0P] = &DD[DIR_P0M * size_Mat];
   //   D.f[DIR_P0M] = &DD[DIR_M0P * size_Mat];
   //   D.f[DIR_0MM] = &DD[DIR_0PP * size_Mat];
   //   D.f[DIR_0PP] = &DD[DIR_0MM * size_Mat];
   //   D.f[DIR_0MP] = &DD[DIR_0PM * size_Mat];
   //   D.f[DIR_0PM] = &DD[DIR_0MP * size_Mat];
   //   D.f[DIR_000] = &DD[DIR_000 * size_Mat];
   //   D.f[DIR_PPP] = &DD[DIR_MMM * size_Mat];
   //   D.f[DIR_MMP] = &DD[DIR_PPM * size_Mat];
   //   D.f[DIR_PMP] = &DD[DIR_MPM * size_Mat];
   //   D.f[DIR_MPP] = &DD[DIR_PMM * size_Mat];
   //   D.f[DIR_PPM] = &DD[DIR_MMP * size_Mat];
   //   D.f[DIR_MMM] = &DD[DIR_PPP * size_Mat];
   //   D.f[DIR_PMM] = &DD[DIR_MPP * size_Mat];
   //   D.f[DIR_MPM] = &DD[DIR_PMP * size_Mat];
   //}

   Distributions27 D27;
   if (isEvenTimestep==true)
   {
      D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
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
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
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
      //real f_W    = (D.f[DIR_P00])[ke   ];
      //real f_E    = (D.f[DIR_M00])[kw   ];
      //real f_S    = (D.f[DIR_0P0])[kn   ];
      //real f_N    = (D.f[DIR_0M0])[ks   ];
      //real f_B    = (D.f[DIR_00P])[kt   ];
      //real f_T    = (D.f[DIR_00M])[kb   ];
      //real f_SW   = (D.f[DIR_PP0])[kne  ];
      //real f_NE   = (D.f[DIR_MM0])[ksw  ];
      //real f_NW   = (D.f[DIR_PM0])[kse  ];
      //real f_SE   = (D.f[DIR_MP0])[knw  ];
      //real f_BW   = (D.f[DIR_P0P])[kte  ];
      //real f_TE   = (D.f[DIR_M0M])[kbw  ];
      //real f_TW   = (D.f[DIR_P0M])[kbe  ];
      //real f_BE   = (D.f[DIR_M0P])[ktw  ];
      //real f_BS   = (D.f[DIR_0PP])[ktn  ];
      //real f_TN   = (D.f[DIR_0MM])[kbs  ];
      //real f_TS   = (D.f[DIR_0PM])[kbn  ];
      //real f_BN   = (D.f[DIR_0MP])[kts  ];
      //real f_ZERO = (D.f[DIR_000])[kzero];
      //real f_BSW  = (D.f[DIR_PPP])[ktne ];
      //real f_BNE  = (D.f[DIR_MMP])[ktsw ];
      //real f_BNW  = (D.f[DIR_PMP])[ktse ];
      //real f_BSE  = (D.f[DIR_MPP])[ktnw ];
      //real f_TSW  = (D.f[DIR_PPM])[kbne ];
      //real f_TNE  = (D.f[DIR_MMM])[kbsw ];
      //real f_TNW  = (D.f[DIR_PMM])[kbse ];
      //real f_TSE  = (D.f[DIR_MPM])[kbnw ];
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
      real f27_W    = (D27.f[DIR_P00])[ke   ];
      real f27_E    = (D27.f[DIR_M00])[kw   ];
      real f27_S    = (D27.f[DIR_0P0])[kn   ];
      real f27_N    = (D27.f[DIR_0M0])[ks   ];
      real f27_B    = (D27.f[DIR_00P])[kt   ];
      real f27_T    = (D27.f[DIR_00M])[kb   ];
      real f27_SW   = (D27.f[DIR_PP0])[kne  ];
      real f27_NE   = (D27.f[DIR_MM0])[ksw  ];
      real f27_NW   = (D27.f[DIR_PM0])[kse  ];
      real f27_SE   = (D27.f[DIR_MP0])[knw  ];
      real f27_BW   = (D27.f[DIR_P0P])[kte  ];
      real f27_TE   = (D27.f[DIR_M0M])[kbw  ];
      real f27_TW   = (D27.f[DIR_P0M])[kbe  ];
      real f27_BE   = (D27.f[DIR_M0P])[ktw  ];
      real f27_BS   = (D27.f[DIR_0PP])[ktn  ];
      real f27_TN   = (D27.f[DIR_0MM])[kbs  ];
      real f27_TS   = (D27.f[DIR_0PM])[kbn  ];
      real f27_BN   = (D27.f[DIR_0MP])[kts  ];
      //real f27_ZERO = (D27.f[DIR_000])[kzero];
      real f27_BSW  = (D27.f[DIR_PPP])[ktne ];
      real f27_BNE  = (D27.f[DIR_MMP])[ktsw ];
      real f27_BNW  = (D27.f[DIR_PMP])[ktse ];
      real f27_BSE  = (D27.f[DIR_MPP])[ktnw ];
      real f27_TSW  = (D27.f[DIR_PPM])[kbne ];
      real f27_TNE  = (D27.f[DIR_MMM])[kbsw ];
      real f27_TNW  = (D27.f[DIR_PMM])[kbse ];
      real f27_TSE  = (D27.f[DIR_MPM])[kbnw ];
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
         D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[DIR_000])[k]=0.1f;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M00])[kw  ]=f27_E  ;
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P00])[ke  ]=f27_W  ;
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0M0])[ks  ]=f27_N  ;
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0P0])[kn  ]=f27_S  ;
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00M])[kb  ]=f27_T  ;
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00P])[kt  ]=f27_B  ;
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MM0])[ksw ]=f27_NE ;
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PP0])[kne ]=f27_SW ;
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MP0])[knw ]=f27_SE ;
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PM0])[kse ]=f27_NW ;
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0M])[kbw ]=f27_TE ;
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0P])[kte ]=f27_BW ;
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0P])[ktw ]=f27_BE ;
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0M])[kbe ]=f27_TW ;
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MM])[kbs ]=f27_TN ;
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PP])[ktn ]=f27_BS ;
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MP])[kts ]=f27_BN ;
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PM])[kbn ]=f27_TS ;
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMM])[kbsw]=f27_TNE;
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPP])[ktne]=f27_BSW;
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMP])[ktsw]=f27_BNE;
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPM])[kbne]=f27_TSW;
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPM])[kbnw]=f27_TSE;
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMP])[ktse]=f27_BNW;
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPP])[ktnw]=f27_BSE;
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMM])[kbse]=f27_TNW;
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
__global__ void QNoSlipADincomp7(
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
											 unsigned long long numberOfLBnodes, 
											 bool isEvenTimestep)
{
   //Distributions27 D;
   //if (isEvenTimestep==true)
   //{
   //   D.f[DIR_P00] = &DD[DIR_P00 * size_Mat];
   //   D.f[DIR_M00] = &DD[DIR_M00 * size_Mat];
   //   D.f[DIR_0P0] = &DD[DIR_0P0 * size_Mat];
   //   D.f[DIR_0M0] = &DD[DIR_0M0 * size_Mat];
   //   D.f[DIR_00P] = &DD[DIR_00P * size_Mat];
   //   D.f[DIR_00M] = &DD[DIR_00M * size_Mat];
   //   D.f[DIR_PP0] = &DD[DIR_PP0 * size_Mat];
   //   D.f[DIR_MM0] = &DD[DIR_MM0 * size_Mat];
   //   D.f[DIR_PM0] = &DD[DIR_PM0 * size_Mat];
   //   D.f[DIR_MP0] = &DD[DIR_MP0 * size_Mat];
   //   D.f[DIR_P0P] = &DD[DIR_P0P * size_Mat];
   //   D.f[DIR_M0M] = &DD[DIR_M0M * size_Mat];
   //   D.f[DIR_P0M] = &DD[DIR_P0M * size_Mat];
   //   D.f[DIR_M0P] = &DD[DIR_M0P * size_Mat];
   //   D.f[DIR_0PP] = &DD[DIR_0PP * size_Mat];
   //   D.f[DIR_0MM] = &DD[DIR_0MM * size_Mat];
   //   D.f[DIR_0PM] = &DD[DIR_0PM * size_Mat];
   //   D.f[DIR_0MP] = &DD[DIR_0MP * size_Mat];
   //   D.f[DIR_000] = &DD[DIR_000 * size_Mat];
   //   D.f[DIR_PPP] = &DD[DIR_PPP * size_Mat];
   //   D.f[DIR_MMP] = &DD[DIR_MMP * size_Mat];
   //   D.f[DIR_PMP] = &DD[DIR_PMP * size_Mat];
   //   D.f[DIR_MPP] = &DD[DIR_MPP * size_Mat];
   //   D.f[DIR_PPM] = &DD[DIR_PPM * size_Mat];
   //   D.f[DIR_MMM] = &DD[DIR_MMM * size_Mat];
   //   D.f[DIR_PMM] = &DD[DIR_PMM * size_Mat];
   //   D.f[DIR_MPM] = &DD[DIR_MPM * size_Mat];
   //} 
   //else
   //{
   //   D.f[DIR_M00] = &DD[DIR_P00 * size_Mat];
   //   D.f[DIR_P00] = &DD[DIR_M00 * size_Mat];
   //   D.f[DIR_0M0] = &DD[DIR_0P0 * size_Mat];
   //   D.f[DIR_0P0] = &DD[DIR_0M0 * size_Mat];
   //   D.f[DIR_00M] = &DD[DIR_00P * size_Mat];
   //   D.f[DIR_00P] = &DD[DIR_00M * size_Mat];
   //   D.f[DIR_MM0] = &DD[DIR_PP0 * size_Mat];
   //   D.f[DIR_PP0] = &DD[DIR_MM0 * size_Mat];
   //   D.f[DIR_MP0] = &DD[DIR_PM0 * size_Mat];
   //   D.f[DIR_PM0] = &DD[DIR_MP0 * size_Mat];
   //   D.f[DIR_M0M] = &DD[DIR_P0P * size_Mat];
   //   D.f[DIR_P0P] = &DD[DIR_M0M * size_Mat];
   //   D.f[DIR_M0P] = &DD[DIR_P0M * size_Mat];
   //   D.f[DIR_P0M] = &DD[DIR_M0P * size_Mat];
   //   D.f[DIR_0MM] = &DD[DIR_0PP * size_Mat];
   //   D.f[DIR_0PP] = &DD[DIR_0MM * size_Mat];
   //   D.f[DIR_0MP] = &DD[DIR_0PM * size_Mat];
   //   D.f[DIR_0PM] = &DD[DIR_0MP * size_Mat];
   //   D.f[DIR_000] = &DD[DIR_000 * size_Mat];
   //   D.f[DIR_PPP] = &DD[DIR_MMM * size_Mat];
   //   D.f[DIR_MMP] = &DD[DIR_PPM * size_Mat];
   //   D.f[DIR_PMP] = &DD[DIR_MPM * size_Mat];
   //   D.f[DIR_MPP] = &DD[DIR_PMM * size_Mat];
   //   D.f[DIR_PPM] = &DD[DIR_MMP * size_Mat];
   //   D.f[DIR_MMM] = &DD[DIR_PPP * size_Mat];
   //   D.f[DIR_PMM] = &DD[DIR_MPP * size_Mat];
   //   D.f[DIR_MPM] = &DD[DIR_PMP * size_Mat];
   //}

   Distributions7 D7;
   if (isEvenTimestep==true)
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[1] = &DD7[1*numberOfLBnodes];
      D7.f[2] = &DD7[2*numberOfLBnodes];
      D7.f[3] = &DD7[3*numberOfLBnodes];
      D7.f[4] = &DD7[4*numberOfLBnodes];
      D7.f[5] = &DD7[5*numberOfLBnodes];
      D7.f[6] = &DD7[6*numberOfLBnodes];
   }
   else
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[2] = &DD7[1*numberOfLBnodes];
      D7.f[1] = &DD7[2*numberOfLBnodes];
      D7.f[4] = &DD7[3*numberOfLBnodes];
      D7.f[3] = &DD7[4*numberOfLBnodes];
      D7.f[6] = &DD7[5*numberOfLBnodes];
      D7.f[5] = &DD7[6*numberOfLBnodes];
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

      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
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
      //real f_W    = (D.f[DIR_P00])[ke   ];
      //real f_E    = (D.f[DIR_M00])[kw   ];
      //real f_S    = (D.f[DIR_0P0])[kn   ];
      //real f_N    = (D.f[DIR_0M0])[ks   ];
      //real f_B    = (D.f[DIR_00P])[kt   ];
      //real f_T    = (D.f[DIR_00M])[kb   ];
      //real f_SW   = (D.f[DIR_PP0])[kne  ];
      //real f_NE   = (D.f[DIR_MM0])[ksw  ];
      //real f_NW   = (D.f[DIR_PM0])[kse  ];
      //real f_SE   = (D.f[DIR_MP0])[knw  ];
      //real f_BW   = (D.f[DIR_P0P])[kte  ];
      //real f_TE   = (D.f[DIR_M0M])[kbw  ];
      //real f_TW   = (D.f[DIR_P0M])[kbe  ];
      //real f_BE   = (D.f[DIR_M0P])[ktw  ];
      //real f_BS   = (D.f[DIR_0PP])[ktn  ];
      //real f_TN   = (D.f[DIR_0MM])[kbs  ];
      //real f_TS   = (D.f[DIR_0PM])[kbn  ];
      //real f_BN   = (D.f[DIR_0MP])[kts  ];
      //real f_BSW  = (D.f[DIR_PPP])[ktne ];
      //real f_BNE  = (D.f[DIR_MMP])[ktsw ];
      //real f_BNW  = (D.f[DIR_PMP])[ktse ];
      //real f_BSE  = (D.f[DIR_MPP])[ktnw ];
      //real f_TSW  = (D.f[DIR_PPM])[kbne ];
      //real f_TNE  = (D.f[DIR_MMM])[kbsw ];
      //real f_TNW  = (D.f[DIR_PMM])[kbse ];
      //real f_TSE  = (D.f[DIR_MPM])[kbnw ];
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
         D7.f[0] = &DD7[0*numberOfLBnodes];
         D7.f[1] = &DD7[1*numberOfLBnodes];
         D7.f[2] = &DD7[2*numberOfLBnodes];
         D7.f[3] = &DD7[3*numberOfLBnodes];
         D7.f[4] = &DD7[4*numberOfLBnodes];
         D7.f[5] = &DD7[5*numberOfLBnodes];
         D7.f[6] = &DD7[6*numberOfLBnodes];
      }
      else
      {
         D7.f[0] = &DD7[0*numberOfLBnodes];
         D7.f[2] = &DD7[1*numberOfLBnodes];
         D7.f[1] = &DD7[2*numberOfLBnodes];
         D7.f[4] = &DD7[3*numberOfLBnodes];
         D7.f[3] = &DD7[4*numberOfLBnodes];
         D7.f[6] = &DD7[5*numberOfLBnodes];
         D7.f[5] = &DD7[6*numberOfLBnodes];
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
__global__ void QNoSlipADincomp27(
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
											 unsigned long long numberOfLBnodes, 
											 bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
   }

   Distributions27 D27;
   if (isEvenTimestep==true)
   {
      D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
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
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
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
      real f_W    = (D.f[DIR_P00])[ke   ];
      real f_E    = (D.f[DIR_M00])[kw   ];
      real f_S    = (D.f[DIR_0P0])[kn   ];
      real f_N    = (D.f[DIR_0M0])[ks   ];
      real f_B    = (D.f[DIR_00P])[kt   ];
      real f_T    = (D.f[DIR_00M])[kb   ];
      real f_SW   = (D.f[DIR_PP0])[kne  ];
      real f_NE   = (D.f[DIR_MM0])[ksw  ];
      real f_NW   = (D.f[DIR_PM0])[kse  ];
      real f_SE   = (D.f[DIR_MP0])[knw  ];
      real f_BW   = (D.f[DIR_P0P])[kte  ];
      real f_TE   = (D.f[DIR_M0M])[kbw  ];
      real f_TW   = (D.f[DIR_P0M])[kbe  ];
      real f_BE   = (D.f[DIR_M0P])[ktw  ];
      real f_BS   = (D.f[DIR_0PP])[ktn  ];
      real f_TN   = (D.f[DIR_0MM])[kbs  ];
      real f_TS   = (D.f[DIR_0PM])[kbn  ];
      real f_BN   = (D.f[DIR_0MP])[kts  ];
      //real f_ZERO = (D.f[DIR_000])[kzero];
      real f_BSW  = (D.f[DIR_PPP])[ktne ];
      real f_BNE  = (D.f[DIR_MMP])[ktsw ];
      real f_BNW  = (D.f[DIR_PMP])[ktse ];
      real f_BSE  = (D.f[DIR_MPP])[ktnw ];
      real f_TSW  = (D.f[DIR_PPM])[kbne ];
      real f_TNE  = (D.f[DIR_MMM])[kbsw ];
      real f_TNW  = (D.f[DIR_PMM])[kbse ];
      real f_TSE  = (D.f[DIR_MPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1 =  ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      real vx2 =  ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      real vx3 =  ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      real f27_W    = (D27.f[DIR_P00])[ke   ];
      real f27_E    = (D27.f[DIR_M00])[kw   ];
      real f27_S    = (D27.f[DIR_0P0])[kn   ];
      real f27_N    = (D27.f[DIR_0M0])[ks   ];
      real f27_B    = (D27.f[DIR_00P])[kt   ];
      real f27_T    = (D27.f[DIR_00M])[kb   ];
      real f27_SW   = (D27.f[DIR_PP0])[kne  ];
      real f27_NE   = (D27.f[DIR_MM0])[ksw  ];
      real f27_NW   = (D27.f[DIR_PM0])[kse  ];
      real f27_SE   = (D27.f[DIR_MP0])[knw  ];
      real f27_BW   = (D27.f[DIR_P0P])[kte  ];
      real f27_TE   = (D27.f[DIR_M0M])[kbw  ];
      real f27_TW   = (D27.f[DIR_P0M])[kbe  ];
      real f27_BE   = (D27.f[DIR_M0P])[ktw  ];
      real f27_BS   = (D27.f[DIR_0PP])[ktn  ];
      real f27_TN   = (D27.f[DIR_0MM])[kbs  ];
      real f27_TS   = (D27.f[DIR_0PM])[kbn  ];
      real f27_BN   = (D27.f[DIR_0MP])[kts  ];
      real f27_ZERO = (D27.f[DIR_000])[kzero];
      real f27_BSW  = (D27.f[DIR_PPP])[ktne ];
      real f27_BNE  = (D27.f[DIR_MMP])[ktsw ];
      real f27_BNW  = (D27.f[DIR_PMP])[ktse ];
      real f27_BSE  = (D27.f[DIR_MPP])[ktnw ];
      real f27_TSW  = (D27.f[DIR_PPM])[kbne ];
      real f27_TNE  = (D27.f[DIR_MMM])[kbsw ];
      real f27_TNW  = (D27.f[DIR_PMM])[kbse ];
      real f27_TSE  = (D27.f[DIR_MPM])[kbnw ];
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
         D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[DIR_000])[k]=0.1f;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real q;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M00])[kw  ]=(c2o1*feqW27_W  -(f27_E  *(q*omegaD-c1o1)-omegaD*feq27_E  *(q-c1o1))/(omegaD-c1o1)+f27_W  *q)/(q+c1o1);
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P00])[ke  ]=(c2o1*feqW27_E  -(f27_W  *(q*omegaD-c1o1)-omegaD*feq27_W  *(q-c1o1))/(omegaD-c1o1)+f27_E  *q)/(q+c1o1);
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0M0])[ks  ]=(c2o1*feqW27_S  -(f27_N  *(q*omegaD-c1o1)-omegaD*feq27_N  *(q-c1o1))/(omegaD-c1o1)+f27_S  *q)/(q+c1o1);
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0P0])[kn  ]=(c2o1*feqW27_N  -(f27_S  *(q*omegaD-c1o1)-omegaD*feq27_S  *(q-c1o1))/(omegaD-c1o1)+f27_N  *q)/(q+c1o1);
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00M])[kb  ]=(c2o1*feqW27_B  -(f27_T  *(q*omegaD-c1o1)-omegaD*feq27_T  *(q-c1o1))/(omegaD-c1o1)+f27_B  *q)/(q+c1o1);
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00P])[kt  ]=(c2o1*feqW27_T  -(f27_B  *(q*omegaD-c1o1)-omegaD*feq27_B  *(q-c1o1))/(omegaD-c1o1)+f27_T  *q)/(q+c1o1);
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MM0])[ksw ]=(c2o1*feqW27_SW -(f27_NE *(q*omegaD-c1o1)-omegaD*feq27_NE *(q-c1o1))/(omegaD-c1o1)+f27_SW *q)/(q+c1o1);
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PP0])[kne ]=(c2o1*feqW27_NE -(f27_SW *(q*omegaD-c1o1)-omegaD*feq27_SW *(q-c1o1))/(omegaD-c1o1)+f27_NE *q)/(q+c1o1);
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MP0])[knw ]=(c2o1*feqW27_NW -(f27_SE *(q*omegaD-c1o1)-omegaD*feq27_SE *(q-c1o1))/(omegaD-c1o1)+f27_NW *q)/(q+c1o1);
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PM0])[kse ]=(c2o1*feqW27_SE -(f27_NW *(q*omegaD-c1o1)-omegaD*feq27_NW *(q-c1o1))/(omegaD-c1o1)+f27_SE *q)/(q+c1o1);
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0M])[kbw ]=(c2o1*feqW27_BW -(f27_TE *(q*omegaD-c1o1)-omegaD*feq27_TE *(q-c1o1))/(omegaD-c1o1)+f27_BW *q)/(q+c1o1);
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0P])[kte ]=(c2o1*feqW27_TE -(f27_BW *(q*omegaD-c1o1)-omegaD*feq27_BW *(q-c1o1))/(omegaD-c1o1)+f27_TE *q)/(q+c1o1);
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0P])[ktw ]=(c2o1*feqW27_TW -(f27_BE *(q*omegaD-c1o1)-omegaD*feq27_BE *(q-c1o1))/(omegaD-c1o1)+f27_TW *q)/(q+c1o1);
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0M])[kbe ]=(c2o1*feqW27_BE -(f27_TW *(q*omegaD-c1o1)-omegaD*feq27_TW *(q-c1o1))/(omegaD-c1o1)+f27_BE *q)/(q+c1o1);
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MM])[kbs ]=(c2o1*feqW27_BS -(f27_TN *(q*omegaD-c1o1)-omegaD*feq27_TN *(q-c1o1))/(omegaD-c1o1)+f27_BS *q)/(q+c1o1);
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PP])[ktn ]=(c2o1*feqW27_TN -(f27_BS *(q*omegaD-c1o1)-omegaD*feq27_BS *(q-c1o1))/(omegaD-c1o1)+f27_TN *q)/(q+c1o1);
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MP])[kts ]=(c2o1*feqW27_TS -(f27_BN *(q*omegaD-c1o1)-omegaD*feq27_BN *(q-c1o1))/(omegaD-c1o1)+f27_TS *q)/(q+c1o1);
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PM])[kbn ]=(c2o1*feqW27_BN -(f27_TS *(q*omegaD-c1o1)-omegaD*feq27_TS *(q-c1o1))/(omegaD-c1o1)+f27_BN *q)/(q+c1o1);
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMM])[kbsw]=(c2o1*feqW27_BSW-(f27_TNE*(q*omegaD-c1o1)-omegaD*feq27_TNE*(q-c1o1))/(omegaD-c1o1)+f27_BSW*q)/(q+c1o1);
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPP])[ktne]=(c2o1*feqW27_TNE-(f27_BSW*(q*omegaD-c1o1)-omegaD*feq27_BSW*(q-c1o1))/(omegaD-c1o1)+f27_TNE*q)/(q+c1o1);
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMP])[ktsw]=(c2o1*feqW27_TSW-(f27_BNE*(q*omegaD-c1o1)-omegaD*feq27_BNE*(q-c1o1))/(omegaD-c1o1)+f27_TSW*q)/(q+c1o1);
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPM])[kbne]=(c2o1*feqW27_BNE-(f27_TSW*(q*omegaD-c1o1)-omegaD*feq27_TSW*(q-c1o1))/(omegaD-c1o1)+f27_BNE*q)/(q+c1o1);
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPM])[kbnw]=(c2o1*feqW27_BNW-(f27_TSE*(q*omegaD-c1o1)-omegaD*feq27_TSE*(q-c1o1))/(omegaD-c1o1)+f27_BNW*q)/(q+c1o1);
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMP])[ktse]=(c2o1*feqW27_TSE-(f27_BNW*(q*omegaD-c1o1)-omegaD*feq27_BNW*(q-c1o1))/(omegaD-c1o1)+f27_TSE*q)/(q+c1o1);
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPP])[ktnw]=(c2o1*feqW27_TNW-(f27_BSE*(q*omegaD-c1o1)-omegaD*feq27_BSE*(q-c1o1))/(omegaD-c1o1)+f27_TNW*q)/(q+c1o1);
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMM])[kbse]=(c2o1*feqW27_BSE-(f27_TNW*(q*omegaD-c1o1)-omegaD*feq27_TNW*(q-c1o1))/(omegaD-c1o1)+f27_BSE*q)/(q+c1o1);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QADVeloIncomp7(
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
											unsigned long long numberOfLBnodes, 
											bool isEvenTimestep)
{
   //Distributions27 D;
   //if (isEvenTimestep==true)
   //{
   //   D.f[DIR_P00] = &DD[DIR_P00 * size_Mat];
   //   D.f[DIR_M00] = &DD[DIR_M00 * size_Mat];
   //   D.f[DIR_0P0] = &DD[DIR_0P0 * size_Mat];
   //   D.f[DIR_0M0] = &DD[DIR_0M0 * size_Mat];
   //   D.f[DIR_00P] = &DD[DIR_00P * size_Mat];
   //   D.f[DIR_00M] = &DD[DIR_00M * size_Mat];
   //   D.f[DIR_PP0] = &DD[DIR_PP0 * size_Mat];
   //   D.f[DIR_MM0] = &DD[DIR_MM0 * size_Mat];
   //   D.f[DIR_PM0] = &DD[DIR_PM0 * size_Mat];
   //   D.f[DIR_MP0] = &DD[DIR_MP0 * size_Mat];
   //   D.f[DIR_P0P] = &DD[DIR_P0P * size_Mat];
   //   D.f[DIR_M0M] = &DD[DIR_M0M * size_Mat];
   //   D.f[DIR_P0M] = &DD[DIR_P0M * size_Mat];
   //   D.f[DIR_M0P] = &DD[DIR_M0P * size_Mat];
   //   D.f[DIR_0PP] = &DD[DIR_0PP * size_Mat];
   //   D.f[DIR_0MM] = &DD[DIR_0MM * size_Mat];
   //   D.f[DIR_0PM] = &DD[DIR_0PM * size_Mat];
   //   D.f[DIR_0MP] = &DD[DIR_0MP * size_Mat];
   //   D.f[DIR_000] = &DD[DIR_000 * size_Mat];
   //   D.f[DIR_PPP] = &DD[DIR_PPP * size_Mat];
   //   D.f[DIR_MMP] = &DD[DIR_MMP * size_Mat];
   //   D.f[DIR_PMP] = &DD[DIR_PMP * size_Mat];
   //   D.f[DIR_MPP] = &DD[DIR_MPP * size_Mat];
   //   D.f[DIR_PPM] = &DD[DIR_PPM * size_Mat];
   //   D.f[DIR_MMM] = &DD[DIR_MMM * size_Mat];
   //   D.f[DIR_PMM] = &DD[DIR_PMM * size_Mat];
   //   D.f[DIR_MPM] = &DD[DIR_MPM * size_Mat];
   //} 
   //else
   //{
   //   D.f[DIR_M00] = &DD[DIR_P00 * size_Mat];
   //   D.f[DIR_P00] = &DD[DIR_M00 * size_Mat];
   //   D.f[DIR_0M0] = &DD[DIR_0P0 * size_Mat];
   //   D.f[DIR_0P0] = &DD[DIR_0M0 * size_Mat];
   //   D.f[DIR_00M] = &DD[DIR_00P * size_Mat];
   //   D.f[DIR_00P] = &DD[DIR_00M * size_Mat];
   //   D.f[DIR_MM0] = &DD[DIR_PP0 * size_Mat];
   //   D.f[DIR_PP0] = &DD[DIR_MM0 * size_Mat];
   //   D.f[DIR_MP0] = &DD[DIR_PM0 * size_Mat];
   //   D.f[DIR_PM0] = &DD[DIR_MP0 * size_Mat];
   //   D.f[DIR_M0M] = &DD[DIR_P0P * size_Mat];
   //   D.f[DIR_P0P] = &DD[DIR_M0M * size_Mat];
   //   D.f[DIR_M0P] = &DD[DIR_P0M * size_Mat];
   //   D.f[DIR_P0M] = &DD[DIR_M0P * size_Mat];
   //   D.f[DIR_0MM] = &DD[DIR_0PP * size_Mat];
   //   D.f[DIR_0PP] = &DD[DIR_0MM * size_Mat];
   //   D.f[DIR_0MP] = &DD[DIR_0PM * size_Mat];
   //   D.f[DIR_0PM] = &DD[DIR_0MP * size_Mat];
   //   D.f[DIR_000] = &DD[DIR_000 * size_Mat];
   //   D.f[DIR_PPP] = &DD[DIR_MMM * size_Mat];
   //   D.f[DIR_MMP] = &DD[DIR_PPM * size_Mat];
   //   D.f[DIR_PMP] = &DD[DIR_MPM * size_Mat];
   //   D.f[DIR_MPP] = &DD[DIR_PMM * size_Mat];
   //   D.f[DIR_PPM] = &DD[DIR_MMP * size_Mat];
   //   D.f[DIR_MMM] = &DD[DIR_PPP * size_Mat];
   //   D.f[DIR_PMM] = &DD[DIR_MPP * size_Mat];
   //   D.f[DIR_MPM] = &DD[DIR_PMP * size_Mat];
   //}

   Distributions7 D7;
   if (isEvenTimestep==true)
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[1] = &DD7[1*numberOfLBnodes];
      D7.f[2] = &DD7[2*numberOfLBnodes];
      D7.f[3] = &DD7[3*numberOfLBnodes];
      D7.f[4] = &DD7[4*numberOfLBnodes];
      D7.f[5] = &DD7[5*numberOfLBnodes];
      D7.f[6] = &DD7[6*numberOfLBnodes];
   }
   else
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[2] = &DD7[1*numberOfLBnodes];
      D7.f[1] = &DD7[2*numberOfLBnodes];
      D7.f[4] = &DD7[3*numberOfLBnodes];
      D7.f[3] = &DD7[4*numberOfLBnodes];
      D7.f[6] = &DD7[5*numberOfLBnodes];
      D7.f[5] = &DD7[6*numberOfLBnodes];
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

      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
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
      //real f_W    = (D.f[DIR_P00])[ke   ];
      //real f_E    = (D.f[DIR_M00])[kw   ];
      //real f_S    = (D.f[DIR_0P0])[kn   ];
      //real f_N    = (D.f[DIR_0M0])[ks   ];
      //real f_B    = (D.f[DIR_00P])[kt   ];
      //real f_T    = (D.f[DIR_00M])[kb   ];
      //real f_SW   = (D.f[DIR_PP0])[kne  ];
      //real f_NE   = (D.f[DIR_MM0])[ksw  ];
      //real f_NW   = (D.f[DIR_PM0])[kse  ];
      //real f_SE   = (D.f[DIR_MP0])[knw  ];
      //real f_BW   = (D.f[DIR_P0P])[kte  ];
      //real f_TE   = (D.f[DIR_M0M])[kbw  ];
      //real f_TW   = (D.f[DIR_P0M])[kbe  ];
      //real f_BE   = (D.f[DIR_M0P])[ktw  ];
      //real f_BS   = (D.f[DIR_0PP])[ktn  ];
      //real f_TN   = (D.f[DIR_0MM])[kbs  ];
      //real f_TS   = (D.f[DIR_0PM])[kbn  ];
      //real f_BN   = (D.f[DIR_0MP])[kts  ];
      //real f_BSW  = (D.f[DIR_PPP])[ktne ];
      //real f_BNE  = (D.f[DIR_MMP])[ktsw ];
      //real f_BNW  = (D.f[DIR_PMP])[ktse ];
      //real f_BSE  = (D.f[DIR_MPP])[ktnw ];
      //real f_TSW  = (D.f[DIR_PPM])[kbne ];
      //real f_TNE  = (D.f[DIR_MMM])[kbsw ];
      //real f_TNW  = (D.f[DIR_PMM])[kbse ];
      //real f_TSE  = (D.f[DIR_MPM])[kbnw ];
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
         D7.f[0] = &DD7[0*numberOfLBnodes];
         D7.f[1] = &DD7[1*numberOfLBnodes];
         D7.f[2] = &DD7[2*numberOfLBnodes];
         D7.f[3] = &DD7[3*numberOfLBnodes];
         D7.f[4] = &DD7[4*numberOfLBnodes];
         D7.f[5] = &DD7[5*numberOfLBnodes];
         D7.f[6] = &DD7[6*numberOfLBnodes];
      }
      else
      {
         D7.f[0] = &DD7[0*numberOfLBnodes];
         D7.f[2] = &DD7[1*numberOfLBnodes];
         D7.f[1] = &DD7[2*numberOfLBnodes];
         D7.f[4] = &DD7[3*numberOfLBnodes];
         D7.f[3] = &DD7[4*numberOfLBnodes];
         D7.f[6] = &DD7[5*numberOfLBnodes];
         D7.f[5] = &DD7[6*numberOfLBnodes];
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
__global__ void QADVeloIncomp27(
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
											unsigned long long numberOfLBnodes, 
											bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
   }

   Distributions27 D27;
   if (isEvenTimestep==true)
   {
      D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
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
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
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
      real f_W    = (D.f[DIR_P00])[ke   ];
      real f_E    = (D.f[DIR_M00])[kw   ];
      real f_S    = (D.f[DIR_0P0])[kn   ];
      real f_N    = (D.f[DIR_0M0])[ks   ];
      real f_B    = (D.f[DIR_00P])[kt   ];
      real f_T    = (D.f[DIR_00M])[kb   ];
      real f_SW   = (D.f[DIR_PP0])[kne  ];
      real f_NE   = (D.f[DIR_MM0])[ksw  ];
      real f_NW   = (D.f[DIR_PM0])[kse  ];
      real f_SE   = (D.f[DIR_MP0])[knw  ];
      real f_BW   = (D.f[DIR_P0P])[kte  ];
      real f_TE   = (D.f[DIR_M0M])[kbw  ];
      real f_TW   = (D.f[DIR_P0M])[kbe  ];
      real f_BE   = (D.f[DIR_M0P])[ktw  ];
      real f_BS   = (D.f[DIR_0PP])[ktn  ];
      real f_TN   = (D.f[DIR_0MM])[kbs  ];
      real f_TS   = (D.f[DIR_0PM])[kbn  ];
      real f_BN   = (D.f[DIR_0MP])[kts  ];
      //real f_ZERO = (D.f[DIR_000])[kzero];
      real f_BSW  = (D.f[DIR_PPP])[ktne ];
      real f_BNE  = (D.f[DIR_MMP])[ktsw ];
      real f_BNW  = (D.f[DIR_PMP])[ktse ];
      real f_BSE  = (D.f[DIR_MPP])[ktnw ];
      real f_TSW  = (D.f[DIR_PPM])[kbne ];
      real f_TNE  = (D.f[DIR_MMM])[kbsw ];
      real f_TNW  = (D.f[DIR_PMM])[kbse ];
      real f_TSE  = (D.f[DIR_MPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      real vx2 = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      real vx3 = ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      //real f27_W    = (D27.f[DIR_P00])[ke   ];
      //real f27_E    = (D27.f[DIR_M00])[kw   ];
      //real f27_S    = (D27.f[DIR_0P0])[kn   ];
      //real f27_N    = (D27.f[DIR_0M0])[ks   ];
      //real f27_B    = (D27.f[DIR_00P])[kt   ];
      //real f27_T    = (D27.f[DIR_00M])[kb   ];
      //real f27_SW   = (D27.f[DIR_PP0])[kne  ];
      //real f27_NE   = (D27.f[DIR_MM0])[ksw  ];
      //real f27_NW   = (D27.f[DIR_PM0])[kse  ];
      //real f27_SE   = (D27.f[DIR_MP0])[knw  ];
      //real f27_BW   = (D27.f[DIR_P0P])[kte  ];
      //real f27_TE   = (D27.f[DIR_M0M])[kbw  ];
      //real f27_TW   = (D27.f[DIR_P0M])[kbe  ];
      //real f27_BE   = (D27.f[DIR_M0P])[ktw  ];
      //real f27_BS   = (D27.f[DIR_0PP])[ktn  ];
      //real f27_TN   = (D27.f[DIR_0MM])[kbs  ];
      //real f27_TS   = (D27.f[DIR_0PM])[kbn  ];
      //real f27_BN   = (D27.f[DIR_0MP])[kts  ];
      //real f27_ZERO = (D27.f[DIR_000])[kzero];
      //real f27_BSW  = (D27.f[DIR_PPP])[ktne ];
      //real f27_BNE  = (D27.f[DIR_MMP])[ktsw ];
      //real f27_BNW  = (D27.f[DIR_PMP])[ktse ];
      //real f27_BSE  = (D27.f[DIR_MPP])[ktnw ];
      //real f27_TSW  = (D27.f[DIR_PPM])[kbne ];
      //real f27_TNE  = (D27.f[DIR_MMM])[kbsw ];
      //real f27_TNW  = (D27.f[DIR_PMM])[kbse ];
      //real f27_TSE  = (D27.f[DIR_MPM])[kbnw ];
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
         D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[DIR_000])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real q;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M00])[kw  ]= -feqW27_W  + c2o1 * c2o27  * TempD;
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P00])[ke  ]= -feqW27_E  + c2o1 * c2o27  * TempD;
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0M0])[ks  ]= -feqW27_S  + c2o1 * c2o27  * TempD;
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0P0])[kn  ]= -feqW27_N  + c2o1 * c2o27  * TempD;
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00M])[kb  ]= -feqW27_B  + c2o1 * c2o27  * TempD;
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00P])[kt  ]= -feqW27_T  + c2o1 * c2o27  * TempD;
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MM0])[ksw ]= -feqW27_SW + c2o1 * c1o54  * TempD;
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PP0])[kne ]= -feqW27_NE + c2o1 * c1o54  * TempD;
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MP0])[knw ]= -feqW27_NW + c2o1 * c1o54  * TempD;
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PM0])[kse ]= -feqW27_SE + c2o1 * c1o54  * TempD;
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0M])[kbw ]= -feqW27_BW + c2o1 * c1o54  * TempD;
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0P])[kte ]= -feqW27_TE + c2o1 * c1o54  * TempD;
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0P])[ktw ]= -feqW27_TW + c2o1 * c1o54  * TempD;
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0M])[kbe ]= -feqW27_BE + c2o1 * c1o54  * TempD;
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MM])[kbs ]= -feqW27_BS + c2o1 * c1o54  * TempD;
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PP])[ktn ]= -feqW27_TN + c2o1 * c1o54  * TempD;
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MP])[kts ]= -feqW27_TS + c2o1 * c1o54  * TempD;
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PM])[kbn ]= -feqW27_BN + c2o1 * c1o54  * TempD;
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMM])[kbsw]= -feqW27_BSW+ c2o1 * c1o216 * TempD;
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPP])[ktne]= -feqW27_TNE+ c2o1 * c1o216 * TempD;
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMP])[ktsw]= -feqW27_TSW+ c2o1 * c1o216 * TempD;
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPM])[kbne]= -feqW27_BNE+ c2o1 * c1o216 * TempD;
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPM])[kbnw]= -feqW27_BNW+ c2o1 * c1o216 * TempD;
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMP])[ktse]= -feqW27_TSE+ c2o1 * c1o216 * TempD;
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPP])[ktnw]= -feqW27_TNW+ c2o1 * c1o216 * TempD;
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMM])[kbse]= -feqW27_BSE+ c2o1 * c1o216 * TempD;
      //q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[DIR_M00])[kw  ]=(two*feqW27_W  -(f27_E  *(q*omegaD-one)-omegaD*feq27_E  *(q-one))/(omegaD-one)+f27_W  *q)/(q+one);
      //q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[DIR_P00])[ke  ]=(two*feqW27_E  -(f27_W  *(q*omegaD-one)-omegaD*feq27_W  *(q-one))/(omegaD-one)+f27_E  *q)/(q+one);
      //q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[DIR_0M0])[ks  ]=(two*feqW27_S  -(f27_N  *(q*omegaD-one)-omegaD*feq27_N  *(q-one))/(omegaD-one)+f27_S  *q)/(q+one);
      //q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[DIR_0P0])[kn  ]=(two*feqW27_N  -(f27_S  *(q*omegaD-one)-omegaD*feq27_S  *(q-one))/(omegaD-one)+f27_N  *q)/(q+one);
      //q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[DIR_00M])[kb  ]=(two*feqW27_B  -(f27_T  *(q*omegaD-one)-omegaD*feq27_T  *(q-one))/(omegaD-one)+f27_B  *q)/(q+one);
      //q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[DIR_00P])[kt  ]=(two*feqW27_T  -(f27_B  *(q*omegaD-one)-omegaD*feq27_B  *(q-one))/(omegaD-one)+f27_T  *q)/(q+one);
      //q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[DIR_MM0])[ksw ]=(two*feqW27_SW -(f27_NE *(q*omegaD-one)-omegaD*feq27_NE *(q-one))/(omegaD-one)+f27_SW *q)/(q+one);
      //q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[DIR_PP0])[kne ]=(two*feqW27_NE -(f27_SW *(q*omegaD-one)-omegaD*feq27_SW *(q-one))/(omegaD-one)+f27_NE *q)/(q+one);
      //q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[DIR_MP0])[knw ]=(two*feqW27_NW -(f27_SE *(q*omegaD-one)-omegaD*feq27_SE *(q-one))/(omegaD-one)+f27_NW *q)/(q+one);
      //q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[DIR_PM0])[kse ]=(two*feqW27_SE -(f27_NW *(q*omegaD-one)-omegaD*feq27_NW *(q-one))/(omegaD-one)+f27_SE *q)/(q+one);
      //q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[DIR_M0M])[kbw ]=(two*feqW27_BW -(f27_TE *(q*omegaD-one)-omegaD*feq27_TE *(q-one))/(omegaD-one)+f27_BW *q)/(q+one);
      //q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[DIR_P0P])[kte ]=(two*feqW27_TE -(f27_BW *(q*omegaD-one)-omegaD*feq27_BW *(q-one))/(omegaD-one)+f27_TE *q)/(q+one);
      //q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[DIR_M0P])[ktw ]=(two*feqW27_TW -(f27_BE *(q*omegaD-one)-omegaD*feq27_BE *(q-one))/(omegaD-one)+f27_TW *q)/(q+one);
      //q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[DIR_P0M])[kbe ]=(two*feqW27_BE -(f27_TW *(q*omegaD-one)-omegaD*feq27_TW *(q-one))/(omegaD-one)+f27_BE *q)/(q+one);
      //q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[DIR_0MM])[kbs ]=(two*feqW27_BS -(f27_TN *(q*omegaD-one)-omegaD*feq27_TN *(q-one))/(omegaD-one)+f27_BS *q)/(q+one);
      //q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[DIR_0PP])[ktn ]=(two*feqW27_TN -(f27_BS *(q*omegaD-one)-omegaD*feq27_BS *(q-one))/(omegaD-one)+f27_TN *q)/(q+one);
      //q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[DIR_0MP])[kts ]=(two*feqW27_TS -(f27_BN *(q*omegaD-one)-omegaD*feq27_BN *(q-one))/(omegaD-one)+f27_TS *q)/(q+one);
      //q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[DIR_0PM])[kbn ]=(two*feqW27_BN -(f27_TS *(q*omegaD-one)-omegaD*feq27_TS *(q-one))/(omegaD-one)+f27_BN *q)/(q+one);
      //q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[DIR_MMM])[kbsw]=(two*feqW27_BSW-(f27_TNE*(q*omegaD-one)-omegaD*feq27_TNE*(q-one))/(omegaD-one)+f27_BSW*q)/(q+one);
      //q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[DIR_PPP])[ktne]=(two*feqW27_TNE-(f27_BSW*(q*omegaD-one)-omegaD*feq27_BSW*(q-one))/(omegaD-one)+f27_TNE*q)/(q+one);
      //q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[DIR_MMP])[ktsw]=(two*feqW27_TSW-(f27_BNE*(q*omegaD-one)-omegaD*feq27_BNE*(q-one))/(omegaD-one)+f27_TSW*q)/(q+one);
      //q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[DIR_PPM])[kbne]=(two*feqW27_BNE-(f27_TSW*(q*omegaD-one)-omegaD*feq27_TSW*(q-one))/(omegaD-one)+f27_BNE*q)/(q+one);
      //q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[DIR_MPM])[kbnw]=(two*feqW27_BNW-(f27_TSE*(q*omegaD-one)-omegaD*feq27_TSE*(q-one))/(omegaD-one)+f27_BNW*q)/(q+one);
      //q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[DIR_PMP])[ktse]=(two*feqW27_TSE-(f27_BNW*(q*omegaD-one)-omegaD*feq27_BNW*(q-one))/(omegaD-one)+f27_TSE*q)/(q+one);
      //q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[DIR_MPP])[ktnw]=(two*feqW27_TNW-(f27_BSE*(q*omegaD-one)-omegaD*feq27_BSE*(q-one))/(omegaD-one)+f27_TNW*q)/(q+one);
      //q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[DIR_PMM])[kbse]=(two*feqW27_BSE-(f27_TNW*(q*omegaD-one)-omegaD*feq27_TNW*(q-one))/(omegaD-one)+f27_BSE*q)/(q+one);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QADPressIncomp7( real* DD, 
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
										   unsigned long long numberOfLBnodes, 
										   bool isEvenTimestep)
{
  /* Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &DD[DIR_P00 * size_Mat];
      D.f[DIR_M00] = &DD[DIR_M00 * size_Mat];
      D.f[DIR_0P0] = &DD[DIR_0P0 * size_Mat];
      D.f[DIR_0M0] = &DD[DIR_0M0 * size_Mat];
      D.f[DIR_00P] = &DD[DIR_00P * size_Mat];
      D.f[DIR_00M] = &DD[DIR_00M * size_Mat];
      D.f[DIR_PP0] = &DD[DIR_PP0 * size_Mat];
      D.f[DIR_MM0] = &DD[DIR_MM0 * size_Mat];
      D.f[DIR_PM0] = &DD[DIR_PM0 * size_Mat];
      D.f[DIR_MP0] = &DD[DIR_MP0 * size_Mat];
      D.f[DIR_P0P] = &DD[DIR_P0P * size_Mat];
      D.f[DIR_M0M] = &DD[DIR_M0M * size_Mat];
      D.f[DIR_P0M] = &DD[DIR_P0M * size_Mat];
      D.f[DIR_M0P] = &DD[DIR_M0P * size_Mat];
      D.f[DIR_0PP] = &DD[DIR_0PP * size_Mat];
      D.f[DIR_0MM] = &DD[DIR_0MM * size_Mat];
      D.f[DIR_0PM] = &DD[DIR_0PM * size_Mat];
      D.f[DIR_0MP] = &DD[DIR_0MP * size_Mat];
      D.f[DIR_000] = &DD[DIR_000 * size_Mat];
      D.f[DIR_PPP] = &DD[DIR_PPP * size_Mat];
      D.f[DIR_MMP] = &DD[DIR_MMP * size_Mat];
      D.f[DIR_PMP] = &DD[DIR_PMP * size_Mat];
      D.f[DIR_MPP] = &DD[DIR_MPP * size_Mat];
      D.f[DIR_PPM] = &DD[DIR_PPM * size_Mat];
      D.f[DIR_MMM] = &DD[DIR_MMM * size_Mat];
      D.f[DIR_PMM] = &DD[DIR_PMM * size_Mat];
      D.f[DIR_MPM] = &DD[DIR_MPM * size_Mat];
   } 
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * size_Mat];
      D.f[DIR_P00] = &DD[DIR_M00 * size_Mat];
      D.f[DIR_0M0] = &DD[DIR_0P0 * size_Mat];
      D.f[DIR_0P0] = &DD[DIR_0M0 * size_Mat];
      D.f[DIR_00M] = &DD[DIR_00P * size_Mat];
      D.f[DIR_00P] = &DD[DIR_00M * size_Mat];
      D.f[DIR_MM0] = &DD[DIR_PP0 * size_Mat];
      D.f[DIR_PP0] = &DD[DIR_MM0 * size_Mat];
      D.f[DIR_MP0] = &DD[DIR_PM0 * size_Mat];
      D.f[DIR_PM0] = &DD[DIR_MP0 * size_Mat];
      D.f[DIR_M0M] = &DD[DIR_P0P * size_Mat];
      D.f[DIR_P0P] = &DD[DIR_M0M * size_Mat];
      D.f[DIR_M0P] = &DD[DIR_P0M * size_Mat];
      D.f[DIR_P0M] = &DD[DIR_M0P * size_Mat];
      D.f[DIR_0MM] = &DD[DIR_0PP * size_Mat];
      D.f[DIR_0PP] = &DD[DIR_0MM * size_Mat];
      D.f[DIR_0MP] = &DD[DIR_0PM * size_Mat];
      D.f[DIR_0PM] = &DD[DIR_0MP * size_Mat];
      D.f[DIR_000] = &DD[DIR_000 * size_Mat];
      D.f[DIR_PPP] = &DD[DIR_MMM * size_Mat];
      D.f[DIR_MMP] = &DD[DIR_PPM * size_Mat];
      D.f[DIR_PMP] = &DD[DIR_MPM * size_Mat];
      D.f[DIR_MPP] = &DD[DIR_PMM * size_Mat];
      D.f[DIR_PPM] = &DD[DIR_MMP * size_Mat];
      D.f[DIR_MMM] = &DD[DIR_PPP * size_Mat];
      D.f[DIR_PMM] = &DD[DIR_MPP * size_Mat];
      D.f[DIR_MPM] = &DD[DIR_PMP * size_Mat];
   }*/

   Distributions7 D7;
   if (isEvenTimestep==true)
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[1] = &DD7[1*numberOfLBnodes];
      D7.f[2] = &DD7[2*numberOfLBnodes];
      D7.f[3] = &DD7[3*numberOfLBnodes];
      D7.f[4] = &DD7[4*numberOfLBnodes];
      D7.f[5] = &DD7[5*numberOfLBnodes];
      D7.f[6] = &DD7[6*numberOfLBnodes];
   }
   else
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[2] = &DD7[1*numberOfLBnodes];
      D7.f[1] = &DD7[2*numberOfLBnodes];
      D7.f[4] = &DD7[3*numberOfLBnodes];
      D7.f[3] = &DD7[4*numberOfLBnodes];
      D7.f[6] = &DD7[5*numberOfLBnodes];
      D7.f[5] = &DD7[6*numberOfLBnodes];
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

      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
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

      f_W    = (D.f[DIR_P00])[ke   ];
      f_E    = (D.f[DIR_M00])[kw   ];
      f_S    = (D.f[DIR_0P0])[kn   ];
      f_N    = (D.f[DIR_0M0])[ks   ];
      f_B    = (D.f[DIR_00P])[kt   ];
      f_T    = (D.f[DIR_00M])[kb   ];
      f_SW   = (D.f[DIR_PP0])[kne  ];
      f_NE   = (D.f[DIR_MM0])[ksw  ];
      f_NW   = (D.f[DIR_PM0])[kse  ];
      f_SE   = (D.f[DIR_MP0])[knw  ];
      f_BW   = (D.f[DIR_P0P])[kte  ];
      f_TE   = (D.f[DIR_M0M])[kbw  ];
      f_TW   = (D.f[DIR_P0M])[kbe  ];
      f_BE   = (D.f[DIR_M0P])[ktw  ];
      f_BS   = (D.f[DIR_0PP])[ktn  ];
      f_TN   = (D.f[DIR_0MM])[kbs  ];
      f_TS   = (D.f[DIR_0PM])[kbn  ];
      f_BN   = (D.f[DIR_0MP])[kts  ];
      f_BSW  = (D.f[DIR_PPP])[ktne ];
      f_BNE  = (D.f[DIR_MMP])[ktsw ];
      f_BNW  = (D.f[DIR_PMP])[ktse ];
      f_BSE  = (D.f[DIR_MPP])[ktnw ];
      f_TSW  = (D.f[DIR_PPM])[kbne ];
      f_TNE  = (D.f[DIR_MMM])[kbsw ];
      f_TNW  = (D.f[DIR_PMM])[kbse ];
      f_TSE  = (D.f[DIR_MPM])[kbnw ];*/
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
         D7.f[0] = &DD7[0*numberOfLBnodes];
         D7.f[1] = &DD7[1*numberOfLBnodes];
         D7.f[2] = &DD7[2*numberOfLBnodes];
         D7.f[3] = &DD7[3*numberOfLBnodes];
         D7.f[4] = &DD7[4*numberOfLBnodes];
         D7.f[5] = &DD7[5*numberOfLBnodes];
         D7.f[6] = &DD7[6*numberOfLBnodes];
      }
      else
      {
         D7.f[0] = &DD7[0*numberOfLBnodes];
         D7.f[2] = &DD7[1*numberOfLBnodes];
         D7.f[1] = &DD7[2*numberOfLBnodes];
         D7.f[4] = &DD7[3*numberOfLBnodes];
         D7.f[3] = &DD7[4*numberOfLBnodes];
         D7.f[6] = &DD7[5*numberOfLBnodes];
         D7.f[5] = &DD7[6*numberOfLBnodes];
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
__global__ void QADPressIncomp27(
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
											   unsigned long long numberOfLBnodes, 
											   bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
   }

   Distributions27 D27;
   if (isEvenTimestep==true)
   {
      D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
      D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
      D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
      D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
      D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
      D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
      D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
      D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
      D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
      D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
      D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
      D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
      D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
      D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
      D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
      D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
      D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
      D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
      D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
      D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
      D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
      D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
      D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
      D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
      D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
      D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
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
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
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
      real f_W    = (D.f[DIR_P00])[ke   ];
      real f_E    = (D.f[DIR_M00])[kw   ];
      real f_S    = (D.f[DIR_0P0])[kn   ];
      real f_N    = (D.f[DIR_0M0])[ks   ];
      real f_B    = (D.f[DIR_00P])[kt   ];
      real f_T    = (D.f[DIR_00M])[kb   ];
      real f_SW   = (D.f[DIR_PP0])[kne  ];
      real f_NE   = (D.f[DIR_MM0])[ksw  ];
      real f_NW   = (D.f[DIR_PM0])[kse  ];
      real f_SE   = (D.f[DIR_MP0])[knw  ];
      real f_BW   = (D.f[DIR_P0P])[kte  ];
      real f_TE   = (D.f[DIR_M0M])[kbw  ];
      real f_TW   = (D.f[DIR_P0M])[kbe  ];
      real f_BE   = (D.f[DIR_M0P])[ktw  ];
      real f_BS   = (D.f[DIR_0PP])[ktn  ];
      real f_TN   = (D.f[DIR_0MM])[kbs  ];
      real f_TS   = (D.f[DIR_0PM])[kbn  ];
      real f_BN   = (D.f[DIR_0MP])[kts  ];
      //real f_ZERO = (D.f[DIR_000])[kzero];
      real f_BSW  = (D.f[DIR_PPP])[ktne ];
      real f_BNE  = (D.f[DIR_MMP])[ktsw ];
      real f_BNW  = (D.f[DIR_PMP])[ktse ];
      real f_BSE  = (D.f[DIR_MPP])[ktnw ];
      real f_TSW  = (D.f[DIR_PPM])[kbne ];
      real f_TNE  = (D.f[DIR_MMM])[kbsw ];
      real f_TNW  = (D.f[DIR_PMM])[kbse ];
      real f_TSE  = (D.f[DIR_MPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1      = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      real vx2      = ((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      real vx3      = ((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
      ////////////////////////////////////////////////////////////////////////////////
      //real f27_W    = (D27.f[DIR_P00])[ke   ];
      //real f27_E    = (D27.f[DIR_M00])[kw   ];
      //real f27_S    = (D27.f[DIR_0P0])[kn   ];
      //real f27_N    = (D27.f[DIR_0M0])[ks   ];
      //real f27_B    = (D27.f[DIR_00P])[kt   ];
      //real f27_T    = (D27.f[DIR_00M])[kb   ];
      //real f27_SW   = (D27.f[DIR_PP0])[kne  ];
      //real f27_NE   = (D27.f[DIR_MM0])[ksw  ];
      //real f27_NW   = (D27.f[DIR_PM0])[kse  ];
      //real f27_SE   = (D27.f[DIR_MP0])[knw  ];
      //real f27_BW   = (D27.f[DIR_P0P])[kte  ];
      //real f27_TE   = (D27.f[DIR_M0M])[kbw  ];
      //real f27_TW   = (D27.f[DIR_P0M])[kbe  ];
      //real f27_BE   = (D27.f[DIR_M0P])[ktw  ];
      //real f27_BS   = (D27.f[DIR_0PP])[ktn  ];
      //real f27_TN   = (D27.f[DIR_0MM])[kbs  ];
      //real f27_TS   = (D27.f[DIR_0PM])[kbn  ];
      //real f27_BN   = (D27.f[DIR_0MP])[kts  ];
      //real f27_ZERO = (D27.f[DIR_000])[kzero];
      //real f27_BSW  = (D27.f[DIR_PPP])[ktne ];
      //real f27_BNE  = (D27.f[DIR_MMP])[ktsw ];
      //real f27_BNW  = (D27.f[DIR_PMP])[ktse ];
      //real f27_BSE  = (D27.f[DIR_MPP])[ktnw ];
      //real f27_TSW  = (D27.f[DIR_PPM])[kbne ];
      //real f27_TNE  = (D27.f[DIR_MMM])[kbsw ];
      //real f27_TNW  = (D27.f[DIR_PMM])[kbse ];
      //real f27_TSE  = (D27.f[DIR_MPM])[kbnw ];
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
         D27.f[DIR_P00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_M00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_PMP * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         D27.f[DIR_M00] = &DD27[DIR_P00 * numberOfLBnodes];
         D27.f[DIR_P00] = &DD27[DIR_M00 * numberOfLBnodes];
         D27.f[DIR_0M0] = &DD27[DIR_0P0 * numberOfLBnodes];
         D27.f[DIR_0P0] = &DD27[DIR_0M0 * numberOfLBnodes];
         D27.f[DIR_00M] = &DD27[DIR_00P * numberOfLBnodes];
         D27.f[DIR_00P] = &DD27[DIR_00M * numberOfLBnodes];
         D27.f[DIR_MM0] = &DD27[DIR_PP0 * numberOfLBnodes];
         D27.f[DIR_PP0] = &DD27[DIR_MM0 * numberOfLBnodes];
         D27.f[DIR_MP0] = &DD27[DIR_PM0 * numberOfLBnodes];
         D27.f[DIR_PM0] = &DD27[DIR_MP0 * numberOfLBnodes];
         D27.f[DIR_M0M] = &DD27[DIR_P0P * numberOfLBnodes];
         D27.f[DIR_P0P] = &DD27[DIR_M0M * numberOfLBnodes];
         D27.f[DIR_M0P] = &DD27[DIR_P0M * numberOfLBnodes];
         D27.f[DIR_P0M] = &DD27[DIR_M0P * numberOfLBnodes];
         D27.f[DIR_0MM] = &DD27[DIR_0PP * numberOfLBnodes];
         D27.f[DIR_0PP] = &DD27[DIR_0MM * numberOfLBnodes];
         D27.f[DIR_0MP] = &DD27[DIR_0PM * numberOfLBnodes];
         D27.f[DIR_0PM] = &DD27[DIR_0MP * numberOfLBnodes];
         D27.f[DIR_000] = &DD27[DIR_000 * numberOfLBnodes];
         D27.f[DIR_PPP] = &DD27[DIR_MMM * numberOfLBnodes];
         D27.f[DIR_MMP] = &DD27[DIR_PPM * numberOfLBnodes];
         D27.f[DIR_PMP] = &DD27[DIR_MPM * numberOfLBnodes];
         D27.f[DIR_MPP] = &DD27[DIR_PMM * numberOfLBnodes];
         D27.f[DIR_PPM] = &DD27[DIR_MMP * numberOfLBnodes];
         D27.f[DIR_MMM] = &DD27[DIR_PPP * numberOfLBnodes];
         D27.f[DIR_PMM] = &DD27[DIR_MPP * numberOfLBnodes];
         D27.f[DIR_MPM] = &DD27[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[DIR_000])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real q;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M00])[kw  ]= -feqW27_W  + c2o1 * c2o27  * TempD;
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P00])[ke  ]= -feqW27_E  + c2o1 * c2o27  * TempD;
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0M0])[ks  ]= -feqW27_S  + c2o1 * c2o27  * TempD;
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0P0])[kn  ]= -feqW27_N  + c2o1 * c2o27  * TempD;
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00M])[kb  ]= -feqW27_B  + c2o1 * c2o27  * TempD;
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[DIR_00P])[kt  ]= -feqW27_T  + c2o1 * c2o27  * TempD;
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MM0])[ksw ]= -feqW27_SW + c2o1 * c1o54  * TempD;
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PP0])[kne ]= -feqW27_NE + c2o1 * c1o54  * TempD;
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MP0])[knw ]= -feqW27_NW + c2o1 * c1o54  * TempD;
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PM0])[kse ]= -feqW27_SE + c2o1 * c1o54  * TempD;
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0M])[kbw ]= -feqW27_BW + c2o1 * c1o54  * TempD;
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0P])[kte ]= -feqW27_TE + c2o1 * c1o54  * TempD;
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_M0P])[ktw ]= -feqW27_TW + c2o1 * c1o54  * TempD;
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_P0M])[kbe ]= -feqW27_BE + c2o1 * c1o54  * TempD;
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MM])[kbs ]= -feqW27_BS + c2o1 * c1o54  * TempD;
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PP])[ktn ]= -feqW27_TN + c2o1 * c1o54  * TempD;
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0MP])[kts ]= -feqW27_TS + c2o1 * c1o54  * TempD;
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[DIR_0PM])[kbn ]= -feqW27_BN + c2o1 * c1o54  * TempD;
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMM])[kbsw]= -feqW27_BSW+ c2o1 * c1o216 * TempD;
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPP])[ktne]= -feqW27_TNE+ c2o1 * c1o216 * TempD;
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MMP])[ktsw]= -feqW27_TSW+ c2o1 * c1o216 * TempD;
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PPM])[kbne]= -feqW27_BNE+ c2o1 * c1o216 * TempD;
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPM])[kbnw]= -feqW27_BNW+ c2o1 * c1o216 * TempD;
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMP])[ktse]= -feqW27_TSE+ c2o1 * c1o216 * TempD;
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_MPP])[ktnw]= -feqW27_TNW+ c2o1 * c1o216 * TempD;
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[DIR_PMM])[kbse]= -feqW27_BSE+ c2o1 * c1o216 * TempD;
      //q = q_dirE[k];   if (q>=zero && q<=one) (D27.f[DIR_M00])[kw  ]=(two*feqW27_W  -(f27_E  *(q*omegaD-one)-omegaD*feq27_E  *(q-one))/(omegaD-one)+f27_W  *q)/(q+one);
      //q = q_dirW[k];   if (q>=zero && q<=one) (D27.f[DIR_P00])[ke  ]=(two*feqW27_E  -(f27_W  *(q*omegaD-one)-omegaD*feq27_W  *(q-one))/(omegaD-one)+f27_E  *q)/(q+one);
      //q = q_dirN[k];   if (q>=zero && q<=one) (D27.f[DIR_0M0])[ks  ]=(two*feqW27_S  -(f27_N  *(q*omegaD-one)-omegaD*feq27_N  *(q-one))/(omegaD-one)+f27_S  *q)/(q+one);
      //q = q_dirS[k];   if (q>=zero && q<=one) (D27.f[DIR_0P0])[kn  ]=(two*feqW27_N  -(f27_S  *(q*omegaD-one)-omegaD*feq27_S  *(q-one))/(omegaD-one)+f27_N  *q)/(q+one);
      //q = q_dirT[k];   if (q>=zero && q<=one) (D27.f[DIR_00M])[kb  ]=(two*feqW27_B  -(f27_T  *(q*omegaD-one)-omegaD*feq27_T  *(q-one))/(omegaD-one)+f27_B  *q)/(q+one);
      //q = q_dirB[k];   if (q>=zero && q<=one) (D27.f[DIR_00P])[kt  ]=(two*feqW27_T  -(f27_B  *(q*omegaD-one)-omegaD*feq27_B  *(q-one))/(omegaD-one)+f27_T  *q)/(q+one);
      //q = q_dirNE[k];  if (q>=zero && q<=one) (D27.f[DIR_MM0])[ksw ]=(two*feqW27_SW -(f27_NE *(q*omegaD-one)-omegaD*feq27_NE *(q-one))/(omegaD-one)+f27_SW *q)/(q+one);
      //q = q_dirSW[k];  if (q>=zero && q<=one) (D27.f[DIR_PP0])[kne ]=(two*feqW27_NE -(f27_SW *(q*omegaD-one)-omegaD*feq27_SW *(q-one))/(omegaD-one)+f27_NE *q)/(q+one);
      //q = q_dirSE[k];  if (q>=zero && q<=one) (D27.f[DIR_MP0])[knw ]=(two*feqW27_NW -(f27_SE *(q*omegaD-one)-omegaD*feq27_SE *(q-one))/(omegaD-one)+f27_NW *q)/(q+one);
      //q = q_dirNW[k];  if (q>=zero && q<=one) (D27.f[DIR_PM0])[kse ]=(two*feqW27_SE -(f27_NW *(q*omegaD-one)-omegaD*feq27_NW *(q-one))/(omegaD-one)+f27_SE *q)/(q+one);
      //q = q_dirTE[k];  if (q>=zero && q<=one) (D27.f[DIR_M0M])[kbw ]=(two*feqW27_BW -(f27_TE *(q*omegaD-one)-omegaD*feq27_TE *(q-one))/(omegaD-one)+f27_BW *q)/(q+one);
      //q = q_dirBW[k];  if (q>=zero && q<=one) (D27.f[DIR_P0P])[kte ]=(two*feqW27_TE -(f27_BW *(q*omegaD-one)-omegaD*feq27_BW *(q-one))/(omegaD-one)+f27_TE *q)/(q+one);
      //q = q_dirBE[k];  if (q>=zero && q<=one) (D27.f[DIR_M0P])[ktw ]=(two*feqW27_TW -(f27_BE *(q*omegaD-one)-omegaD*feq27_BE *(q-one))/(omegaD-one)+f27_TW *q)/(q+one);
      //q = q_dirTW[k];  if (q>=zero && q<=one) (D27.f[DIR_P0M])[kbe ]=(two*feqW27_BE -(f27_TW *(q*omegaD-one)-omegaD*feq27_TW *(q-one))/(omegaD-one)+f27_BE *q)/(q+one);
      //q = q_dirTN[k];  if (q>=zero && q<=one) (D27.f[DIR_0MM])[kbs ]=(two*feqW27_BS -(f27_TN *(q*omegaD-one)-omegaD*feq27_TN *(q-one))/(omegaD-one)+f27_BS *q)/(q+one);
      //q = q_dirBS[k];  if (q>=zero && q<=one) (D27.f[DIR_0PP])[ktn ]=(two*feqW27_TN -(f27_BS *(q*omegaD-one)-omegaD*feq27_BS *(q-one))/(omegaD-one)+f27_TN *q)/(q+one);
      //q = q_dirBN[k];  if (q>=zero && q<=one) (D27.f[DIR_0MP])[kts ]=(two*feqW27_TS -(f27_BN *(q*omegaD-one)-omegaD*feq27_BN *(q-one))/(omegaD-one)+f27_TS *q)/(q+one);
      //q = q_dirTS[k];  if (q>=zero && q<=one) (D27.f[DIR_0PM])[kbn ]=(two*feqW27_BN -(f27_TS *(q*omegaD-one)-omegaD*feq27_TS *(q-one))/(omegaD-one)+f27_BN *q)/(q+one);
      //q = q_dirTNE[k]; if (q>=zero && q<=one) (D27.f[DIR_MMM])[kbsw]=(two*feqW27_BSW-(f27_TNE*(q*omegaD-one)-omegaD*feq27_TNE*(q-one))/(omegaD-one)+f27_BSW*q)/(q+one);
      //q = q_dirBSW[k]; if (q>=zero && q<=one) (D27.f[DIR_PPP])[ktne]=(two*feqW27_TNE-(f27_BSW*(q*omegaD-one)-omegaD*feq27_BSW*(q-one))/(omegaD-one)+f27_TNE*q)/(q+one);
      //q = q_dirBNE[k]; if (q>=zero && q<=one) (D27.f[DIR_MMP])[ktsw]=(two*feqW27_TSW-(f27_BNE*(q*omegaD-one)-omegaD*feq27_BNE*(q-one))/(omegaD-one)+f27_TSW*q)/(q+one);
      //q = q_dirTSW[k]; if (q>=zero && q<=one) (D27.f[DIR_PPM])[kbne]=(two*feqW27_BNE-(f27_TSW*(q*omegaD-one)-omegaD*feq27_TSW*(q-one))/(omegaD-one)+f27_BNE*q)/(q+one);
      //q = q_dirTSE[k]; if (q>=zero && q<=one) (D27.f[DIR_MPM])[kbnw]=(two*feqW27_BNW-(f27_TSE*(q*omegaD-one)-omegaD*feq27_TSE*(q-one))/(omegaD-one)+f27_BNW*q)/(q+one);
      //q = q_dirBNW[k]; if (q>=zero && q<=one) (D27.f[DIR_PMP])[ktse]=(two*feqW27_TSE-(f27_BNW*(q*omegaD-one)-omegaD*feq27_BNW*(q-one))/(omegaD-one)+f27_TSE*q)/(q+one);
      //q = q_dirBSE[k]; if (q>=zero && q<=one) (D27.f[DIR_MPP])[ktnw]=(two*feqW27_TNW-(f27_BSE*(q*omegaD-one)-omegaD*feq27_BSE*(q-one))/(omegaD-one)+f27_TNW*q)/(q+one);
      //q = q_dirTNW[k]; if (q>=zero && q<=one) (D27.f[DIR_PMM])[kbse]=(two*feqW27_BSE-(f27_TNW*(q*omegaD-one)-omegaD*feq27_TNW*(q-one))/(omegaD-one)+f27_BSE*q)/(q+one);
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
__global__ void AD_SlipVelDeviceComp(
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
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    Distributions27 D;
    if (isEvenTimestep)
    {
        D.f[DIR_P00] = &distributions[DIR_P00 * numberOfLBnodes];
        D.f[DIR_M00] = &distributions[DIR_M00 * numberOfLBnodes];
        D.f[DIR_0P0] = &distributions[DIR_0P0 * numberOfLBnodes];
        D.f[DIR_0M0] = &distributions[DIR_0M0 * numberOfLBnodes];
        D.f[DIR_00P] = &distributions[DIR_00P * numberOfLBnodes];
        D.f[DIR_00M] = &distributions[DIR_00M * numberOfLBnodes];
        D.f[DIR_PP0] = &distributions[DIR_PP0 * numberOfLBnodes];
        D.f[DIR_MM0] = &distributions[DIR_MM0 * numberOfLBnodes];
        D.f[DIR_PM0] = &distributions[DIR_PM0 * numberOfLBnodes];
        D.f[DIR_MP0] = &distributions[DIR_MP0 * numberOfLBnodes];
        D.f[DIR_P0P] = &distributions[DIR_P0P * numberOfLBnodes];
        D.f[DIR_M0M] = &distributions[DIR_M0M * numberOfLBnodes];
        D.f[DIR_P0M] = &distributions[DIR_P0M * numberOfLBnodes];
        D.f[DIR_M0P] = &distributions[DIR_M0P * numberOfLBnodes];
        D.f[DIR_0PP] = &distributions[DIR_0PP * numberOfLBnodes];
        D.f[DIR_0MM] = &distributions[DIR_0MM * numberOfLBnodes];
        D.f[DIR_0PM] = &distributions[DIR_0PM * numberOfLBnodes];
        D.f[DIR_0MP] = &distributions[DIR_0MP * numberOfLBnodes];
        D.f[DIR_000] = &distributions[DIR_000 * numberOfLBnodes];
        D.f[DIR_PPP] = &distributions[DIR_PPP * numberOfLBnodes];
        D.f[DIR_MMP] = &distributions[DIR_MMP * numberOfLBnodes];
        D.f[DIR_PMP] = &distributions[DIR_PMP * numberOfLBnodes];
        D.f[DIR_MPP] = &distributions[DIR_MPP * numberOfLBnodes];
        D.f[DIR_PPM] = &distributions[DIR_PPM * numberOfLBnodes];
        D.f[DIR_MMM] = &distributions[DIR_MMM * numberOfLBnodes];
        D.f[DIR_PMM] = &distributions[DIR_PMM * numberOfLBnodes];
        D.f[DIR_MPM] = &distributions[DIR_MPM * numberOfLBnodes];
    }
    else
    {
        D.f[DIR_M00] = &distributions[DIR_P00 * numberOfLBnodes];
        D.f[DIR_P00] = &distributions[DIR_M00 * numberOfLBnodes];
        D.f[DIR_0M0] = &distributions[DIR_0P0 * numberOfLBnodes];
        D.f[DIR_0P0] = &distributions[DIR_0M0 * numberOfLBnodes];
        D.f[DIR_00M] = &distributions[DIR_00P * numberOfLBnodes];
        D.f[DIR_00P] = &distributions[DIR_00M * numberOfLBnodes];
        D.f[DIR_MM0] = &distributions[DIR_PP0 * numberOfLBnodes];
        D.f[DIR_PP0] = &distributions[DIR_MM0 * numberOfLBnodes];
        D.f[DIR_MP0] = &distributions[DIR_PM0 * numberOfLBnodes];
        D.f[DIR_PM0] = &distributions[DIR_MP0 * numberOfLBnodes];
        D.f[DIR_M0M] = &distributions[DIR_P0P * numberOfLBnodes];
        D.f[DIR_P0P] = &distributions[DIR_M0M * numberOfLBnodes];
        D.f[DIR_M0P] = &distributions[DIR_P0M * numberOfLBnodes];
        D.f[DIR_P0M] = &distributions[DIR_M0P * numberOfLBnodes];
        D.f[DIR_0MM] = &distributions[DIR_0PP * numberOfLBnodes];
        D.f[DIR_0PP] = &distributions[DIR_0MM * numberOfLBnodes];
        D.f[DIR_0MP] = &distributions[DIR_0PM * numberOfLBnodes];
        D.f[DIR_0PM] = &distributions[DIR_0MP * numberOfLBnodes];
        D.f[DIR_000] = &distributions[DIR_000 * numberOfLBnodes];
        D.f[DIR_PPP] = &distributions[DIR_MMM * numberOfLBnodes];
        D.f[DIR_MMP] = &distributions[DIR_PPM * numberOfLBnodes];
        D.f[DIR_PMP] = &distributions[DIR_MPM * numberOfLBnodes];
        D.f[DIR_MPP] = &distributions[DIR_PMM * numberOfLBnodes];
        D.f[DIR_PPM] = &distributions[DIR_MMP * numberOfLBnodes];
        D.f[DIR_MMM] = &distributions[DIR_PPP * numberOfLBnodes];
        D.f[DIR_PMM] = &distributions[DIR_MPP * numberOfLBnodes];
        D.f[DIR_MPM] = &distributions[DIR_PMP * numberOfLBnodes];
    }
    ////////////////////////////////////////////////////////////////////////////////
    Distributions27 DAD;
    if (isEvenTimestep)
    {
        DAD.f[DIR_P00] = &distributionsAD[DIR_P00 * numberOfLBnodes];
        DAD.f[DIR_M00] = &distributionsAD[DIR_M00 * numberOfLBnodes];
        DAD.f[DIR_0P0] = &distributionsAD[DIR_0P0 * numberOfLBnodes];
        DAD.f[DIR_0M0] = &distributionsAD[DIR_0M0 * numberOfLBnodes];
        DAD.f[DIR_00P] = &distributionsAD[DIR_00P * numberOfLBnodes];
        DAD.f[DIR_00M] = &distributionsAD[DIR_00M * numberOfLBnodes];
        DAD.f[DIR_PP0] = &distributionsAD[DIR_PP0 * numberOfLBnodes];
        DAD.f[DIR_MM0] = &distributionsAD[DIR_MM0 * numberOfLBnodes];
        DAD.f[DIR_PM0] = &distributionsAD[DIR_PM0 * numberOfLBnodes];
        DAD.f[DIR_MP0] = &distributionsAD[DIR_MP0 * numberOfLBnodes];
        DAD.f[DIR_P0P] = &distributionsAD[DIR_P0P * numberOfLBnodes];
        DAD.f[DIR_M0M] = &distributionsAD[DIR_M0M * numberOfLBnodes];
        DAD.f[DIR_P0M] = &distributionsAD[DIR_P0M * numberOfLBnodes];
        DAD.f[DIR_M0P] = &distributionsAD[DIR_M0P * numberOfLBnodes];
        DAD.f[DIR_0PP] = &distributionsAD[DIR_0PP * numberOfLBnodes];
        DAD.f[DIR_0MM] = &distributionsAD[DIR_0MM * numberOfLBnodes];
        DAD.f[DIR_0PM] = &distributionsAD[DIR_0PM * numberOfLBnodes];
        DAD.f[DIR_0MP] = &distributionsAD[DIR_0MP * numberOfLBnodes];
        DAD.f[DIR_000] = &distributionsAD[DIR_000 * numberOfLBnodes];
        DAD.f[DIR_PPP] = &distributionsAD[DIR_PPP * numberOfLBnodes];
        DAD.f[DIR_MMP] = &distributionsAD[DIR_MMP * numberOfLBnodes];
        DAD.f[DIR_PMP] = &distributionsAD[DIR_PMP * numberOfLBnodes];
        DAD.f[DIR_MPP] = &distributionsAD[DIR_MPP * numberOfLBnodes];
        DAD.f[DIR_PPM] = &distributionsAD[DIR_PPM * numberOfLBnodes];
        DAD.f[DIR_MMM] = &distributionsAD[DIR_MMM * numberOfLBnodes];
        DAD.f[DIR_PMM] = &distributionsAD[DIR_PMM * numberOfLBnodes];
        DAD.f[DIR_MPM] = &distributionsAD[DIR_MPM * numberOfLBnodes];
    }
    else
    {
        DAD.f[DIR_M00] = &distributionsAD[DIR_P00 * numberOfLBnodes];
        DAD.f[DIR_P00] = &distributionsAD[DIR_M00 * numberOfLBnodes];
        DAD.f[DIR_0M0] = &distributionsAD[DIR_0P0 * numberOfLBnodes];
        DAD.f[DIR_0P0] = &distributionsAD[DIR_0M0 * numberOfLBnodes];
        DAD.f[DIR_00M] = &distributionsAD[DIR_00P * numberOfLBnodes];
        DAD.f[DIR_00P] = &distributionsAD[DIR_00M * numberOfLBnodes];
        DAD.f[DIR_MM0] = &distributionsAD[DIR_PP0 * numberOfLBnodes];
        DAD.f[DIR_PP0] = &distributionsAD[DIR_MM0 * numberOfLBnodes];
        DAD.f[DIR_MP0] = &distributionsAD[DIR_PM0 * numberOfLBnodes];
        DAD.f[DIR_PM0] = &distributionsAD[DIR_MP0 * numberOfLBnodes];
        DAD.f[DIR_M0M] = &distributionsAD[DIR_P0P * numberOfLBnodes];
        DAD.f[DIR_P0P] = &distributionsAD[DIR_M0M * numberOfLBnodes];
        DAD.f[DIR_M0P] = &distributionsAD[DIR_P0M * numberOfLBnodes];
        DAD.f[DIR_P0M] = &distributionsAD[DIR_M0P * numberOfLBnodes];
        DAD.f[DIR_0MM] = &distributionsAD[DIR_0PP * numberOfLBnodes];
        DAD.f[DIR_0PP] = &distributionsAD[DIR_0MM * numberOfLBnodes];
        DAD.f[DIR_0MP] = &distributionsAD[DIR_0PM * numberOfLBnodes];
        DAD.f[DIR_0PM] = &distributionsAD[DIR_0MP * numberOfLBnodes];
        DAD.f[DIR_000] = &distributionsAD[DIR_000 * numberOfLBnodes];
        DAD.f[DIR_PPP] = &distributionsAD[DIR_MMM * numberOfLBnodes];
        DAD.f[DIR_MMP] = &distributionsAD[DIR_PPM * numberOfLBnodes];
        DAD.f[DIR_PMP] = &distributionsAD[DIR_MPM * numberOfLBnodes];
        DAD.f[DIR_MPP] = &distributionsAD[DIR_PMM * numberOfLBnodes];
        DAD.f[DIR_PPM] = &distributionsAD[DIR_MMP * numberOfLBnodes];
        DAD.f[DIR_MMM] = &distributionsAD[DIR_PPP * numberOfLBnodes];
        DAD.f[DIR_PMM] = &distributionsAD[DIR_MPP * numberOfLBnodes];
        DAD.f[DIR_MPM] = &distributionsAD[DIR_PMP * numberOfLBnodes];
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
        q_dirE   = &Qarrays[DIR_P00 * numberOfBCnodes];
        q_dirW   = &Qarrays[DIR_M00 * numberOfBCnodes];
        q_dirN   = &Qarrays[DIR_0P0 * numberOfBCnodes];
        q_dirS   = &Qarrays[DIR_0M0 * numberOfBCnodes];
        q_dirT   = &Qarrays[DIR_00P * numberOfBCnodes];
        q_dirB   = &Qarrays[DIR_00M * numberOfBCnodes];
        q_dirNE  = &Qarrays[DIR_PP0 * numberOfBCnodes];
        q_dirSW  = &Qarrays[DIR_MM0 * numberOfBCnodes];
        q_dirSE  = &Qarrays[DIR_PM0 * numberOfBCnodes];
        q_dirNW  = &Qarrays[DIR_MP0 * numberOfBCnodes];
        q_dirTE  = &Qarrays[DIR_P0P * numberOfBCnodes];
        q_dirBW  = &Qarrays[DIR_M0M * numberOfBCnodes];
        q_dirBE  = &Qarrays[DIR_P0M * numberOfBCnodes];
        q_dirTW  = &Qarrays[DIR_M0P * numberOfBCnodes];
        q_dirTN  = &Qarrays[DIR_0PP * numberOfBCnodes];
        q_dirBS  = &Qarrays[DIR_0MM * numberOfBCnodes];
        q_dirBN  = &Qarrays[DIR_0PM * numberOfBCnodes];
        q_dirTS  = &Qarrays[DIR_0MP * numberOfBCnodes];
        q_dirTNE = &Qarrays[DIR_PPP * numberOfBCnodes];
        q_dirTSW = &Qarrays[DIR_MMP * numberOfBCnodes];
        q_dirTSE = &Qarrays[DIR_PMP * numberOfBCnodes];
        q_dirTNW = &Qarrays[DIR_MPP * numberOfBCnodes];
        q_dirBNE = &Qarrays[DIR_PPM * numberOfBCnodes];
        q_dirBSW = &Qarrays[DIR_MMM * numberOfBCnodes];
        q_dirBSE = &Qarrays[DIR_PMM * numberOfBCnodes];
        q_dirBNW = &Qarrays[DIR_MPM * numberOfBCnodes];
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

        f_W   = (D.f[DIR_P00])[ke];
        f_E   = (D.f[DIR_M00])[kw];
        f_S   = (D.f[DIR_0P0])[kn];
        f_N   = (D.f[DIR_0M0])[ks];
        f_B   = (D.f[DIR_00P])[kt];
        f_T   = (D.f[DIR_00M])[kb];
        f_SW  = (D.f[DIR_PP0])[kne];
        f_NE  = (D.f[DIR_MM0])[ksw];
        f_NW  = (D.f[DIR_PM0])[kse];
        f_SE  = (D.f[DIR_MP0])[knw];
        f_BW  = (D.f[DIR_P0P])[kte];
        f_TE  = (D.f[DIR_M0M])[kbw];
        f_TW  = (D.f[DIR_P0M])[kbe];
        f_BE  = (D.f[DIR_M0P])[ktw];
        f_BS  = (D.f[DIR_0PP])[ktn];
        f_TN  = (D.f[DIR_0MM])[kbs];
        f_TS  = (D.f[DIR_0PM])[kbn];
        f_BN  = (D.f[DIR_0MP])[kts];
        f_BSW = (D.f[DIR_PPP])[ktne];
        f_BNE = (D.f[DIR_MMP])[ktsw];
        f_BNW = (D.f[DIR_PMP])[ktse];
        f_BSE = (D.f[DIR_MPP])[ktnw];
        f_TSW = (D.f[DIR_PPM])[kbne];
        f_TNE = (D.f[DIR_MMM])[kbsw];
        f_TNW = (D.f[DIR_PMM])[kbse];
        f_TSE = (D.f[DIR_MPM])[kbnw];
        ////////////////////////////////////////////////////////////////////////////////
        real vx1, vx2, vx3, drho, q;
        drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]);

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
        f_W   = (DAD.f[DIR_P00])[ke];
        f_E   = (DAD.f[DIR_M00])[kw];
        f_S   = (DAD.f[DIR_0P0])[kn];
        f_N   = (DAD.f[DIR_0M0])[ks];
        f_B   = (DAD.f[DIR_00P])[kt];
        f_T   = (DAD.f[DIR_00M])[kb];
        f_SW  = (DAD.f[DIR_PP0])[kne];
        f_NE  = (DAD.f[DIR_MM0])[ksw];
        f_NW  = (DAD.f[DIR_PM0])[kse];
        f_SE  = (DAD.f[DIR_MP0])[knw];
        f_BW  = (DAD.f[DIR_P0P])[kte];
        f_TE  = (DAD.f[DIR_M0M])[kbw];
        f_TW  = (DAD.f[DIR_P0M])[kbe];
        f_BE  = (DAD.f[DIR_M0P])[ktw];
        f_BS  = (DAD.f[DIR_0PP])[ktn];
        f_TN  = (DAD.f[DIR_0MM])[kbs];
        f_TS  = (DAD.f[DIR_0PM])[kbn];
        f_BN  = (DAD.f[DIR_0MP])[kts];
        f_BSW = (DAD.f[DIR_PPP])[ktne];
        f_BNE = (DAD.f[DIR_MMP])[ktsw];
        f_BNW = (DAD.f[DIR_PMP])[ktse];
        f_BSE = (DAD.f[DIR_MPP])[ktnw];
        f_TSW = (DAD.f[DIR_PPM])[kbne];
        f_TNE = (DAD.f[DIR_MMM])[kbsw];
        f_TNW = (DAD.f[DIR_PMM])[kbse];
        f_TSE = (DAD.f[DIR_MPM])[kbnw];
        //////////////////////////////////////////////////////////////////////////
        if (!isEvenTimestep)
        {
            DAD.f[DIR_P00] = &distributionsAD[DIR_P00 * numberOfLBnodes];
            DAD.f[DIR_M00] = &distributionsAD[DIR_M00 * numberOfLBnodes];
            DAD.f[DIR_0P0] = &distributionsAD[DIR_0P0 * numberOfLBnodes];
            DAD.f[DIR_0M0] = &distributionsAD[DIR_0M0 * numberOfLBnodes];
            DAD.f[DIR_00P] = &distributionsAD[DIR_00P * numberOfLBnodes];
            DAD.f[DIR_00M] = &distributionsAD[DIR_00M * numberOfLBnodes];
            DAD.f[DIR_PP0] = &distributionsAD[DIR_PP0 * numberOfLBnodes];
            DAD.f[DIR_MM0] = &distributionsAD[DIR_MM0 * numberOfLBnodes];
            DAD.f[DIR_PM0] = &distributionsAD[DIR_PM0 * numberOfLBnodes];
            DAD.f[DIR_MP0] = &distributionsAD[DIR_MP0 * numberOfLBnodes];
            DAD.f[DIR_P0P] = &distributionsAD[DIR_P0P * numberOfLBnodes];
            DAD.f[DIR_M0M] = &distributionsAD[DIR_M0M * numberOfLBnodes];
            DAD.f[DIR_P0M] = &distributionsAD[DIR_P0M * numberOfLBnodes];
            DAD.f[DIR_M0P] = &distributionsAD[DIR_M0P * numberOfLBnodes];
            DAD.f[DIR_0PP] = &distributionsAD[DIR_0PP * numberOfLBnodes];
            DAD.f[DIR_0MM] = &distributionsAD[DIR_0MM * numberOfLBnodes];
            DAD.f[DIR_0PM] = &distributionsAD[DIR_0PM * numberOfLBnodes];
            DAD.f[DIR_0MP] = &distributionsAD[DIR_0MP * numberOfLBnodes];
            DAD.f[DIR_000] = &distributionsAD[DIR_000 * numberOfLBnodes];
            DAD.f[DIR_PPP] = &distributionsAD[DIR_PPP * numberOfLBnodes];
            DAD.f[DIR_MMP] = &distributionsAD[DIR_MMP * numberOfLBnodes];
            DAD.f[DIR_PMP] = &distributionsAD[DIR_PMP * numberOfLBnodes];
            DAD.f[DIR_MPP] = &distributionsAD[DIR_MPP * numberOfLBnodes];
            DAD.f[DIR_PPM] = &distributionsAD[DIR_PPM * numberOfLBnodes];
            DAD.f[DIR_MMM] = &distributionsAD[DIR_MMM * numberOfLBnodes];
            DAD.f[DIR_PMM] = &distributionsAD[DIR_PMM * numberOfLBnodes];
            DAD.f[DIR_MPM] = &distributionsAD[DIR_MPM * numberOfLBnodes];
        }
        else
        {
            DAD.f[DIR_M00] = &distributionsAD[DIR_P00 * numberOfLBnodes];
            DAD.f[DIR_P00] = &distributionsAD[DIR_M00 * numberOfLBnodes];
            DAD.f[DIR_0M0] = &distributionsAD[DIR_0P0 * numberOfLBnodes];
            DAD.f[DIR_0P0] = &distributionsAD[DIR_0M0 * numberOfLBnodes];
            DAD.f[DIR_00M] = &distributionsAD[DIR_00P * numberOfLBnodes];
            DAD.f[DIR_00P] = &distributionsAD[DIR_00M * numberOfLBnodes];
            DAD.f[DIR_MM0] = &distributionsAD[DIR_PP0 * numberOfLBnodes];
            DAD.f[DIR_PP0] = &distributionsAD[DIR_MM0 * numberOfLBnodes];
            DAD.f[DIR_MP0] = &distributionsAD[DIR_PM0 * numberOfLBnodes];
            DAD.f[DIR_PM0] = &distributionsAD[DIR_MP0 * numberOfLBnodes];
            DAD.f[DIR_M0M] = &distributionsAD[DIR_P0P * numberOfLBnodes];
            DAD.f[DIR_P0P] = &distributionsAD[DIR_M0M * numberOfLBnodes];
            DAD.f[DIR_M0P] = &distributionsAD[DIR_P0M * numberOfLBnodes];
            DAD.f[DIR_P0M] = &distributionsAD[DIR_M0P * numberOfLBnodes];
            DAD.f[DIR_0MM] = &distributionsAD[DIR_0PP * numberOfLBnodes];
            DAD.f[DIR_0PP] = &distributionsAD[DIR_0MM * numberOfLBnodes];
            DAD.f[DIR_0MP] = &distributionsAD[DIR_0PM * numberOfLBnodes];
            DAD.f[DIR_0PM] = &distributionsAD[DIR_0MP * numberOfLBnodes];
            DAD.f[DIR_000] = &distributionsAD[DIR_000 * numberOfLBnodes];
            DAD.f[DIR_PPP] = &distributionsAD[DIR_MMM * numberOfLBnodes];
            DAD.f[DIR_MMP] = &distributionsAD[DIR_PPM * numberOfLBnodes];
            DAD.f[DIR_PMP] = &distributionsAD[DIR_MPM * numberOfLBnodes];
            DAD.f[DIR_MPP] = &distributionsAD[DIR_PMM * numberOfLBnodes];
            DAD.f[DIR_PPM] = &distributionsAD[DIR_MMP * numberOfLBnodes];
            DAD.f[DIR_MMM] = &distributionsAD[DIR_PPP * numberOfLBnodes];
            DAD.f[DIR_PMM] = &distributionsAD[DIR_MPP * numberOfLBnodes];
            DAD.f[DIR_MPM] = &distributionsAD[DIR_PMP * numberOfLBnodes];
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        real concentration =
            f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]);

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
        q = q_dirE[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_M00])[kw  ] = calcDistributionBC_AD(q, c2o27,   vx1,         cu_sq, f_E,   f_W,   omegaDiffusivity,        jTan1,       concentration); }
        q = q_dirW[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_P00])[ke  ] = calcDistributionBC_AD(q, c2o27,  -vx1,         cu_sq, f_W,   f_E,   omegaDiffusivity,       -jTan1,       concentration); }
        q = q_dirN[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_0M0])[ks  ] = calcDistributionBC_AD(q, c2o27,   vx2,         cu_sq, f_N,   f_S,   omegaDiffusivity,        jTan2,       concentration); }
        q = q_dirS[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_0P0])[kn  ] = calcDistributionBC_AD(q, c2o27,  -vx2,         cu_sq, f_S,   f_N,   omegaDiffusivity,       -jTan2,       concentration); }
        q = q_dirT[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_00M])[kb  ] = calcDistributionBC_AD(q, c2o27,   vx3,         cu_sq, f_T,   f_B,   omegaDiffusivity,        jTan3,       concentration); }
        q = q_dirB[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_00P])[kt  ] = calcDistributionBC_AD(q, c2o27,  -vx3,         cu_sq, f_B,   f_T,   omegaDiffusivity,       -jTan3,       concentration); }
        q = q_dirNE[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_MM0])[ksw ] = calcDistributionBC_AD(q, c1o54,   vx1+vx2,     cu_sq, f_NE,  f_SW,  omegaDiffusivity,  jTan1+jTan2,       concentration); }
        q = q_dirSW[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_PP0])[kne ] = calcDistributionBC_AD(q, c1o54,  -vx1-vx2,     cu_sq, f_SW,  f_NE,  omegaDiffusivity, -jTan1-jTan2,       concentration); }
        q = q_dirSE[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_MP0])[knw ] = calcDistributionBC_AD(q, c1o54,   vx1-vx2,     cu_sq, f_SE,  f_NW,  omegaDiffusivity,  jTan1-jTan2,       concentration); }
        q = q_dirNW[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_PM0])[kse ] = calcDistributionBC_AD(q, c1o54,  -vx1+vx2,     cu_sq, f_NW,  f_SE,  omegaDiffusivity, -jTan1+jTan2,       concentration); }
        q = q_dirTE[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_M0M])[kbw ] = calcDistributionBC_AD(q, c1o54,   vx1    +vx3, cu_sq, f_TE,  f_BW,  omegaDiffusivity,  jTan1      +jTan3, concentration); }
        q = q_dirBW[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_P0P])[kte ] = calcDistributionBC_AD(q, c1o54,  -vx1    -vx3, cu_sq, f_BW,  f_TE,  omegaDiffusivity, -jTan1      -jTan3, concentration); }
        q = q_dirBE[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_M0P])[ktw ] = calcDistributionBC_AD(q, c1o54,   vx1    -vx3, cu_sq, f_BE,  f_TW,  omegaDiffusivity,  jTan1      -jTan3, concentration); }
        q = q_dirTW[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_P0M])[kbe ] = calcDistributionBC_AD(q, c1o54,  -vx1    +vx3, cu_sq, f_TW,  f_BE,  omegaDiffusivity, -jTan1      +jTan3, concentration); }
        q = q_dirTN[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_0MM])[kbs ] = calcDistributionBC_AD(q, c1o54,       vx2+vx3, cu_sq, f_TN,  f_BS,  omegaDiffusivity,        jTan2+jTan3, concentration); }
        q = q_dirBS[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_0PP])[ktn ] = calcDistributionBC_AD(q, c1o54,      -vx2-vx3, cu_sq, f_BS,  f_TN,  omegaDiffusivity,       -jTan2-jTan3, concentration); }
        q = q_dirBN[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_0MP])[kts ] = calcDistributionBC_AD(q, c1o54,       vx2-vx3, cu_sq, f_BN,  f_TS,  omegaDiffusivity,        jTan2-jTan3, concentration); }
        q = q_dirTS[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_0PM])[kbn ] = calcDistributionBC_AD(q, c1o54,      -vx2+vx3, cu_sq, f_TS,  f_BN,  omegaDiffusivity,       -jTan2+jTan3, concentration); }
        q = q_dirTNE[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_MMM])[kbsw] = calcDistributionBC_AD(q, c1o216,  vx1+vx2+vx3, cu_sq, f_TNE, f_BSW, omegaDiffusivity,  jTan1+jTan2+jTan3, concentration); }
        q = q_dirBSW[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_PPP])[ktne] = calcDistributionBC_AD(q, c1o216, -vx1-vx2-vx3, cu_sq, f_BSW, f_TNE, omegaDiffusivity, -jTan1-jTan2-jTan3, concentration); }
        q = q_dirBNE[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_MMP])[ktsw] = calcDistributionBC_AD(q, c1o216,  vx1+vx2-vx3, cu_sq, f_BNE, f_TSW, omegaDiffusivity,  jTan1+jTan2-jTan3, concentration); }
        q = q_dirTSW[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_PPM])[kbne] = calcDistributionBC_AD(q, c1o216, -vx1-vx2+vx3, cu_sq, f_TSW, f_BNE, omegaDiffusivity, -jTan1-jTan2+jTan3, concentration); }
        q = q_dirTSE[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_MPM])[kbnw] = calcDistributionBC_AD(q, c1o216,  vx1-vx2+vx3, cu_sq, f_TSE, f_BNW, omegaDiffusivity,  jTan1-jTan2+jTan3, concentration); }
        q = q_dirBNW[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_PMP])[ktse] = calcDistributionBC_AD(q, c1o216, -vx1+vx2-vx3, cu_sq, f_BNW, f_TSE, omegaDiffusivity, -jTan1+jTan2-jTan3, concentration); }
        q = q_dirBSE[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_MPP])[ktnw] = calcDistributionBC_AD(q, c1o216,  vx1-vx2-vx3, cu_sq, f_BSE, f_TNW, omegaDiffusivity,  jTan1-jTan2-jTan3, concentration); }
        q = q_dirTNW[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[DIR_PMM])[kbse] = calcDistributionBC_AD(q, c1o216, -vx1+vx2+vx3, cu_sq, f_TNW, f_BSE, omegaDiffusivity, -jTan1+jTan2+jTan3, concentration); }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
