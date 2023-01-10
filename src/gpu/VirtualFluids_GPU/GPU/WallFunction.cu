/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;


//////////////////////////////////////////////////////////////////////////////
__global__ void WallFunction27(
										  real* vx,
										  real* vy,
										  real* vz,
										  real* DD, 
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
      real VeloX = vx[k];
      //real VeloY = vy[k];
      //real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
      ////////////////////////////////////////////////////////////////////////////////
      //real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
      //      *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
      //      *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
      //      *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
      //      *q_dirBSE, *q_dirBNW; 
      //q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      //q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      //q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      //q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      //q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      //q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
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
      // real vx2, vx3, feq, q;
      real vx1, drho;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]); 

       vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                 (f_E - f_W)) / (c1o1 + drho); 
         

    //   vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
    //              ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
    //              (f_N - f_S)) / (c1o1 + drho); 

    //   vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
    //              (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
    //              (f_T - f_B)) / (c1o1 + drho); 

      //real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (one + drho);

	  real nu = c1o3 * (c1o1 / om1 - c1o2);
	  real qw = c1o1;
	  real uTau = sqrt(nu * (vx1 - VeloX) / qw);

	  if (abs(uTau)/nu>11){
	  uTau = vx1 * 0.41 / (log10(9.8 * uTau * qw / nu));
	  

	  
	  vx[k] = vx1 - uTau * uTau * qw / nu;
	  vx[k] = (vx[k]> 0.05) ? 0.05 : ((vx[k]< -0.05) ? -0.05 : vx[k] );  
	  }
	  else{ vx[k]=c0o1; }
	  //vy[k] = 0.01;							//Test...muss wieder raus
	  //vz[k] = 0.01;							//Test...muss wieder raus

   //   //////////////////////////////////////////////////////////////////////////
   //   if (isEvenTimestep==false)
   //   {
   //      D.f[DIR_P00] = &DD[DIR_P00 * size_Mat];
   //      D.f[DIR_M00] = &DD[DIR_M00 * size_Mat];
   //      D.f[DIR_0P0] = &DD[DIR_0P0 * size_Mat];
   //      D.f[DIR_0M0] = &DD[DIR_0M0 * size_Mat];
   //      D.f[DIR_00P] = &DD[DIR_00P * size_Mat];
   //      D.f[DIR_00M] = &DD[DIR_00M * size_Mat];
   //      D.f[DIR_PP0] = &DD[DIR_PP0 * size_Mat];
   //      D.f[DIR_MM0] = &DD[DIR_MM0 * size_Mat];
   //      D.f[DIR_PM0] = &DD[DIR_PM0 * size_Mat];
   //      D.f[DIR_MP0] = &DD[DIR_MP0 * size_Mat];
   //      D.f[DIR_P0P] = &DD[DIR_P0P * size_Mat];
   //      D.f[DIR_M0M] = &DD[DIR_M0M * size_Mat];
   //      D.f[DIR_P0M] = &DD[DIR_P0M * size_Mat];
   //      D.f[DIR_M0P] = &DD[DIR_M0P * size_Mat];
   //      D.f[DIR_0PP] = &DD[DIR_0PP * size_Mat];
   //      D.f[DIR_0MM] = &DD[DIR_0MM * size_Mat];
   //      D.f[DIR_0PM] = &DD[DIR_0PM * size_Mat];
   //      D.f[DIR_0MP] = &DD[DIR_0MP * size_Mat];
   //      D.f[DIR_000] = &DD[DIR_000 * size_Mat];
   //      D.f[DIR_PPP] = &DD[DIR_PPP * size_Mat];
   //      D.f[DIR_MMP] = &DD[DIR_MMP * size_Mat];
   //      D.f[DIR_PMP] = &DD[DIR_PMP * size_Mat];
   //      D.f[DIR_MPP] = &DD[DIR_MPP * size_Mat];
   //      D.f[DIR_PPM] = &DD[DIR_PPM * size_Mat];
   //      D.f[DIR_MMM] = &DD[DIR_MMM * size_Mat];
   //      D.f[DIR_PMM] = &DD[DIR_PMM * size_Mat];
   //      D.f[DIR_MPM] = &DD[DIR_MPM * size_Mat];
   //   } 
   //   else
   //   {
   //      D.f[DIR_M00] = &DD[DIR_P00 * size_Mat];
   //      D.f[DIR_P00] = &DD[DIR_M00 * size_Mat];
   //      D.f[DIR_0M0] = &DD[DIR_0P0 * size_Mat];
   //      D.f[DIR_0P0] = &DD[DIR_0M0 * size_Mat];
   //      D.f[DIR_00M] = &DD[DIR_00P * size_Mat];
   //      D.f[DIR_00P] = &DD[DIR_00M * size_Mat];
   //      D.f[DIR_MM0] = &DD[DIR_PP0 * size_Mat];
   //      D.f[DIR_PP0] = &DD[DIR_MM0 * size_Mat];
   //      D.f[DIR_MP0] = &DD[DIR_PM0 * size_Mat];
   //      D.f[DIR_PM0] = &DD[DIR_MP0 * size_Mat];
   //      D.f[DIR_M0M] = &DD[DIR_P0P * size_Mat];
   //      D.f[DIR_P0P] = &DD[DIR_M0M * size_Mat];
   //      D.f[DIR_M0P] = &DD[DIR_P0M * size_Mat];
   //      D.f[DIR_P0M] = &DD[DIR_M0P * size_Mat];
   //      D.f[DIR_0MM] = &DD[DIR_0PP * size_Mat];
   //      D.f[DIR_0PP] = &DD[DIR_0MM * size_Mat];
   //      D.f[DIR_0MP] = &DD[DIR_0PM * size_Mat];
   //      D.f[DIR_0PM] = &DD[DIR_0MP * size_Mat];
   //      D.f[DIR_000] = &DD[DIR_000 * size_Mat];
   //      D.f[DIR_PPP] = &DD[DIR_MMM * size_Mat];
   //      D.f[DIR_MMP] = &DD[DIR_PPM * size_Mat];
   //      D.f[DIR_PMP] = &DD[DIR_MPM * size_Mat];
   //      D.f[DIR_MPP] = &DD[DIR_PMM * size_Mat];
   //      D.f[DIR_PPM] = &DD[DIR_MMP * size_Mat];
   //      D.f[DIR_MMM] = &DD[DIR_PPP * size_Mat];
   //      D.f[DIR_PMM] = &DD[DIR_MPP * size_Mat];
   //      D.f[DIR_MPM] = &DD[DIR_PMP * size_Mat];
   //   }
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   //Test
   //   //(D.f[DIR_000])[k]=c1o10;
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  ////ToDo anders Klammern

   //   q = q_dirE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
   //      (D.f[DIR_M00])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[DIR_M00])[kw]=zero;
   //   }

   //   q = q_dirW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
   //      (D.f[DIR_P00])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[DIR_P00])[ke]=zero;
   //   }

   //   q = q_dirN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[DIR_0M0])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[DIR_0M0])[ks]=zero;
   //   }

   //   q = q_dirS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[DIR_0P0])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[DIR_0P0])[kn]=zero;
   //   }

   //   q = q_dirT[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_00M])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[DIR_00M])[kb]=one;
   //   }

   //   q = q_dirB[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_00P])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[DIR_00P])[kt]=zero;
   //   }

   //   q = q_dirNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[DIR_MM0])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);// - c1over54 * drho;
   //      //(D.f[DIR_MM0])[ksw]=zero;
   //   }

   //   q = q_dirSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[DIR_PP0])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);// - c1over54 * drho;
   //      //(D.f[DIR_PP0])[kne]=zero;
   //   }

   //   q = q_dirSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[DIR_MP0])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);// - c1over54 * drho;
   //      //(D.f[DIR_MP0])[knw]=zero;
   //   }

   //   q = q_dirNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[DIR_PM0])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);// - c1over54 * drho;
   //      //(D.f[DIR_PM0])[kse]=zero;
   //   }

   //   q = q_dirTE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_M0M])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[DIR_M0M])[kbw]=zero;
   //   }

   //   q = q_dirBW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_P0P])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[DIR_P0P])[kte]=zero;
   //   }

   //   q = q_dirBE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_M0P])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[DIR_M0P])[ktw]=zero;
   //   }

   //   q = q_dirTW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_P0M])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[DIR_P0M])[kbe]=zero;
   //   }

   //   q = q_dirTN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_0MM])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[DIR_0MM])[kbs]=zero;
   //   }

   //   q = q_dirBS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_0PP])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[DIR_0PP])[ktn]=zero;
   //   }

   //   q = q_dirBN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_0MP])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[DIR_0MP])[kts]=zero;
   //   }

   //   q = q_dirTS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_0PM])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[DIR_0PM])[kbn]=zero;
   //   }

   //   q = q_dirTNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_MMM])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[DIR_MMM])[kbsw]=zero;
   //   }

   //   q = q_dirBSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_PPP])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[DIR_PPP])[ktne]=zero;
   //   }

   //   q = q_dirBNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_MMP])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[DIR_MMP])[ktsw]=zero;
   //   }

   //   q = q_dirTSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_PPM])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[DIR_PPM])[kbne]=zero;
   //   }

   //   q = q_dirTSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_MPM])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[DIR_MPM])[kbnw]=zero;
   //   }

   //   q = q_dirBNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_PMP])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[DIR_PMP])[ktse]=zero;
   //   }

   //   q = q_dirBSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_MPP])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[DIR_MPP])[ktnw]=zero;
   //   }

   //   q = q_dirTNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[DIR_PMM])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[DIR_PMM])[kbse]=zero;
   //   }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////









