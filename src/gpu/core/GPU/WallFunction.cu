/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
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
      D.f[dP00] = &DD[dP00 * numberOfLBnodes];
      D.f[dM00] = &DD[dM00 * numberOfLBnodes];
      D.f[d0P0] = &DD[d0P0 * numberOfLBnodes];
      D.f[d0M0] = &DD[d0M0 * numberOfLBnodes];
      D.f[d00P] = &DD[d00P * numberOfLBnodes];
      D.f[d00M] = &DD[d00M * numberOfLBnodes];
      D.f[dPP0] = &DD[dPP0 * numberOfLBnodes];
      D.f[dMM0] = &DD[dMM0 * numberOfLBnodes];
      D.f[dPM0] = &DD[dPM0 * numberOfLBnodes];
      D.f[dMP0] = &DD[dMP0 * numberOfLBnodes];
      D.f[dP0P] = &DD[dP0P * numberOfLBnodes];
      D.f[dM0M] = &DD[dM0M * numberOfLBnodes];
      D.f[dP0M] = &DD[dP0M * numberOfLBnodes];
      D.f[dM0P] = &DD[dM0P * numberOfLBnodes];
      D.f[d0PP] = &DD[d0PP * numberOfLBnodes];
      D.f[d0MM] = &DD[d0MM * numberOfLBnodes];
      D.f[d0PM] = &DD[d0PM * numberOfLBnodes];
      D.f[d0MP] = &DD[d0MP * numberOfLBnodes];
      D.f[d000] = &DD[d000 * numberOfLBnodes];
      D.f[dPPP] = &DD[dPPP * numberOfLBnodes];
      D.f[dMMP] = &DD[dMMP * numberOfLBnodes];
      D.f[dPMP] = &DD[dPMP * numberOfLBnodes];
      D.f[dMPP] = &DD[dMPP * numberOfLBnodes];
      D.f[dPPM] = &DD[dPPM * numberOfLBnodes];
      D.f[dMMM] = &DD[dMMM * numberOfLBnodes];
      D.f[dPMM] = &DD[dPMM * numberOfLBnodes];
      D.f[dMPM] = &DD[dMPM * numberOfLBnodes];
   } 
   else
   {
      D.f[dM00] = &DD[dP00 * numberOfLBnodes];
      D.f[dP00] = &DD[dM00 * numberOfLBnodes];
      D.f[d0M0] = &DD[d0P0 * numberOfLBnodes];
      D.f[d0P0] = &DD[d0M0 * numberOfLBnodes];
      D.f[d00M] = &DD[d00P * numberOfLBnodes];
      D.f[d00P] = &DD[d00M * numberOfLBnodes];
      D.f[dMM0] = &DD[dPP0 * numberOfLBnodes];
      D.f[dPP0] = &DD[dMM0 * numberOfLBnodes];
      D.f[dMP0] = &DD[dPM0 * numberOfLBnodes];
      D.f[dPM0] = &DD[dMP0 * numberOfLBnodes];
      D.f[dM0M] = &DD[dP0P * numberOfLBnodes];
      D.f[dP0P] = &DD[dM0M * numberOfLBnodes];
      D.f[dM0P] = &DD[dP0M * numberOfLBnodes];
      D.f[dP0M] = &DD[dM0P * numberOfLBnodes];
      D.f[d0MM] = &DD[d0PP * numberOfLBnodes];
      D.f[d0PP] = &DD[d0MM * numberOfLBnodes];
      D.f[d0MP] = &DD[d0PM * numberOfLBnodes];
      D.f[d0PM] = &DD[d0MP * numberOfLBnodes];
      D.f[d000] = &DD[d000 * numberOfLBnodes];
      D.f[dPPP] = &DD[dMMM * numberOfLBnodes];
      D.f[dMMP] = &DD[dPPM * numberOfLBnodes];
      D.f[dPMP] = &DD[dMPM * numberOfLBnodes];
      D.f[dMPP] = &DD[dPMM * numberOfLBnodes];
      D.f[dPPM] = &DD[dMMP * numberOfLBnodes];
      D.f[dMMM] = &DD[dPPP * numberOfLBnodes];
      D.f[dPMM] = &DD[dMPP * numberOfLBnodes];
      D.f[dMPM] = &DD[dPMP * numberOfLBnodes];
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
      //q_dirE   = &QQ[dP00 * numberOfBCnodes];
      //q_dirW   = &QQ[dM00 * numberOfBCnodes];
      //q_dirN   = &QQ[d0P0 * numberOfBCnodes];
      //q_dirS   = &QQ[d0M0 * numberOfBCnodes];
      //q_dirT   = &QQ[d00P * numberOfBCnodes];
      //q_dirB   = &QQ[d00M * numberOfBCnodes];
      //q_dirNE  = &QQ[dPP0 * numberOfBCnodes];
      //q_dirSW  = &QQ[dMM0 * numberOfBCnodes];
      //q_dirSE  = &QQ[dPM0 * numberOfBCnodes];
      //q_dirNW  = &QQ[dMP0 * numberOfBCnodes];
      //q_dirTE  = &QQ[dP0P * numberOfBCnodes];
      //q_dirBW  = &QQ[dM0M * numberOfBCnodes];
      //q_dirBE  = &QQ[dP0M * numberOfBCnodes];
      //q_dirTW  = &QQ[dM0P * numberOfBCnodes];
      //q_dirTN  = &QQ[d0PP * numberOfBCnodes];
      //q_dirBS  = &QQ[d0MM * numberOfBCnodes];
      //q_dirBN  = &QQ[d0PM * numberOfBCnodes];
      //q_dirTS  = &QQ[d0MP * numberOfBCnodes];
      //q_dirTNE = &QQ[dPPP * numberOfBCnodes];
      //q_dirTSW = &QQ[dMMP * numberOfBCnodes];
      //q_dirTSE = &QQ[dPMP * numberOfBCnodes];
      //q_dirTNW = &QQ[dMPP * numberOfBCnodes];
      //q_dirBNE = &QQ[dPPM * numberOfBCnodes];
      //q_dirBSW = &QQ[dMMM * numberOfBCnodes];
      //q_dirBSE = &QQ[dPMM * numberOfBCnodes];
      //q_dirBNW = &QQ[dMPM * numberOfBCnodes];
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

      f_W    = (D.f[dP00])[ke   ];
      f_E    = (D.f[dM00])[kw   ];
      f_S    = (D.f[d0P0])[kn   ];
      f_N    = (D.f[d0M0])[ks   ];
      f_B    = (D.f[d00P])[kt   ];
      f_T    = (D.f[d00M])[kb   ];
      f_SW   = (D.f[dPP0])[kne  ];
      f_NE   = (D.f[dMM0])[ksw  ];
      f_NW   = (D.f[dPM0])[kse  ];
      f_SE   = (D.f[dMP0])[knw  ];
      f_BW   = (D.f[dP0P])[kte  ];
      f_TE   = (D.f[dM0M])[kbw  ];
      f_TW   = (D.f[dP0M])[kbe  ];
      f_BE   = (D.f[dM0P])[ktw  ];
      f_BS   = (D.f[d0PP])[ktn  ];
      f_TN   = (D.f[d0MM])[kbs  ];
      f_TS   = (D.f[d0PM])[kbn  ];
      f_BN   = (D.f[d0MP])[kts  ];
      f_BSW  = (D.f[dPPP])[ktne ];
      f_BNE  = (D.f[dMMP])[ktsw ];
      f_BNW  = (D.f[dPMP])[ktse ];
      f_BSE  = (D.f[dMPP])[ktnw ];
      f_TSW  = (D.f[dPPM])[kbne ];
      f_TNE  = (D.f[dMMM])[kbsw ];
      f_TNW  = (D.f[dPMM])[kbse ];
      f_TSE  = (D.f[dMPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      // real vx2, vx3, feq, q;
      real vx1, drho;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]); 

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
   //      D.f[dP00] = &DD[dP00 * size_Mat];
   //      D.f[dM00] = &DD[dM00 * size_Mat];
   //      D.f[d0P0] = &DD[d0P0 * size_Mat];
   //      D.f[d0M0] = &DD[d0M0 * size_Mat];
   //      D.f[d00P] = &DD[d00P * size_Mat];
   //      D.f[d00M] = &DD[d00M * size_Mat];
   //      D.f[dPP0] = &DD[dPP0 * size_Mat];
   //      D.f[dMM0] = &DD[dMM0 * size_Mat];
   //      D.f[dPM0] = &DD[dPM0 * size_Mat];
   //      D.f[dMP0] = &DD[dMP0 * size_Mat];
   //      D.f[dP0P] = &DD[dP0P * size_Mat];
   //      D.f[dM0M] = &DD[dM0M * size_Mat];
   //      D.f[dP0M] = &DD[dP0M * size_Mat];
   //      D.f[dM0P] = &DD[dM0P * size_Mat];
   //      D.f[d0PP] = &DD[d0PP * size_Mat];
   //      D.f[d0MM] = &DD[d0MM * size_Mat];
   //      D.f[d0PM] = &DD[d0PM * size_Mat];
   //      D.f[d0MP] = &DD[d0MP * size_Mat];
   //      D.f[d000] = &DD[d000 * size_Mat];
   //      D.f[dPPP] = &DD[dPPP * size_Mat];
   //      D.f[dMMP] = &DD[dMMP * size_Mat];
   //      D.f[dPMP] = &DD[dPMP * size_Mat];
   //      D.f[dMPP] = &DD[dMPP * size_Mat];
   //      D.f[dPPM] = &DD[dPPM * size_Mat];
   //      D.f[dMMM] = &DD[dMMM * size_Mat];
   //      D.f[dPMM] = &DD[dPMM * size_Mat];
   //      D.f[dMPM] = &DD[dMPM * size_Mat];
   //   } 
   //   else
   //   {
   //      D.f[dM00] = &DD[dP00 * size_Mat];
   //      D.f[dP00] = &DD[dM00 * size_Mat];
   //      D.f[d0M0] = &DD[d0P0 * size_Mat];
   //      D.f[d0P0] = &DD[d0M0 * size_Mat];
   //      D.f[d00M] = &DD[d00P * size_Mat];
   //      D.f[d00P] = &DD[d00M * size_Mat];
   //      D.f[dMM0] = &DD[dPP0 * size_Mat];
   //      D.f[dPP0] = &DD[dMM0 * size_Mat];
   //      D.f[dMP0] = &DD[dPM0 * size_Mat];
   //      D.f[dPM0] = &DD[dMP0 * size_Mat];
   //      D.f[dM0M] = &DD[dP0P * size_Mat];
   //      D.f[dP0P] = &DD[dM0M * size_Mat];
   //      D.f[dM0P] = &DD[dP0M * size_Mat];
   //      D.f[dP0M] = &DD[dM0P * size_Mat];
   //      D.f[d0MM] = &DD[d0PP * size_Mat];
   //      D.f[d0PP] = &DD[d0MM * size_Mat];
   //      D.f[d0MP] = &DD[d0PM * size_Mat];
   //      D.f[d0PM] = &DD[d0MP * size_Mat];
   //      D.f[d000] = &DD[d000 * size_Mat];
   //      D.f[dPPP] = &DD[dMMM * size_Mat];
   //      D.f[dMMP] = &DD[dPPM * size_Mat];
   //      D.f[dPMP] = &DD[dMPM * size_Mat];
   //      D.f[dMPP] = &DD[dPMM * size_Mat];
   //      D.f[dPPM] = &DD[dMMP * size_Mat];
   //      D.f[dMMM] = &DD[dPPP * size_Mat];
   //      D.f[dPMM] = &DD[dMPP * size_Mat];
   //      D.f[dMPM] = &DD[dPMP * size_Mat];
   //   }
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   //Test
   //   //(D.f[d000])[k]=c1o10;
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  ////ToDo anders Klammern

   //   q = q_dirE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
   //      (D.f[dM00])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[dM00])[kw]=zero;
   //   }

   //   q = q_dirW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
   //      (D.f[dP00])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[dP00])[ke]=zero;
   //   }

   //   q = q_dirN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[d0M0])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[d0M0])[ks]=zero;
   //   }

   //   q = q_dirS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[d0P0])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[d0P0])[kn]=zero;
   //   }

   //   q = q_dirT[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
   //      (D.f[d00M])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[d00M])[kb]=one;
   //   }

   //   q = q_dirB[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
   //      (D.f[d00P])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[d00P])[kt]=zero;
   //   }

   //   q = q_dirNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[dMM0])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);// - c1over54 * drho;
   //      //(D.f[dMM0])[ksw]=zero;
   //   }

   //   q = q_dirSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[dPP0])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);// - c1over54 * drho;
   //      //(D.f[dPP0])[kne]=zero;
   //   }

   //   q = q_dirSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[dMP0])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);// - c1over54 * drho;
   //      //(D.f[dMP0])[knw]=zero;
   //   }

   //   q = q_dirNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[dPM0])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);// - c1over54 * drho;
   //      //(D.f[dPM0])[kse]=zero;
   //   }

   //   q = q_dirTE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
   //      (D.f[dM0M])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[dM0M])[kbw]=zero;
   //   }

   //   q = q_dirBW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
   //      (D.f[dP0P])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[dP0P])[kte]=zero;
   //   }

   //   q = q_dirBE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
   //      (D.f[dM0P])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[dM0P])[ktw]=zero;
   //   }

   //   q = q_dirTW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
   //      (D.f[dP0M])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[dP0M])[kbe]=zero;
   //   }

   //   q = q_dirTN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[d0MM])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[d0MM])[kbs]=zero;
   //   }

   //   q = q_dirBS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[d0PP])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[d0PP])[ktn]=zero;
   //   }

   //   q = q_dirBN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[d0MP])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[d0MP])[kts]=zero;
   //   }

   //   q = q_dirTS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[d0PM])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[d0PM])[kbn]=zero;
   //   }

   //   q = q_dirTNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[dMMM])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[dMMM])[kbsw]=zero;
   //   }

   //   q = q_dirBSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[dPPP])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[dPPP])[ktne]=zero;
   //   }

   //   q = q_dirBNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[dMMP])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[dMMP])[ktsw]=zero;
   //   }

   //   q = q_dirTSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[dPPM])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[dPPM])[kbne]=zero;
   //   }

   //   q = q_dirTSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[dMPM])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[dMPM])[kbnw]=zero;
   //   }

   //   q = q_dirBNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[dPMP])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[dPMP])[ktse]=zero;
   //   }

   //   q = q_dirBSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[dMPP])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[dMPP])[ktnw]=zero;
   //   }

   //   q = q_dirTNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[dPMM])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[dPMM])[kbse]=zero;
   //   }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////









