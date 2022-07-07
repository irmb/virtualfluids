/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;


//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void WallFunction27(
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      real VeloX = vx[k];
      //real VeloY = vy[k];
      //real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
      ////////////////////////////////////////////////////////////////////////////////
      //real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
      //      *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
      //      *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
      //      *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
      //      *q_dirBSE, *q_dirBNW; 
      //q_dirE   = &QQ[E   * numberOfBCnodes];
      //q_dirW   = &QQ[W   * numberOfBCnodes];
      //q_dirN   = &QQ[N   * numberOfBCnodes];
      //q_dirS   = &QQ[S   * numberOfBCnodes];
      //q_dirT   = &QQ[T   * numberOfBCnodes];
      //q_dirB   = &QQ[B   * numberOfBCnodes];
      //q_dirNE  = &QQ[NE  * numberOfBCnodes];
      //q_dirSW  = &QQ[SW  * numberOfBCnodes];
      //q_dirSE  = &QQ[SE  * numberOfBCnodes];
      //q_dirNW  = &QQ[NW  * numberOfBCnodes];
      //q_dirTE  = &QQ[TE  * numberOfBCnodes];
      //q_dirBW  = &QQ[BW  * numberOfBCnodes];
      //q_dirBE  = &QQ[BE  * numberOfBCnodes];
      //q_dirTW  = &QQ[TW  * numberOfBCnodes];
      //q_dirTN  = &QQ[TN  * numberOfBCnodes];
      //q_dirBS  = &QQ[BS  * numberOfBCnodes];
      //q_dirBN  = &QQ[BN  * numberOfBCnodes];
      //q_dirTS  = &QQ[TS  * numberOfBCnodes];
      //q_dirTNE = &QQ[TNE * numberOfBCnodes];
      //q_dirTSW = &QQ[TSW * numberOfBCnodes];
      //q_dirTSE = &QQ[TSE * numberOfBCnodes];
      //q_dirTNW = &QQ[TNW * numberOfBCnodes];
      //q_dirBNE = &QQ[BNE * numberOfBCnodes];
      //q_dirBSW = &QQ[BSW * numberOfBCnodes];
      //q_dirBSE = &QQ[BSE * numberOfBCnodes];
      //q_dirBNW = &QQ[BNW * numberOfBCnodes];
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

      f_W    = (D.f[E   ])[ke   ];
      f_E    = (D.f[W   ])[kw   ];
      f_S    = (D.f[N   ])[kn   ];
      f_N    = (D.f[S   ])[ks   ];
      f_B    = (D.f[T   ])[kt   ];
      f_T    = (D.f[B   ])[kb   ];
      f_SW   = (D.f[NE  ])[kne  ];
      f_NE   = (D.f[SW  ])[ksw  ];
      f_NW   = (D.f[SE  ])[kse  ];
      f_SE   = (D.f[NW  ])[knw  ];
      f_BW   = (D.f[TE  ])[kte  ];
      f_TE   = (D.f[BW  ])[kbw  ];
      f_TW   = (D.f[BE  ])[kbe  ];
      f_BE   = (D.f[TW  ])[ktw  ];
      f_BS   = (D.f[TN  ])[ktn  ];
      f_TN   = (D.f[BS  ])[kbs  ];
      f_TS   = (D.f[BN  ])[kbn  ];
      f_BN   = (D.f[TS  ])[kts  ];
      f_BSW  = (D.f[TNE ])[ktne ];
      f_BNE  = (D.f[TSW ])[ktsw ];
      f_BNW  = (D.f[TSE ])[ktse ];
      f_BSE  = (D.f[TNW ])[ktnw ];
      f_TSW  = (D.f[BNE ])[kbne ];
      f_TNE  = (D.f[BSW ])[kbsw ];
      f_TNW  = (D.f[BSE ])[kbse ];
      f_TSE  = (D.f[BNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      // real vx2, vx3, feq, q;
      real vx1, drho;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]); 

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
   //      D.f[E   ] = &DD[E   *size_Mat];
   //      D.f[W   ] = &DD[W   *size_Mat];
   //      D.f[N   ] = &DD[N   *size_Mat];
   //      D.f[S   ] = &DD[S   *size_Mat];
   //      D.f[T   ] = &DD[T   *size_Mat];
   //      D.f[B   ] = &DD[B   *size_Mat];
   //      D.f[NE  ] = &DD[NE  *size_Mat];
   //      D.f[SW  ] = &DD[SW  *size_Mat];
   //      D.f[SE  ] = &DD[SE  *size_Mat];
   //      D.f[NW  ] = &DD[NW  *size_Mat];
   //      D.f[TE  ] = &DD[TE  *size_Mat];
   //      D.f[BW  ] = &DD[BW  *size_Mat];
   //      D.f[BE  ] = &DD[BE  *size_Mat];
   //      D.f[TW  ] = &DD[TW  *size_Mat];
   //      D.f[TN  ] = &DD[TN  *size_Mat];
   //      D.f[BS  ] = &DD[BS  *size_Mat];
   //      D.f[BN  ] = &DD[BN  *size_Mat];
   //      D.f[TS  ] = &DD[TS  *size_Mat];
   //      D.f[dirREST] = &DD[dirREST*size_Mat];
   //      D.f[TNE ] = &DD[TNE *size_Mat];
   //      D.f[TSW ] = &DD[TSW *size_Mat];
   //      D.f[TSE ] = &DD[TSE *size_Mat];
   //      D.f[TNW ] = &DD[TNW *size_Mat];
   //      D.f[BNE ] = &DD[BNE *size_Mat];
   //      D.f[BSW ] = &DD[BSW *size_Mat];
   //      D.f[BSE ] = &DD[BSE *size_Mat];
   //      D.f[BNW ] = &DD[BNW *size_Mat];
   //   } 
   //   else
   //   {
   //      D.f[W   ] = &DD[E   *size_Mat];
   //      D.f[E   ] = &DD[W   *size_Mat];
   //      D.f[S   ] = &DD[N   *size_Mat];
   //      D.f[N   ] = &DD[S   *size_Mat];
   //      D.f[B   ] = &DD[T   *size_Mat];
   //      D.f[T   ] = &DD[B   *size_Mat];
   //      D.f[SW  ] = &DD[NE  *size_Mat];
   //      D.f[NE  ] = &DD[SW  *size_Mat];
   //      D.f[NW  ] = &DD[SE  *size_Mat];
   //      D.f[SE  ] = &DD[NW  *size_Mat];
   //      D.f[BW  ] = &DD[TE  *size_Mat];
   //      D.f[TE  ] = &DD[BW  *size_Mat];
   //      D.f[TW  ] = &DD[BE  *size_Mat];
   //      D.f[BE  ] = &DD[TW  *size_Mat];
   //      D.f[BS  ] = &DD[TN  *size_Mat];
   //      D.f[TN  ] = &DD[BS  *size_Mat];
   //      D.f[TS  ] = &DD[BN  *size_Mat];
   //      D.f[BN  ] = &DD[TS  *size_Mat];
   //      D.f[dirREST] = &DD[dirREST*size_Mat];
   //      D.f[TNE ] = &DD[BSW *size_Mat];
   //      D.f[TSW ] = &DD[BNE *size_Mat];
   //      D.f[TSE ] = &DD[BNW *size_Mat];
   //      D.f[TNW ] = &DD[BSE *size_Mat];
   //      D.f[BNE ] = &DD[TSW *size_Mat];
   //      D.f[BSW ] = &DD[TNE *size_Mat];
   //      D.f[BSE ] = &DD[TNW *size_Mat];
   //      D.f[BNW ] = &DD[TSE *size_Mat];
   //   }
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   //Test
   //   //(D.f[dirREST])[k]=c1o10;
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  ////ToDo anders Klammern

   //   q = q_dirE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
   //      (D.f[W])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[W])[kw]=zero;
   //   }

   //   q = q_dirW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
   //      (D.f[E])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[E])[ke]=zero;
   //   }

   //   q = q_dirN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[S])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[S])[ks]=zero;
   //   }

   //   q = q_dirS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[N])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[N])[kn]=zero;
   //   }

   //   q = q_dirT[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
   //      (D.f[B])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[B])[kb]=one;
   //   }

   //   q = q_dirB[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
   //      (D.f[T])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);// - c2over27 * drho;
   //      //(D.f[T])[kt]=zero;
   //   }

   //   q = q_dirNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[SW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);// - c1over54 * drho;
   //      //(D.f[SW])[ksw]=zero;
   //   }

   //   q = q_dirSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[NE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);// - c1over54 * drho;
   //      //(D.f[NE])[kne]=zero;
   //   }

   //   q = q_dirSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[NW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);// - c1over54 * drho;
   //      //(D.f[NW])[knw]=zero;
   //   }

   //   q = q_dirNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
   //      (D.f[SE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);// - c1over54 * drho;
   //      //(D.f[SE])[kse]=zero;
   //   }

   //   q = q_dirTE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
   //      (D.f[BW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[BW])[kbw]=zero;
   //   }

   //   q = q_dirBW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
   //      (D.f[TE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[TE])[kte]=zero;
   //   }

   //   q = q_dirBE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
   //      (D.f[TW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[TW])[ktw]=zero;
   //   }

   //   q = q_dirTW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
   //      (D.f[BE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[BE])[kbe]=zero;
   //   }

   //   q = q_dirTN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[BS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[BS])[kbs]=zero;
   //   }

   //   q = q_dirBS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[TN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[TN])[ktn]=zero;
   //   }

   //   q = q_dirBN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[TS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[TS])[kts]=zero;
   //   }

   //   q = q_dirTS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[BN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);// - c1over54 * drho;
   //      //(D.f[BN])[kbn]=zero;
   //   }

   //   q = q_dirTNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[BSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[BSW])[kbsw]=zero;
   //   }

   //   q = q_dirBSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[TNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[TNE])[ktne]=zero;
   //   }

   //   q = q_dirBNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[TSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[TSW])[ktsw]=zero;
   //   }

   //   q = q_dirTSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[BNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[BNE])[kbne]=zero;
   //   }

   //   q = q_dirTSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[BNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[BNW])[kbnw]=zero;
   //   }

   //   q = q_dirBNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[TSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[TSE])[ktse]=zero;
   //   }

   //   q = q_dirBSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
   //      (D.f[TNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[TNW])[ktnw]=zero;
   //   }

   //   q = q_dirTNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
   //      (D.f[BSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);// - c1over216 * drho;
   //      //(D.f[BSE])[kbse]=zero;
   //   }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////









