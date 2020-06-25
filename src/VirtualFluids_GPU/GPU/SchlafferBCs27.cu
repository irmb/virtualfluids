/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include "Core/RealConstants.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void PressSchlaff27(real* rhoBC,
                                          real* DD,
                                          real* vx0,
                                          real* vy0,
                                          real* vz0,
                                          real* deltaVz0,
                                          int* k_Q, 
                                          int* k_N, 
                                          int kQ, 
                                          real om1, 
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_E    = (D.f[dirE   ])[ke   ];
      f1_W    = (D.f[dirW   ])[kw   ];
      f1_N    = (D.f[dirN   ])[kn   ];
      f1_S    = (D.f[dirS   ])[ks   ];
      f1_T    = (D.f[dirT   ])[kt   ];
      f1_B    = (D.f[dirB   ])[kb   ];
      f1_NE   = (D.f[dirNE  ])[kne  ];
      f1_SW   = (D.f[dirSW  ])[ksw  ];
      f1_SE   = (D.f[dirSE  ])[kse  ];
      f1_NW   = (D.f[dirNW  ])[knw  ];
      f1_TE   = (D.f[dirTE  ])[kte  ];
      f1_BW   = (D.f[dirBW  ])[kbw  ];
      f1_BE   = (D.f[dirBE  ])[kbe  ];
      f1_TW   = (D.f[dirTW  ])[ktw  ];
      f1_TN   = (D.f[dirTN  ])[ktn  ];
      f1_BS   = (D.f[dirBS  ])[kbs  ];
      f1_BN   = (D.f[dirBN  ])[kbn  ];
      f1_TS   = (D.f[dirTS  ])[kts  ];
      f1_ZERO = (D.f[dirZERO])[kzero];
      f1_TNE  = (D.f[dirTNE ])[ktne ];
      f1_TSW  = (D.f[dirTSW ])[ktsw ];
      f1_TSE  = (D.f[dirTSE ])[ktse ];
      f1_TNW  = (D.f[dirTNW ])[ktnw ];
      f1_BNE  = (D.f[dirBNE ])[kbne ];
      f1_BSW  = (D.f[dirBSW ])[kbsw ];
      f1_BSE  = (D.f[dirBSE ])[kbse ];
      f1_BNW  = (D.f[dirBNW ])[kbnw ];
      //////////////////////////////////////////////////////////////////////////
      real cs       = c1o1/sqrt(c3o1);
      real csp1     = cs + c1o1;
      real csp1Sq  = (c1o1 + cs)*(c1o1 + cs);
      real relFac   = c21o20; // 0.9...1.0
      //////////////////////////////////////////////////////////////////////////
      // For adaption:
      //     Pressure limits with rho0 = 1:
      //      2.2e-10 ~  0.94 dB
      //      6.2e-10 ~   9.9 dB
      //      6.2e-9  ~  29.9 dB
      //      2.0e-7  ~  60.1 dB   /Vel
      //      2.0e-5  ~ 100.1 dB   /press
      const double dPlimit  = Op0000002;
      const double dRlimit  = dPlimit * c3o1;// three = c1oCs2;
      const double uSlimit  = dRlimit * c1o1;// one = c1oRho0;
      //////////////////////////////////////////////////////////////////////////
      real VX = vx0[k];
      real VY = vy0[k];
      real VZ = vz0[k];
      //////////////////////////////////////////////////////////////////////////

      real temp = c2o1*(f1_TNE + f1_TSE + f1_TSW + f1_TNW) + c2o1*(f1_TE + f1_TW + f1_TN + f1_TS) + f1_NE + f1_SW + f1_SE + f1_NW + c2o1*f1_T + f1_E + f1_W + f1_N + f1_S + f1_ZERO;

      real vs_z = relFac * (VZ+cs) * ( csp1 - sqrt(csp1Sq + c2o1*VZ - c2o1*temp) );    //old =  relFac * cs * ( csp1 - sqrt(csp1Sq + two*VZ - two*temp) );

      // 3. Compute density of compensated velocity:		
      real tempDeltaV = deltaVz0[k];
      real rholoc = temp - c1o1 * (VZ + tempDeltaV + vs_z);

      // 4. Compute density deviation:
      real drho = rholoc - rhoBC[k];

      // 5. Adapt Speed:
      real dv = tempDeltaV + vs_z;

      if( drho > dRlimit) {
         VZ += dv + uSlimit;
         tempDeltaV += uSlimit;
      }
      else if( drho < -dRlimit) {
         VZ += dv - uSlimit;
         tempDeltaV -= uSlimit;
      }
      else {
         VZ += dv + drho;
         tempDeltaV += drho;
      }

      //VZ = vz0[k] + vs_z;
      // 6. Set unknown distributions:
      f1_B   = f1_T   - c4o9  * VZ;
      f1_BW  = f1_TE  - c1o9  * (VX + VZ);
      f1_BE  = f1_TW  + c1o9  * (VX - VZ);
      f1_BS  = f1_TN  - c1o9  * (VY + VZ);
      f1_BN  = f1_TS  + c1o9  * (VY - VZ);
      f1_BSW = f1_TNE - c1o36 * (VX + VY + VZ);
      f1_BNW = f1_TSE - c1o36 * (VX - VY + VZ);
      f1_BNE = f1_TSW + c1o36 * (VX + VY - VZ);
      f1_BSE = f1_TNW + c1o36 * (VX - VY - VZ);

      deltaVz0[k] = tempDeltaV;

      (D.f[dirB   ])[kb   ] = f1_B   ;
      (D.f[dirBW  ])[kbw  ] = f1_BW  ;
      (D.f[dirBE  ])[kbe  ] = f1_BE  ;
      (D.f[dirBS  ])[kbs  ] = f1_BS  ;
      (D.f[dirBN  ])[kbn  ] = f1_BN  ;
      (D.f[dirBNE ])[kbne ] = f1_BNE ;
      (D.f[dirBSW ])[kbsw ] = f1_BSW ;
      (D.f[dirBSE ])[kbse ] = f1_BSE ;
      (D.f[dirBNW ])[kbnw ] = f1_BNW ;       
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void VelSchlaff27(  int t,
                                          real* DD,
                                          real* vz0,
                                          real* deltaVz0,
                                          int* k_Q, 
                                          int* k_N, 
                                          int kQ, 
                                          real om1, 
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                     f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_E    = (D.f[dirE   ])[ke   ];
      f1_W    = (D.f[dirW   ])[kw   ];
      f1_N    = (D.f[dirN   ])[kn   ];
      f1_S    = (D.f[dirS   ])[ks   ];
      f1_T    = (D.f[dirT   ])[kt   ];
      f1_B    = (D.f[dirB   ])[kb   ];
      f1_NE   = (D.f[dirNE  ])[kne  ];
      f1_SW   = (D.f[dirSW  ])[ksw  ];
      f1_SE   = (D.f[dirSE  ])[kse  ];
      f1_NW   = (D.f[dirNW  ])[knw  ];
      f1_TE   = (D.f[dirTE  ])[kte  ];
      f1_BW   = (D.f[dirBW  ])[kbw  ];
      f1_BE   = (D.f[dirBE  ])[kbe  ];
      f1_TW   = (D.f[dirTW  ])[ktw  ];
      f1_TN   = (D.f[dirTN  ])[ktn  ];
      f1_BS   = (D.f[dirBS  ])[kbs  ];
      f1_BN   = (D.f[dirBN  ])[kbn  ];
      f1_TS   = (D.f[dirTS  ])[kts  ];
      f1_ZERO = (D.f[dirZERO])[kzero];
      f1_TNE  = (D.f[dirTNE ])[ktne ];
      f1_TSW  = (D.f[dirTSW ])[ktsw ];
      f1_TSE  = (D.f[dirTSE ])[ktse ];
      f1_TNW  = (D.f[dirTNW ])[ktnw ];
      f1_BNE  = (D.f[dirBNE ])[kbne ];
      f1_BSW  = (D.f[dirBSW ])[kbsw ];
      f1_BSE  = (D.f[dirBSE ])[kbse ];
      f1_BNW  = (D.f[dirBNW ])[kbnw ];
      //f1_W    = (D.f[dirE   ])[ke   ];
      //f1_E    = (D.f[dirW   ])[kw   ];
      //f1_S    = (D.f[dirN   ])[kn   ];
      //f1_N    = (D.f[dirS   ])[ks   ];
      //f1_B    = (D.f[dirT   ])[kt   ];
      //f1_T    = (D.f[dirB   ])[kb   ];
      //f1_SW   = (D.f[dirNE  ])[kne  ];
      //f1_NE   = (D.f[dirSW  ])[ksw  ];
      //f1_NW   = (D.f[dirSE  ])[kse  ];
      //f1_SE   = (D.f[dirNW  ])[knw  ];
      //f1_BW   = (D.f[dirTE  ])[kte  ];
      //f1_TE   = (D.f[dirBW  ])[kbw  ];
      //f1_TW   = (D.f[dirBE  ])[kbe  ];
      //f1_BE   = (D.f[dirTW  ])[ktw  ];
      //f1_BS   = (D.f[dirTN  ])[ktn  ];
      //f1_TN   = (D.f[dirBS  ])[kbs  ];
      //f1_TS   = (D.f[dirBN  ])[kbn  ];
      //f1_BN   = (D.f[dirTS  ])[kts  ];
      //f1_ZERO = (D.f[dirZERO])[kzero];
      //f1_BSW  = (D.f[dirTNE ])[ktne ];
      //f1_BNE  = (D.f[dirTSW ])[ktsw ];
      //f1_BNW  = (D.f[dirTSE ])[ktse ];
      //f1_BSE  = (D.f[dirTNW ])[ktnw ];
      //f1_TSW  = (D.f[dirBNE ])[kbne ];
      //f1_TNE  = (D.f[dirBSW ])[kbsw ];
      //f1_TNW  = (D.f[dirBSE ])[kbse ];
      //f1_TSE  = (D.f[dirBNW ])[kbnw ];
      //////////////////////////////////////////////////////////////////////////
      real cs       = c1o1/sqrt(c3o1);
      real csp1     = cs + c1o1;
      real csp1Sq  = (c1o1 + cs)*(c1o1 + cs);
      real relFac   = c19o20; // 0.9...1.0
      //////////////////////////////////////////////////////////////////////////
      // For adaption:
      //     Pressure limits with rho0 = 1:
      //      2.2e-10 ~  0.94 dB
      //      6.2e-10 ~   9.9 dB
      //      6.2e-9  ~  29.9 dB
      //      2.0e-7  ~  60.1 dB   /Vel
      //      2.0e-5  ~ 100.1 dB   /press
      real uSlimit  = Op0000002; 
      //////////////////////////////////////////////////////////////////////////
      real VX = c0o1;
      real VY = c0o1;
      real VZ = vz0[k];
      //////////////////////////////////////////////////////////////////////////
      real temp = f1_ZERO + f1_E + f1_W + f1_N + f1_S + f1_NE + f1_SW + f1_SE + f1_NW + c2o1*(f1_B + f1_BE + f1_BW + f1_BN + f1_BS + f1_BNE + f1_BSE + f1_BSW + f1_BNW);
      //real temp = f1_ZERO + f1_E + f1_W + f1_N + f1_S + f1_NE + f1_SW + f1_SE + f1_NW + two*(f1_T + f1_TE + f1_TW + f1_TN + f1_TS + f1_TNE + f1_TSE + f1_TSW + f1_TNW);
      ////real temp2= c1mcsSq + two*VZ - two*temp;
      real vs_z;
      //if (t < 5)
      //{
      //   vs_z = zero;
      //}
      //else
      //{
         vs_z = relFac * (cs-VZ) * ( sqrt(csp1Sq - c2o1*VZ - c2o1*temp) - csp1 );         //old = relFac * cs * ( sqrt(csp1Sq - two*VZ - two*temp) - csp1 ); 
      //}
      
      // 3. Adapt Speed:
      real tempDeltaV = deltaVz0[k];
      real dv = tempDeltaV + vs_z;

      if( dv > uSlimit) {
         VZ  += dv - uSlimit;
         tempDeltaV -= uSlimit;
      }
      else if( dv < -uSlimit) {
         VZ  += dv + uSlimit;
         tempDeltaV += uSlimit;
      }
      else {
         tempDeltaV = -vs_z;
      }

      //VZ = vz0[k]+vs_z;
      // 4. Set unknown distributions:
      //f1_B   = f1_T   - c4o9  * VZ;
      //f1_BW  = f1_TE  - c1o9  * (VX + VZ);
      //f1_BE  = f1_TW  + c1o9  * (VX - VZ);
      //f1_BS  = f1_TN  - c1o9  * (VY + VZ);
      //f1_BN  = f1_TS  + c1o9  * (VY - VZ);
      //f1_BSW = f1_TNE - c1o36 * (VX + VY + VZ);
      //f1_BNW = f1_TSE - c1o36 * (VX - VY + VZ);
      //f1_BNE = f1_TSW + c1o36 * (VX + VY - VZ);
      //f1_BSE = f1_TNW + c1o36 * (VX - VY - VZ);

      f1_T   = f1_B   + c4o9  * VZ;
      f1_TE  = f1_BW  + c1o9  * (VX + VZ);
      f1_TW  = f1_BE  - c1o9  * (VX - VZ);
      f1_TN  = f1_BS  + c1o9  * (VY + VZ);
      f1_TS  = f1_BN  - c1o9  * (VY - VZ);
      f1_TNE = f1_BSW + c1o36 * (VX + VY + VZ);
      f1_TSE = f1_BNW + c1o36 * (VX - VY + VZ);
      f1_TSW = f1_BNE - c1o36 * (VX + VY - VZ);
      f1_TNW = f1_BSE - c1o36 * (VX - VY - VZ);

      deltaVz0[k] = tempDeltaV;
      (D.f[dirT   ])[kt   ] = f1_T  ;
      (D.f[dirTE  ])[kte  ] = f1_TE ;
      (D.f[dirTW  ])[ktw  ] = f1_TW ;
      (D.f[dirTN  ])[ktn  ] = f1_TN ;
      (D.f[dirTS  ])[kts  ] = f1_TS ;
      (D.f[dirTNE ])[ktne ] = f1_TNE;
      (D.f[dirTSW ])[ktsw ] = f1_TSW;
      (D.f[dirTSE ])[ktse ] = f1_TSE;
      (D.f[dirTNW ])[ktnw ] = f1_TNW;

      //(D.f[dirB   ])[kb   ] = f1_B   ;
      //(D.f[dirBW  ])[kbw  ] = f1_BW  ;
      //(D.f[dirBE  ])[kbe  ] = f1_BE  ;
      //(D.f[dirBS  ])[kbs  ] = f1_BS  ;
      //(D.f[dirBN  ])[kbn  ] = f1_BN  ;
      //(D.f[dirBNE ])[kbne ] = f1_BNE ;
      //(D.f[dirBSW ])[kbsw ] = f1_BSW ;
      //(D.f[dirBSE ])[kbse ] = f1_BSE ;
      //(D.f[dirBNW ])[kbnw ] = f1_BNW ;       


      //(D.f[dirT   ])[kt   ] = f1_B  ;
      //(D.f[dirTE  ])[kte  ] = f1_BW ;
      //(D.f[dirTW  ])[ktw  ] = f1_BE ;
      //(D.f[dirTN  ])[ktn  ] = f1_BS ;
      //(D.f[dirTS  ])[kts  ] = f1_BN ;
      //(D.f[dirTNE ])[ktne ] = f1_BSW;
      //(D.f[dirTSW ])[ktsw ] = f1_BNE;
      //(D.f[dirTSE ])[ktse ] = f1_BNW;
      //(D.f[dirTNW ])[ktnw ] = f1_BSE;

      //(D.f[dirE   ])[ke   ] = f1_W   -c2over27*drho1;
      //(D.f[dirW   ])[kw   ] = f1_E   -c2over27*drho1;
      //(D.f[dirN   ])[kn   ] = f1_S   -c2over27*drho1;
      //(D.f[dirS   ])[ks   ] = f1_N   -c2over27*drho1;
      //(D.f[dirT   ])[kt   ] = f1_B   -c2over27*drho1;
      //(D.f[dirB   ])[kb   ] = f1_T   -c2over27*drho1;
      //(D.f[dirNE  ])[kne  ] = f1_SW  -c1over54*drho1;
      //(D.f[dirSW  ])[ksw  ] = f1_NE  -c1over54*drho1;
      //(D.f[dirSE  ])[kse  ] = f1_NW  -c1over54*drho1;
      //(D.f[dirNW  ])[knw  ] = f1_SE  -c1over54*drho1;
      //(D.f[dirTE  ])[kte  ] = f1_BW  -c1over54*drho1;
      //(D.f[dirBW  ])[kbw  ] = f1_TE  -c1over54*drho1;
      //(D.f[dirBE  ])[kbe  ] = f1_TW  -c1over54*drho1;
      //(D.f[dirTW  ])[ktw  ] = f1_BE  -c1over54*drho1;
      //(D.f[dirTN  ])[ktn  ] = f1_BS  -c1over54*drho1;
      //(D.f[dirBS  ])[kbs  ] = f1_TN  -c1over54*drho1;
      //(D.f[dirBN  ])[kbn  ] = f1_TS  -c1over54*drho1;
      //(D.f[dirTS  ])[kts  ] = f1_BN  -c1over54*drho1;
      //(D.f[dirZERO])[kzero] = f1_ZERO-c8over27*drho1;
      //(D.f[dirTNE ])[ktne ] = f1_BSW -c1over216*drho1;
      //(D.f[dirTSW ])[ktsw ] = f1_BNE -c1over216*drho1;
      //(D.f[dirTSE ])[ktse ] = f1_BNW -c1over216*drho1;
      //(D.f[dirTNW ])[ktnw ] = f1_BSE -c1over216*drho1;
      //(D.f[dirBNE ])[kbne ] = f1_TSW -c1over216*drho1;
      //(D.f[dirBSW ])[kbsw ] = f1_TNE -c1over216*drho1;
      //(D.f[dirBSE ])[kbse ] = f1_TNW -c1over216*drho1;
      //(D.f[dirBNW ])[kbnw ] = f1_TSE -c1over216*drho1;       
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





