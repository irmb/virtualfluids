/* Device code */
#include "LBM/LB.h"
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;

// TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void PressSchlaff27(real* rhoBC,
                                          real* DD,
                                          real* vx0,
                                          real* vy0,
                                          real* vz0,
                                          real* deltaVz0,
                                          int* k_Q,
                                          int* k_N,
                                          int numberOfBCnodes,
                                          real om1,
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
                                          unsigned int size_Mat,
                                          bool isEvenTimestep)
{
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
         D.f[REST] = &DD[REST*size_Mat];
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
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_E    = (D.f[E   ])[ke   ];
      f1_W    = (D.f[W   ])[kw   ];
      f1_N    = (D.f[N   ])[kn   ];
      f1_S    = (D.f[S   ])[ks   ];
      f1_T    = (D.f[T   ])[kt   ];
      f1_B    = (D.f[B   ])[kb   ];
      f1_NE   = (D.f[NE  ])[kne  ];
      f1_SW   = (D.f[SW  ])[ksw  ];
      f1_SE   = (D.f[SE  ])[kse  ];
      f1_NW   = (D.f[NW  ])[knw  ];
      f1_TE   = (D.f[TE  ])[kte  ];
      f1_BW   = (D.f[BW  ])[kbw  ];
      f1_BE   = (D.f[BE  ])[kbe  ];
      f1_TW   = (D.f[TW  ])[ktw  ];
      f1_TN   = (D.f[TN  ])[ktn  ];
      f1_BS   = (D.f[BS  ])[kbs  ];
      f1_BN   = (D.f[BN  ])[kbn  ];
      f1_TS   = (D.f[TS  ])[kts  ];
      f1_ZERO = (D.f[REST])[kzero];
      f1_TNE  = (D.f[TNE ])[ktne ];
      f1_TSW  = (D.f[TSW ])[ktsw ];
      f1_TSE  = (D.f[TSE ])[ktse ];
      f1_TNW  = (D.f[TNW ])[ktnw ];
      f1_BNE  = (D.f[BNE ])[kbne ];
      f1_BSW  = (D.f[BSW ])[kbsw ];
      f1_BSE  = (D.f[BSE ])[kbse ];
      f1_BNW  = (D.f[BNW ])[kbnw ];
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

      (D.f[B   ])[kb   ] = f1_B   ;
      (D.f[BW  ])[kbw  ] = f1_BW  ;
      (D.f[BE  ])[kbe  ] = f1_BE  ;
      (D.f[BS  ])[kbs  ] = f1_BS  ;
      (D.f[BN  ])[kbn  ] = f1_BN  ;
      (D.f[BNE ])[kbne ] = f1_BNE ;
      (D.f[BSW ])[kbsw ] = f1_BSW ;
      (D.f[BSE ])[kbse ] = f1_BSE ;
      (D.f[BNW ])[kbnw ] = f1_BNW ;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






































// TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void VelSchlaff27(  int t,
                                          real* DD,
                                          real* vz0,
                                          real* deltaVz0,
                                          int* k_Q,
                                          int* k_N,
                                          int numberOfBCnodes,
                                          real om1,
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
                                          unsigned int size_Mat,
                                          bool isEvenTimestep)
{
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
         D.f[REST] = &DD[REST*size_Mat];
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
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                     f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_E    = (D.f[E   ])[ke   ];
      f1_W    = (D.f[W   ])[kw   ];
      f1_N    = (D.f[N   ])[kn   ];
      f1_S    = (D.f[S   ])[ks   ];
      f1_T    = (D.f[T   ])[kt   ];
      f1_B    = (D.f[B   ])[kb   ];
      f1_NE   = (D.f[NE  ])[kne  ];
      f1_SW   = (D.f[SW  ])[ksw  ];
      f1_SE   = (D.f[SE  ])[kse  ];
      f1_NW   = (D.f[NW  ])[knw  ];
      f1_TE   = (D.f[TE  ])[kte  ];
      f1_BW   = (D.f[BW  ])[kbw  ];
      f1_BE   = (D.f[BE  ])[kbe  ];
      f1_TW   = (D.f[TW  ])[ktw  ];
      f1_TN   = (D.f[TN  ])[ktn  ];
      f1_BS   = (D.f[BS  ])[kbs  ];
      f1_BN   = (D.f[BN  ])[kbn  ];
      f1_TS   = (D.f[TS  ])[kts  ];
      f1_ZERO = (D.f[REST])[kzero];
      f1_TNE  = (D.f[TNE ])[ktne ];
      f1_TSW  = (D.f[TSW ])[ktsw ];
      f1_TSE  = (D.f[TSE ])[ktse ];
      f1_TNW  = (D.f[TNW ])[ktnw ];
      f1_BNE  = (D.f[BNE ])[kbne ];
      f1_BSW  = (D.f[BSW ])[kbsw ];
      f1_BSE  = (D.f[BSE ])[kbse ];
      f1_BNW  = (D.f[BNW ])[kbnw ];
      //f1_W    = (D.f[E   ])[ke   ];
      //f1_E    = (D.f[W   ])[kw   ];
      //f1_S    = (D.f[N   ])[kn   ];
      //f1_N    = (D.f[S   ])[ks   ];
      //f1_B    = (D.f[T   ])[kt   ];
      //f1_T    = (D.f[B   ])[kb   ];
      //f1_SW   = (D.f[NE  ])[kne  ];
      //f1_NE   = (D.f[SW  ])[ksw  ];
      //f1_NW   = (D.f[SE  ])[kse  ];
      //f1_SE   = (D.f[NW  ])[knw  ];
      //f1_BW   = (D.f[TE  ])[kte  ];
      //f1_TE   = (D.f[BW  ])[kbw  ];
      //f1_TW   = (D.f[BE  ])[kbe  ];
      //f1_BE   = (D.f[TW  ])[ktw  ];
      //f1_BS   = (D.f[TN  ])[ktn  ];
      //f1_TN   = (D.f[BS  ])[kbs  ];
      //f1_TS   = (D.f[BN  ])[kbn  ];
      //f1_BN   = (D.f[TS  ])[kts  ];
      //f1_ZERO = (D.f[REST])[kzero];
      //f1_BSW  = (D.f[TNE ])[ktne ];
      //f1_BNE  = (D.f[TSW ])[ktsw ];
      //f1_BNW  = (D.f[TSE ])[ktse ];
      //f1_BSE  = (D.f[TNW ])[ktnw ];
      //f1_TSW  = (D.f[BNE ])[kbne ];
      //f1_TNE  = (D.f[BSW ])[kbsw ];
      //f1_TNW  = (D.f[BSE ])[kbse ];
      //f1_TSE  = (D.f[BNW ])[kbnw ];
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
      (D.f[T   ])[kt   ] = f1_T  ;
      (D.f[TE  ])[kte  ] = f1_TE ;
      (D.f[TW  ])[ktw  ] = f1_TW ;
      (D.f[TN  ])[ktn  ] = f1_TN ;
      (D.f[TS  ])[kts  ] = f1_TS ;
      (D.f[TNE ])[ktne ] = f1_TNE;
      (D.f[TSW ])[ktsw ] = f1_TSW;
      (D.f[TSE ])[ktse ] = f1_TSE;
      (D.f[TNW ])[ktnw ] = f1_TNW;

      //(D.f[B   ])[kb   ] = f1_B   ;
      //(D.f[BW  ])[kbw  ] = f1_BW  ;
      //(D.f[BE  ])[kbe  ] = f1_BE  ;
      //(D.f[BS  ])[kbs  ] = f1_BS  ;
      //(D.f[BN  ])[kbn  ] = f1_BN  ;
      //(D.f[BNE ])[kbne ] = f1_BNE ;
      //(D.f[BSW ])[kbsw ] = f1_BSW ;
      //(D.f[BSE ])[kbse ] = f1_BSE ;
      //(D.f[BNW ])[kbnw ] = f1_BNW ;


      //(D.f[T   ])[kt   ] = f1_B  ;
      //(D.f[TE  ])[kte  ] = f1_BW ;
      //(D.f[TW  ])[ktw  ] = f1_BE ;
      //(D.f[TN  ])[ktn  ] = f1_BS ;
      //(D.f[TS  ])[kts  ] = f1_BN ;
      //(D.f[TNE ])[ktne ] = f1_BSW;
      //(D.f[TSW ])[ktsw ] = f1_BNE;
      //(D.f[TSE ])[ktse ] = f1_BNW;
      //(D.f[TNW ])[ktnw ] = f1_BSE;

      //(D.f[E   ])[ke   ] = f1_W   -c2over27*drho1;
      //(D.f[W   ])[kw   ] = f1_E   -c2over27*drho1;
      //(D.f[N   ])[kn   ] = f1_S   -c2over27*drho1;
      //(D.f[S   ])[ks   ] = f1_N   -c2over27*drho1;
      //(D.f[T   ])[kt   ] = f1_B   -c2over27*drho1;
      //(D.f[B   ])[kb   ] = f1_T   -c2over27*drho1;
      //(D.f[NE  ])[kne  ] = f1_SW  -c1over54*drho1;
      //(D.f[SW  ])[ksw  ] = f1_NE  -c1over54*drho1;
      //(D.f[SE  ])[kse  ] = f1_NW  -c1over54*drho1;
      //(D.f[NW  ])[knw  ] = f1_SE  -c1over54*drho1;
      //(D.f[TE  ])[kte  ] = f1_BW  -c1over54*drho1;
      //(D.f[BW  ])[kbw  ] = f1_TE  -c1over54*drho1;
      //(D.f[BE  ])[kbe  ] = f1_TW  -c1over54*drho1;
      //(D.f[TW  ])[ktw  ] = f1_BE  -c1over54*drho1;
      //(D.f[TN  ])[ktn  ] = f1_BS  -c1over54*drho1;
      //(D.f[BS  ])[kbs  ] = f1_TN  -c1over54*drho1;
      //(D.f[BN  ])[kbn  ] = f1_TS  -c1over54*drho1;
      //(D.f[TS  ])[kts  ] = f1_BN  -c1over54*drho1;
      //(D.f[REST])[kzero] = f1_ZERO-c8over27*drho1;
      //(D.f[TNE ])[ktne ] = f1_BSW -c1over216*drho1;
      //(D.f[TSW ])[ktsw ] = f1_BNE -c1over216*drho1;
      //(D.f[TSE ])[ktse ] = f1_BNW -c1over216*drho1;
      //(D.f[TNW ])[ktnw ] = f1_BSE -c1over216*drho1;
      //(D.f[BNE ])[kbne ] = f1_TSW -c1over216*drho1;
      //(D.f[BSW ])[kbsw ] = f1_TNE -c1over216*drho1;
      //(D.f[BSE ])[kbse ] = f1_TNW -c1over216*drho1;
      //(D.f[BNW ])[kbnw ] = f1_TSE -c1over216*drho1;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





