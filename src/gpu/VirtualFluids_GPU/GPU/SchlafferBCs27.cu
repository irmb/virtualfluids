/* Device code */
#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

// TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void PressSchlaff27(real* rhoBC,
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
                                          unsigned long long numberOfLBnodes,
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_E    = (D.f[dP00])[ke   ];
      f1_W    = (D.f[dM00])[kw   ];
      f1_N    = (D.f[d0P0])[kn   ];
      f1_S    = (D.f[d0M0])[ks   ];
      f1_T    = (D.f[d00P])[kt   ];
      f1_B    = (D.f[d00M])[kb   ];
      f1_NE   = (D.f[dPP0])[kne  ];
      f1_SW   = (D.f[dMM0])[ksw  ];
      f1_SE   = (D.f[dPM0])[kse  ];
      f1_NW   = (D.f[dMP0])[knw  ];
      f1_TE   = (D.f[dP0P])[kte  ];
      f1_BW   = (D.f[dM0M])[kbw  ];
      f1_BE   = (D.f[dP0M])[kbe  ];
      f1_TW   = (D.f[dM0P])[ktw  ];
      f1_TN   = (D.f[d0PP])[ktn  ];
      f1_BS   = (D.f[d0MM])[kbs  ];
      f1_BN   = (D.f[d0PM])[kbn  ];
      f1_TS   = (D.f[d0MP])[kts  ];
      f1_ZERO = (D.f[d000])[kzero];
      f1_TNE  = (D.f[dPPP])[ktne ];
      f1_TSW  = (D.f[dMMP])[ktsw ];
      f1_TSE  = (D.f[dPMP])[ktse ];
      f1_TNW  = (D.f[dMPP])[ktnw ];
      f1_BNE  = (D.f[dPPM])[kbne ];
      f1_BSW  = (D.f[dMMM])[kbsw ];
      f1_BSE  = (D.f[dPMM])[kbse ];
      f1_BNW  = (D.f[dMPM])[kbnw ];
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

      (D.f[d00M])[kb   ] = f1_B   ;
      (D.f[dM0M])[kbw  ] = f1_BW  ;
      (D.f[dP0M])[kbe  ] = f1_BE  ;
      (D.f[d0MM])[kbs  ] = f1_BS  ;
      (D.f[d0PM])[kbn  ] = f1_BN  ;
      (D.f[dPPM])[kbne ] = f1_BNE ;
      (D.f[dMMM])[kbsw ] = f1_BSW ;
      (D.f[dPMM])[kbse ] = f1_BSE ;
      (D.f[dMPM])[kbnw ] = f1_BNW ;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






































// TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void VelSchlaff27(  int t,
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
                                          unsigned long long numberOfLBnodes,
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                     f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_E    = (D.f[dP00])[ke   ];
      f1_W    = (D.f[dM00])[kw   ];
      f1_N    = (D.f[d0P0])[kn   ];
      f1_S    = (D.f[d0M0])[ks   ];
      f1_T    = (D.f[d00P])[kt   ];
      f1_B    = (D.f[d00M])[kb   ];
      f1_NE   = (D.f[dPP0])[kne  ];
      f1_SW   = (D.f[dMM0])[ksw  ];
      f1_SE   = (D.f[dPM0])[kse  ];
      f1_NW   = (D.f[dMP0])[knw  ];
      f1_TE   = (D.f[dP0P])[kte  ];
      f1_BW   = (D.f[dM0M])[kbw  ];
      f1_BE   = (D.f[dP0M])[kbe  ];
      f1_TW   = (D.f[dM0P])[ktw  ];
      f1_TN   = (D.f[d0PP])[ktn  ];
      f1_BS   = (D.f[d0MM])[kbs  ];
      f1_BN   = (D.f[d0PM])[kbn  ];
      f1_TS   = (D.f[d0MP])[kts  ];
      f1_ZERO = (D.f[d000])[kzero];
      f1_TNE  = (D.f[dPPP])[ktne ];
      f1_TSW  = (D.f[dMMP])[ktsw ];
      f1_TSE  = (D.f[dPMP])[ktse ];
      f1_TNW  = (D.f[dMPP])[ktnw ];
      f1_BNE  = (D.f[dPPM])[kbne ];
      f1_BSW  = (D.f[dMMM])[kbsw ];
      f1_BSE  = (D.f[dPMM])[kbse ];
      f1_BNW  = (D.f[dMPM])[kbnw ];
      //f1_W    = (D.f[dP00])[ke   ];
      //f1_E    = (D.f[dM00])[kw   ];
      //f1_S    = (D.f[d0P0])[kn   ];
      //f1_N    = (D.f[d0M0])[ks   ];
      //f1_B    = (D.f[d00P])[kt   ];
      //f1_T    = (D.f[d00M])[kb   ];
      //f1_SW   = (D.f[dPP0])[kne  ];
      //f1_NE   = (D.f[dMM0])[ksw  ];
      //f1_NW   = (D.f[dPM0])[kse  ];
      //f1_SE   = (D.f[dMP0])[knw  ];
      //f1_BW   = (D.f[dP0P])[kte  ];
      //f1_TE   = (D.f[dM0M])[kbw  ];
      //f1_TW   = (D.f[dP0M])[kbe  ];
      //f1_BE   = (D.f[dM0P])[ktw  ];
      //f1_BS   = (D.f[d0PP])[ktn  ];
      //f1_TN   = (D.f[d0MM])[kbs  ];
      //f1_TS   = (D.f[d0PM])[kbn  ];
      //f1_BN   = (D.f[d0MP])[kts  ];
      //f1_ZERO = (D.f[d000])[kzero];
      //f1_BSW  = (D.f[dPPP])[ktne ];
      //f1_BNE  = (D.f[dMMP])[ktsw ];
      //f1_BNW  = (D.f[dPMP])[ktse ];
      //f1_BSE  = (D.f[dMPP])[ktnw ];
      //f1_TSW  = (D.f[dPPM])[kbne ];
      //f1_TNE  = (D.f[dMMM])[kbsw ];
      //f1_TNW  = (D.f[dPMM])[kbse ];
      //f1_TSE  = (D.f[dMPM])[kbnw ];
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
      (D.f[d00P])[kt   ] = f1_T  ;
      (D.f[dP0P])[kte  ] = f1_TE ;
      (D.f[dM0P])[ktw  ] = f1_TW ;
      (D.f[d0PP])[ktn  ] = f1_TN ;
      (D.f[d0MP])[kts  ] = f1_TS ;
      (D.f[dPPP])[ktne ] = f1_TNE;
      (D.f[dMMP])[ktsw ] = f1_TSW;
      (D.f[dPMP])[ktse ] = f1_TSE;
      (D.f[dMPP])[ktnw ] = f1_TNW;

      //(D.f[d00M])[kb   ] = f1_B   ;
      //(D.f[dM0M])[kbw  ] = f1_BW  ;
      //(D.f[dP0M])[kbe  ] = f1_BE  ;
      //(D.f[d0MM])[kbs  ] = f1_BS  ;
      //(D.f[d0PM])[kbn  ] = f1_BN  ;
      //(D.f[dPPM])[kbne ] = f1_BNE ;
      //(D.f[dMMM])[kbsw ] = f1_BSW ;
      //(D.f[dPMM])[kbse ] = f1_BSE ;
      //(D.f[dMPM])[kbnw ] = f1_BNW ;


      //(D.f[d00P])[kt   ] = f1_B  ;
      //(D.f[dP0P])[kte  ] = f1_BW ;
      //(D.f[dM0P])[ktw  ] = f1_BE ;
      //(D.f[d0PP])[ktn  ] = f1_BS ;
      //(D.f[d0MP])[kts  ] = f1_BN ;
      //(D.f[dPPP])[ktne ] = f1_BSW;
      //(D.f[dMMP])[ktsw ] = f1_BNE;
      //(D.f[dPMP])[ktse ] = f1_BNW;
      //(D.f[dMPP])[ktnw ] = f1_BSE;

      //(D.f[dP00])[ke   ] = f1_W   -c2over27*drho1;
      //(D.f[dM00])[kw   ] = f1_E   -c2over27*drho1;
      //(D.f[d0P0])[kn   ] = f1_S   -c2over27*drho1;
      //(D.f[d0M0])[ks   ] = f1_N   -c2over27*drho1;
      //(D.f[d00P])[kt   ] = f1_B   -c2over27*drho1;
      //(D.f[d00M])[kb   ] = f1_T   -c2over27*drho1;
      //(D.f[dPP0])[kne  ] = f1_SW  -c1over54*drho1;
      //(D.f[dMM0])[ksw  ] = f1_NE  -c1over54*drho1;
      //(D.f[dPM0])[kse  ] = f1_NW  -c1over54*drho1;
      //(D.f[dMP0])[knw  ] = f1_SE  -c1over54*drho1;
      //(D.f[dP0P])[kte  ] = f1_BW  -c1over54*drho1;
      //(D.f[dM0M])[kbw  ] = f1_TE  -c1over54*drho1;
      //(D.f[dP0M])[kbe  ] = f1_TW  -c1over54*drho1;
      //(D.f[dM0P])[ktw  ] = f1_BE  -c1over54*drho1;
      //(D.f[d0PP])[ktn  ] = f1_BS  -c1over54*drho1;
      //(D.f[d0MM])[kbs  ] = f1_TN  -c1over54*drho1;
      //(D.f[d0PM])[kbn  ] = f1_TS  -c1over54*drho1;
      //(D.f[d0MP])[kts  ] = f1_BN  -c1over54*drho1;
      //(D.f[d000])[kzero] = f1_ZERO-c8over27*drho1;
      //(D.f[dPPP])[ktne ] = f1_BSW -c1over216*drho1;
      //(D.f[dMMP])[ktsw ] = f1_BNE -c1over216*drho1;
      //(D.f[dPMP])[ktse ] = f1_BNW -c1over216*drho1;
      //(D.f[dMPP])[ktnw ] = f1_BSE -c1over216*drho1;
      //(D.f[dPPM])[kbne ] = f1_TSW -c1over216*drho1;
      //(D.f[dMMM])[kbsw ] = f1_TNE -c1over216*drho1;
      //(D.f[dPMM])[kbse ] = f1_TNW -c1over216*drho1;
      //(D.f[dMPM])[kbnw ] = f1_TSE -c1over216*drho1;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





