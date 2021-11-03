//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////

/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;

//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void CalcTurbulenceIntensity(
   real* vxx,
   real* vyy,
   real* vzz,
   real* vx_mean,
   real* vy_mean,
   real* vz_mean,
   real* DD, 
   uint* typeOfGridNode, 
   real om1, 
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

   if (k < size_Mat && (typeOfGridNode[k] == GEO_FLUID)) {
       ////////////////////////////////////////////////////////////////////////////////
       // index
       unsigned int KQK   = k;
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
       real f_E, f_W, f_N, f_S, f_T, f_B, f_NE, f_SW, f_SE, f_NW, f_TE, f_BW, f_BE, f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE,
           f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

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
       real vx, vy, vz, drho;
       drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW + f_BN + f_TS + f_TN + f_BS + f_BE + f_TW +
              f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]);

       vx = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
              ((f_BE - f_TW) + (f_TE - f_BW)) + ((f_SE - f_NW) + (f_NE - f_SW)) + (f_E - f_W)) /
             (c1o1 + drho);

       vy = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
              ((f_BN - f_TS) + (f_TN - f_BS)) + (-(f_SE - f_NW) + (f_NE - f_SW)) + (f_N - f_S)) /
             (c1o1 + drho);

       vz = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
              (-(f_BN - f_TS) + (f_TN - f_BS)) + ((f_TE - f_BW) - (f_BE - f_TW)) + (f_T - f_B)) /
             (c1o1 + drho);

       // compute subtotals:
       // fluctuations
       vxx[k] = vxx[k] + vx * vx;
       vyy[k] = vyy[k] + vy * vy;
       vzz[k] = vzz[k] + vz * vz;

       // velocity (for mean velocity)
       vx_mean[k] = vx_mean[k] + vx;
       vy_mean[k] = vy_mean[k] + vy;
       vz_mean[k] = vz_mean[k] + vz;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
