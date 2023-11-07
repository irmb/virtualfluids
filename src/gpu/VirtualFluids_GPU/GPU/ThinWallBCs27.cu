//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////

/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

/////////////////////////////////////////////////////////////////////////
__global__ void QVelDeviceCompThinWallsPartOne27(
	real* vx,
	real* vy,
	real* vz,
	real* DD, 
	int* k_Q, 
	real* QQ,
	uint numberOfBCnodes, 
	real om1, 
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	unsigned long long numberOfLBnodes, 
	bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[dP00] = &DD[dP00 * numberOfLBnodes];
      D.f[dM00] = &DD[dM00 * numberOfLBnodes];
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
      D.f[d000] = &DD[d000 * numberOfLBnodes];
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
      D.f[dM00] = &DD[dP00 * numberOfLBnodes];
      D.f[dP00] = &DD[dM00 * numberOfLBnodes];
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
      D.f[d000] = &DD[d000 * numberOfLBnodes];
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
      real VeloY = vy[k];
      real VeloZ = vz[k]; 
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dP00 * numberOfBCnodes];
      q_dirW   = &QQ[dM00 * numberOfBCnodes];
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
      uint KQK  = k_Q[k];
      uint kzero= KQK;
      uint ke   = KQK;
      uint kw   = neighborX[KQK];
      uint kn   = KQK;
      uint ks   = neighborY[KQK];
      uint kt   = KQK;
      uint kb   = neighborZ[KQK];
      uint ksw  = neighborY[kw];
      uint kne  = KQK;
      uint kse  = ks;
      uint knw  = kw;
      uint kbw  = neighborZ[kw];
      uint kte  = KQK;
      uint kbe  = kb;
      uint ktw  = kw;
      uint kbs  = neighborZ[ks];
      uint ktn  = KQK;
      uint kbn  = kb;
      uint kts  = ks;
      uint ktse = ks;
      uint kbnw = kbw;
      uint ktnw = kw;
      uint kbse = kbs;
      uint ktsw = ksw;
      uint kbne = kb;
      uint ktne = KQK;
      uint kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[dP00])[ke   ];
      f_E    = (D.f[dM00])[kw   ];
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
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (c1o1 + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (c1o1 + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (c1o1 + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

      //////////////////////////////////////////////////////////////////////////

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho + c9o2 * ( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq);
		 (D.f[dM00])[kw] = (c1o1 - q) / (c1o1 + q)*(f_E - f_W + (f_E + f_W - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_E + f_W) - c6o1*c2o27*(VeloX)) / (c1o1 + q);
	  }

	  q = q_dirW[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (-vx1)*(-vx1) * (c1o1 + drho) - cu_sq);
		  (D.f[dP00])[ke] = (c1o1 - q) / (c1o1 + q)*(f_W - f_E + (f_W + f_E - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_W + f_E) - c6o1*c2o27*(-VeloX)) / (c1o1 + q);
	  }

	  q = q_dirN[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (vx2)*(vx2) * (c1o1 + drho) - cu_sq);
		  (D.f[DIR_0M0])[ks] = (c1o1 - q) / (c1o1 + q)*(f_N - f_S + (f_N + f_S - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_N + f_S) - c6o1*c2o27*(VeloY)) / (c1o1 + q);
	  }

	  q = q_dirS[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (-vx2)*(-vx2) * (c1o1 + drho) - cu_sq);
		  (D.f[DIR_0P0])[kn] = (c1o1 - q) / (c1o1 + q)*(f_S - f_N + (f_S + f_N - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_S + f_N) - c6o1*c2o27*(-VeloY)) / (c1o1 + q);
	  }

	  q = q_dirT[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (vx3)*(vx3) * (c1o1 + drho) - cu_sq);
		  (D.f[DIR_00M])[kb] = (c1o1 - q) / (c1o1 + q)*(f_T - f_B + (f_T + f_B - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_T + f_B) - c6o1*c2o27*(VeloZ)) / (c1o1 + q);
	  }

	  q = q_dirB[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (-vx3)*(-vx3) * (c1o1 + drho) - cu_sq);
		  (D.f[DIR_00P])[kt] = (c1o1 - q) / (c1o1 + q)*(f_B - f_T + (f_B + f_T - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_B + f_T) - c6o1*c2o27*(-VeloZ)) / (c1o1 + q);
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * ( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq);
         (D.f[DIR_MM0])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*( VeloX+VeloY))/(c1o1+q);
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq);
         (D.f[DIR_PP0])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * ( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq);
         (D.f[DIR_MP0])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq);
         (D.f[DIR_PM0])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * ( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq);
         (D.f[DIR_M0M])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq);
         (D.f[DIR_P0P])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * ( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq);
         (D.f[DIR_M0P])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq);
         (D.f[DIR_P0M])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq);
         (D.f[DIR_0MM])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq);
         (D.f[DIR_0PP])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*(-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq);
         (D.f[DIR_0MP])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq);
         (D.f[DIR_0PM])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*(-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * ( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MMM])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * (-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PPP])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * ( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MMP])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * (-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PPM])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * ( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MPM])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * (-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PMP])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * ( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MPP])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * (-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PMM])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QDeviceCompThinWallsPartOne27(
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
	if (isEvenTimestep == true)
	{
		D.f[dP00] = &DD[dP00 * numberOfLBnodes];
		D.f[dM00] = &DD[dM00 * numberOfLBnodes];
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
		D.f[d000] = &DD[d000 * numberOfLBnodes];
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
		D.f[dM00] = &DD[dP00 * numberOfLBnodes];
		D.f[dP00] = &DD[dM00 * numberOfLBnodes];
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
		D.f[d000] = &DD[d000 * numberOfLBnodes];
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

	if (k < numberOfBCnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		real *q_dirE, *q_dirW, *q_dirN, *q_dirS, *q_dirT, *q_dirB,
			*q_dirNE, *q_dirSW, *q_dirSE, *q_dirNW, *q_dirTE, *q_dirBW,
			*q_dirBE, *q_dirTW, *q_dirTN, *q_dirBS, *q_dirBN, *q_dirTS,
			*q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			*q_dirBSE, *q_dirBNW;
		q_dirE = &QQ[dP00 * numberOfBCnodes];
		q_dirW = &QQ[dM00 * numberOfBCnodes];
		q_dirN = &QQ[DIR_0P0 * numberOfBCnodes];
		q_dirS = &QQ[DIR_0M0 * numberOfBCnodes];
		q_dirT = &QQ[DIR_00P * numberOfBCnodes];
		q_dirB = &QQ[DIR_00M * numberOfBCnodes];
		q_dirNE = &QQ[DIR_PP0 * numberOfBCnodes];
		q_dirSW = &QQ[DIR_MM0 * numberOfBCnodes];
		q_dirSE = &QQ[DIR_PM0 * numberOfBCnodes];
		q_dirNW = &QQ[DIR_MP0 * numberOfBCnodes];
		q_dirTE = &QQ[DIR_P0P * numberOfBCnodes];
		q_dirBW = &QQ[DIR_M0M * numberOfBCnodes];
		q_dirBE = &QQ[DIR_P0M * numberOfBCnodes];
		q_dirTW = &QQ[DIR_M0P * numberOfBCnodes];
		q_dirTN = &QQ[DIR_0PP * numberOfBCnodes];
		q_dirBS = &QQ[DIR_0MM * numberOfBCnodes];
		q_dirBN = &QQ[DIR_0PM * numberOfBCnodes];
		q_dirTS = &QQ[DIR_0MP * numberOfBCnodes];
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
		real f_E, f_W, f_N, f_S, f_T, f_B, f_NE, f_SW, f_SE, f_NW, f_TE, f_BW, f_BE,
			f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

		f_W = (D.f[dP00])[ke];
		f_E = (D.f[dM00])[kw];
		f_S = (D.f[DIR_0P0])[kn];
		f_N = (D.f[DIR_0M0])[ks];
		f_B = (D.f[DIR_00P])[kt];
		f_T = (D.f[DIR_00M])[kb];
		f_SW = (D.f[DIR_PP0])[kne];
		f_NE = (D.f[DIR_MM0])[ksw];
		f_NW = (D.f[DIR_PM0])[kse];
		f_SE = (D.f[DIR_MP0])[knw];
		f_BW = (D.f[DIR_P0P])[kte];
		f_TE = (D.f[DIR_M0M])[kbw];
		f_TW = (D.f[DIR_P0M])[kbe];
		f_BE = (D.f[DIR_M0P])[ktw];
		f_BS = (D.f[DIR_0PP])[ktn];
		f_TN = (D.f[DIR_0MM])[kbs];
		f_TS = (D.f[DIR_0PM])[kbn];
		f_BN = (D.f[DIR_0MP])[kts];
		f_BSW = (D.f[DIR_PPP])[ktne];
		f_BNE = (D.f[DIR_MMP])[ktsw];
		f_BNW = (D.f[DIR_PMP])[ktse];
		f_BSE = (D.f[DIR_MPP])[ktnw];
		f_TSW = (D.f[DIR_PPM])[kbne];
		f_TNE = (D.f[DIR_MMM])[kbsw];
		f_TNW = (D.f[DIR_PMM])[kbse];
		f_TSE = (D.f[DIR_MPM])[kbnw];
		////////////////////////////////////////////////////////////////////////////////
		real vx1, vx2, vx3, drho, feq, q;
		drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
			f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
			f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]);

		vx1 = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
			((f_BE - f_TW) + (f_TE - f_BW)) + ((f_SE - f_NW) + (f_NE - f_SW)) +
			(f_E - f_W)) / (c1o1 + drho);


		vx2 = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
			((f_BN - f_TS) + (f_TN - f_BS)) + (-(f_SE - f_NW) + (f_NE - f_SW)) +
			(f_N - f_S)) / (c1o1 + drho);

		vx3 = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
			(-(f_BN - f_TS) + (f_TN - f_BS)) + ((f_TE - f_BW) - (f_BE - f_TW)) +
			(f_T - f_B)) / (c1o1 + drho);

		////////////////////////////////////////////////////////////////////////////////
		real cu_sq = c3o2*(vx1*vx1 + vx2*vx2 + vx3*vx3) * (c1o1 + drho);
		////////////////////////////////////////////////////////////////////////////////

		q = q_dirE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(vx1)*(vx1) * (c1o1 + drho) - cu_sq);
			(D.f[dM00])[kw] = (c1o1 - q) / (c1o1 + q)*(f_E - f_W + (f_E + f_W - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_E + f_W)) / (c1o1 + q);
		}

		q = q_dirW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(-vx1)*(-vx1) * (c1o1 + drho) - cu_sq);
			(D.f[dP00])[ke] = (c1o1 - q) / (c1o1 + q)*(f_W - f_E + (f_W + f_E - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_W + f_E)) / (c1o1 + q);
		}

		q = q_dirN[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(vx2)*(vx2) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_0M0])[ks] = (c1o1 - q) / (c1o1 + q)*(f_N - f_S + (f_N + f_S - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_N + f_S)) / (c1o1 + q);
		}

		q = q_dirS[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(-vx2)*(-vx2) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_0P0])[kn] = (c1o1 - q) / (c1o1 + q)*(f_S - f_N + (f_S + f_N - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_S + f_N)) / (c1o1 + q);
		}

		q = q_dirT[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(vx3)*(vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_00M])[kb] = (c1o1 - q) / (c1o1 + q)*(f_T - f_B + (f_T + f_B - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_T + f_B)) / (c1o1 + q);
		}

		q = q_dirB[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(-vx3)*(-vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_00P])[kt] = (c1o1 - q) / (c1o1 + q)*(f_B - f_T + (f_B + f_T - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_B + f_T)) / (c1o1 + q);
		}

		q = q_dirNE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx1 + vx2)*(vx1 + vx2) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_MM0])[ksw] = (c1o1 - q) / (c1o1 + q)*(f_NE - f_SW + (f_NE + f_SW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_NE + f_SW)) / (c1o1 + q);
		}

		q = q_dirSW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx1 - vx2)*(-vx1 - vx2) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_PP0])[kne] = (c1o1 - q) / (c1o1 + q)*(f_SW - f_NE + (f_SW + f_NE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_SW + f_NE)) / (c1o1 + q);
		}

		q = q_dirSE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx1 - vx2)*(vx1 - vx2) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_MP0])[knw] = (c1o1 - q) / (c1o1 + q)*(f_SE - f_NW + (f_SE + f_NW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_SE + f_NW)) / (c1o1 + q);
		}

		q = q_dirNW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx1 + vx2)*(-vx1 + vx2) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_PM0])[kse] = (c1o1 - q) / (c1o1 + q)*(f_NW - f_SE + (f_NW + f_SE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_NW + f_SE)) / (c1o1 + q);
		}

		q = q_dirTE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx1 + vx3)*(vx1 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_M0M])[kbw] = (c1o1 - q) / (c1o1 + q)*(f_TE - f_BW + (f_TE + f_BW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TE + f_BW)) / (c1o1 + q);
		}

		q = q_dirBW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx1 - vx3)*(-vx1 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_P0P])[kte] = (c1o1 - q) / (c1o1 + q)*(f_BW - f_TE + (f_BW + f_TE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BW + f_TE)) / (c1o1 + q);
		}

		q = q_dirBE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx1 - vx3)*(vx1 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_M0P])[ktw] = (c1o1 - q) / (c1o1 + q)*(f_BE - f_TW + (f_BE + f_TW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BE + f_TW)) / (c1o1 + q);
		}

		q = q_dirTW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx1 + vx3)*(-vx1 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_P0M])[kbe] = (c1o1 - q) / (c1o1 + q)*(f_TW - f_BE + (f_TW + f_BE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TW + f_BE)) / (c1o1 + q);
		}

		q = q_dirTN[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx2 + vx3)*(vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_0MM])[kbs] = (c1o1 - q) / (c1o1 + q)*(f_TN - f_BS + (f_TN + f_BS - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TN + f_BS)) / (c1o1 + q);
		}

		q = q_dirBS[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx2 - vx3)*(-vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_0PP])[ktn] = (c1o1 - q) / (c1o1 + q)*(f_BS - f_TN + (f_BS + f_TN - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BS + f_TN)) / (c1o1 + q);
		}

		q = q_dirBN[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx2 - vx3)*(vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_0MP])[kts] = (c1o1 - q) / (c1o1 + q)*(f_BN - f_TS + (f_BN + f_TS - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BN + f_TS)) / (c1o1 + q);
		}

		q = q_dirTS[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx2 + vx3)*(-vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_0PM])[kbn] = (c1o1 - q) / (c1o1 + q)*(f_TS - f_BN + (f_TS + f_BN - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TS + f_BN)) / (c1o1 + q);
		}

		q = q_dirTNE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(vx1 + vx2 + vx3)*(vx1 + vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_MMM])[kbsw] = (c1o1 - q) / (c1o1 + q)*(f_TNE - f_BSW + (f_TNE + f_BSW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TNE + f_BSW)) / (c1o1 + q);
		}

		q = q_dirBSW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(-vx1 - vx2 - vx3)*(-vx1 - vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_PPP])[ktne] = (c1o1 - q) / (c1o1 + q)*(f_BSW - f_TNE + (f_BSW + f_TNE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BSW + f_TNE)) / (c1o1 + q);
		}

		q = q_dirBNE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(vx1 + vx2 - vx3)*(vx1 + vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_MMP])[ktsw] = (c1o1 - q) / (c1o1 + q)*(f_BNE - f_TSW + (f_BNE + f_TSW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BNE + f_TSW)) / (c1o1 + q);
		}

		q = q_dirTSW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(-vx1 - vx2 + vx3)*(-vx1 - vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_PPM])[kbne] = (c1o1 - q) / (c1o1 + q)*(f_TSW - f_BNE + (f_TSW + f_BNE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TSW + f_BNE)) / (c1o1 + q);
		}

		q = q_dirTSE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(vx1 - vx2 + vx3)*(vx1 - vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_MPM])[kbnw] = (c1o1 - q) / (c1o1 + q)*(f_TSE - f_BNW + (f_TSE + f_BNW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TSE + f_BNW)) / (c1o1 + q);
		}

		q = q_dirBNW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(-vx1 + vx2 - vx3)*(-vx1 + vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_PMP])[ktse] = (c1o1 - q) / (c1o1 + q)*(f_BNW - f_TSE + (f_BNW + f_TSE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BNW + f_TSE)) / (c1o1 + q);
		}

		q = q_dirBSE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(vx1 - vx2 - vx3)*(vx1 - vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_MPP])[ktnw] = (c1o1 - q) / (c1o1 + q)*(f_BSE - f_TNW + (f_BSE + f_TNW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BSE + f_TNW)) / (c1o1 + q);
		}

		q = q_dirTNW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(-vx1 + vx2 + vx3)*(-vx1 + vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[DIR_PMM])[kbse] = (c1o1 - q) / (c1o1 + q)*(f_TNW - f_BSE + (f_TNW + f_BSE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TNW + f_BSE)) / (c1o1 + q);
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QThinWallsPartTwo27(
	real* DD, 
	int* k_Q, 
	real* QQ,
	uint numberOfBCnodes, 
	uint* geom,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint* neighborWSB,
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
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dP00 * numberOfBCnodes];
      q_dirW   = &QQ[dM00 * numberOfBCnodes];
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
      uint KQK  = k_Q[k];
      //uint kzero= KQK;
      uint ke   = KQK;
      uint kw   = neighborX[KQK];
      uint kn   = KQK;
      uint ks   = neighborY[KQK];
      uint kt   = KQK;
      uint kb   = neighborZ[KQK];
      uint ksw  = neighborY[kw];
      uint kne  = KQK;
      uint kse  = ks;
      uint knw  = kw;
      uint kbw  = neighborZ[kw];
      uint kte  = KQK;
      uint kbe  = kb;
      uint ktw  = kw;
      uint kbs  = neighborZ[ks];
      uint ktn  = KQK;
      uint kbn  = kb;
      uint kts  = ks;
      uint ktse = ks;
      uint kbnw = kbw;
      uint ktnw = kw;
      uint kbse = kbs;
      uint ktsw = ksw;
      uint kbne = kb;
      uint ktne = KQK;
      uint kbsw = neighborZ[ksw];
	  ////////////////////////////////////////////////////////////////////////////////
	  //anti ET intermediate steps
	  uint kmmm = neighborWSB[KQK]; // -1 -1 -1
	  uint k0mm = neighborX[kmmm];  //  0 -1 -1
	  uint km0m = neighborY[kmmm];  // -1  0 -1
	  uint kmm0 = neighborZ[kmmm];  // -1 -1  0
	  uint k0m0 = neighborX[kmm0];  //  0 -1  0
	  uint km00 = neighborY[kmm0];  // -1  0  0
	  /////////////////////////////////////////////////
	  //final indices for anti ET
	  uint kpmm = neighborX[k0mm];  //  1 -1 -1
	  uint kmpm = neighborY[km0m];  // -1  1 -1
	  uint kmmp = neighborZ[kmm0];  // -1 -1  1
	  uint kmp0 = neighborY[km00];  // -1  1  0
	  uint km0p = neighborZ[km00];  // -1  0  1
	  uint k0mp = neighborZ[k0m0];  //  0 -1  1
	  ////////////////////////////////////////////////////////////////////////////////
	  Distributions27 D, DN;
	  if (isEvenTimestep == true)
	  {
		  D.f[dP00] = &DD[dP00 * numberOfLBnodes];
		  D.f[dM00] = &DD[dM00 * numberOfLBnodes];
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
		  D.f[d000] = &DD[d000 * numberOfLBnodes];
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
		  D.f[dM00] = &DD[dP00 * numberOfLBnodes];
		  D.f[dP00] = &DD[dM00 * numberOfLBnodes];
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
		  D.f[d000] = &DD[d000 * numberOfLBnodes];
		  D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
		  D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
		  D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
		  D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
		  D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
		  D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
		  D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
		  D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
	  }
	  if (isEvenTimestep==false)
      {
         DN.f[dP00] = &DD[dP00 * numberOfLBnodes];
         DN.f[dM00] = &DD[dM00 * numberOfLBnodes];
         DN.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         DN.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         DN.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         DN.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         DN.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         DN.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         DN.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         DN.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         DN.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         DN.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         DN.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         DN.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         DN.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         DN.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         DN.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         DN.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         DN.f[d000] = &DD[d000 * numberOfLBnodes];
         DN.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         DN.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         DN.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         DN.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         DN.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         DN.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         DN.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         DN.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         DN.f[dM00] = &DD[dP00 * numberOfLBnodes];
         DN.f[dP00] = &DD[dM00 * numberOfLBnodes];
         DN.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         DN.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         DN.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         DN.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         DN.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         DN.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         DN.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         DN.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         DN.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         DN.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         DN.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         DN.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         DN.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         DN.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         DN.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         DN.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         DN.f[d000] = &DD[d000 * numberOfLBnodes];
         DN.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         DN.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         DN.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         DN.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         DN.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         DN.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         DN.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         DN.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //directions allways exchange
	  //(-1 -1 -1) (-1  0  0) ( 0 -1  0) ( 0  0 -1) (-1 -1  0) (-1  0 -1) ( 0 -1 -1) ( 1  1 -1) ( 1 -1  1) (-1  1  1) ( 1 -1  0) ( 1  0 -1) ( 0  1 -1)
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //directions exchange if solid neighbor
	  //( 1  1  1) ( 1  0  0) ( 0  1  0) ( 0  0  1) ( 1  1  0) ( 1  0  1) ( 0  1  1) (-1 -1  1) (-1  1 -1) ( 1 -1 -1) (-1  1  0) (-1  0  1) ( 0 -1  1)
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real q, tmp;
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1){ if (geom[kw  ] < GEO_FLUID){tmp = (DN.f[dM00])[kw  ]; (DN.f[dM00])[kw  ]=(D.f[dM00])[kw  ]; (D.f[dM00])[kw  ]=tmp;}}
	  q = q_dirW[k];   if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dP00])[ke  ]; (DN.f[dP00])[ke  ]=(D.f[dP00])[ke  ]; (D.f[dP00])[ke  ]=tmp;}}
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1){ if (geom[ks  ] < GEO_FLUID){tmp = (DN.f[DIR_0M0])[ks  ]; (DN.f[DIR_0M0])[ks  ]=(D.f[DIR_0M0])[ks  ]; (D.f[DIR_0M0])[ks  ]=tmp;}}
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[DIR_0P0])[kn  ]; (DN.f[DIR_0P0])[kn  ]=(D.f[DIR_0P0])[kn  ]; (D.f[DIR_0P0])[kn  ]=tmp;}}
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1){ if (geom[kb  ] < GEO_FLUID){tmp = (DN.f[DIR_00M])[kb  ]; (DN.f[DIR_00M])[kb  ]=(D.f[DIR_00M])[kb  ]; (D.f[DIR_00M])[kb  ]=tmp;}}
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[DIR_00P])[kt  ]; (DN.f[DIR_00P])[kt  ]=(D.f[DIR_00P])[kt  ]; (D.f[DIR_00P])[kt  ]=tmp;}}
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1){ if (geom[ksw ] < GEO_FLUID){tmp = (DN.f[DIR_MM0])[ksw ]; (DN.f[DIR_MM0])[ksw ]=(D.f[DIR_MM0])[ksw ]; (D.f[DIR_MM0])[ksw ]=tmp;}}
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[DIR_PP0])[kne ]; (DN.f[DIR_PP0])[kne ]=(D.f[DIR_PP0])[kne ]; (D.f[DIR_PP0])[kne ]=tmp;}}
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[DIR_MP0])[knw ]; (DN.f[DIR_MP0])[knw ]=(D.f[DIR_MP0])[knw ]; (D.f[DIR_MP0])[knw ]=tmp;}}
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1){ if (geom[kmp0] < GEO_FLUID){tmp = (DN.f[DIR_PM0])[kse ]; (DN.f[DIR_PM0])[kse ]=(D.f[DIR_PM0])[kse ]; (D.f[DIR_PM0])[kse ]=tmp;}}
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1){ if (geom[kbw ] < GEO_FLUID){tmp = (DN.f[DIR_M0M])[kbw ]; (DN.f[DIR_M0M])[kbw ]=(D.f[DIR_M0M])[kbw ]; (D.f[DIR_M0M])[kbw ]=tmp;}}
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[DIR_P0P])[kte ]; (DN.f[DIR_P0P])[kte ]=(D.f[DIR_P0P])[kte ]; (D.f[DIR_P0P])[kte ]=tmp;}}
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[DIR_M0P])[ktw ]; (DN.f[DIR_M0P])[ktw ]=(D.f[DIR_M0P])[ktw ]; (D.f[DIR_M0P])[ktw ]=tmp;}}
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1){ if (geom[km0p] < GEO_FLUID){tmp = (DN.f[DIR_P0M])[kbe ]; (DN.f[DIR_P0M])[kbe ]=(D.f[DIR_P0M])[kbe ]; (D.f[DIR_P0M])[kbe ]=tmp;}}
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1){ if (geom[kbs ] < GEO_FLUID){tmp = (DN.f[DIR_0MM])[kbs ]; (DN.f[DIR_0MM])[kbs ]=(D.f[DIR_0MM])[kbs ]; (D.f[DIR_0MM])[kbs ]=tmp;}}
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[DIR_0PP])[ktn ]; (DN.f[DIR_0PP])[ktn ]=(D.f[DIR_0PP])[ktn ]; (D.f[DIR_0PP])[ktn ]=tmp;}}
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[DIR_0MP])[kts ]; (DN.f[DIR_0MP])[kts ]=(D.f[DIR_0MP])[kts ]; (D.f[DIR_0MP])[kts ]=tmp;}}
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1){ if (geom[k0mp] < GEO_FLUID){tmp = (DN.f[DIR_0PM])[kbn ]; (DN.f[DIR_0PM])[kbn ]=(D.f[DIR_0PM])[kbn ]; (D.f[DIR_0PM])[kbn ]=tmp;}}
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1){ if (geom[kbsw] < GEO_FLUID){tmp = (DN.f[DIR_MMM])[kbsw]; (DN.f[DIR_MMM])[kbsw]=(D.f[DIR_MMM])[kbsw]; (D.f[DIR_MMM])[kbsw]=tmp;}}
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[DIR_PPP])[ktne]; (DN.f[DIR_PPP])[ktne]=(D.f[DIR_PPP])[ktne]; (D.f[DIR_PPP])[ktne]=tmp;}}
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[DIR_MMP])[ktsw]; (DN.f[DIR_MMP])[ktsw]=(D.f[DIR_MMP])[ktsw]; (D.f[DIR_MMP])[ktsw]=tmp;}}
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1){ if (geom[kmmp] < GEO_FLUID){tmp = (DN.f[DIR_PPM])[kbne]; (DN.f[DIR_PPM])[kbne]=(D.f[DIR_PPM])[kbne]; (D.f[DIR_PPM])[kbne]=tmp;}}
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[DIR_MPM])[kbnw]; (DN.f[DIR_MPM])[kbnw]=(D.f[DIR_MPM])[kbnw]; (D.f[DIR_MPM])[kbnw]=tmp;}}
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1){ if (geom[kmpm] < GEO_FLUID){tmp = (DN.f[DIR_PMP])[ktse]; (DN.f[DIR_PMP])[ktse]=(D.f[DIR_PMP])[ktse]; (D.f[DIR_PMP])[ktse]=tmp;}}
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1){ if (geom[kpmm] < GEO_FLUID){tmp = (DN.f[DIR_MPP])[ktnw]; (DN.f[DIR_MPP])[ktnw]=(D.f[DIR_MPP])[ktnw]; (D.f[DIR_MPP])[ktnw]=tmp;}}
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[DIR_PMM])[kbse]; (DN.f[DIR_PMM])[kbse]=(D.f[DIR_PMM])[kbse]; (D.f[DIR_PMM])[kbse]=tmp;}}
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
