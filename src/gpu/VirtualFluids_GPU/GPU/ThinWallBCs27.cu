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
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

/////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceCompThinWallsPartOne27(
	real* vx,
	real* vy,
	real* vz,
	real* DD, 
	int* k_Q, 
	real* QQ,
	int numberOfBCnodes, 
	real om1, 
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint size_Mat, 
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
      q_dirE   = &QQ[E   * numberOfBCnodes];
      q_dirW   = &QQ[W   * numberOfBCnodes];
      q_dirN   = &QQ[N   * numberOfBCnodes];
      q_dirS   = &QQ[S   * numberOfBCnodes];
      q_dirT   = &QQ[T   * numberOfBCnodes];
      q_dirB   = &QQ[B   * numberOfBCnodes];
      q_dirNE  = &QQ[NE  * numberOfBCnodes];
      q_dirSW  = &QQ[SW  * numberOfBCnodes];
      q_dirSE  = &QQ[SE  * numberOfBCnodes];
      q_dirNW  = &QQ[NW  * numberOfBCnodes];
      q_dirTE  = &QQ[TE  * numberOfBCnodes];
      q_dirBW  = &QQ[BW  * numberOfBCnodes];
      q_dirBE  = &QQ[BE  * numberOfBCnodes];
      q_dirTW  = &QQ[TW  * numberOfBCnodes];
      q_dirTN  = &QQ[TN  * numberOfBCnodes];
      q_dirBS  = &QQ[BS  * numberOfBCnodes];
      q_dirBN  = &QQ[BN  * numberOfBCnodes];
      q_dirTS  = &QQ[TS  * numberOfBCnodes];
      q_dirTNE = &QQ[TNE * numberOfBCnodes];
      q_dirTSW = &QQ[TSW * numberOfBCnodes];
      q_dirTSE = &QQ[TSE * numberOfBCnodes];
      q_dirTNW = &QQ[TNW * numberOfBCnodes];
      q_dirBNE = &QQ[BNE * numberOfBCnodes];
      q_dirBSW = &QQ[BSW * numberOfBCnodes];
      q_dirBSE = &QQ[BSE * numberOfBCnodes];
      q_dirBNW = &QQ[BNW * numberOfBCnodes];
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
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[REST])[kzero]); 

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
		 (D.f[W])[kw] = (c1o1 - q) / (c1o1 + q)*(f_E - f_W + (f_E + f_W - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_E + f_W) - c6o1*c2o27*(VeloX)) / (c1o1 + q);
	  }

	  q = q_dirW[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (-vx1)*(-vx1) * (c1o1 + drho) - cu_sq);
		  (D.f[E])[ke] = (c1o1 - q) / (c1o1 + q)*(f_W - f_E + (f_W + f_E - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_W + f_E) - c6o1*c2o27*(-VeloX)) / (c1o1 + q);
	  }

	  q = q_dirN[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (vx2)*(vx2) * (c1o1 + drho) - cu_sq);
		  (D.f[S])[ks] = (c1o1 - q) / (c1o1 + q)*(f_N - f_S + (f_N + f_S - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_N + f_S) - c6o1*c2o27*(VeloY)) / (c1o1 + q);
	  }

	  q = q_dirS[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (-vx2)*(-vx2) * (c1o1 + drho) - cu_sq);
		  (D.f[N])[kn] = (c1o1 - q) / (c1o1 + q)*(f_S - f_N + (f_S + f_N - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_S + f_N) - c6o1*c2o27*(-VeloY)) / (c1o1 + q);
	  }

	  q = q_dirT[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (vx3)*(vx3) * (c1o1 + drho) - cu_sq);
		  (D.f[B])[kb] = (c1o1 - q) / (c1o1 + q)*(f_T - f_B + (f_T + f_B - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_T + f_B) - c6o1*c2o27*(VeloZ)) / (c1o1 + q);
	  }

	  q = q_dirB[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (-vx3)*(-vx3) * (c1o1 + drho) - cu_sq);
		  (D.f[T])[kt] = (c1o1 - q) / (c1o1 + q)*(f_B - f_T + (f_B + f_T - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_B + f_T) - c6o1*c2o27*(-VeloZ)) / (c1o1 + q);
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * ( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq);
         (D.f[SW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*( VeloX+VeloY))/(c1o1+q);
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq);
         (D.f[NE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * ( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq);
         (D.f[NW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq);
         (D.f[SE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * ( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq);
         (D.f[BW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq);
         (D.f[TE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * ( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq);
         (D.f[TW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq);
         (D.f[BE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq);
         (D.f[BS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq);
         (D.f[TN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*(-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq);
         (D.f[TS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq);
         (D.f[BN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*(-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * ( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * (-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * ( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * (-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * ( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * (-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * ( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * (-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDeviceCompThinWallsPartOne27(
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
	if (isEvenTimestep == true)
	{
		D.f[E] = &DD[E   *size_Mat];
		D.f[W] = &DD[W   *size_Mat];
		D.f[N] = &DD[N   *size_Mat];
		D.f[S] = &DD[S   *size_Mat];
		D.f[T] = &DD[T   *size_Mat];
		D.f[B] = &DD[B   *size_Mat];
		D.f[NE] = &DD[NE  *size_Mat];
		D.f[SW] = &DD[SW  *size_Mat];
		D.f[SE] = &DD[SE  *size_Mat];
		D.f[NW] = &DD[NW  *size_Mat];
		D.f[TE] = &DD[TE  *size_Mat];
		D.f[BW] = &DD[BW  *size_Mat];
		D.f[BE] = &DD[BE  *size_Mat];
		D.f[TW] = &DD[TW  *size_Mat];
		D.f[TN] = &DD[TN  *size_Mat];
		D.f[BS] = &DD[BS  *size_Mat];
		D.f[BN] = &DD[BN  *size_Mat];
		D.f[TS] = &DD[TS  *size_Mat];
		D.f[REST] = &DD[REST*size_Mat];
		D.f[TNE] = &DD[TNE *size_Mat];
		D.f[TSW] = &DD[TSW *size_Mat];
		D.f[TSE] = &DD[TSE *size_Mat];
		D.f[TNW] = &DD[TNW *size_Mat];
		D.f[BNE] = &DD[BNE *size_Mat];
		D.f[BSW] = &DD[BSW *size_Mat];
		D.f[BSE] = &DD[BSE *size_Mat];
		D.f[BNW] = &DD[BNW *size_Mat];
	}
	else
	{
		D.f[W] = &DD[E   *size_Mat];
		D.f[E] = &DD[W   *size_Mat];
		D.f[S] = &DD[N   *size_Mat];
		D.f[N] = &DD[S   *size_Mat];
		D.f[B] = &DD[T   *size_Mat];
		D.f[T] = &DD[B   *size_Mat];
		D.f[SW] = &DD[NE  *size_Mat];
		D.f[NE] = &DD[SW  *size_Mat];
		D.f[NW] = &DD[SE  *size_Mat];
		D.f[SE] = &DD[NW  *size_Mat];
		D.f[BW] = &DD[TE  *size_Mat];
		D.f[TE] = &DD[BW  *size_Mat];
		D.f[TW] = &DD[BE  *size_Mat];
		D.f[BE] = &DD[TW  *size_Mat];
		D.f[BS] = &DD[TN  *size_Mat];
		D.f[TN] = &DD[BS  *size_Mat];
		D.f[TS] = &DD[BN  *size_Mat];
		D.f[BN] = &DD[TS  *size_Mat];
		D.f[REST] = &DD[REST*size_Mat];
		D.f[TNE] = &DD[BSW *size_Mat];
		D.f[TSW] = &DD[BNE *size_Mat];
		D.f[TSE] = &DD[BNW *size_Mat];
		D.f[TNW] = &DD[BSE *size_Mat];
		D.f[BNE] = &DD[TSW *size_Mat];
		D.f[BSW] = &DD[TNE *size_Mat];
		D.f[BSE] = &DD[TNW *size_Mat];
		D.f[BNW] = &DD[TSE *size_Mat];
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
		q_dirE = &QQ[E   * numberOfBCnodes];
		q_dirW = &QQ[W   * numberOfBCnodes];
		q_dirN = &QQ[N   * numberOfBCnodes];
		q_dirS = &QQ[S   * numberOfBCnodes];
		q_dirT = &QQ[T   * numberOfBCnodes];
		q_dirB = &QQ[B   * numberOfBCnodes];
		q_dirNE = &QQ[NE  * numberOfBCnodes];
		q_dirSW = &QQ[SW  * numberOfBCnodes];
		q_dirSE = &QQ[SE  * numberOfBCnodes];
		q_dirNW = &QQ[NW  * numberOfBCnodes];
		q_dirTE = &QQ[TE  * numberOfBCnodes];
		q_dirBW = &QQ[BW  * numberOfBCnodes];
		q_dirBE = &QQ[BE  * numberOfBCnodes];
		q_dirTW = &QQ[TW  * numberOfBCnodes];
		q_dirTN = &QQ[TN  * numberOfBCnodes];
		q_dirBS = &QQ[BS  * numberOfBCnodes];
		q_dirBN = &QQ[BN  * numberOfBCnodes];
		q_dirTS = &QQ[TS  * numberOfBCnodes];
		q_dirTNE = &QQ[TNE * numberOfBCnodes];
		q_dirTSW = &QQ[TSW * numberOfBCnodes];
		q_dirTSE = &QQ[TSE * numberOfBCnodes];
		q_dirTNW = &QQ[TNW * numberOfBCnodes];
		q_dirBNE = &QQ[BNE * numberOfBCnodes];
		q_dirBSW = &QQ[BSW * numberOfBCnodes];
		q_dirBSE = &QQ[BSE * numberOfBCnodes];
		q_dirBNW = &QQ[BNW * numberOfBCnodes];
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

		f_W = (D.f[E])[ke];
		f_E = (D.f[W])[kw];
		f_S = (D.f[N])[kn];
		f_N = (D.f[S])[ks];
		f_B = (D.f[T])[kt];
		f_T = (D.f[B])[kb];
		f_SW = (D.f[NE])[kne];
		f_NE = (D.f[SW])[ksw];
		f_NW = (D.f[SE])[kse];
		f_SE = (D.f[NW])[knw];
		f_BW = (D.f[TE])[kte];
		f_TE = (D.f[BW])[kbw];
		f_TW = (D.f[BE])[kbe];
		f_BE = (D.f[TW])[ktw];
		f_BS = (D.f[TN])[ktn];
		f_TN = (D.f[BS])[kbs];
		f_TS = (D.f[BN])[kbn];
		f_BN = (D.f[TS])[kts];
		f_BSW = (D.f[TNE])[ktne];
		f_BNE = (D.f[TSW])[ktsw];
		f_BNW = (D.f[TSE])[ktse];
		f_BSE = (D.f[TNW])[ktnw];
		f_TSW = (D.f[BNE])[kbne];
		f_TNE = (D.f[BSW])[kbsw];
		f_TNW = (D.f[BSE])[kbse];
		f_TSE = (D.f[BNW])[kbnw];
		////////////////////////////////////////////////////////////////////////////////
		real vx1, vx2, vx3, drho, feq, q;
		drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
			f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
			f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[REST])[kzero]);

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
			(D.f[W])[kw] = (c1o1 - q) / (c1o1 + q)*(f_E - f_W + (f_E + f_W - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_E + f_W)) / (c1o1 + q);
		}

		q = q_dirW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(-vx1)*(-vx1) * (c1o1 + drho) - cu_sq);
			(D.f[E])[ke] = (c1o1 - q) / (c1o1 + q)*(f_W - f_E + (f_W + f_E - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_W + f_E)) / (c1o1 + q);
		}

		q = q_dirN[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(vx2)*(vx2) * (c1o1 + drho) - cu_sq);
			(D.f[S])[ks] = (c1o1 - q) / (c1o1 + q)*(f_N - f_S + (f_N + f_S - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_N + f_S)) / (c1o1 + q);
		}

		q = q_dirS[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(-vx2)*(-vx2) * (c1o1 + drho) - cu_sq);
			(D.f[N])[kn] = (c1o1 - q) / (c1o1 + q)*(f_S - f_N + (f_S + f_N - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_S + f_N)) / (c1o1 + q);
		}

		q = q_dirT[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(vx3)*(vx3) * (c1o1 + drho) - cu_sq);
			(D.f[B])[kb] = (c1o1 - q) / (c1o1 + q)*(f_T - f_B + (f_T + f_B - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_T + f_B)) / (c1o1 + q);
		}

		q = q_dirB[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(-vx3)*(-vx3) * (c1o1 + drho) - cu_sq);
			(D.f[T])[kt] = (c1o1 - q) / (c1o1 + q)*(f_B - f_T + (f_B + f_T - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_B + f_T)) / (c1o1 + q);
		}

		q = q_dirNE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx1 + vx2)*(vx1 + vx2) * (c1o1 + drho) - cu_sq);
			(D.f[SW])[ksw] = (c1o1 - q) / (c1o1 + q)*(f_NE - f_SW + (f_NE + f_SW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_NE + f_SW)) / (c1o1 + q);
		}

		q = q_dirSW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx1 - vx2)*(-vx1 - vx2) * (c1o1 + drho) - cu_sq);
			(D.f[NE])[kne] = (c1o1 - q) / (c1o1 + q)*(f_SW - f_NE + (f_SW + f_NE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_SW + f_NE)) / (c1o1 + q);
		}

		q = q_dirSE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx1 - vx2)*(vx1 - vx2) * (c1o1 + drho) - cu_sq);
			(D.f[NW])[knw] = (c1o1 - q) / (c1o1 + q)*(f_SE - f_NW + (f_SE + f_NW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_SE + f_NW)) / (c1o1 + q);
		}

		q = q_dirNW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx1 + vx2)*(-vx1 + vx2) * (c1o1 + drho) - cu_sq);
			(D.f[SE])[kse] = (c1o1 - q) / (c1o1 + q)*(f_NW - f_SE + (f_NW + f_SE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_NW + f_SE)) / (c1o1 + q);
		}

		q = q_dirTE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx1 + vx3)*(vx1 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[BW])[kbw] = (c1o1 - q) / (c1o1 + q)*(f_TE - f_BW + (f_TE + f_BW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TE + f_BW)) / (c1o1 + q);
		}

		q = q_dirBW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx1 - vx3)*(-vx1 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[TE])[kte] = (c1o1 - q) / (c1o1 + q)*(f_BW - f_TE + (f_BW + f_TE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BW + f_TE)) / (c1o1 + q);
		}

		q = q_dirBE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx1 - vx3)*(vx1 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[TW])[ktw] = (c1o1 - q) / (c1o1 + q)*(f_BE - f_TW + (f_BE + f_TW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BE + f_TW)) / (c1o1 + q);
		}

		q = q_dirTW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx1 + vx3)*(-vx1 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[BE])[kbe] = (c1o1 - q) / (c1o1 + q)*(f_TW - f_BE + (f_TW + f_BE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TW + f_BE)) / (c1o1 + q);
		}

		q = q_dirTN[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx2 + vx3)*(vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[BS])[kbs] = (c1o1 - q) / (c1o1 + q)*(f_TN - f_BS + (f_TN + f_BS - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TN + f_BS)) / (c1o1 + q);
		}

		q = q_dirBS[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx2 - vx3)*(-vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[TN])[ktn] = (c1o1 - q) / (c1o1 + q)*(f_BS - f_TN + (f_BS + f_TN - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BS + f_TN)) / (c1o1 + q);
		}

		q = q_dirBN[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx2 - vx3)*(vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[TS])[kts] = (c1o1 - q) / (c1o1 + q)*(f_BN - f_TS + (f_BN + f_TS - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BN + f_TS)) / (c1o1 + q);
		}

		q = q_dirTS[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx2 + vx3)*(-vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[BN])[kbn] = (c1o1 - q) / (c1o1 + q)*(f_TS - f_BN + (f_TS + f_BN - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TS + f_BN)) / (c1o1 + q);
		}

		q = q_dirTNE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(vx1 + vx2 + vx3)*(vx1 + vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[BSW])[kbsw] = (c1o1 - q) / (c1o1 + q)*(f_TNE - f_BSW + (f_TNE + f_BSW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TNE + f_BSW)) / (c1o1 + q);
		}

		q = q_dirBSW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(-vx1 - vx2 - vx3)*(-vx1 - vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[TNE])[ktne] = (c1o1 - q) / (c1o1 + q)*(f_BSW - f_TNE + (f_BSW + f_TNE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BSW + f_TNE)) / (c1o1 + q);
		}

		q = q_dirBNE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(vx1 + vx2 - vx3)*(vx1 + vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[TSW])[ktsw] = (c1o1 - q) / (c1o1 + q)*(f_BNE - f_TSW + (f_BNE + f_TSW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BNE + f_TSW)) / (c1o1 + q);
		}

		q = q_dirTSW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(-vx1 - vx2 + vx3)*(-vx1 - vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[BNE])[kbne] = (c1o1 - q) / (c1o1 + q)*(f_TSW - f_BNE + (f_TSW + f_BNE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TSW + f_BNE)) / (c1o1 + q);
		}

		q = q_dirTSE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(vx1 - vx2 + vx3)*(vx1 - vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[BNW])[kbnw] = (c1o1 - q) / (c1o1 + q)*(f_TSE - f_BNW + (f_TSE + f_BNW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TSE + f_BNW)) / (c1o1 + q);
		}

		q = q_dirBNW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(-vx1 + vx2 - vx3)*(-vx1 + vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[TSE])[ktse] = (c1o1 - q) / (c1o1 + q)*(f_BNW - f_TSE + (f_BNW + f_TSE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BNW + f_TSE)) / (c1o1 + q);
		}

		q = q_dirBSE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(vx1 - vx2 - vx3)*(vx1 - vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[TNW])[ktnw] = (c1o1 - q) / (c1o1 + q)*(f_BSE - f_TNW + (f_BSE + f_TNW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BSE + f_TNW)) / (c1o1 + q);
		}

		q = q_dirTNW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(-vx1 + vx2 + vx3)*(-vx1 + vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[BSE])[kbse] = (c1o1 - q) / (c1o1 + q)*(f_TNW - f_BSE + (f_TNW + f_BSE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TNW + f_BSE)) / (c1o1 + q);
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QThinWallsPartTwo27(
	real* DD, 
	int* k_Q, 
	real* QQ,
	uint numberOfBCnodes, 
	uint* geom,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint* neighborWSB,
	uint size_Mat, 
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
      q_dirE   = &QQ[E   * numberOfBCnodes];
      q_dirW   = &QQ[W   * numberOfBCnodes];
      q_dirN   = &QQ[N   * numberOfBCnodes];
      q_dirS   = &QQ[S   * numberOfBCnodes];
      q_dirT   = &QQ[T   * numberOfBCnodes];
      q_dirB   = &QQ[B   * numberOfBCnodes];
      q_dirNE  = &QQ[NE  * numberOfBCnodes];
      q_dirSW  = &QQ[SW  * numberOfBCnodes];
      q_dirSE  = &QQ[SE  * numberOfBCnodes];
      q_dirNW  = &QQ[NW  * numberOfBCnodes];
      q_dirTE  = &QQ[TE  * numberOfBCnodes];
      q_dirBW  = &QQ[BW  * numberOfBCnodes];
      q_dirBE  = &QQ[BE  * numberOfBCnodes];
      q_dirTW  = &QQ[TW  * numberOfBCnodes];
      q_dirTN  = &QQ[TN  * numberOfBCnodes];
      q_dirBS  = &QQ[BS  * numberOfBCnodes];
      q_dirBN  = &QQ[BN  * numberOfBCnodes];
      q_dirTS  = &QQ[TS  * numberOfBCnodes];
      q_dirTNE = &QQ[TNE * numberOfBCnodes];
      q_dirTSW = &QQ[TSW * numberOfBCnodes];
      q_dirTSE = &QQ[TSE * numberOfBCnodes];
      q_dirTNW = &QQ[TNW * numberOfBCnodes];
      q_dirBNE = &QQ[BNE * numberOfBCnodes];
      q_dirBSW = &QQ[BSW * numberOfBCnodes];
      q_dirBSE = &QQ[BSE * numberOfBCnodes];
      q_dirBNW = &QQ[BNW * numberOfBCnodes];
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
		  D.f[E] = &DD[E   *size_Mat];
		  D.f[W] = &DD[W   *size_Mat];
		  D.f[N] = &DD[N   *size_Mat];
		  D.f[S] = &DD[S   *size_Mat];
		  D.f[T] = &DD[T   *size_Mat];
		  D.f[B] = &DD[B   *size_Mat];
		  D.f[NE] = &DD[NE  *size_Mat];
		  D.f[SW] = &DD[SW  *size_Mat];
		  D.f[SE] = &DD[SE  *size_Mat];
		  D.f[NW] = &DD[NW  *size_Mat];
		  D.f[TE] = &DD[TE  *size_Mat];
		  D.f[BW] = &DD[BW  *size_Mat];
		  D.f[BE] = &DD[BE  *size_Mat];
		  D.f[TW] = &DD[TW  *size_Mat];
		  D.f[TN] = &DD[TN  *size_Mat];
		  D.f[BS] = &DD[BS  *size_Mat];
		  D.f[BN] = &DD[BN  *size_Mat];
		  D.f[TS] = &DD[TS  *size_Mat];
		  D.f[REST] = &DD[REST*size_Mat];
		  D.f[TNE] = &DD[TNE *size_Mat];
		  D.f[TSW] = &DD[TSW *size_Mat];
		  D.f[TSE] = &DD[TSE *size_Mat];
		  D.f[TNW] = &DD[TNW *size_Mat];
		  D.f[BNE] = &DD[BNE *size_Mat];
		  D.f[BSW] = &DD[BSW *size_Mat];
		  D.f[BSE] = &DD[BSE *size_Mat];
		  D.f[BNW] = &DD[BNW *size_Mat];
	  }
	  else
	  {
		  D.f[W] = &DD[E   *size_Mat];
		  D.f[E] = &DD[W   *size_Mat];
		  D.f[S] = &DD[N   *size_Mat];
		  D.f[N] = &DD[S   *size_Mat];
		  D.f[B] = &DD[T   *size_Mat];
		  D.f[T] = &DD[B   *size_Mat];
		  D.f[SW] = &DD[NE  *size_Mat];
		  D.f[NE] = &DD[SW  *size_Mat];
		  D.f[NW] = &DD[SE  *size_Mat];
		  D.f[SE] = &DD[NW  *size_Mat];
		  D.f[BW] = &DD[TE  *size_Mat];
		  D.f[TE] = &DD[BW  *size_Mat];
		  D.f[TW] = &DD[BE  *size_Mat];
		  D.f[BE] = &DD[TW  *size_Mat];
		  D.f[BS] = &DD[TN  *size_Mat];
		  D.f[TN] = &DD[BS  *size_Mat];
		  D.f[TS] = &DD[BN  *size_Mat];
		  D.f[BN] = &DD[TS  *size_Mat];
		  D.f[REST] = &DD[REST*size_Mat];
		  D.f[TNE] = &DD[BSW *size_Mat];
		  D.f[TSW] = &DD[BNE *size_Mat];
		  D.f[TSE] = &DD[BNW *size_Mat];
		  D.f[TNW] = &DD[BSE *size_Mat];
		  D.f[BNE] = &DD[TSW *size_Mat];
		  D.f[BSW] = &DD[TNE *size_Mat];
		  D.f[BSE] = &DD[TNW *size_Mat];
		  D.f[BNW] = &DD[TSE *size_Mat];
	  }
	  if (isEvenTimestep==false)
      {
         DN.f[E   ] = &DD[E   *size_Mat];
         DN.f[W   ] = &DD[W   *size_Mat];
         DN.f[N   ] = &DD[N   *size_Mat];
         DN.f[S   ] = &DD[S   *size_Mat];
         DN.f[T   ] = &DD[T   *size_Mat];
         DN.f[B   ] = &DD[B   *size_Mat];
         DN.f[NE  ] = &DD[NE  *size_Mat];
         DN.f[SW  ] = &DD[SW  *size_Mat];
         DN.f[SE  ] = &DD[SE  *size_Mat];
         DN.f[NW  ] = &DD[NW  *size_Mat];
         DN.f[TE  ] = &DD[TE  *size_Mat];
         DN.f[BW  ] = &DD[BW  *size_Mat];
         DN.f[BE  ] = &DD[BE  *size_Mat];
         DN.f[TW  ] = &DD[TW  *size_Mat];
         DN.f[TN  ] = &DD[TN  *size_Mat];
         DN.f[BS  ] = &DD[BS  *size_Mat];
         DN.f[BN  ] = &DD[BN  *size_Mat];
         DN.f[TS  ] = &DD[TS  *size_Mat];
         DN.f[REST] = &DD[REST*size_Mat];
         DN.f[TNE ] = &DD[TNE *size_Mat];
         DN.f[TSW ] = &DD[TSW *size_Mat];
         DN.f[TSE ] = &DD[TSE *size_Mat];
         DN.f[TNW ] = &DD[TNW *size_Mat];
         DN.f[BNE ] = &DD[BNE *size_Mat];
         DN.f[BSW ] = &DD[BSW *size_Mat];
         DN.f[BSE ] = &DD[BSE *size_Mat];
         DN.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         DN.f[W   ] = &DD[E   *size_Mat];
         DN.f[E   ] = &DD[W   *size_Mat];
         DN.f[S   ] = &DD[N   *size_Mat];
         DN.f[N   ] = &DD[S   *size_Mat];
         DN.f[B   ] = &DD[T   *size_Mat];
         DN.f[T   ] = &DD[B   *size_Mat];
         DN.f[SW  ] = &DD[NE  *size_Mat];
         DN.f[NE  ] = &DD[SW  *size_Mat];
         DN.f[NW  ] = &DD[SE  *size_Mat];
         DN.f[SE  ] = &DD[NW  *size_Mat];
         DN.f[BW  ] = &DD[TE  *size_Mat];
         DN.f[TE  ] = &DD[BW  *size_Mat];
         DN.f[TW  ] = &DD[BE  *size_Mat];
         DN.f[BE  ] = &DD[TW  *size_Mat];
         DN.f[BS  ] = &DD[TN  *size_Mat];
         DN.f[TN  ] = &DD[BS  *size_Mat];
         DN.f[TS  ] = &DD[BN  *size_Mat];
         DN.f[BN  ] = &DD[TS  *size_Mat];
         DN.f[REST] = &DD[REST*size_Mat];
         DN.f[TNE ] = &DD[BSW *size_Mat];
         DN.f[TSW ] = &DD[BNE *size_Mat];
         DN.f[TSE ] = &DD[BNW *size_Mat];
         DN.f[TNW ] = &DD[BSE *size_Mat];
         DN.f[BNE ] = &DD[TSW *size_Mat];
         DN.f[BSW ] = &DD[TNE *size_Mat];
         DN.f[BSE ] = &DD[TNW *size_Mat];
         DN.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //directions allways exchange
	  //(-1 -1 -1) (-1  0  0) ( 0 -1  0) ( 0  0 -1) (-1 -1  0) (-1  0 -1) ( 0 -1 -1) ( 1  1 -1) ( 1 -1  1) (-1  1  1) ( 1 -1  0) ( 1  0 -1) ( 0  1 -1)
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //directions exchange if solid neighbor
	  //( 1  1  1) ( 1  0  0) ( 0  1  0) ( 0  0  1) ( 1  1  0) ( 1  0  1) ( 0  1  1) (-1 -1  1) (-1  1 -1) ( 1 -1 -1) (-1  1  0) (-1  0  1) ( 0 -1  1)
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real q, tmp;
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1){ if (geom[kw  ] < GEO_FLUID){tmp = (DN.f[W  ])[kw  ]; (DN.f[W  ])[kw  ]=(D.f[W  ])[kw  ]; (D.f[W  ])[kw  ]=tmp;}}
	  q = q_dirW[k];   if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[E  ])[ke  ]; (DN.f[E  ])[ke  ]=(D.f[E  ])[ke  ]; (D.f[E  ])[ke  ]=tmp;}}
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1){ if (geom[ks  ] < GEO_FLUID){tmp = (DN.f[S  ])[ks  ]; (DN.f[S  ])[ks  ]=(D.f[S  ])[ks  ]; (D.f[S  ])[ks  ]=tmp;}}
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[N  ])[kn  ]; (DN.f[N  ])[kn  ]=(D.f[N  ])[kn  ]; (D.f[N  ])[kn  ]=tmp;}}
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1){ if (geom[kb  ] < GEO_FLUID){tmp = (DN.f[B  ])[kb  ]; (DN.f[B  ])[kb  ]=(D.f[B  ])[kb  ]; (D.f[B  ])[kb  ]=tmp;}}
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[T  ])[kt  ]; (DN.f[T  ])[kt  ]=(D.f[T  ])[kt  ]; (D.f[T  ])[kt  ]=tmp;}}
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1){ if (geom[ksw ] < GEO_FLUID){tmp = (DN.f[SW ])[ksw ]; (DN.f[SW ])[ksw ]=(D.f[SW ])[ksw ]; (D.f[SW ])[ksw ]=tmp;}}
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[NE ])[kne ]; (DN.f[NE ])[kne ]=(D.f[NE ])[kne ]; (D.f[NE ])[kne ]=tmp;}}
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[NW ])[knw ]; (DN.f[NW ])[knw ]=(D.f[NW ])[knw ]; (D.f[NW ])[knw ]=tmp;}}
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1){ if (geom[kmp0] < GEO_FLUID){tmp = (DN.f[SE ])[kse ]; (DN.f[SE ])[kse ]=(D.f[SE ])[kse ]; (D.f[SE ])[kse ]=tmp;}}
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1){ if (geom[kbw ] < GEO_FLUID){tmp = (DN.f[BW ])[kbw ]; (DN.f[BW ])[kbw ]=(D.f[BW ])[kbw ]; (D.f[BW ])[kbw ]=tmp;}}
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[TE ])[kte ]; (DN.f[TE ])[kte ]=(D.f[TE ])[kte ]; (D.f[TE ])[kte ]=tmp;}}
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[TW ])[ktw ]; (DN.f[TW ])[ktw ]=(D.f[TW ])[ktw ]; (D.f[TW ])[ktw ]=tmp;}}
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1){ if (geom[km0p] < GEO_FLUID){tmp = (DN.f[BE ])[kbe ]; (DN.f[BE ])[kbe ]=(D.f[BE ])[kbe ]; (D.f[BE ])[kbe ]=tmp;}}
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1){ if (geom[kbs ] < GEO_FLUID){tmp = (DN.f[BS ])[kbs ]; (DN.f[BS ])[kbs ]=(D.f[BS ])[kbs ]; (D.f[BS ])[kbs ]=tmp;}}
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[TN ])[ktn ]; (DN.f[TN ])[ktn ]=(D.f[TN ])[ktn ]; (D.f[TN ])[ktn ]=tmp;}}
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[TS ])[kts ]; (DN.f[TS ])[kts ]=(D.f[TS ])[kts ]; (D.f[TS ])[kts ]=tmp;}}
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1){ if (geom[k0mp] < GEO_FLUID){tmp = (DN.f[BN ])[kbn ]; (DN.f[BN ])[kbn ]=(D.f[BN ])[kbn ]; (D.f[BN ])[kbn ]=tmp;}}
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1){ if (geom[kbsw] < GEO_FLUID){tmp = (DN.f[BSW])[kbsw]; (DN.f[BSW])[kbsw]=(D.f[BSW])[kbsw]; (D.f[BSW])[kbsw]=tmp;}}
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[TNE])[ktne]; (DN.f[TNE])[ktne]=(D.f[TNE])[ktne]; (D.f[TNE])[ktne]=tmp;}}
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[TSW])[ktsw]; (DN.f[TSW])[ktsw]=(D.f[TSW])[ktsw]; (D.f[TSW])[ktsw]=tmp;}}
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1){ if (geom[kmmp] < GEO_FLUID){tmp = (DN.f[BNE])[kbne]; (DN.f[BNE])[kbne]=(D.f[BNE])[kbne]; (D.f[BNE])[kbne]=tmp;}}
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[BNW])[kbnw]; (DN.f[BNW])[kbnw]=(D.f[BNW])[kbnw]; (D.f[BNW])[kbnw]=tmp;}}
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1){ if (geom[kmpm] < GEO_FLUID){tmp = (DN.f[TSE])[ktse]; (DN.f[TSE])[ktse]=(D.f[TSE])[ktse]; (D.f[TSE])[ktse]=tmp;}}
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1){ if (geom[kpmm] < GEO_FLUID){tmp = (DN.f[TNW])[ktnw]; (DN.f[TNW])[ktnw]=(D.f[TNW])[ktnw]; (D.f[TNW])[ktnw]=tmp;}}
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[BSE])[kbse]; (DN.f[BSE])[kbse]=(D.f[BSE])[kbse]; (D.f[BSE])[kbse]=tmp;}}
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
