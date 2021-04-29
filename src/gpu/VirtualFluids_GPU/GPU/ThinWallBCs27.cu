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

/////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceCompThinWallsPartOne27(
	real* vx,
	real* vy,
	real* vz,
	real* DD, 
	int* k_Q, 
	real* QQ,
	uint sizeQ,
	int kQ, 
	real om1, 
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint size_Mat, 
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

   if(k<kQ)
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
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
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

      f_W    = (D.f[dirE   ])[ke   ];
      f_E    = (D.f[dirW   ])[kw   ];
      f_S    = (D.f[dirN   ])[kn   ];
      f_N    = (D.f[dirS   ])[ks   ];
      f_B    = (D.f[dirT   ])[kt   ];
      f_T    = (D.f[dirB   ])[kb   ];
      f_SW   = (D.f[dirNE  ])[kne  ];
      f_NE   = (D.f[dirSW  ])[ksw  ];
      f_NW   = (D.f[dirSE  ])[kse  ];
      f_SE   = (D.f[dirNW  ])[knw  ];
      f_BW   = (D.f[dirTE  ])[kte  ];
      f_TE   = (D.f[dirBW  ])[kbw  ];
      f_TW   = (D.f[dirBE  ])[kbe  ];
      f_BE   = (D.f[dirTW  ])[ktw  ];
      f_BS   = (D.f[dirTN  ])[ktn  ];
      f_TN   = (D.f[dirBS  ])[kbs  ];
      f_TS   = (D.f[dirBN  ])[kbn  ];
      f_BN   = (D.f[dirTS  ])[kts  ];
      f_BSW  = (D.f[dirTNE ])[ktne ];
      f_BNE  = (D.f[dirTSW ])[ktsw ];
      f_BNW  = (D.f[dirTSE ])[ktse ];
      f_BSE  = (D.f[dirTNW ])[ktnw ];
      f_TSW  = (D.f[dirBNE ])[kbne ];
      f_TNE  = (D.f[dirBSW ])[kbsw ];
      f_TNW  = (D.f[dirBSE ])[kbse ];
      f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

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
		 (D.f[dirW])[kw] = (c1o1 - q) / (c1o1 + q)*(f_E - f_W + (f_E + f_W - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_E + f_W) - c6o1*c2o27*(VeloX)) / (c1o1 + q);
	  }

	  q = q_dirW[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (-vx1)*(-vx1) * (c1o1 + drho) - cu_sq);
		  (D.f[dirE])[ke] = (c1o1 - q) / (c1o1 + q)*(f_W - f_E + (f_W + f_E - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_W + f_E) - c6o1*c2o27*(-VeloX)) / (c1o1 + q);
	  }

	  q = q_dirN[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (vx2)*(vx2) * (c1o1 + drho) - cu_sq);
		  (D.f[dirS])[ks] = (c1o1 - q) / (c1o1 + q)*(f_N - f_S + (f_N + f_S - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_N + f_S) - c6o1*c2o27*(VeloY)) / (c1o1 + q);
	  }

	  q = q_dirS[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (-vx2)*(-vx2) * (c1o1 + drho) - cu_sq);
		  (D.f[dirN])[kn] = (c1o1 - q) / (c1o1 + q)*(f_S - f_N + (f_S + f_N - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_S + f_N) - c6o1*c2o27*(-VeloY)) / (c1o1 + q);
	  }

	  q = q_dirT[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (vx3)*(vx3) * (c1o1 + drho) - cu_sq);
		  (D.f[dirB])[kb] = (c1o1 - q) / (c1o1 + q)*(f_T - f_B + (f_T + f_B - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_T + f_B) - c6o1*c2o27*(VeloZ)) / (c1o1 + q);
	  }

	  q = q_dirB[k];
	  if (q >= c0o1 && q <= c1o1)
	  {
		  feq = c2o27* (drho + c9o2 * (-vx3)*(-vx3) * (c1o1 + drho) - cu_sq);
		  (D.f[dirT])[kt] = (c1o1 - q) / (c1o1 + q)*(f_B - f_T + (f_B + f_T - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_B + f_T) - c6o1*c2o27*(-VeloZ)) / (c1o1 + q);
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * ( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq);
         (D.f[dirSW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*( VeloX+VeloY))/(c1o1+q);
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq);
         (D.f[dirNE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * ( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq);
         (D.f[dirNW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq);
         (D.f[dirSE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * ( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq);
         (D.f[dirBW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq);
         (D.f[dirTE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * ( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq);
         (D.f[dirTW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq);
         (D.f[dirBE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq);
         (D.f[dirBS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq);
         (D.f[dirTN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*(-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq);
         (D.f[dirTS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho + c9o2 * (    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq);
         (D.f[dirBN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*(-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * ( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * (-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * ( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * (-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * ( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * (-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * ( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho + c9o2 * (-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDeviceCompThinWallsPartOne27(
	real* DD,
	int* k_Q,
	real* QQ,
	unsigned int sizeQ,
	int kQ,
	real om1,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int size_Mat,
	bool evenOrOdd)
{
	Distributions27 D;
	if (evenOrOdd == true)
	{
		D.f[dirE] = &DD[dirE   *size_Mat];
		D.f[dirW] = &DD[dirW   *size_Mat];
		D.f[dirN] = &DD[dirN   *size_Mat];
		D.f[dirS] = &DD[dirS   *size_Mat];
		D.f[dirT] = &DD[dirT   *size_Mat];
		D.f[dirB] = &DD[dirB   *size_Mat];
		D.f[dirNE] = &DD[dirNE  *size_Mat];
		D.f[dirSW] = &DD[dirSW  *size_Mat];
		D.f[dirSE] = &DD[dirSE  *size_Mat];
		D.f[dirNW] = &DD[dirNW  *size_Mat];
		D.f[dirTE] = &DD[dirTE  *size_Mat];
		D.f[dirBW] = &DD[dirBW  *size_Mat];
		D.f[dirBE] = &DD[dirBE  *size_Mat];
		D.f[dirTW] = &DD[dirTW  *size_Mat];
		D.f[dirTN] = &DD[dirTN  *size_Mat];
		D.f[dirBS] = &DD[dirBS  *size_Mat];
		D.f[dirBN] = &DD[dirBN  *size_Mat];
		D.f[dirTS] = &DD[dirTS  *size_Mat];
		D.f[dirZERO] = &DD[dirZERO*size_Mat];
		D.f[dirTNE] = &DD[dirTNE *size_Mat];
		D.f[dirTSW] = &DD[dirTSW *size_Mat];
		D.f[dirTSE] = &DD[dirTSE *size_Mat];
		D.f[dirTNW] = &DD[dirTNW *size_Mat];
		D.f[dirBNE] = &DD[dirBNE *size_Mat];
		D.f[dirBSW] = &DD[dirBSW *size_Mat];
		D.f[dirBSE] = &DD[dirBSE *size_Mat];
		D.f[dirBNW] = &DD[dirBNW *size_Mat];
	}
	else
	{
		D.f[dirW] = &DD[dirE   *size_Mat];
		D.f[dirE] = &DD[dirW   *size_Mat];
		D.f[dirS] = &DD[dirN   *size_Mat];
		D.f[dirN] = &DD[dirS   *size_Mat];
		D.f[dirB] = &DD[dirT   *size_Mat];
		D.f[dirT] = &DD[dirB   *size_Mat];
		D.f[dirSW] = &DD[dirNE  *size_Mat];
		D.f[dirNE] = &DD[dirSW  *size_Mat];
		D.f[dirNW] = &DD[dirSE  *size_Mat];
		D.f[dirSE] = &DD[dirNW  *size_Mat];
		D.f[dirBW] = &DD[dirTE  *size_Mat];
		D.f[dirTE] = &DD[dirBW  *size_Mat];
		D.f[dirTW] = &DD[dirBE  *size_Mat];
		D.f[dirBE] = &DD[dirTW  *size_Mat];
		D.f[dirBS] = &DD[dirTN  *size_Mat];
		D.f[dirTN] = &DD[dirBS  *size_Mat];
		D.f[dirTS] = &DD[dirBN  *size_Mat];
		D.f[dirBN] = &DD[dirTS  *size_Mat];
		D.f[dirZERO] = &DD[dirZERO*size_Mat];
		D.f[dirTNE] = &DD[dirBSW *size_Mat];
		D.f[dirTSW] = &DD[dirBNE *size_Mat];
		D.f[dirTSE] = &DD[dirBNW *size_Mat];
		D.f[dirTNW] = &DD[dirBSE *size_Mat];
		D.f[dirBNE] = &DD[dirTSW *size_Mat];
		D.f[dirBSW] = &DD[dirTNE *size_Mat];
		D.f[dirBSE] = &DD[dirTNW *size_Mat];
		D.f[dirBNW] = &DD[dirTSE *size_Mat];
	}
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if (k < kQ)
	{
		////////////////////////////////////////////////////////////////////////////////
		real *q_dirE, *q_dirW, *q_dirN, *q_dirS, *q_dirT, *q_dirB,
			*q_dirNE, *q_dirSW, *q_dirSE, *q_dirNW, *q_dirTE, *q_dirBW,
			*q_dirBE, *q_dirTW, *q_dirTN, *q_dirBS, *q_dirBN, *q_dirTS,
			*q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			*q_dirBSE, *q_dirBNW;
		q_dirE = &QQ[dirE   *sizeQ];
		q_dirW = &QQ[dirW   *sizeQ];
		q_dirN = &QQ[dirN   *sizeQ];
		q_dirS = &QQ[dirS   *sizeQ];
		q_dirT = &QQ[dirT   *sizeQ];
		q_dirB = &QQ[dirB   *sizeQ];
		q_dirNE = &QQ[dirNE  *sizeQ];
		q_dirSW = &QQ[dirSW  *sizeQ];
		q_dirSE = &QQ[dirSE  *sizeQ];
		q_dirNW = &QQ[dirNW  *sizeQ];
		q_dirTE = &QQ[dirTE  *sizeQ];
		q_dirBW = &QQ[dirBW  *sizeQ];
		q_dirBE = &QQ[dirBE  *sizeQ];
		q_dirTW = &QQ[dirTW  *sizeQ];
		q_dirTN = &QQ[dirTN  *sizeQ];
		q_dirBS = &QQ[dirBS  *sizeQ];
		q_dirBN = &QQ[dirBN  *sizeQ];
		q_dirTS = &QQ[dirTS  *sizeQ];
		q_dirTNE = &QQ[dirTNE *sizeQ];
		q_dirTSW = &QQ[dirTSW *sizeQ];
		q_dirTSE = &QQ[dirTSE *sizeQ];
		q_dirTNW = &QQ[dirTNW *sizeQ];
		q_dirBNE = &QQ[dirBNE *sizeQ];
		q_dirBSW = &QQ[dirBSW *sizeQ];
		q_dirBSE = &QQ[dirBSE *sizeQ];
		q_dirBNW = &QQ[dirBNW *sizeQ];
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

		f_W = (D.f[dirE])[ke];
		f_E = (D.f[dirW])[kw];
		f_S = (D.f[dirN])[kn];
		f_N = (D.f[dirS])[ks];
		f_B = (D.f[dirT])[kt];
		f_T = (D.f[dirB])[kb];
		f_SW = (D.f[dirNE])[kne];
		f_NE = (D.f[dirSW])[ksw];
		f_NW = (D.f[dirSE])[kse];
		f_SE = (D.f[dirNW])[knw];
		f_BW = (D.f[dirTE])[kte];
		f_TE = (D.f[dirBW])[kbw];
		f_TW = (D.f[dirBE])[kbe];
		f_BE = (D.f[dirTW])[ktw];
		f_BS = (D.f[dirTN])[ktn];
		f_TN = (D.f[dirBS])[kbs];
		f_TS = (D.f[dirBN])[kbn];
		f_BN = (D.f[dirTS])[kts];
		f_BSW = (D.f[dirTNE])[ktne];
		f_BNE = (D.f[dirTSW])[ktsw];
		f_BNW = (D.f[dirTSE])[ktse];
		f_BSE = (D.f[dirTNW])[ktnw];
		f_TSW = (D.f[dirBNE])[kbne];
		f_TNE = (D.f[dirBSW])[kbsw];
		f_TNW = (D.f[dirBSE])[kbse];
		f_TSE = (D.f[dirBNW])[kbnw];
		////////////////////////////////////////////////////////////////////////////////
		real vx1, vx2, vx3, drho, feq, q;
		drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
			f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
			f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]);

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
			(D.f[dirW])[kw] = (c1o1 - q) / (c1o1 + q)*(f_E - f_W + (f_E + f_W - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_E + f_W)) / (c1o1 + q);
		}

		q = q_dirW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(-vx1)*(-vx1) * (c1o1 + drho) - cu_sq);
			(D.f[dirE])[ke] = (c1o1 - q) / (c1o1 + q)*(f_W - f_E + (f_W + f_E - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_W + f_E)) / (c1o1 + q);
		}

		q = q_dirN[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(vx2)*(vx2) * (c1o1 + drho) - cu_sq);
			(D.f[dirS])[ks] = (c1o1 - q) / (c1o1 + q)*(f_N - f_S + (f_N + f_S - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_N + f_S)) / (c1o1 + q);
		}

		q = q_dirS[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(-vx2)*(-vx2) * (c1o1 + drho) - cu_sq);
			(D.f[dirN])[kn] = (c1o1 - q) / (c1o1 + q)*(f_S - f_N + (f_S + f_N - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_S + f_N)) / (c1o1 + q);
		}

		q = q_dirT[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(vx3)*(vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirB])[kb] = (c1o1 - q) / (c1o1 + q)*(f_T - f_B + (f_T + f_B - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_T + f_B)) / (c1o1 + q);
		}

		q = q_dirB[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c2o27* (drho + c9o2*(-vx3)*(-vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirT])[kt] = (c1o1 - q) / (c1o1 + q)*(f_B - f_T + (f_B + f_T - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_B + f_T)) / (c1o1 + q);
		}

		q = q_dirNE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx1 + vx2)*(vx1 + vx2) * (c1o1 + drho) - cu_sq);
			(D.f[dirSW])[ksw] = (c1o1 - q) / (c1o1 + q)*(f_NE - f_SW + (f_NE + f_SW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_NE + f_SW)) / (c1o1 + q);
		}

		q = q_dirSW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx1 - vx2)*(-vx1 - vx2) * (c1o1 + drho) - cu_sq);
			(D.f[dirNE])[kne] = (c1o1 - q) / (c1o1 + q)*(f_SW - f_NE + (f_SW + f_NE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_SW + f_NE)) / (c1o1 + q);
		}

		q = q_dirSE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx1 - vx2)*(vx1 - vx2) * (c1o1 + drho) - cu_sq);
			(D.f[dirNW])[knw] = (c1o1 - q) / (c1o1 + q)*(f_SE - f_NW + (f_SE + f_NW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_SE + f_NW)) / (c1o1 + q);
		}

		q = q_dirNW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx1 + vx2)*(-vx1 + vx2) * (c1o1 + drho) - cu_sq);
			(D.f[dirSE])[kse] = (c1o1 - q) / (c1o1 + q)*(f_NW - f_SE + (f_NW + f_SE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_NW + f_SE)) / (c1o1 + q);
		}

		q = q_dirTE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx1 + vx3)*(vx1 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirBW])[kbw] = (c1o1 - q) / (c1o1 + q)*(f_TE - f_BW + (f_TE + f_BW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TE + f_BW)) / (c1o1 + q);
		}

		q = q_dirBW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx1 - vx3)*(-vx1 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirTE])[kte] = (c1o1 - q) / (c1o1 + q)*(f_BW - f_TE + (f_BW + f_TE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BW + f_TE)) / (c1o1 + q);
		}

		q = q_dirBE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx1 - vx3)*(vx1 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirTW])[ktw] = (c1o1 - q) / (c1o1 + q)*(f_BE - f_TW + (f_BE + f_TW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BE + f_TW)) / (c1o1 + q);
		}

		q = q_dirTW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx1 + vx3)*(-vx1 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirBE])[kbe] = (c1o1 - q) / (c1o1 + q)*(f_TW - f_BE + (f_TW + f_BE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TW + f_BE)) / (c1o1 + q);
		}

		q = q_dirTN[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx2 + vx3)*(vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirBS])[kbs] = (c1o1 - q) / (c1o1 + q)*(f_TN - f_BS + (f_TN + f_BS - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TN + f_BS)) / (c1o1 + q);
		}

		q = q_dirBS[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx2 - vx3)*(-vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirTN])[ktn] = (c1o1 - q) / (c1o1 + q)*(f_BS - f_TN + (f_BS + f_TN - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BS + f_TN)) / (c1o1 + q);
		}

		q = q_dirBN[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(vx2 - vx3)*(vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirTS])[kts] = (c1o1 - q) / (c1o1 + q)*(f_BN - f_TS + (f_BN + f_TS - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BN + f_TS)) / (c1o1 + q);
		}

		q = q_dirTS[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o54* (drho + c9o2*(-vx2 + vx3)*(-vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirBN])[kbn] = (c1o1 - q) / (c1o1 + q)*(f_TS - f_BN + (f_TS + f_BN - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TS + f_BN)) / (c1o1 + q);
		}

		q = q_dirTNE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(vx1 + vx2 + vx3)*(vx1 + vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirBSW])[kbsw] = (c1o1 - q) / (c1o1 + q)*(f_TNE - f_BSW + (f_TNE + f_BSW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TNE + f_BSW)) / (c1o1 + q);
		}

		q = q_dirBSW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(-vx1 - vx2 - vx3)*(-vx1 - vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirTNE])[ktne] = (c1o1 - q) / (c1o1 + q)*(f_BSW - f_TNE + (f_BSW + f_TNE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BSW + f_TNE)) / (c1o1 + q);
		}

		q = q_dirBNE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(vx1 + vx2 - vx3)*(vx1 + vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirTSW])[ktsw] = (c1o1 - q) / (c1o1 + q)*(f_BNE - f_TSW + (f_BNE + f_TSW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BNE + f_TSW)) / (c1o1 + q);
		}

		q = q_dirTSW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(-vx1 - vx2 + vx3)*(-vx1 - vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirBNE])[kbne] = (c1o1 - q) / (c1o1 + q)*(f_TSW - f_BNE + (f_TSW + f_BNE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TSW + f_BNE)) / (c1o1 + q);
		}

		q = q_dirTSE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(vx1 - vx2 + vx3)*(vx1 - vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirBNW])[kbnw] = (c1o1 - q) / (c1o1 + q)*(f_TSE - f_BNW + (f_TSE + f_BNW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TSE + f_BNW)) / (c1o1 + q);
		}

		q = q_dirBNW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(-vx1 + vx2 - vx3)*(-vx1 + vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirTSE])[ktse] = (c1o1 - q) / (c1o1 + q)*(f_BNW - f_TSE + (f_BNW + f_TSE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BNW + f_TSE)) / (c1o1 + q);
		}

		q = q_dirBSE[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(vx1 - vx2 - vx3)*(vx1 - vx2 - vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirTNW])[ktnw] = (c1o1 - q) / (c1o1 + q)*(f_BSE - f_TNW + (f_BSE + f_TNW - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_BSE + f_TNW)) / (c1o1 + q);
		}

		q = q_dirTNW[k];
		if (q >= c0o1 && q <= c1o1)
		{
			feq = c1o216*(drho + c9o2*(-vx1 + vx2 + vx3)*(-vx1 + vx2 + vx3) * (c1o1 + drho) - cu_sq);
			(D.f[dirBSE])[kbse] = (c1o1 - q) / (c1o1 + q)*(f_TNW - f_BSE + (f_TNW + f_BSE - c2o1*feq*om1) / (c1o1 - om1))*c1o2 + (q*(f_TNW + f_BSE)) / (c1o1 + q);
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QThinWallsPartTwo27(
	real* DD, 
	int* k_Q, 
	real* QQ,
	uint sizeQ,
	int kQ, 
	uint* geom,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint* neighborWSB,
	uint size_Mat, 
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
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
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
	  if (evenOrOdd == true)
	  {
		  D.f[dirE] = &DD[dirE   *size_Mat];
		  D.f[dirW] = &DD[dirW   *size_Mat];
		  D.f[dirN] = &DD[dirN   *size_Mat];
		  D.f[dirS] = &DD[dirS   *size_Mat];
		  D.f[dirT] = &DD[dirT   *size_Mat];
		  D.f[dirB] = &DD[dirB   *size_Mat];
		  D.f[dirNE] = &DD[dirNE  *size_Mat];
		  D.f[dirSW] = &DD[dirSW  *size_Mat];
		  D.f[dirSE] = &DD[dirSE  *size_Mat];
		  D.f[dirNW] = &DD[dirNW  *size_Mat];
		  D.f[dirTE] = &DD[dirTE  *size_Mat];
		  D.f[dirBW] = &DD[dirBW  *size_Mat];
		  D.f[dirBE] = &DD[dirBE  *size_Mat];
		  D.f[dirTW] = &DD[dirTW  *size_Mat];
		  D.f[dirTN] = &DD[dirTN  *size_Mat];
		  D.f[dirBS] = &DD[dirBS  *size_Mat];
		  D.f[dirBN] = &DD[dirBN  *size_Mat];
		  D.f[dirTS] = &DD[dirTS  *size_Mat];
		  D.f[dirZERO] = &DD[dirZERO*size_Mat];
		  D.f[dirTNE] = &DD[dirTNE *size_Mat];
		  D.f[dirTSW] = &DD[dirTSW *size_Mat];
		  D.f[dirTSE] = &DD[dirTSE *size_Mat];
		  D.f[dirTNW] = &DD[dirTNW *size_Mat];
		  D.f[dirBNE] = &DD[dirBNE *size_Mat];
		  D.f[dirBSW] = &DD[dirBSW *size_Mat];
		  D.f[dirBSE] = &DD[dirBSE *size_Mat];
		  D.f[dirBNW] = &DD[dirBNW *size_Mat];
	  }
	  else
	  {
		  D.f[dirW] = &DD[dirE   *size_Mat];
		  D.f[dirE] = &DD[dirW   *size_Mat];
		  D.f[dirS] = &DD[dirN   *size_Mat];
		  D.f[dirN] = &DD[dirS   *size_Mat];
		  D.f[dirB] = &DD[dirT   *size_Mat];
		  D.f[dirT] = &DD[dirB   *size_Mat];
		  D.f[dirSW] = &DD[dirNE  *size_Mat];
		  D.f[dirNE] = &DD[dirSW  *size_Mat];
		  D.f[dirNW] = &DD[dirSE  *size_Mat];
		  D.f[dirSE] = &DD[dirNW  *size_Mat];
		  D.f[dirBW] = &DD[dirTE  *size_Mat];
		  D.f[dirTE] = &DD[dirBW  *size_Mat];
		  D.f[dirTW] = &DD[dirBE  *size_Mat];
		  D.f[dirBE] = &DD[dirTW  *size_Mat];
		  D.f[dirBS] = &DD[dirTN  *size_Mat];
		  D.f[dirTN] = &DD[dirBS  *size_Mat];
		  D.f[dirTS] = &DD[dirBN  *size_Mat];
		  D.f[dirBN] = &DD[dirTS  *size_Mat];
		  D.f[dirZERO] = &DD[dirZERO*size_Mat];
		  D.f[dirTNE] = &DD[dirBSW *size_Mat];
		  D.f[dirTSW] = &DD[dirBNE *size_Mat];
		  D.f[dirTSE] = &DD[dirBNW *size_Mat];
		  D.f[dirTNW] = &DD[dirBSE *size_Mat];
		  D.f[dirBNE] = &DD[dirTSW *size_Mat];
		  D.f[dirBSW] = &DD[dirTNE *size_Mat];
		  D.f[dirBSE] = &DD[dirTNW *size_Mat];
		  D.f[dirBNW] = &DD[dirTSE *size_Mat];
	  }
	  if (evenOrOdd==false)
      {
         DN.f[dirE   ] = &DD[dirE   *size_Mat];
         DN.f[dirW   ] = &DD[dirW   *size_Mat];
         DN.f[dirN   ] = &DD[dirN   *size_Mat];
         DN.f[dirS   ] = &DD[dirS   *size_Mat];
         DN.f[dirT   ] = &DD[dirT   *size_Mat];
         DN.f[dirB   ] = &DD[dirB   *size_Mat];
         DN.f[dirNE  ] = &DD[dirNE  *size_Mat];
         DN.f[dirSW  ] = &DD[dirSW  *size_Mat];
         DN.f[dirSE  ] = &DD[dirSE  *size_Mat];
         DN.f[dirNW  ] = &DD[dirNW  *size_Mat];
         DN.f[dirTE  ] = &DD[dirTE  *size_Mat];
         DN.f[dirBW  ] = &DD[dirBW  *size_Mat];
         DN.f[dirBE  ] = &DD[dirBE  *size_Mat];
         DN.f[dirTW  ] = &DD[dirTW  *size_Mat];
         DN.f[dirTN  ] = &DD[dirTN  *size_Mat];
         DN.f[dirBS  ] = &DD[dirBS  *size_Mat];
         DN.f[dirBN  ] = &DD[dirBN  *size_Mat];
         DN.f[dirTS  ] = &DD[dirTS  *size_Mat];
         DN.f[dirZERO] = &DD[dirZERO*size_Mat];
         DN.f[dirTNE ] = &DD[dirTNE *size_Mat];
         DN.f[dirTSW ] = &DD[dirTSW *size_Mat];
         DN.f[dirTSE ] = &DD[dirTSE *size_Mat];
         DN.f[dirTNW ] = &DD[dirTNW *size_Mat];
         DN.f[dirBNE ] = &DD[dirBNE *size_Mat];
         DN.f[dirBSW ] = &DD[dirBSW *size_Mat];
         DN.f[dirBSE ] = &DD[dirBSE *size_Mat];
         DN.f[dirBNW ] = &DD[dirBNW *size_Mat];
      } 
      else
      {
         DN.f[dirW   ] = &DD[dirE   *size_Mat];
         DN.f[dirE   ] = &DD[dirW   *size_Mat];
         DN.f[dirS   ] = &DD[dirN   *size_Mat];
         DN.f[dirN   ] = &DD[dirS   *size_Mat];
         DN.f[dirB   ] = &DD[dirT   *size_Mat];
         DN.f[dirT   ] = &DD[dirB   *size_Mat];
         DN.f[dirSW  ] = &DD[dirNE  *size_Mat];
         DN.f[dirNE  ] = &DD[dirSW  *size_Mat];
         DN.f[dirNW  ] = &DD[dirSE  *size_Mat];
         DN.f[dirSE  ] = &DD[dirNW  *size_Mat];
         DN.f[dirBW  ] = &DD[dirTE  *size_Mat];
         DN.f[dirTE  ] = &DD[dirBW  *size_Mat];
         DN.f[dirTW  ] = &DD[dirBE  *size_Mat];
         DN.f[dirBE  ] = &DD[dirTW  *size_Mat];
         DN.f[dirBS  ] = &DD[dirTN  *size_Mat];
         DN.f[dirTN  ] = &DD[dirBS  *size_Mat];
         DN.f[dirTS  ] = &DD[dirBN  *size_Mat];
         DN.f[dirBN  ] = &DD[dirTS  *size_Mat];
         DN.f[dirZERO] = &DD[dirZERO*size_Mat];
         DN.f[dirTNE ] = &DD[dirBSW *size_Mat];
         DN.f[dirTSW ] = &DD[dirBNE *size_Mat];
         DN.f[dirTSE ] = &DD[dirBNW *size_Mat];
         DN.f[dirTNW ] = &DD[dirBSE *size_Mat];
         DN.f[dirBNE ] = &DD[dirTSW *size_Mat];
         DN.f[dirBSW ] = &DD[dirTNE *size_Mat];
         DN.f[dirBSE ] = &DD[dirTNW *size_Mat];
         DN.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //directions allways exchange
	  //(-1 -1 -1) (-1  0  0) ( 0 -1  0) ( 0  0 -1) (-1 -1  0) (-1  0 -1) ( 0 -1 -1) ( 1  1 -1) ( 1 -1  1) (-1  1  1) ( 1 -1  0) ( 1  0 -1) ( 0  1 -1)
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //directions exchange if solid neighbor
	  //( 1  1  1) ( 1  0  0) ( 0  1  0) ( 0  0  1) ( 1  1  0) ( 1  0  1) ( 0  1  1) (-1 -1  1) (-1  1 -1) ( 1 -1 -1) (-1  1  0) (-1  0  1) ( 0 -1  1)
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real q, tmp;
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1){ if (geom[kw  ] < GEO_FLUID){tmp = (DN.f[dirW  ])[kw  ]; (DN.f[dirW  ])[kw  ]=(D.f[dirW  ])[kw  ]; (D.f[dirW  ])[kw  ]=tmp;}}
	  q = q_dirW[k];   if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirE  ])[ke  ]; (DN.f[dirE  ])[ke  ]=(D.f[dirE  ])[ke  ]; (D.f[dirE  ])[ke  ]=tmp;}}
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1){ if (geom[ks  ] < GEO_FLUID){tmp = (DN.f[dirS  ])[ks  ]; (DN.f[dirS  ])[ks  ]=(D.f[dirS  ])[ks  ]; (D.f[dirS  ])[ks  ]=tmp;}}
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirN  ])[kn  ]; (DN.f[dirN  ])[kn  ]=(D.f[dirN  ])[kn  ]; (D.f[dirN  ])[kn  ]=tmp;}}
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1){ if (geom[kb  ] < GEO_FLUID){tmp = (DN.f[dirB  ])[kb  ]; (DN.f[dirB  ])[kb  ]=(D.f[dirB  ])[kb  ]; (D.f[dirB  ])[kb  ]=tmp;}}
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirT  ])[kt  ]; (DN.f[dirT  ])[kt  ]=(D.f[dirT  ])[kt  ]; (D.f[dirT  ])[kt  ]=tmp;}}
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1){ if (geom[ksw ] < GEO_FLUID){tmp = (DN.f[dirSW ])[ksw ]; (DN.f[dirSW ])[ksw ]=(D.f[dirSW ])[ksw ]; (D.f[dirSW ])[ksw ]=tmp;}}
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirNE ])[kne ]; (DN.f[dirNE ])[kne ]=(D.f[dirNE ])[kne ]; (D.f[dirNE ])[kne ]=tmp;}}
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirNW ])[knw ]; (DN.f[dirNW ])[knw ]=(D.f[dirNW ])[knw ]; (D.f[dirNW ])[knw ]=tmp;}}
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1){ if (geom[kmp0] < GEO_FLUID){tmp = (DN.f[dirSE ])[kse ]; (DN.f[dirSE ])[kse ]=(D.f[dirSE ])[kse ]; (D.f[dirSE ])[kse ]=tmp;}}
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1){ if (geom[kbw ] < GEO_FLUID){tmp = (DN.f[dirBW ])[kbw ]; (DN.f[dirBW ])[kbw ]=(D.f[dirBW ])[kbw ]; (D.f[dirBW ])[kbw ]=tmp;}}
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirTE ])[kte ]; (DN.f[dirTE ])[kte ]=(D.f[dirTE ])[kte ]; (D.f[dirTE ])[kte ]=tmp;}}
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirTW ])[ktw ]; (DN.f[dirTW ])[ktw ]=(D.f[dirTW ])[ktw ]; (D.f[dirTW ])[ktw ]=tmp;}}
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1){ if (geom[km0p] < GEO_FLUID){tmp = (DN.f[dirBE ])[kbe ]; (DN.f[dirBE ])[kbe ]=(D.f[dirBE ])[kbe ]; (D.f[dirBE ])[kbe ]=tmp;}}
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1){ if (geom[kbs ] < GEO_FLUID){tmp = (DN.f[dirBS ])[kbs ]; (DN.f[dirBS ])[kbs ]=(D.f[dirBS ])[kbs ]; (D.f[dirBS ])[kbs ]=tmp;}}
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirTN ])[ktn ]; (DN.f[dirTN ])[ktn ]=(D.f[dirTN ])[ktn ]; (D.f[dirTN ])[ktn ]=tmp;}}
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirTS ])[kts ]; (DN.f[dirTS ])[kts ]=(D.f[dirTS ])[kts ]; (D.f[dirTS ])[kts ]=tmp;}}
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1){ if (geom[k0mp] < GEO_FLUID){tmp = (DN.f[dirBN ])[kbn ]; (DN.f[dirBN ])[kbn ]=(D.f[dirBN ])[kbn ]; (D.f[dirBN ])[kbn ]=tmp;}}
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1){ if (geom[kbsw] < GEO_FLUID){tmp = (DN.f[dirBSW])[kbsw]; (DN.f[dirBSW])[kbsw]=(D.f[dirBSW])[kbsw]; (D.f[dirBSW])[kbsw]=tmp;}}
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirTNE])[ktne]; (DN.f[dirTNE])[ktne]=(D.f[dirTNE])[ktne]; (D.f[dirTNE])[ktne]=tmp;}}
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirTSW])[ktsw]; (DN.f[dirTSW])[ktsw]=(D.f[dirTSW])[ktsw]; (D.f[dirTSW])[ktsw]=tmp;}}
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1){ if (geom[kmmp] < GEO_FLUID){tmp = (DN.f[dirBNE])[kbne]; (DN.f[dirBNE])[kbne]=(D.f[dirBNE])[kbne]; (D.f[dirBNE])[kbne]=tmp;}}
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirBNW])[kbnw]; (DN.f[dirBNW])[kbnw]=(D.f[dirBNW])[kbnw]; (D.f[dirBNW])[kbnw]=tmp;}}
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1){ if (geom[kmpm] < GEO_FLUID){tmp = (DN.f[dirTSE])[ktse]; (DN.f[dirTSE])[ktse]=(D.f[dirTSE])[ktse]; (D.f[dirTSE])[ktse]=tmp;}}
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1){ if (geom[kpmm] < GEO_FLUID){tmp = (DN.f[dirTNW])[ktnw]; (DN.f[dirTNW])[ktnw]=(D.f[dirTNW])[ktnw]; (D.f[dirTNW])[ktnw]=tmp;}}
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1){                            {tmp = (DN.f[dirBSE])[kbse]; (DN.f[dirBSE])[kbse]=(D.f[dirBSE])[kbse]; (D.f[dirBSE])[kbse]=tmp;}}
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
