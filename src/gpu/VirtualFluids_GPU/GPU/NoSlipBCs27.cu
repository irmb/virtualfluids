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
extern "C" __global__ void QDevice3rdMomentsComp27(  int inx,
													 int iny,
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

   if(k < numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      q_dirNE  = &QQ[dirNE  * numberOfBCnodes];
      q_dirSW  = &QQ[dirSW  * numberOfBCnodes];
      q_dirSE  = &QQ[dirSE  * numberOfBCnodes];
      q_dirNW  = &QQ[dirNW  * numberOfBCnodes];
      q_dirTE  = &QQ[dirTE  * numberOfBCnodes];
      q_dirBW  = &QQ[dirBW  * numberOfBCnodes];
      q_dirBE  = &QQ[dirBE  * numberOfBCnodes];
      q_dirTW  = &QQ[dirTW  * numberOfBCnodes];
      q_dirTN  = &QQ[dirTN  * numberOfBCnodes];
      q_dirBS  = &QQ[dirBS  * numberOfBCnodes];
      q_dirBN  = &QQ[dirBN  * numberOfBCnodes];
      q_dirTS  = &QQ[dirTS  * numberOfBCnodes];
      q_dirTNE = &QQ[dirTNE * numberOfBCnodes];
      q_dirTSW = &QQ[dirTSW * numberOfBCnodes];
      q_dirTSE = &QQ[dirTSE * numberOfBCnodes];
      q_dirTNW = &QQ[dirTNW * numberOfBCnodes];
      q_dirBNE = &QQ[dirBNE * numberOfBCnodes];
      q_dirBSW = &QQ[dirBSW * numberOfBCnodes];
      q_dirBSE = &QQ[dirBSE * numberOfBCnodes];
      q_dirBNW = &QQ[dirBNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int numberOfNodesK  = k_Q[k];
      unsigned int kzero= numberOfNodesK;
      unsigned int ke   = numberOfNodesK;
      unsigned int kw   = neighborX[numberOfNodesK];
      unsigned int kn   = numberOfNodesK;
      unsigned int ks   = neighborY[numberOfNodesK];
      unsigned int kt   = numberOfNodesK;
      unsigned int kb   = neighborZ[numberOfNodesK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = numberOfNodesK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = numberOfNodesK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = numberOfNodesK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = numberOfNodesK;
      unsigned int kbsw = neighborZ[ksw];
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
      real vx1, vx2, vx3, drho, feq, q, m3;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (c1o1 + drho); 


      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S)) / (c1o1 + drho); 

      vx3    =    (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                  (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                  (f_T - f_B)) / (c1o1 + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      //Test
         //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_E - f_W - c2o1 * drho * c2o27 * (c3o1*( vx1        ));
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirW])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W-m3+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_W - f_E - c2o1 * drho * c2o27 * (c3o1*(-vx1        ));
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirE])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E-m3+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_N - f_S - c2o1 * drho * c2o27 * (c3o1*( vx2        ));
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirS])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S-m3+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_S - f_N - c2o1 * drho * c2o27 * (c3o1*(   -vx2     ));
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirN])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N-m3+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_T - f_B - c2o1 * drho * c2o27 * (c3o1*(         vx3));
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirB])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B-m3+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_B - f_T - c2o1 * drho * c2o27 * (c3o1*(        -vx3));
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirT])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T-m3+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_NE - f_SW - c2o1 * drho * c1o54 * (c3o1*( vx1+vx2    ));
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW-m3+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_SW - f_NE - c2o1 * drho * c1o54 * (c3o1*(-vx1-vx2    ));
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirNE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE-m3+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_SE - f_NW - c2o1 * drho * c1o54 * (c3o1*( vx1-vx2    ));
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirNW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW-m3+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_NW - f_SE - c2o1 * drho * c1o54 * (c3o1*(-vx1+vx2    ));
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirSE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE-m3+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TE - f_BW - c2o1 * drho * c1o54 * (c3o1*( vx1    +vx3));
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW-m3+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BW - f_TE - c2o1 * drho * c1o54 * (c3o1*(-vx1    -vx3));
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE-m3+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BE - f_TW - c2o1 * drho * c1o54 * (c3o1*( vx1    -vx3));
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW-m3+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TW - f_BE - c2o1 * drho * c1o54 * (c3o1*(-vx1    +vx3));
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE-m3+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TN - f_BS - c2o1 * drho * c1o54 * (c3o1*(     vx2+vx3));
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS-m3+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BS - f_TN - c2o1 * drho * c1o54 * (c3o1*(    -vx2-vx3));
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN-m3+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BN - f_TS - c2o1 * drho * c1o54 * (c3o1*(     vx2-vx3));
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS-m3+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TS - f_BN - c2o1 * drho * c1o54 * (c3o1*(    -vx2+vx3));
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN-m3+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TNE - f_BSW - c2o1 * drho * c1o216 * (c3o1*( vx1+vx2+vx3));
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW-m3+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BSW - f_TNE - c2o1 * drho * c1o216 * (c3o1*(-vx1-vx2-vx3));
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE-m3+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BNE - f_TSW - c2o1 * drho * c1o216 * (c3o1*( vx1+vx2-vx3));
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW-m3+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TSW - f_BNE - c2o1 * drho * c1o216 * (c3o1*(-vx1-vx2+vx3));
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE-m3+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TSE - f_BNW - c2o1 * drho * c1o216 * (c3o1*( vx1-vx2+vx3));
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW-m3+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BNW - f_TSE - c2o1 * drho * c1o216 * (c3o1*(-vx1+vx2-vx3));
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE-m3+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BSE - f_TNW - c2o1 * drho * c1o216 * (c3o1*( vx1-vx2-vx3));
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW-m3+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TNW - f_BSE - c2o1 * drho * c1o216 * (c3o1*(-vx1+vx2+vx3));
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE-m3+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE))/(c1o1+q)+(m3*c1o2);
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDeviceIncompHighNu27(int inx,
												 int iny,
												 real* DD, 
												 int* k_Q, 
												 real* QQ,
												 unsigned int sizeQ,
												 int numberOfNodes, 
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

   if(k<numberOfNodes)
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
      unsigned int numberOfNodesK  = k_Q[k];
      unsigned int kzero= numberOfNodesK;
      unsigned int ke   = numberOfNodesK;
      unsigned int kw   = neighborX[numberOfNodesK];
      unsigned int kn   = numberOfNodesK;
      unsigned int ks   = neighborY[numberOfNodesK];
      unsigned int kt   = numberOfNodesK;
      unsigned int kb   = neighborZ[numberOfNodesK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = numberOfNodesK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = numberOfNodesK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = numberOfNodesK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = numberOfNodesK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
            f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_E   = (D.f[dirE   ])[ke   ];
      f_W   = (D.f[dirW   ])[kw   ];
      f_N   = (D.f[dirN   ])[kn   ];
      f_S   = (D.f[dirS   ])[ks   ];
      f_T   = (D.f[dirT   ])[kt   ];
      f_B   = (D.f[dirB   ])[kb   ];
      f_NE  = (D.f[dirNE  ])[kne  ];
      f_SW  = (D.f[dirSW  ])[ksw  ];
      f_SE  = (D.f[dirSE  ])[kse  ];
      f_NW  = (D.f[dirNW  ])[knw  ];
      f_TE  = (D.f[dirTE  ])[kte  ];
      f_BW  = (D.f[dirBW  ])[kbw  ];
      f_BE  = (D.f[dirBE  ])[kbe  ];
      f_TW  = (D.f[dirTW  ])[ktw  ];
      f_TN  = (D.f[dirTN  ])[ktn  ];
      f_BS  = (D.f[dirBS  ])[kbs  ];
      f_BN  = (D.f[dirBN  ])[kbn  ];
      f_TS  = (D.f[dirTS  ])[kts  ];
      f_TNE = (D.f[dirTNE ])[ktne ];
      f_TSW = (D.f[dirTSW ])[ktsw ];
      f_TSE = (D.f[dirTSE ])[ktse ];
      f_TNW = (D.f[dirTNW ])[ktnw ];
      f_BNE = (D.f[dirBNE ])[kbne ];
      f_BSW = (D.f[dirBSW ])[kbsw ];
      f_BSE = (D.f[dirBSE ])[kbse ];
      f_BNW = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W));// / (one + drho); 


      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S));// / (one + drho); 

      vx3    =    (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                  (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                  (f_T - f_B));// / (one + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);// * (one + drho);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      //Test
         //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real VeloX = c0o1;
      real VeloY = c0o1;
      real VeloZ = c0o1;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) /** (one + drho)*/-cu_sq); 
         (D.f[dirW])[kw]=((c1o1 - q) * f_E + q * ((f_E + f_W) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloX     )) / (q + c1o1) ;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) /** (one + drho)*/-cu_sq); 
         (D.f[dirE])[ke]=((c1o1 - q) * f_W + q * ((f_W + f_E) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloX     )) / (q + c1o1) ;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirS])[ks]=((c1o1 - q) * f_N + q * ((f_N + f_S) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloY     )) / (q + c1o1) ;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirN])[kn]=((c1o1 - q) * f_S + q * ((f_S + f_N) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloY     )) / (q + c1o1) ;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirB])[kb]=((c1o1 - q) * f_T + q * ((f_T + f_B) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloZ     )) / (q + c1o1) ;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirT])[kt]=((c1o1 - q) * f_B + q * ((f_B + f_T) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloZ     )) / (q + c1o1) ;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirSW])[ksw]=((c1o1 - q) * f_NE + q * ((f_NE + f_SW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloY)) / (q + c1o1) ;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirNE])[kne]=((c1o1 - q) * f_SW + q * ((f_SW + f_NE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloY)) / (q + c1o1) ;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirNW])[knw]=((c1o1 - q) * f_SE + q * ((f_SE + f_NW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloY)) / (q + c1o1) ;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirSE])[kse]=((c1o1 - q) * f_NW + q * ((f_NW + f_SE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloY)) / (q + c1o1) ;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBW])[kbw]=((c1o1 - q) * f_TE + q * ((f_TE + f_BW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTE])[kte]=((c1o1 - q) * f_BW + q * ((f_BW + f_TE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTW])[ktw]=((c1o1 - q) * f_BE + q * ((f_BE + f_TW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBE])[kbe]=((c1o1 - q) * f_TW + q * ((f_TW + f_BE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBS])[kbs]=((c1o1 - q) * f_TN + q * ((f_TN + f_BS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTN])[ktn]=((c1o1 - q) * f_BS + q * ((f_BS + f_TN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTS])[kts]=((c1o1 - q) * f_BN + q * ((f_BN + f_TS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBN])[kbn]=((c1o1 - q) * f_TS + q * ((f_TS + f_BN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBSW])[kbsw]=((c1o1 - q) * f_TNE + q * ((f_TNE + f_BSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTNE])[ktne]=((c1o1 - q) * f_BSW + q * ((f_BSW + f_TNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTSW])[ktsw]=((c1o1 - q) * f_BNE + q * ((f_BNE + f_TSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBNE])[kbne]=((c1o1 - q) * f_TSW + q * ((f_TSW + f_BNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBNW])[kbnw]=((c1o1 - q) * f_TSE + q * ((f_TSE + f_BNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTSE])[ktse]=((c1o1 - q) * f_BNW + q * ((f_BNW + f_TSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTNW])[ktnw]=((c1o1 - q) * f_BSE + q * ((f_BSE + f_TNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBSE])[kbse]=((c1o1 - q) * f_TNW + q * ((f_TNW + f_BSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY+VeloZ)) / (q + c1o1) ;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDeviceCompHighNu27(  int inx,
												 int iny,
												 real* DD, 
												 int* k_Q, 
												 real* QQ,
												 unsigned int sizeQ,
												 int numberOfNodes, 
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

   if(k<numberOfNodes)
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
      unsigned int numberOfNodesK  = k_Q[k];
      unsigned int kzero= numberOfNodesK;
      unsigned int ke   = numberOfNodesK;
      unsigned int kw   = neighborX[numberOfNodesK];
      unsigned int kn   = numberOfNodesK;
      unsigned int ks   = neighborY[numberOfNodesK];
      unsigned int kt   = numberOfNodesK;
      unsigned int kb   = neighborZ[numberOfNodesK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = numberOfNodesK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = numberOfNodesK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = numberOfNodesK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = numberOfNodesK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
            f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_E   = (D.f[dirE   ])[ke   ];
      f_W   = (D.f[dirW   ])[kw   ];
      f_N   = (D.f[dirN   ])[kn   ];
      f_S   = (D.f[dirS   ])[ks   ];
      f_T   = (D.f[dirT   ])[kt   ];
      f_B   = (D.f[dirB   ])[kb   ];
      f_NE  = (D.f[dirNE  ])[kne  ];
      f_SW  = (D.f[dirSW  ])[ksw  ];
      f_SE  = (D.f[dirSE  ])[kse  ];
      f_NW  = (D.f[dirNW  ])[knw  ];
      f_TE  = (D.f[dirTE  ])[kte  ];
      f_BW  = (D.f[dirBW  ])[kbw  ];
      f_BE  = (D.f[dirBE  ])[kbe  ];
      f_TW  = (D.f[dirTW  ])[ktw  ];
      f_TN  = (D.f[dirTN  ])[ktn  ];
      f_BS  = (D.f[dirBS  ])[kbs  ];
      f_BN  = (D.f[dirBN  ])[kbn  ];
      f_TS  = (D.f[dirTS  ])[kts  ];
      f_TNE = (D.f[dirTNE ])[ktne ];
      f_TSW = (D.f[dirTSW ])[ktsw ];
      f_TSE = (D.f[dirTSE ])[ktse ];
      f_TNW = (D.f[dirTNW ])[ktnw ];
      f_BNE = (D.f[dirBNE ])[kbne ];
      f_BSW = (D.f[dirBSW ])[kbsw ];
      f_BSE = (D.f[dirBSE ])[kbse ];
      f_BNW = (D.f[dirBNW ])[kbnw ];
      //f_W    = (D.f[dirE   ])[ke   ];
      //f_E    = (D.f[dirW   ])[kw   ];
      //f_S    = (D.f[dirN   ])[kn   ];
      //f_N    = (D.f[dirS   ])[ks   ];
      //f_B    = (D.f[dirT   ])[kt   ];
      //f_T    = (D.f[dirB   ])[kb   ];
      //f_SW   = (D.f[dirNE  ])[kne  ];
      //f_NE   = (D.f[dirSW  ])[ksw  ];
      //f_NW   = (D.f[dirSE  ])[kse  ];
      //f_SE   = (D.f[dirNW  ])[knw  ];
      //f_BW   = (D.f[dirTE  ])[kte  ];
      //f_TE   = (D.f[dirBW  ])[kbw  ];
      //f_TW   = (D.f[dirBE  ])[kbe  ];
      //f_BE   = (D.f[dirTW  ])[ktw  ];
      //f_BS   = (D.f[dirTN  ])[ktn  ];
      //f_TN   = (D.f[dirBS  ])[kbs  ];
      //f_TS   = (D.f[dirBN  ])[kbn  ];
      //f_BN   = (D.f[dirTS  ])[kts  ];
      //f_BSW  = (D.f[dirTNE ])[ktne ];
      //f_BNE  = (D.f[dirTSW ])[ktsw ];
      //f_BNW  = (D.f[dirTSE ])[ktse ];
      //f_BSE  = (D.f[dirTNW ])[ktnw ];
      //f_TSW  = (D.f[dirBNE ])[kbne ];
      //f_TNE  = (D.f[dirBSW ])[kbsw ];
      //f_TNW  = (D.f[dirBSE ])[kbse ];
      //f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (c1o1 + drho); 


      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S)) / (c1o1 + drho); 

      vx3    =    (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                  (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                  (f_T - f_B)) / (c1o1 + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      //Test
         //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real VeloX = c0o1;
      real VeloY = c0o1;
      real VeloZ = c0o1;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirW])[kw]=((c1o1 - q) * f_E + q * ((f_E + f_W) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloX     )) / (q + c1o1) ;
         //(D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirE])[ke]=((c1o1 - q) * f_W + q * ((f_W + f_E) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloX     )) / (q + c1o1) ;
         //(D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirS])[ks]=((c1o1 - q) * f_N + q * ((f_N + f_S) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloY     )) / (q + c1o1) ;
         //(D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirN])[kn]=((c1o1 - q) * f_S + q * ((f_S + f_N) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloY     )) / (q + c1o1) ;
         //(D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirB])[kb]=((c1o1 - q) * f_T + q * ((f_T + f_B) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloZ     )) / (q + c1o1) ;
         //(D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirT])[kt]=((c1o1 - q) * f_B + q * ((f_B + f_T) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloZ     )) / (q + c1o1) ;
         //(D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=((c1o1 - q) * f_NE + q * ((f_NE + f_SW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloY)) / (q + c1o1) ;
         //(D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirNE])[kne]=((c1o1 - q) * f_SW + q * ((f_SW + f_NE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloY)) / (q + c1o1) ;
         //(D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirNW])[knw]=((c1o1 - q) * f_SE + q * ((f_SE + f_NW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloY)) / (q + c1o1) ;
         //(D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirSE])[kse]=((c1o1 - q) * f_NW + q * ((f_NW + f_SE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloY)) / (q + c1o1) ;
         //(D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=((c1o1 - q) * f_TE + q * ((f_TE + f_BW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloZ)) / (q + c1o1) ;
         //(D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTE])[kte]=((c1o1 - q) * f_BW + q * ((f_BW + f_TE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloZ)) / (q + c1o1) ;
         //(D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=((c1o1 - q) * f_BE + q * ((f_BE + f_TW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloZ)) / (q + c1o1) ;
         //(D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=((c1o1 - q) * f_TW + q * ((f_TW + f_BE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloZ)) / (q + c1o1) ;
         //(D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=((c1o1 - q) * f_TN + q * ((f_TN + f_BS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=((c1o1 - q) * f_BS + q * ((f_BS + f_TN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTS])[kts]=((c1o1 - q) * f_BN + q * ((f_BN + f_TS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=((c1o1 - q) * f_TS + q * ((f_TS + f_BN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=((c1o1 - q) * f_TNE + q * ((f_TNE + f_BSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=((c1o1 - q) * f_BSW + q * ((f_BSW + f_TNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=((c1o1 - q) * f_BNE + q * ((f_BNE + f_TSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=((c1o1 - q) * f_TSW + q * ((f_TSW + f_BNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=((c1o1 - q) * f_TSE + q * ((f_TSE + f_BNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=((c1o1 - q) * f_BNW + q * ((f_BNW + f_TSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=((c1o1 - q) * f_BSE + q * ((f_BSE + f_TNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=((c1o1 - q) * f_TNW + q * ((f_TNW + f_BSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDeviceComp27(int inx,
										 int iny,
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

   if(k < numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   * numberOfBCnodes];
      q_dirW   = &QQ[dirW   * numberOfBCnodes];
      q_dirN   = &QQ[dirN   * numberOfBCnodes];
      q_dirS   = &QQ[dirS   * numberOfBCnodes];
      q_dirT   = &QQ[dirT   * numberOfBCnodes];
      q_dirB   = &QQ[dirB   * numberOfBCnodes];
      q_dirNE  = &QQ[dirNE  * numberOfBCnodes];
      q_dirSW  = &QQ[dirSW  * numberOfBCnodes];
      q_dirSE  = &QQ[dirSE  * numberOfBCnodes];
      q_dirNW  = &QQ[dirNW  * numberOfBCnodes];
      q_dirTE  = &QQ[dirTE  * numberOfBCnodes];
      q_dirBW  = &QQ[dirBW  * numberOfBCnodes];
      q_dirBE  = &QQ[dirBE  * numberOfBCnodes];
      q_dirTW  = &QQ[dirTW  * numberOfBCnodes];
      q_dirTN  = &QQ[dirTN  * numberOfBCnodes];
      q_dirBS  = &QQ[dirBS  * numberOfBCnodes];
      q_dirBN  = &QQ[dirBN  * numberOfBCnodes];
      q_dirTS  = &QQ[dirTS  * numberOfBCnodes];
      q_dirTNE = &QQ[dirTNE * numberOfBCnodes];
      q_dirTSW = &QQ[dirTSW * numberOfBCnodes];
      q_dirTSE = &QQ[dirTSE * numberOfBCnodes];
      q_dirTNW = &QQ[dirTNW * numberOfBCnodes];
      q_dirBNE = &QQ[dirBNE * numberOfBCnodes];
      q_dirBSW = &QQ[dirBSW * numberOfBCnodes];
      q_dirBSE = &QQ[dirBSE * numberOfBCnodes];
      q_dirBNW = &QQ[dirBNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int numberOfNodesK  = k_Q[k];
      unsigned int kzero= numberOfNodesK;
      unsigned int ke   = numberOfNodesK;
      unsigned int kw   = neighborX[numberOfNodesK];
      unsigned int kn   = numberOfNodesK;
      unsigned int ks   = neighborY[numberOfNodesK];
      unsigned int kt   = numberOfNodesK;
      unsigned int kb   = neighborZ[numberOfNodesK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = numberOfNodesK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = numberOfNodesK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = numberOfNodesK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = numberOfNodesK;
      unsigned int kbsw = neighborZ[ksw];
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


      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S)) / (c1o1 + drho); 

      vx3    =    (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                  (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                  (f_T - f_B)) / (c1o1 + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      //Test
         //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirW])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W))/(c1o1+q);
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirE])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E))/(c1o1+q);
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirS])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S))/(c1o1+q);
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirN])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N))/(c1o1+q);
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirB])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B))/(c1o1+q);
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirT])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T))/(c1o1+q);
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW))/(c1o1+q);
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirNE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE))/(c1o1+q);
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirNW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW))/(c1o1+q);
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dirSE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE))/(c1o1+q);
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW))/(c1o1+q);
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE))/(c1o1+q);
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW))/(c1o1+q);
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE))/(c1o1+q);
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS))/(c1o1+q);
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN))/(c1o1+q);
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS))/(c1o1+q);
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN))/(c1o1+q);
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW))/(c1o1+q);
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE))/(c1o1+q);
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW))/(c1o1+q);
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE))/(c1o1+q);
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW))/(c1o1+q);
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE))/(c1o1+q);
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW))/(c1o1+q);
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE))/(c1o1+q);
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDevice27(int inx,
                                     int iny,
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

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *numberOfBCnodes];
      q_dirW   = &QQ[dirW   *numberOfBCnodes];
      q_dirN   = &QQ[dirN   *numberOfBCnodes];
      q_dirS   = &QQ[dirS   *numberOfBCnodes];
      q_dirT   = &QQ[dirT   *numberOfBCnodes];
      q_dirB   = &QQ[dirB   *numberOfBCnodes];
      q_dirNE  = &QQ[dirNE  *numberOfBCnodes];
      q_dirSW  = &QQ[dirSW  *numberOfBCnodes];
      q_dirSE  = &QQ[dirSE  *numberOfBCnodes];
      q_dirNW  = &QQ[dirNW  *numberOfBCnodes];
      q_dirTE  = &QQ[dirTE  *numberOfBCnodes];
      q_dirBW  = &QQ[dirBW  *numberOfBCnodes];
      q_dirBE  = &QQ[dirBE  *numberOfBCnodes];
      q_dirTW  = &QQ[dirTW  *numberOfBCnodes];
      q_dirTN  = &QQ[dirTN  *numberOfBCnodes];
      q_dirBS  = &QQ[dirBS  *numberOfBCnodes];
      q_dirBN  = &QQ[dirBN  *numberOfBCnodes];
      q_dirTS  = &QQ[dirTS  *numberOfBCnodes];
      q_dirTNE = &QQ[dirTNE *numberOfBCnodes];
      q_dirTSW = &QQ[dirTSW *numberOfBCnodes];
      q_dirTSE = &QQ[dirTSE *numberOfBCnodes];
      q_dirTNW = &QQ[dirTNW *numberOfBCnodes];
      q_dirBNE = &QQ[dirBNE *numberOfBCnodes];
      q_dirBSW = &QQ[dirBSW *numberOfBCnodes];
      q_dirBSE = &QQ[dirBSE *numberOfBCnodes];
      q_dirBNW = &QQ[dirBNW *numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int numberOfNodesK  = k_Q[k];
      unsigned int kzero= numberOfNodesK;
      unsigned int ke   = numberOfNodesK;
      unsigned int kw   = neighborX[numberOfNodesK];
      unsigned int kn   = numberOfNodesK;
      unsigned int ks   = neighborY[numberOfNodesK];
      unsigned int kt   = numberOfNodesK;
      unsigned int kb   = neighborZ[numberOfNodesK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = numberOfNodesK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = numberOfNodesK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = numberOfNodesK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = numberOfNodesK;
      unsigned int kbsw = neighborZ[ksw];
      //unsigned int nxny = nx*ny;
      //unsigned int kzero= numberOfNodesK;
      //unsigned int ke   = numberOfNodesK;
      //unsigned int kw   = numberOfNodesK + 1;
      //unsigned int kn   = numberOfNodesK;
      //unsigned int ks   = numberOfNodesK + nx;
      //unsigned int kt   = numberOfNodesK;
      //unsigned int kb   = numberOfNodesK + nxny;
      //unsigned int ksw  = numberOfNodesK + nx + 1;
      //unsigned int kne  = numberOfNodesK;
      //unsigned int kse  = numberOfNodesK + nx;
      //unsigned int knw  = numberOfNodesK + 1;
      //unsigned int kbw  = numberOfNodesK + nxny + 1;
      //unsigned int kte  = numberOfNodesK;
      //unsigned int kbe  = numberOfNodesK + nxny;
      //unsigned int ktw  = numberOfNodesK + 1;
      //unsigned int kbs  = numberOfNodesK + nxny + nx;
      //unsigned int ktn  = numberOfNodesK;
      //unsigned int kbn  = numberOfNodesK + nxny;
      //unsigned int kts  = numberOfNodesK + nx;
      //unsigned int ktse = numberOfNodesK + nx;
      //unsigned int kbnw = numberOfNodesK + nxny + 1;
      //unsigned int ktnw = numberOfNodesK + 1;
      //unsigned int kbse = numberOfNodesK + nxny + nx;
      //unsigned int ktsw = numberOfNodesK + nx + 1;
      //unsigned int kbne = numberOfNodesK + nxny;
      //unsigned int ktne = numberOfNodesK;
      //unsigned int kbsw = numberOfNodesK + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
      //real vx1, vx2, vx3, drho, feq, q;
      //drho    =   (D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
      //            (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
      //            (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
      //            (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
      //            (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
      //            (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
      //            (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
      //            (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
      //            (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
      //            (D.f[dirZERO])[kzero]+ 
      //            (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
      //            (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
      //            (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
      //            (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw];

      //vx1    =    (D.f[dirE   ])[ke  ]- (D.f[dirW   ])[kw  ]+ 
      //            (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]+
      //            (D.f[dirSE  ])[kse ]- (D.f[dirNW  ])[knw ]+
      //            (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]+
      //            (D.f[dirBE  ])[kbe ]- (D.f[dirTW  ])[ktw ]+
      //            (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]+ 
      //            (D.f[dirTSE ])[ktse]- (D.f[dirTNW ])[ktnw]+ 
      //            (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]+ 
      //            (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw];

      //vx2    =    (D.f[dirN   ])[kn  ]- (D.f[dirS   ])[ks  ]+
      //            (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]-
      //            (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
      //            (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]+
      //            (D.f[dirBN  ])[kbn ]- (D.f[dirTS  ])[kts ]+
      //            (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]- 
      //            (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
      //            (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
      //            (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw];

      //vx3    =    (D.f[dirT   ])[kt  ]- (D.f[dirB   ])[kb  ]+
      //            (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]-
      //            (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
      //            (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]-
      //            (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
      //            (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
      //            (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]- 
      //            (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
      //            (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw];

      //real cu_sq=1.5f*(vx1*vx1+vx2*vx2+vx3*vx3);
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

      vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W); 


      vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE))  + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S); 

      vx3    =    ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                  (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                  (f_T - f_B); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

	  //b�ser lecktest
	  //q = q_dirE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirW])[kw]=999.f;
   //   }

   //   q = q_dirW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirE])[ke]=999.f;
   //   }

   //   q = q_dirN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirS])[ks]=999.f;
   //   }

   //   q = q_dirS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirN])[kn]=999.f;
   //   }

   //   q = q_dirT[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirB])[kb]=999.f;
   //   }

   //   q = q_dirB[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirT])[kt]=999.f;
   //   }

   //   q = q_dirNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirSW])[ksw]=999.f;
   //   }

   //   q = q_dirSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirNE])[kne]=999.f;
   //   }

   //   q = q_dirSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirNW])[knw]=999.f;
   //   }

   //   q = q_dirNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirSE])[kse]=999.f;
   //   }

   //   q = q_dirTE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBW])[kbw]=999.f;
   //   }

   //   q = q_dirBW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTE])[kte]=999.f;
   //   }

   //   q = q_dirBE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTW])[ktw]=999.f;
   //   }

   //   q = q_dirTW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBE])[kbe]=999.f;
   //   }

   //   q = q_dirTN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBS])[kbs]=999.f;
   //   }

   //   q = q_dirBS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTN])[ktn]=999.f;
   //   }

   //   q = q_dirBN[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTS])[kts]=999.f;
   //   }

   //   q = q_dirTS[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBN])[kbn]=999.f;
   //   }

   //   q = q_dirTNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBSW])[kbsw]=999.f;
   //   }

   //   q = q_dirBSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTNE])[ktne]=999.f;
   //   }

   //   q = q_dirBNE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTSW])[ktsw]=999.f;
   //   }

   //   q = q_dirTSW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBNE])[kbne]=999.f;
   //   }

   //   q = q_dirTSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBNW])[kbnw]=999.f;
   //   }

   //   q = q_dirBNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTSE])[ktse]=999.f;
   //   }

   //   q = q_dirBSE[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirTNW])[ktnw]=999.f;
   //   }

   //   q = q_dirTNW[k];
   //   if (q>=zero && q<=one)
   //   {
   //      (D.f[dirBSE])[kbse]=999.f;
   //   }

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      //Test
         //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        )-cu_sq); 
         (D.f[dirW])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W))/(c1o1+q);
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        )-cu_sq); 
         (D.f[dirE])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E))/(c1o1+q);
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    )-cu_sq); 
         (D.f[dirS])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S))/(c1o1+q);
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    )-cu_sq); 
         (D.f[dirN])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N))/(c1o1+q);
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3)-cu_sq); 
         (D.f[dirB])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B))/(c1o1+q);
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3)-cu_sq); 
         (D.f[dirT])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T))/(c1o1+q);
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         (D.f[dirSW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW))/(c1o1+q);
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         (D.f[dirNE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE))/(c1o1+q);
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         (D.f[dirNW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW))/(c1o1+q);
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         (D.f[dirSE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE))/(c1o1+q);
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         (D.f[dirBW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW))/(c1o1+q);
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         (D.f[dirTE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE))/(c1o1+q);
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         (D.f[dirTW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW))/(c1o1+q);
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         (D.f[dirBE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE))/(c1o1+q);
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         (D.f[dirBS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS))/(c1o1+q);
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         (D.f[dirTN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN))/(c1o1+q);
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         (D.f[dirTS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS))/(c1o1+q);
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         (D.f[dirBN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN))/(c1o1+q);
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW))/(c1o1+q);
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         (D.f[dirTNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE))/(c1o1+q);
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW))/(c1o1+q);
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         (D.f[dirBNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE))/(c1o1+q);
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW))/(c1o1+q);
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         (D.f[dirTSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE))/(c1o1+q);
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW))/(c1o1+q);
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         (D.f[dirBSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE))/(c1o1+q);
         //(D.f[dirBSE])[kbse]=zero;
      }

	 // q = q_dirE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq); 
  //       (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-feq*om1)/(one-om1)+q/(one+q)*(f_E+f_W);
		//// (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);
  //    }

  //    q = q_dirW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
  //       (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-feq*om1)/(one-om1)+q/(one+q)*(f_W+f_E);
  //    }

  //    q = q_dirN[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
  //       (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-feq*om1)/(one-om1)+q/(one+q)*(f_N+f_S);
  //    }

  //    q = q_dirS[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
  //       (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-feq*om1)/(one-om1)+q/(one+q)*(f_S+f_N);
  //    }

  //    q = q_dirT[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq); 
  //       (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-feq*om1)/(one-om1)+q/(one+q)*(f_T+f_B);
  //    }

  //    q = q_dirB[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
  //       (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-feq*om1)/(one-om1)+q/(one+q)*(f_B+f_T);
  //    }

  //    q = q_dirNE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
  //       (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-feq*om1)/(one-om1)+q/(one+q)*(f_NE+f_SW);
  //    }

  //    q = q_dirSW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
  //       (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-feq*om1)/(one-om1)+q/(one+q)*(f_SW+f_NE);
  //    }

  //    q = q_dirSE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
  //       (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-feq*om1)/(one-om1)+q/(one+q)*(f_SE+f_NW);
  //    }

  //    q = q_dirNW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
  //       (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-feq*om1)/(one-om1)+q/(one+q)*(f_NW+f_SE);
  //    }

  //    q = q_dirTE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
  //       (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-feq*om1)/(one-om1)+q/(one+q)*(f_TE+f_BW);
  //    }

  //    q = q_dirBW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
  //       (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-feq*om1)/(one-om1)+q/(one+q)*(f_BW+f_TE);
  //    }

  //    q = q_dirBE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
  //       (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-feq*om1)/(one-om1)+q/(one+q)*(f_BE+f_TW);
  //    }

  //    q = q_dirTW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
  //       (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-feq*om1)/(one-om1)+q/(one+q)*(f_TW+f_BE);
  //    }

  //    q = q_dirTN[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
  //       (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-feq*om1)/(one-om1)+q/(one+q)*(f_TN+f_BS);
  //    }

  //    q = q_dirBS[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
  //       (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-feq*om1)/(one-om1)+q/(one+q)*(f_BS+f_TN);
  //    }

  //    q = q_dirBN[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
  //       (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-feq*om1)/(one-om1)+q/(one+q)*(f_BN+f_TS);
  //    }

  //    q = q_dirTS[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
  //       (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-feq*om1)/(one-om1)+q/(one+q)*(f_TS+f_BN);
  //    }

  //    q = q_dirTNE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
  //       (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-feq*om1)/(one-om1)+q/(one+q)*(f_TNE+f_BSW);
  //    }

  //    q = q_dirBSW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
  //       (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-feq*om1)/(one-om1)+q/(one+q)*(f_BSW+f_TNE);
  //    }

  //    q = q_dirBNE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
  //       (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-feq*om1)/(one-om1)+q/(one+q)*(f_BNE+f_TSW);
  //    }

  //    q = q_dirTSW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
  //       (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-feq*om1)/(one-om1)+q/(one+q)*(f_TSW+f_BNE);
  //    }

  //    q = q_dirTSE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
  //       (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-feq*om1)/(one-om1)+q/(one+q)*(f_TSE+f_BNW);
  //    }

  //    q = q_dirBNW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
  //       (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-feq*om1)/(one-om1)+q/(one+q)*(f_BNW+f_TSE);
  //    }

  //    q = q_dirBSE[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
  //       (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-feq*om1)/(one-om1)+q/(one+q)*(f_BSE+f_TNW);
  //    }

  //    q = q_dirTNW[k];
  //    if (q>=zero && q<=one)
  //    {
  //       feq=c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
  //       (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-feq*om1)/(one-om1)+q/(one+q)*(f_TNW+f_BSE);
  //    }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void BBDevice27(int inx,
                                     int iny,
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

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *numberOfBCnodes];
      q_dirW   = &QQ[dirW   *numberOfBCnodes];
      q_dirN   = &QQ[dirN   *numberOfBCnodes];
      q_dirS   = &QQ[dirS   *numberOfBCnodes];
      q_dirT   = &QQ[dirT   *numberOfBCnodes];
      q_dirB   = &QQ[dirB   *numberOfBCnodes];
      q_dirNE  = &QQ[dirNE  *numberOfBCnodes];
      q_dirSW  = &QQ[dirSW  *numberOfBCnodes];
      q_dirSE  = &QQ[dirSE  *numberOfBCnodes];
      q_dirNW  = &QQ[dirNW  *numberOfBCnodes];
      q_dirTE  = &QQ[dirTE  *numberOfBCnodes];
      q_dirBW  = &QQ[dirBW  *numberOfBCnodes];
      q_dirBE  = &QQ[dirBE  *numberOfBCnodes];
      q_dirTW  = &QQ[dirTW  *numberOfBCnodes];
      q_dirTN  = &QQ[dirTN  *numberOfBCnodes];
      q_dirBS  = &QQ[dirBS  *numberOfBCnodes];
      q_dirBN  = &QQ[dirBN  *numberOfBCnodes];
      q_dirTS  = &QQ[dirTS  *numberOfBCnodes];
      q_dirTNE = &QQ[dirTNE *numberOfBCnodes];
      q_dirTSW = &QQ[dirTSW *numberOfBCnodes];
      q_dirTSE = &QQ[dirTSE *numberOfBCnodes];
      q_dirTNW = &QQ[dirTNW *numberOfBCnodes];
      q_dirBNE = &QQ[dirBNE *numberOfBCnodes];
      q_dirBSW = &QQ[dirBSW *numberOfBCnodes];
      q_dirBSE = &QQ[dirBSE *numberOfBCnodes];
      q_dirBNW = &QQ[dirBNW *numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int numberOfNodesK  = k_Q[k];
      //unsigned int kzero= numberOfNodesK;
      unsigned int ke   = numberOfNodesK;
      unsigned int kw   = neighborX[numberOfNodesK];
      unsigned int kn   = numberOfNodesK;
      unsigned int ks   = neighborY[numberOfNodesK];
      unsigned int kt   = numberOfNodesK;
      unsigned int kb   = neighborZ[numberOfNodesK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = numberOfNodesK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = numberOfNodesK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = numberOfNodesK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = numberOfNodesK;
      unsigned int kbsw = neighborZ[ksw];
      //unsigned int nxny = nx*ny;
      //unsigned int kzero= numberOfNodesK;
      //unsigned int ke   = numberOfNodesK;
      //unsigned int kw   = numberOfNodesK + 1;
      //unsigned int kn   = numberOfNodesK;
      //unsigned int ks   = numberOfNodesK + nx;
      //unsigned int kt   = numberOfNodesK;
      //unsigned int kb   = numberOfNodesK + nxny;
      //unsigned int ksw  = numberOfNodesK + nx + 1;
      //unsigned int kne  = numberOfNodesK;
      //unsigned int kse  = numberOfNodesK + nx;
      //unsigned int knw  = numberOfNodesK + 1;
      //unsigned int kbw  = numberOfNodesK + nxny + 1;
      //unsigned int kte  = numberOfNodesK;
      //unsigned int kbe  = numberOfNodesK + nxny;
      //unsigned int ktw  = numberOfNodesK + 1;
      //unsigned int kbs  = numberOfNodesK + nxny + nx;
      //unsigned int ktn  = numberOfNodesK;
      //unsigned int kbn  = numberOfNodesK + nxny;
      //unsigned int kts  = numberOfNodesK + nx;
      //unsigned int ktse = numberOfNodesK + nx;
      //unsigned int kbnw = numberOfNodesK + nxny + 1;
      //unsigned int ktnw = numberOfNodesK + 1;
      //unsigned int kbse = numberOfNodesK + nxny + nx;
      //unsigned int ktsw = numberOfNodesK + nx + 1;
      //unsigned int kbne = numberOfNodesK + nxny;
      //unsigned int ktne = numberOfNodesK;
      //unsigned int kbsw = numberOfNodesK + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
     
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

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real q;
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirW])[kw]=f_E;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirE])[ke]=f_W;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirS])[ks]=f_N;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirN])[kn]=f_S;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirB])[kb]=f_T;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirT])[kt]=f_B;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirSW])[ksw]=f_NE;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirNE])[kne]=f_SW;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirNW])[knw]=f_SE;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirSE])[kse]=f_NW;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirBW])[kbw]=f_TE;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirTE])[kte]=f_BW;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirTW])[ktw]=f_BE;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirBE])[kbe]=f_TW;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirBS])[kbs]=f_TN;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirTN])[ktn]=f_BS;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirTS])[kts]=f_BN;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirBN])[kbn]=f_TS;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirBSW])[kbsw]=f_TNE;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirTNE])[ktne]=f_BSW;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirTSW])[ktsw]=f_BNE;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirBNE])[kbne]=f_TSW;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirBNW])[kbnw]=f_TSE;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirTSE])[ktse]=f_BNW;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirTNW])[ktnw]=f_BSE;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dirBSE])[kbse]=f_TNW;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

