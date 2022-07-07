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
#include "KernelUtilities.h"

using namespace vf::lbm::constant;

//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDevice3rdMomentsComp27(
													 real* distributions, 
													 int* subgridDistanceIndices, 
													 real* subgridDistances,
													 unsigned int numberOfBCnodes, 
													 real omega, 
													 unsigned int* neighborX,
													 unsigned int* neighborY,
													 unsigned int* neighborZ,
													 unsigned int numberOfLBnodes, 
													 bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[E   ] = &distributions[E   *numberOfLBnodes];
      D.f[W   ] = &distributions[W   *numberOfLBnodes];
      D.f[N   ] = &distributions[N   *numberOfLBnodes];
      D.f[S   ] = &distributions[S   *numberOfLBnodes];
      D.f[T   ] = &distributions[T   *numberOfLBnodes];
      D.f[B   ] = &distributions[B   *numberOfLBnodes];
      D.f[NE  ] = &distributions[NE  *numberOfLBnodes];
      D.f[SW  ] = &distributions[SW  *numberOfLBnodes];
      D.f[SE  ] = &distributions[SE  *numberOfLBnodes];
      D.f[NW  ] = &distributions[NW  *numberOfLBnodes];
      D.f[TE  ] = &distributions[TE  *numberOfLBnodes];
      D.f[BW  ] = &distributions[BW  *numberOfLBnodes];
      D.f[BE  ] = &distributions[BE  *numberOfLBnodes];
      D.f[TW  ] = &distributions[TW  *numberOfLBnodes];
      D.f[TN  ] = &distributions[TN  *numberOfLBnodes];
      D.f[BS  ] = &distributions[BS  *numberOfLBnodes];
      D.f[BN  ] = &distributions[BN  *numberOfLBnodes];
      D.f[TS  ] = &distributions[TS  *numberOfLBnodes];
      D.f[dirREST] = &distributions[dirREST*numberOfLBnodes];
      D.f[TNE ] = &distributions[TNE *numberOfLBnodes];
      D.f[TSW ] = &distributions[TSW *numberOfLBnodes];
      D.f[TSE ] = &distributions[TSE *numberOfLBnodes];
      D.f[TNW ] = &distributions[TNW *numberOfLBnodes];
      D.f[BNE ] = &distributions[BNE *numberOfLBnodes];
      D.f[BSW ] = &distributions[BSW *numberOfLBnodes];
      D.f[BSE ] = &distributions[BSE *numberOfLBnodes];
      D.f[BNW ] = &distributions[BNW *numberOfLBnodes];
   } 
   else
   {
      D.f[W   ] = &distributions[E   *numberOfLBnodes];
      D.f[E   ] = &distributions[W   *numberOfLBnodes];
      D.f[S   ] = &distributions[N   *numberOfLBnodes];
      D.f[N   ] = &distributions[S   *numberOfLBnodes];
      D.f[B   ] = &distributions[T   *numberOfLBnodes];
      D.f[T   ] = &distributions[B   *numberOfLBnodes];
      D.f[SW  ] = &distributions[NE  *numberOfLBnodes];
      D.f[NE  ] = &distributions[SW  *numberOfLBnodes];
      D.f[NW  ] = &distributions[SE  *numberOfLBnodes];
      D.f[SE  ] = &distributions[NW  *numberOfLBnodes];
      D.f[BW  ] = &distributions[TE  *numberOfLBnodes];
      D.f[TE  ] = &distributions[BW  *numberOfLBnodes];
      D.f[TW  ] = &distributions[BE  *numberOfLBnodes];
      D.f[BE  ] = &distributions[TW  *numberOfLBnodes];
      D.f[BS  ] = &distributions[TN  *numberOfLBnodes];
      D.f[TN  ] = &distributions[BS  *numberOfLBnodes];
      D.f[TS  ] = &distributions[BN  *numberOfLBnodes];
      D.f[BN  ] = &distributions[TS  *numberOfLBnodes];
      D.f[dirREST] = &distributions[dirREST*numberOfLBnodes];
      D.f[TNE ] = &distributions[BSW *numberOfLBnodes];
      D.f[TSW ] = &distributions[BNE *numberOfLBnodes];
      D.f[TSE ] = &distributions[BNW *numberOfLBnodes];
      D.f[TNW ] = &distributions[BSE *numberOfLBnodes];
      D.f[BNE ] = &distributions[TSW *numberOfLBnodes];
      D.f[BSW ] = &distributions[TNE *numberOfLBnodes];
      D.f[BSE ] = &distributions[TNW *numberOfLBnodes];
      D.f[BNW ] = &distributions[TSE *numberOfLBnodes];
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
      q_dirE   = &subgridDistances[E   * numberOfBCnodes];
      q_dirW   = &subgridDistances[W   * numberOfBCnodes];
      q_dirN   = &subgridDistances[N   * numberOfBCnodes];
      q_dirS   = &subgridDistances[S   * numberOfBCnodes];
      q_dirT   = &subgridDistances[T   * numberOfBCnodes];
      q_dirB   = &subgridDistances[B   * numberOfBCnodes];
      q_dirNE  = &subgridDistances[NE  * numberOfBCnodes];
      q_dirSW  = &subgridDistances[SW  * numberOfBCnodes];
      q_dirSE  = &subgridDistances[SE  * numberOfBCnodes];
      q_dirNW  = &subgridDistances[NW  * numberOfBCnodes];
      q_dirTE  = &subgridDistances[TE  * numberOfBCnodes];
      q_dirBW  = &subgridDistances[BW  * numberOfBCnodes];
      q_dirBE  = &subgridDistances[BE  * numberOfBCnodes];
      q_dirTW  = &subgridDistances[TW  * numberOfBCnodes];
      q_dirTN  = &subgridDistances[TN  * numberOfBCnodes];
      q_dirBS  = &subgridDistances[BS  * numberOfBCnodes];
      q_dirBN  = &subgridDistances[BN  * numberOfBCnodes];
      q_dirTS  = &subgridDistances[TS  * numberOfBCnodes];
      q_dirTNE = &subgridDistances[TNE * numberOfBCnodes];
      q_dirTSW = &subgridDistances[TSW * numberOfBCnodes];
      q_dirTSE = &subgridDistances[TSE * numberOfBCnodes];
      q_dirTNW = &subgridDistances[TNW * numberOfBCnodes];
      q_dirBNE = &subgridDistances[BNE * numberOfBCnodes];
      q_dirBSW = &subgridDistances[BSW * numberOfBCnodes];
      q_dirBSE = &subgridDistances[BSE * numberOfBCnodes];
      q_dirBNW = &subgridDistances[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int numberOfNodesK  = subgridDistanceIndices[k];
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
      real vx1, vx2, vx3, drho, feq, q, m3;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]); 

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
         D.f[E   ] = &distributions[E   *numberOfLBnodes];
         D.f[W   ] = &distributions[W   *numberOfLBnodes];
         D.f[N   ] = &distributions[N   *numberOfLBnodes];
         D.f[S   ] = &distributions[S   *numberOfLBnodes];
         D.f[T   ] = &distributions[T   *numberOfLBnodes];
         D.f[B   ] = &distributions[B   *numberOfLBnodes];
         D.f[NE  ] = &distributions[NE  *numberOfLBnodes];
         D.f[SW  ] = &distributions[SW  *numberOfLBnodes];
         D.f[SE  ] = &distributions[SE  *numberOfLBnodes];
         D.f[NW  ] = &distributions[NW  *numberOfLBnodes];
         D.f[TE  ] = &distributions[TE  *numberOfLBnodes];
         D.f[BW  ] = &distributions[BW  *numberOfLBnodes];
         D.f[BE  ] = &distributions[BE  *numberOfLBnodes];
         D.f[TW  ] = &distributions[TW  *numberOfLBnodes];
         D.f[TN  ] = &distributions[TN  *numberOfLBnodes];
         D.f[BS  ] = &distributions[BS  *numberOfLBnodes];
         D.f[BN  ] = &distributions[BN  *numberOfLBnodes];
         D.f[TS  ] = &distributions[TS  *numberOfLBnodes];
         D.f[dirREST] = &distributions[dirREST*numberOfLBnodes];
         D.f[TNE ] = &distributions[TNE *numberOfLBnodes];
         D.f[TSW ] = &distributions[TSW *numberOfLBnodes];
         D.f[TSE ] = &distributions[TSE *numberOfLBnodes];
         D.f[TNW ] = &distributions[TNW *numberOfLBnodes];
         D.f[BNE ] = &distributions[BNE *numberOfLBnodes];
         D.f[BSW ] = &distributions[BSW *numberOfLBnodes];
         D.f[BSE ] = &distributions[BSE *numberOfLBnodes];
         D.f[BNW ] = &distributions[BNW *numberOfLBnodes];
      } 
      else
      {
         D.f[W   ] = &distributions[E   *numberOfLBnodes];
         D.f[E   ] = &distributions[W   *numberOfLBnodes];
         D.f[S   ] = &distributions[N   *numberOfLBnodes];
         D.f[N   ] = &distributions[S   *numberOfLBnodes];
         D.f[B   ] = &distributions[T   *numberOfLBnodes];
         D.f[T   ] = &distributions[B   *numberOfLBnodes];
         D.f[SW  ] = &distributions[NE  *numberOfLBnodes];
         D.f[NE  ] = &distributions[SW  *numberOfLBnodes];
         D.f[NW  ] = &distributions[SE  *numberOfLBnodes];
         D.f[SE  ] = &distributions[NW  *numberOfLBnodes];
         D.f[BW  ] = &distributions[TE  *numberOfLBnodes];
         D.f[TE  ] = &distributions[BW  *numberOfLBnodes];
         D.f[TW  ] = &distributions[BE  *numberOfLBnodes];
         D.f[BE  ] = &distributions[TW  *numberOfLBnodes];
         D.f[BS  ] = &distributions[TN  *numberOfLBnodes];
         D.f[TN  ] = &distributions[BS  *numberOfLBnodes];
         D.f[TS  ] = &distributions[BN  *numberOfLBnodes];
         D.f[BN  ] = &distributions[TS  *numberOfLBnodes];
         D.f[dirREST] = &distributions[dirREST*numberOfLBnodes];
         D.f[TNE ] = &distributions[BSW *numberOfLBnodes];
         D.f[TSW ] = &distributions[BNE *numberOfLBnodes];
         D.f[TSE ] = &distributions[BNW *numberOfLBnodes];
         D.f[TNW ] = &distributions[BSE *numberOfLBnodes];
         D.f[BNE ] = &distributions[TSW *numberOfLBnodes];
         D.f[BSW ] = &distributions[TNE *numberOfLBnodes];
         D.f[BSE ] = &distributions[TNW *numberOfLBnodes];
         D.f[BNW ] = &distributions[TSE *numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
         //(D.f[dirREST])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_E - f_W - c2o1 * drho * c2o27 * (c3o1*( vx1        ));
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[W])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W-m3+(f_E+f_W-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_E+f_W))/(c1o1+q)+(m3*c1o2);
         //(D.f[W])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_W - f_E - c2o1 * drho * c2o27 * (c3o1*(-vx1        ));
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[E])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E-m3+(f_W+f_E-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_W+f_E))/(c1o1+q)+(m3*c1o2);
         //(D.f[E])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_N - f_S - c2o1 * drho * c2o27 * (c3o1*( vx2        ));
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[S])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S-m3+(f_N+f_S-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_N+f_S))/(c1o1+q)+(m3*c1o2);
         //(D.f[S])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_S - f_N - c2o1 * drho * c2o27 * (c3o1*(   -vx2     ));
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[N])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N-m3+(f_S+f_N-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_S+f_N))/(c1o1+q)+(m3*c1o2);
         //(D.f[N])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_T - f_B - c2o1 * drho * c2o27 * (c3o1*(         vx3));
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[B])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B-m3+(f_T+f_B-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_T+f_B))/(c1o1+q)+(m3*c1o2);
         //(D.f[B])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_B - f_T - c2o1 * drho * c2o27 * (c3o1*(        -vx3));
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[T])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T-m3+(f_B+f_T-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_B+f_T))/(c1o1+q)+(m3*c1o2);
         //(D.f[T])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_NE - f_SW - c2o1 * drho * c1o54 * (c3o1*( vx1+vx2    ));
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW-m3+(f_NE+f_SW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_NE+f_SW))/(c1o1+q)+(m3*c1o2);
         //(D.f[SW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_SW - f_NE - c2o1 * drho * c1o54 * (c3o1*(-vx1-vx2    ));
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE-m3+(f_SW+f_NE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_SW+f_NE))/(c1o1+q)+(m3*c1o2);
         //(D.f[NE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_SE - f_NW - c2o1 * drho * c1o54 * (c3o1*( vx1-vx2    ));
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW-m3+(f_SE+f_NW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_SE+f_NW))/(c1o1+q)+(m3*c1o2);
         //(D.f[NW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_NW - f_SE - c2o1 * drho * c1o54 * (c3o1*(-vx1+vx2    ));
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE-m3+(f_NW+f_SE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_NW+f_SE))/(c1o1+q)+(m3*c1o2);
         //(D.f[SE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TE - f_BW - c2o1 * drho * c1o54 * (c3o1*( vx1    +vx3));
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW-m3+(f_TE+f_BW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TE+f_BW))/(c1o1+q)+(m3*c1o2);
         //(D.f[BW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BW - f_TE - c2o1 * drho * c1o54 * (c3o1*(-vx1    -vx3));
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE-m3+(f_BW+f_TE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BW+f_TE))/(c1o1+q)+(m3*c1o2);
         //(D.f[TE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BE - f_TW - c2o1 * drho * c1o54 * (c3o1*( vx1    -vx3));
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW-m3+(f_BE+f_TW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BE+f_TW))/(c1o1+q)+(m3*c1o2);
         //(D.f[TW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TW - f_BE - c2o1 * drho * c1o54 * (c3o1*(-vx1    +vx3));
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE-m3+(f_TW+f_BE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TW+f_BE))/(c1o1+q)+(m3*c1o2);
         //(D.f[BE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TN - f_BS - c2o1 * drho * c1o54 * (c3o1*(     vx2+vx3));
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS-m3+(f_TN+f_BS-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TN+f_BS))/(c1o1+q)+(m3*c1o2);
         //(D.f[BS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BS - f_TN - c2o1 * drho * c1o54 * (c3o1*(    -vx2-vx3));
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN-m3+(f_BS+f_TN-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BS+f_TN))/(c1o1+q)+(m3*c1o2);
         //(D.f[TN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BN - f_TS - c2o1 * drho * c1o54 * (c3o1*(     vx2-vx3));
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS-m3+(f_BN+f_TS-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BN+f_TS))/(c1o1+q)+(m3*c1o2);
         //(D.f[TS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TS - f_BN - c2o1 * drho * c1o54 * (c3o1*(    -vx2+vx3));
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN-m3+(f_TS+f_BN-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TS+f_BN))/(c1o1+q)+(m3*c1o2);
         //(D.f[BN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TNE - f_BSW - c2o1 * drho * c1o216 * (c3o1*( vx1+vx2+vx3));
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW-m3+(f_TNE+f_BSW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TNE+f_BSW))/(c1o1+q)+(m3*c1o2);
         //(D.f[BSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BSW - f_TNE - c2o1 * drho * c1o216 * (c3o1*(-vx1-vx2-vx3));
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE-m3+(f_BSW+f_TNE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BSW+f_TNE))/(c1o1+q)+(m3*c1o2);
         //(D.f[TNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BNE - f_TSW - c2o1 * drho * c1o216 * (c3o1*( vx1+vx2-vx3));
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW-m3+(f_BNE+f_TSW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BNE+f_TSW))/(c1o1+q)+(m3*c1o2);
         //(D.f[TSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TSW - f_BNE - c2o1 * drho * c1o216 * (c3o1*(-vx1-vx2+vx3));
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE-m3+(f_TSW+f_BNE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TSW+f_BNE))/(c1o1+q)+(m3*c1o2);
         //(D.f[BNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TSE - f_BNW - c2o1 * drho * c1o216 * (c3o1*( vx1-vx2+vx3));
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW-m3+(f_TSE+f_BNW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TSE+f_BNW))/(c1o1+q)+(m3*c1o2);
         //(D.f[BNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BNW - f_TSE - c2o1 * drho * c1o216 * (c3o1*(-vx1+vx2-vx3));
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE-m3+(f_BNW+f_TSE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BNW+f_TSE))/(c1o1+q)+(m3*c1o2);
         //(D.f[TSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BSE - f_TNW - c2o1 * drho * c1o216 * (c3o1*( vx1-vx2-vx3));
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW-m3+(f_BSE+f_TNW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BSE+f_TNW))/(c1o1+q)+(m3*c1o2);
         //(D.f[TNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TNW - f_BSE - c2o1 * drho * c1o216 * (c3o1*(-vx1+vx2+vx3));
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE-m3+(f_TNW+f_BSE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TNW+f_BSE))/(c1o1+q)+(m3*c1o2);
         //(D.f[BSE])[kbse]=zero;
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
												 unsigned int  numberOfBCnodes,
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

   if(k<numberOfNodes)
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

      f_E   = (D.f[E   ])[ke   ];
      f_W   = (D.f[W   ])[kw   ];
      f_N   = (D.f[N   ])[kn   ];
      f_S   = (D.f[S   ])[ks   ];
      f_T   = (D.f[T   ])[kt   ];
      f_B   = (D.f[B   ])[kb   ];
      f_NE  = (D.f[NE  ])[kne  ];
      f_SW  = (D.f[SW  ])[ksw  ];
      f_SE  = (D.f[SE  ])[kse  ];
      f_NW  = (D.f[NW  ])[knw  ];
      f_TE  = (D.f[TE  ])[kte  ];
      f_BW  = (D.f[BW  ])[kbw  ];
      f_BE  = (D.f[BE  ])[kbe  ];
      f_TW  = (D.f[TW  ])[ktw  ];
      f_TN  = (D.f[TN  ])[ktn  ];
      f_BS  = (D.f[BS  ])[kbs  ];
      f_BN  = (D.f[BN  ])[kbn  ];
      f_TS  = (D.f[TS  ])[kts  ];
      f_TNE = (D.f[TNE ])[ktne ];
      f_TSW = (D.f[TSW ])[ktsw ];
      f_TSE = (D.f[TSE ])[ktse ];
      f_TNW = (D.f[TNW ])[ktnw ];
      f_BNE = (D.f[BNE ])[kbne ];
      f_BSW = (D.f[BSW ])[kbsw ];
      f_BSE = (D.f[BSE ])[kbse ];
      f_BNW = (D.f[BNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]); 

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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
         //(D.f[dirREST])[k]=c1o10;
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
         (D.f[W])[kw]=((c1o1 - q) * f_E + q * ((f_E + f_W) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloX     )) / (q + c1o1) ;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) /** (one + drho)*/-cu_sq); 
         (D.f[E])[ke]=((c1o1 - q) * f_W + q * ((f_W + f_E) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloX     )) / (q + c1o1) ;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[S])[ks]=((c1o1 - q) * f_N + q * ((f_N + f_S) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloY     )) / (q + c1o1) ;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[N])[kn]=((c1o1 - q) * f_S + q * ((f_S + f_N) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloY     )) / (q + c1o1) ;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) /** (one + drho)*/-cu_sq); 
         (D.f[B])[kb]=((c1o1 - q) * f_T + q * ((f_T + f_B) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloZ     )) / (q + c1o1) ;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[T])[kt]=((c1o1 - q) * f_B + q * ((f_B + f_T) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloZ     )) / (q + c1o1) ;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[SW])[ksw]=((c1o1 - q) * f_NE + q * ((f_NE + f_SW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloY)) / (q + c1o1) ;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[NE])[kne]=((c1o1 - q) * f_SW + q * ((f_SW + f_NE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloY)) / (q + c1o1) ;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[NW])[knw]=((c1o1 - q) * f_SE + q * ((f_SE + f_NW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloY)) / (q + c1o1) ;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[SE])[kse]=((c1o1 - q) * f_NW + q * ((f_NW + f_SE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloY)) / (q + c1o1) ;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BW])[kbw]=((c1o1 - q) * f_TE + q * ((f_TE + f_BW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TE])[kte]=((c1o1 - q) * f_BW + q * ((f_BW + f_TE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TW])[ktw]=((c1o1 - q) * f_BE + q * ((f_BE + f_TW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BE])[kbe]=((c1o1 - q) * f_TW + q * ((f_TW + f_BE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BS])[kbs]=((c1o1 - q) * f_TN + q * ((f_TN + f_BS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TN])[ktn]=((c1o1 - q) * f_BS + q * ((f_BS + f_TN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TS])[kts]=((c1o1 - q) * f_BN + q * ((f_BN + f_TS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BN])[kbn]=((c1o1 - q) * f_TS + q * ((f_TS + f_BN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BSW])[kbsw]=((c1o1 - q) * f_TNE + q * ((f_TNE + f_BSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TNE])[ktne]=((c1o1 - q) * f_BSW + q * ((f_BSW + f_TNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TSW])[ktsw]=((c1o1 - q) * f_BNE + q * ((f_BNE + f_TSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BNE])[kbne]=((c1o1 - q) * f_TSW + q * ((f_TSW + f_BNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BNW])[kbnw]=((c1o1 - q) * f_TSE + q * ((f_TSE + f_BNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TSE])[ktse]=((c1o1 - q) * f_BNW + q * ((f_BNW + f_TSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TNW])[ktnw]=((c1o1 - q) * f_BSE + q * ((f_BSE + f_TNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BSE])[kbse]=((c1o1 - q) * f_TNW + q * ((f_TNW + f_BSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY+VeloZ)) / (q + c1o1) ;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDeviceCompHighNu27(
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

      f_E   = (D.f[E   ])[ke   ];
      f_W   = (D.f[W   ])[kw   ];
      f_N   = (D.f[N   ])[kn   ];
      f_S   = (D.f[S   ])[ks   ];
      f_T   = (D.f[T   ])[kt   ];
      f_B   = (D.f[B   ])[kb   ];
      f_NE  = (D.f[NE  ])[kne  ];
      f_SW  = (D.f[SW  ])[ksw  ];
      f_SE  = (D.f[SE  ])[kse  ];
      f_NW  = (D.f[NW  ])[knw  ];
      f_TE  = (D.f[TE  ])[kte  ];
      f_BW  = (D.f[BW  ])[kbw  ];
      f_BE  = (D.f[BE  ])[kbe  ];
      f_TW  = (D.f[TW  ])[ktw  ];
      f_TN  = (D.f[TN  ])[ktn  ];
      f_BS  = (D.f[BS  ])[kbs  ];
      f_BN  = (D.f[BN  ])[kbn  ];
      f_TS  = (D.f[TS  ])[kts  ];
      f_TNE = (D.f[TNE ])[ktne ];
      f_TSW = (D.f[TSW ])[ktsw ];
      f_TSE = (D.f[TSE ])[ktse ];
      f_TNW = (D.f[TNW ])[ktnw ];
      f_BNE = (D.f[BNE ])[kbne ];
      f_BSW = (D.f[BSW ])[kbsw ];
      f_BSE = (D.f[BSE ])[kbse ];
      f_BNW = (D.f[BNW ])[kbnw ];
      //f_W    = (D.f[E   ])[ke   ];
      //f_E    = (D.f[W   ])[kw   ];
      //f_S    = (D.f[N   ])[kn   ];
      //f_N    = (D.f[S   ])[ks   ];
      //f_B    = (D.f[T   ])[kt   ];
      //f_T    = (D.f[B   ])[kb   ];
      //f_SW   = (D.f[NE  ])[kne  ];
      //f_NE   = (D.f[SW  ])[ksw  ];
      //f_NW   = (D.f[SE  ])[kse  ];
      //f_SE   = (D.f[NW  ])[knw  ];
      //f_BW   = (D.f[TE  ])[kte  ];
      //f_TE   = (D.f[BW  ])[kbw  ];
      //f_TW   = (D.f[BE  ])[kbe  ];
      //f_BE   = (D.f[TW  ])[ktw  ];
      //f_BS   = (D.f[TN  ])[ktn  ];
      //f_TN   = (D.f[BS  ])[kbs  ];
      //f_TS   = (D.f[BN  ])[kbn  ];
      //f_BN   = (D.f[TS  ])[kts  ];
      //f_BSW  = (D.f[TNE ])[ktne ];
      //f_BNE  = (D.f[TSW ])[ktsw ];
      //f_BNW  = (D.f[TSE ])[ktse ];
      //f_BSE  = (D.f[TNW ])[ktnw ];
      //f_TSW  = (D.f[BNE ])[kbne ];
      //f_TNE  = (D.f[BSW ])[kbsw ];
      //f_TNW  = (D.f[BSE ])[kbse ];
      //f_TSE  = (D.f[BNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]); 

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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
         //(D.f[dirREST])[k]=c1o10;
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
         (D.f[W])[kw]=((c1o1 - q) * f_E + q * ((f_E + f_W) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloX     )) / (q + c1o1) ;
         //(D.f[W])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[W])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[E])[ke]=((c1o1 - q) * f_W + q * ((f_W + f_E) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloX     )) / (q + c1o1) ;
         //(D.f[E])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[E])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[S])[ks]=((c1o1 - q) * f_N + q * ((f_N + f_S) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloY     )) / (q + c1o1) ;
         //(D.f[S])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[S])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[N])[kn]=((c1o1 - q) * f_S + q * ((f_S + f_N) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloY     )) / (q + c1o1) ;
         //(D.f[N])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[N])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[B])[kb]=((c1o1 - q) * f_T + q * ((f_T + f_B) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloZ     )) / (q + c1o1) ;
         //(D.f[B])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[B])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[T])[kt]=((c1o1 - q) * f_B + q * ((f_B + f_T) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloZ     )) / (q + c1o1) ;
         //(D.f[T])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[T])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SW])[ksw]=((c1o1 - q) * f_NE + q * ((f_NE + f_SW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloY)) / (q + c1o1) ;
         //(D.f[SW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[SW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NE])[kne]=((c1o1 - q) * f_SW + q * ((f_SW + f_NE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloY)) / (q + c1o1) ;
         //(D.f[NE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[NE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NW])[knw]=((c1o1 - q) * f_SE + q * ((f_SE + f_NW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloY)) / (q + c1o1) ;
         //(D.f[NW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[NW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SE])[kse]=((c1o1 - q) * f_NW + q * ((f_NW + f_SE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloY)) / (q + c1o1) ;
         //(D.f[SE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[SE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BW])[kbw]=((c1o1 - q) * f_TE + q * ((f_TE + f_BW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloZ)) / (q + c1o1) ;
         //(D.f[BW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[BW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TE])[kte]=((c1o1 - q) * f_BW + q * ((f_BW + f_TE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloZ)) / (q + c1o1) ;
         //(D.f[TE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[TE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TW])[ktw]=((c1o1 - q) * f_BE + q * ((f_BE + f_TW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloZ)) / (q + c1o1) ;
         //(D.f[TW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[TW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BE])[kbe]=((c1o1 - q) * f_TW + q * ((f_TW + f_BE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloZ)) / (q + c1o1) ;
         //(D.f[BE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[BE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BS])[kbs]=((c1o1 - q) * f_TN + q * ((f_TN + f_BS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[BS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[BS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TN])[ktn]=((c1o1 - q) * f_BS + q * ((f_BS + f_TN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[TN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[TN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TS])[kts]=((c1o1 - q) * f_BN + q * ((f_BN + f_TS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[TS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[TS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BN])[kbn]=((c1o1 - q) * f_TS + q * ((f_TS + f_BN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[BN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[BN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSW])[kbsw]=((c1o1 - q) * f_TNE + q * ((f_TNE + f_BSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[BSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[BSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNE])[ktne]=((c1o1 - q) * f_BSW + q * ((f_BSW + f_TNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[TNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[TNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSW])[ktsw]=((c1o1 - q) * f_BNE + q * ((f_BNE + f_TSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[TSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[TSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNE])[kbne]=((c1o1 - q) * f_TSW + q * ((f_TSW + f_BNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[BNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[BNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNW])[kbnw]=((c1o1 - q) * f_TSE + q * ((f_TSE + f_BNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[BNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[BNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSE])[ktse]=((c1o1 - q) * f_BNW + q * ((f_BNW + f_TSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[TSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[TSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNW])[ktnw]=((c1o1 - q) * f_BSE + q * ((f_BSE + f_TNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[TNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[TNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSE])[kbse]=((c1o1 - q) * f_TNW + q * ((f_TNW + f_BSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[BSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[BSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDeviceComp27(
										 real* distributions, 
										 int* subgridDistanceIndices, 
										 real* subgridDistances,
										 unsigned int numberOfBCnodes, 
										 real omega, 
										 unsigned int* neighborX,
										 unsigned int* neighborY,
										 unsigned int* neighborZ,
										 unsigned int numberOfLBnodes, 
										 bool isEvenTimestep)
{
   //////////////////////////////////////////////////////////////////////////
   //! The velocity boundary condition is executed in the following steps
   //!
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned  x = threadIdx.x;  // global x-index 
   const unsigned  y = blockIdx.x;   // global y-index 
   const unsigned  z = blockIdx.y;   // global z-index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;

   if(k < numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist;
      getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local subgrid distances (q's)
      //!
      SubgridDistances27 subgridD;
      getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);
      
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int indexOfBCnode  = subgridDistanceIndices[k];
      unsigned int kzero= indexOfBCnode;
      unsigned int ke   = indexOfBCnode;
      unsigned int kw   = neighborX[indexOfBCnode];
      unsigned int kn   = indexOfBCnode;
      unsigned int ks   = neighborY[indexOfBCnode];
      unsigned int kt   = indexOfBCnode;
      unsigned int kb   = neighborZ[indexOfBCnode];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = indexOfBCnode;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = indexOfBCnode;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = indexOfBCnode;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = indexOfBCnode;
      unsigned int kbsw = neighborZ[ksw];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local distributions
      //!
      real f_W    = (dist.f[E   ])[ke   ];
      real f_E    = (dist.f[W   ])[kw   ];
      real f_S    = (dist.f[N   ])[kn   ];
      real f_N    = (dist.f[S   ])[ks   ];
      real f_B    = (dist.f[T   ])[kt   ];
      real f_T    = (dist.f[B   ])[kb   ];
      real f_SW   = (dist.f[NE  ])[kne  ];
      real f_NE   = (dist.f[SW  ])[ksw  ];
      real f_NW   = (dist.f[SE  ])[kse  ];
      real f_SE   = (dist.f[NW  ])[knw  ];
      real f_BW   = (dist.f[TE  ])[kte  ];
      real f_TE   = (dist.f[BW  ])[kbw  ];
      real f_TW   = (dist.f[BE  ])[kbe  ];
      real f_BE   = (dist.f[TW  ])[ktw  ];
      real f_BS   = (dist.f[TN  ])[ktn  ];
      real f_TN   = (dist.f[BS  ])[kbs  ];
      real f_TS   = (dist.f[BN  ])[kbn  ];
      real f_BN   = (dist.f[TS  ])[kts  ];
      real f_BSW  = (dist.f[TNE ])[ktne ];
      real f_BNE  = (dist.f[TSW ])[ktsw ];
      real f_BNW  = (dist.f[TSE ])[ktse ];
      real f_BSE  = (dist.f[TNW ])[ktnw ];
      real f_TSW  = (dist.f[BNE ])[kbne ];
      real f_TNE  = (dist.f[BSW ])[kbsw ];
      real f_TNW  = (dist.f[BSE ])[kbse ];
      real f_TSE  = (dist.f[BNW ])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Calculate macroscopic quantities
      //!
      real drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                  f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[dirREST])[kzero]); 

      real vx1  = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                   (f_E - f_W)) / (c1o1 + drho);          

      real vx2  = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                   (f_N - f_S)) / (c1o1 + drho); 

      real vx3  = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                   (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                   (f_T - f_B)) / (c1o1 + drho); 

      real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3) * (c1o1 + drho);

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

       ////////////////////////////////////////////////////////////////////////////////
      //! - Update distributions with subgrid distance (q) between zero and one
      real feq, q, velocityLB;
      q = (subgridD.q[E])[k];
      if (q>=c0o1 && q<=c1o1) // only update distribution for q between zero and one
      {
         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[W])[kw] = getInterpolatedDistributionForNoSlipBC(q, f_E, f_W, feq, omega);
      }

      q = (subgridD.q[W])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1;
         
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[E])[ke] = getInterpolatedDistributionForNoSlipBC(q, f_W, f_E, feq, omega);
      }

      q = (subgridD.q[N])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[S])[ks] = getInterpolatedDistributionForNoSlipBC(q, f_N, f_S, feq, omega);
      }

      q = (subgridD.q[S])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[N])[kn] = getInterpolatedDistributionForNoSlipBC(q, f_S, f_N, feq, omega);
      }

      q = (subgridD.q[T])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[B])[kb] = getInterpolatedDistributionForNoSlipBC(q, f_T, f_B, feq, omega);
      }

      q = (subgridD.q[B])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[T])[kt] = getInterpolatedDistributionForNoSlipBC(q, f_B, f_T, feq, omega);
      }

      q = (subgridD.q[NE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[SW])[ksw] = getInterpolatedDistributionForNoSlipBC(q, f_NE, f_SW, feq, omega);
      }

      q = (subgridD.q[SW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[NE])[kne] = getInterpolatedDistributionForNoSlipBC(q, f_SW, f_NE, feq, omega);
      }

      q = (subgridD.q[SE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[NW])[knw] = getInterpolatedDistributionForNoSlipBC(q, f_SE, f_NW, feq, omega);
      }

      q = (subgridD.q[NW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[SE])[kse] = getInterpolatedDistributionForNoSlipBC(q, f_NW, f_SE, feq, omega);
      }

      q = (subgridD.q[TE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[BW])[kbw] = getInterpolatedDistributionForNoSlipBC(q, f_TE, f_BW, feq, omega);
      }

      q = (subgridD.q[BW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[TE])[kte] = getInterpolatedDistributionForNoSlipBC(q, f_BW, f_TE, feq, omega);
      }

      q = (subgridD.q[BE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[TW])[ktw] = getInterpolatedDistributionForNoSlipBC(q, f_BE, f_TW, feq, omega);
      }

      q = (subgridD.q[TW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[BE])[kbe] = getInterpolatedDistributionForNoSlipBC(q, f_TW, f_BE, feq, omega);
      }

      q = (subgridD.q[TN])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[BS])[kbs] = getInterpolatedDistributionForNoSlipBC(q, f_TN, f_BS, feq, omega);
      }

      q = (subgridD.q[BS])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[TN])[ktn] = getInterpolatedDistributionForNoSlipBC(q, f_BS, f_TN, feq, omega);
      }

      q = (subgridD.q[BN])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[TS])[kts] = getInterpolatedDistributionForNoSlipBC(q, f_BN, f_TS, feq, omega);
      }

      q = (subgridD.q[TS])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[BN])[kbn] = getInterpolatedDistributionForNoSlipBC(q, f_TS, f_BN, feq, omega);
      }

      q = (subgridD.q[TNE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[BSW])[kbsw] = getInterpolatedDistributionForNoSlipBC(q, f_TNE, f_BSW, feq, omega);
      }

      q = (subgridD.q[BSW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[TNE])[ktne] = getInterpolatedDistributionForNoSlipBC(q, f_BSW, f_TNE, feq, omega);
      }

      q = (subgridD.q[BNE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[TSW])[ktsw] = getInterpolatedDistributionForNoSlipBC(q, f_BNE, f_TSW, feq, omega);
      }

      q = (subgridD.q[TSW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[BNE])[kbne] = getInterpolatedDistributionForNoSlipBC(q, f_TSW, f_BNE, feq, omega);
      }

      q = (subgridD.q[TSE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[BNW])[kbnw] = getInterpolatedDistributionForNoSlipBC(q, f_TSE, f_BNW, feq, omega);
      }

      q = (subgridD.q[BNW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[TSE])[ktse] = getInterpolatedDistributionForNoSlipBC(q, f_BNW, f_TSE, feq, omega);
      }

      q = (subgridD.q[BSE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[TNW])[ktnw] = getInterpolatedDistributionForNoSlipBC(q, f_BSE, f_TNW, feq, omega);
      }

      q = (subgridD.q[TNW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[BSE])[kbse] = getInterpolatedDistributionForNoSlipBC(q, f_TNW, f_BSE, feq, omega);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QDevice27(real* distributions, 
                                     int* subgridDistanceIndices, 
                                     real* subgridDistances,
                                     unsigned int numberOfBCnodes, 
                                     real omega, 
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     unsigned int numberOfLBnodes, 
                                     bool isEvenTimestep)
{
   //////////////////////////////////////////////////////////////////////////
   //! The no-slip boundary condition is executed in the following steps
   //!
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned  x = threadIdx.x;  // global x-index 
   const unsigned  y = blockIdx.x;   // global y-index 
   const unsigned  z = blockIdx.y;   // global z-index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;

   //////////////////////////////////////////////////////////////////////////
   //! - Run for all indices in size of boundary condition (numberOfBCnodes)
   //!
   if(k < numberOfBCnodes)
   {

      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist;
      getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local subgrid distances (q's)
      //!
      SubgridDistances27 subgridD;
      getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int indexOfBCnode  = subgridDistanceIndices[k];
      unsigned int kzero= indexOfBCnode;
      unsigned int ke   = indexOfBCnode;
      unsigned int kw   = neighborX[indexOfBCnode];
      unsigned int kn   = indexOfBCnode;
      unsigned int ks   = neighborY[indexOfBCnode];
      unsigned int kt   = indexOfBCnode;
      unsigned int kb   = neighborZ[indexOfBCnode];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = indexOfBCnode;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = indexOfBCnode;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = indexOfBCnode;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = indexOfBCnode;
      unsigned int kbsw = neighborZ[ksw];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local distributions
      //!
      real f_W    = (dist.f[E   ])[ke   ];
      real f_E    = (dist.f[W   ])[kw   ];
      real f_S    = (dist.f[N   ])[kn   ];
      real f_N    = (dist.f[S   ])[ks   ];
      real f_B    = (dist.f[T   ])[kt   ];
      real f_T    = (dist.f[B   ])[kb   ];
      real f_SW   = (dist.f[NE  ])[kne  ];
      real f_NE   = (dist.f[SW  ])[ksw  ];
      real f_NW   = (dist.f[SE  ])[kse  ];
      real f_SE   = (dist.f[NW  ])[knw  ];
      real f_BW   = (dist.f[TE  ])[kte  ];
      real f_TE   = (dist.f[BW  ])[kbw  ];
      real f_TW   = (dist.f[BE  ])[kbe  ];
      real f_BE   = (dist.f[TW  ])[ktw  ];
      real f_BS   = (dist.f[TN  ])[ktn  ];
      real f_TN   = (dist.f[BS  ])[kbs  ];
      real f_TS   = (dist.f[BN  ])[kbn  ];
      real f_BN   = (dist.f[TS  ])[kts  ];
      real f_BSW  = (dist.f[TNE ])[ktne ];
      real f_BNE  = (dist.f[TSW ])[ktsw ];
      real f_BNW  = (dist.f[TSE ])[ktse ];
      real f_BSE  = (dist.f[TNW ])[ktnw ];
      real f_TSW  = (dist.f[BNE ])[kbne ];
      real f_TNE  = (dist.f[BSW ])[kbsw ];
      real f_TNW  = (dist.f[BSE ])[kbse ];
      real f_TSE  = (dist.f[BNW ])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Calculate macroscopic quantities
      //!
      real drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                  f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[dirREST])[kzero]); 

      real vx1  = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                   (f_E - f_W));          

      real vx2  = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                   (f_N - f_S)); 

      real vx3  = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                   (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                   (f_T - f_B)); 

      real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Update distributions with subgrid distance (q) between zero and one
      //!
      real feq, q, velocityLB;
      q = (subgridD.q[E])[k];
      if (q>=c0o1 && q<=c1o1) // only update distribution for q between zero and one
      {
         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[W])[kw] = getInterpolatedDistributionForNoSlipBC(q, f_E, f_W, feq, omega);
      }

      q = (subgridD.q[W])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[E])[ke] = getInterpolatedDistributionForNoSlipBC(q, f_W, f_E, feq, omega);
      }

      q = (subgridD.q[N])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[S])[ks] = getInterpolatedDistributionForNoSlipBC(q, f_N, f_S, feq, omega);
      }

      q = (subgridD.q[S])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[N])[kn] = getInterpolatedDistributionForNoSlipBC(q, f_S, f_N, feq, omega);
      }

      q = (subgridD.q[T])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[B])[kb] = getInterpolatedDistributionForNoSlipBC(q, f_T, f_B, feq, omega);
      }

      q = (subgridD.q[B])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[T])[kt] = getInterpolatedDistributionForNoSlipBC(q, f_B, f_T, feq, omega);
      }

      q = (subgridD.q[NE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[SW])[ksw] = getInterpolatedDistributionForNoSlipBC(q, f_NE, f_SW, feq, omega);
      }

      q = (subgridD.q[SW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[NE])[kne] = getInterpolatedDistributionForNoSlipBC(q, f_SW, f_NE, feq, omega);
      }

      q = (subgridD.q[SE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[NW])[knw] = getInterpolatedDistributionForNoSlipBC(q, f_SE, f_NW, feq, omega);
      }

      q = (subgridD.q[NW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[SE])[kse] = getInterpolatedDistributionForNoSlipBC(q, f_NW, f_SE, feq, omega);
      }

      q = (subgridD.q[TE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[BW])[kbw] = getInterpolatedDistributionForNoSlipBC(q, f_TE, f_BW, feq, omega);
      }

      q = (subgridD.q[BW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[TE])[kte] = getInterpolatedDistributionForNoSlipBC(q, f_BW, f_TE, feq, omega);
      }

      q = (subgridD.q[BE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[TW])[ktw] = getInterpolatedDistributionForNoSlipBC(q, f_BE, f_TW, feq, omega);
      }

      q = (subgridD.q[TW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[BE])[kbe] = getInterpolatedDistributionForNoSlipBC(q, f_TW, f_BE, feq, omega);
      }

      q = (subgridD.q[TN])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[BS])[kbs] = getInterpolatedDistributionForNoSlipBC(q, f_TN, f_BS, feq, omega);
      }

      q = (subgridD.q[BS])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[TN])[ktn] = getInterpolatedDistributionForNoSlipBC(q, f_BS, f_TN, feq, omega);
      }

      q = (subgridD.q[BN])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[TS])[kts] = getInterpolatedDistributionForNoSlipBC(q, f_BN, f_TS, feq, omega);
      }

      q = (subgridD.q[TS])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[BN])[kbn] = getInterpolatedDistributionForNoSlipBC(q, f_TS, f_BN, feq, omega);
      }

      q = (subgridD.q[TNE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[BSW])[kbsw] = getInterpolatedDistributionForNoSlipBC(q, f_TNE, f_BSW, feq, omega);
      }

      q = (subgridD.q[BSW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[TNE])[ktne] = getInterpolatedDistributionForNoSlipBC(q, f_BSW, f_TNE, feq, omega);
      }

      q = (subgridD.q[BNE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[TSW])[ktsw] = getInterpolatedDistributionForNoSlipBC(q, f_BNE, f_TSW, feq, omega);
      }

      q = (subgridD.q[TSW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[BNE])[kbne] = getInterpolatedDistributionForNoSlipBC(q, f_TSW, f_BNE, feq, omega);
      }

      q = (subgridD.q[TSE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[BNW])[kbnw] = getInterpolatedDistributionForNoSlipBC(q, f_TSE, f_BNW, feq, omega);
      }

      q = (subgridD.q[BNW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[TSE])[ktse] = getInterpolatedDistributionForNoSlipBC(q, f_BNW, f_TSE, feq, omega);
      }

      q = (subgridD.q[BSE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[TNW])[ktnw] = getInterpolatedDistributionForNoSlipBC(q, f_BSE, f_TNW, feq, omega);
      }

      q = (subgridD.q[TNW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[BSE])[kbse] = getInterpolatedDistributionForNoSlipBC(q, f_TNW, f_BSE, feq, omega);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void BBDevice27(real* distributions, 
                                     int* subgridDistanceIndices, 
                                     real* subgridDistances,
                                     unsigned int numberOfBCnodes, 
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     unsigned int numberOfLBnodes, 
                                     bool isEvenTimestep)
{
   //////////////////////////////////////////////////////////////////////////
   //! The no-slip boundary condition is executed in the following steps
   //!
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned  x = threadIdx.x;   // global x-index
   const unsigned  y = blockIdx.x;    // global y-index
   const unsigned  z = blockIdx.y;    // global z-index

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;

   //////////////////////////////////////////////////////////////////////////
   // run for all indices in size of boundary condition (numberOfBCnodes)
   if(k < numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist;
      getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local subgrid distances (q's)
      //!
      SubgridDistances27 subgridD;
      getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int indexOfBCnode  = subgridDistanceIndices[k];
      unsigned int ke   = indexOfBCnode;
      unsigned int kw   = neighborX[indexOfBCnode];
      unsigned int kn   = indexOfBCnode;
      unsigned int ks   = neighborY[indexOfBCnode];
      unsigned int kt   = indexOfBCnode;
      unsigned int kb   = neighborZ[indexOfBCnode];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = indexOfBCnode;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = indexOfBCnode;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = indexOfBCnode;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = indexOfBCnode;
      unsigned int kbsw = neighborZ[ksw];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local distributions
      //!
      real f_W    = (dist.f[E   ])[ke   ];
      real f_E    = (dist.f[W   ])[kw   ];
      real f_S    = (dist.f[N   ])[kn   ];
      real f_N    = (dist.f[S   ])[ks   ];
      real f_B    = (dist.f[T   ])[kt   ];
      real f_T    = (dist.f[B   ])[kb   ];
      real f_SW   = (dist.f[NE  ])[kne  ];
      real f_NE   = (dist.f[SW  ])[ksw  ];
      real f_NW   = (dist.f[SE  ])[kse  ];
      real f_SE   = (dist.f[NW  ])[knw  ];
      real f_BW   = (dist.f[TE  ])[kte  ];
      real f_TE   = (dist.f[BW  ])[kbw  ];
      real f_TW   = (dist.f[BE  ])[kbe  ];
      real f_BE   = (dist.f[TW  ])[ktw  ];
      real f_BS   = (dist.f[TN  ])[ktn  ];
      real f_TN   = (dist.f[BS  ])[kbs  ];
      real f_TS   = (dist.f[BN  ])[kbn  ];
      real f_BN   = (dist.f[TS  ])[kts  ];
      real f_BSW  = (dist.f[TNE ])[ktne ];
      real f_BNE  = (dist.f[TSW ])[ktsw ];
      real f_BNW  = (dist.f[TSE ])[ktse ];
      real f_BSE  = (dist.f[TNW ])[ktnw ];
      real f_TSW  = (dist.f[BNE ])[kbne ];
      real f_TNE  = (dist.f[BSW ])[kbsw ];
      real f_TNW  = (dist.f[BSE ])[kbse ];
      real f_TSE  = (dist.f[BNW ])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - rewrite distributions if there is a sub-grid distance (q) in same direction
      real q;
      q = (subgridD.q[E  ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[W  ])[kw  ]=f_E  ;
      q = (subgridD.q[W  ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[E  ])[ke  ]=f_W  ;
      q = (subgridD.q[N  ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[S  ])[ks  ]=f_N  ;
      q = (subgridD.q[S  ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[N  ])[kn  ]=f_S  ;
      q = (subgridD.q[T  ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[B  ])[kb  ]=f_T  ;
      q = (subgridD.q[B  ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[T  ])[kt  ]=f_B  ;
      q = (subgridD.q[NE ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[SW ])[ksw ]=f_NE ;
      q = (subgridD.q[SW ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[NE ])[kne ]=f_SW ;
      q = (subgridD.q[SE ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[NW ])[knw ]=f_SE ;
      q = (subgridD.q[NW ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[SE ])[kse ]=f_NW ;
      q = (subgridD.q[TE ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BW ])[kbw ]=f_TE ;
      q = (subgridD.q[BW ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TE ])[kte ]=f_BW ;
      q = (subgridD.q[BE ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TW ])[ktw ]=f_BE ;
      q = (subgridD.q[TW ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BE ])[kbe ]=f_TW ;
      q = (subgridD.q[TN ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BS ])[kbs ]=f_TN ;
      q = (subgridD.q[BS ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TN ])[ktn ]=f_BS ;
      q = (subgridD.q[BN ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TS ])[kts ]=f_BN ;
      q = (subgridD.q[TS ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BN ])[kbn ]=f_TS ;
      q = (subgridD.q[TNE])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BSW])[kbsw]=f_TNE;
      q = (subgridD.q[BSW])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TNE])[ktne]=f_BSW;
      q = (subgridD.q[BNE])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TSW])[ktsw]=f_BNE;
      q = (subgridD.q[TSW])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BNE])[kbne]=f_TSW;
      q = (subgridD.q[TSE])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BNW])[kbnw]=f_TSE;
      q = (subgridD.q[BNW])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TSE])[ktse]=f_BNW;
      q = (subgridD.q[BSE])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TNW])[ktnw]=f_BSE;
      q = (subgridD.q[TNW])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BSE])[kbse]=f_TNW;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

