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
#include "KernelUtilities.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

//////////////////////////////////////////////////////////////////////////////
__global__ void QDevice3rdMomentsComp27(
													 real* distributions, 
													 int* subgridDistanceIndices, 
													 real* subgridDistances,
													 unsigned int numberOfBCnodes, 
													 real omega, 
													 unsigned int* neighborX,
													 unsigned int* neighborY,
													 unsigned int* neighborZ,
													 unsigned long long numberOfLBnodes, 
													 bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &distributions[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &distributions[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &distributions[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &distributions[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &distributions[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &distributions[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &distributions[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &distributions[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &distributions[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &distributions[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &distributions[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &distributions[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &distributions[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &distributions[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &distributions[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &distributions[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &distributions[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &distributions[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &distributions[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &distributions[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &distributions[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &distributions[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &distributions[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &distributions[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &distributions[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &distributions[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &distributions[DIR_MPM * numberOfLBnodes];
   } 
   else
   {
      D.f[DIR_M00] = &distributions[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &distributions[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &distributions[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &distributions[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &distributions[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &distributions[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &distributions[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &distributions[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &distributions[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &distributions[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &distributions[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &distributions[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &distributions[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &distributions[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &distributions[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &distributions[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &distributions[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &distributions[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &distributions[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &distributions[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &distributions[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &distributions[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &distributions[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &distributions[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &distributions[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &distributions[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &distributions[DIR_PMP * numberOfLBnodes];
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
      q_dirE   = &subgridDistances[DIR_P00 * numberOfBCnodes];
      q_dirW   = &subgridDistances[DIR_M00 * numberOfBCnodes];
      q_dirN   = &subgridDistances[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &subgridDistances[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &subgridDistances[DIR_00P * numberOfBCnodes];
      q_dirB   = &subgridDistances[DIR_00M * numberOfBCnodes];
      q_dirNE  = &subgridDistances[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &subgridDistances[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &subgridDistances[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &subgridDistances[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &subgridDistances[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &subgridDistances[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &subgridDistances[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &subgridDistances[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &subgridDistances[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &subgridDistances[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &subgridDistances[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &subgridDistances[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &subgridDistances[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &subgridDistances[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &subgridDistances[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &subgridDistances[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &subgridDistances[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &subgridDistances[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &subgridDistances[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &subgridDistances[DIR_MPM * numberOfBCnodes];
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

      f_W    = (D.f[DIR_P00])[ke   ];
      f_E    = (D.f[DIR_M00])[kw   ];
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
      real vx1, vx2, vx3, drho, feq, q, m3;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]); 

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
         D.f[DIR_P00] = &distributions[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &distributions[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &distributions[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &distributions[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &distributions[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &distributions[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &distributions[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &distributions[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &distributions[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &distributions[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &distributions[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &distributions[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &distributions[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &distributions[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &distributions[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &distributions[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &distributions[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &distributions[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &distributions[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &distributions[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &distributions[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &distributions[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &distributions[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &distributions[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &distributions[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &distributions[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &distributions[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         D.f[DIR_M00] = &distributions[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &distributions[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &distributions[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &distributions[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &distributions[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &distributions[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &distributions[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &distributions[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &distributions[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &distributions[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &distributions[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &distributions[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &distributions[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &distributions[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &distributions[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &distributions[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &distributions[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &distributions[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &distributions[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &distributions[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &distributions[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &distributions[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &distributions[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &distributions[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &distributions[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &distributions[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &distributions[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
         //(D.f[DIR_000])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  //ToDo anders klammern !!!!!!
	  
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_E - f_W - c2o1 * drho * c2o27 * (c3o1*( vx1        ));
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_M00])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W-m3+(f_E+f_W-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_E+f_W))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_M00])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_W - f_E - c2o1 * drho * c2o27 * (c3o1*(-vx1        ));
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_P00])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E-m3+(f_W+f_E-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_W+f_E))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_P00])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_N - f_S - c2o1 * drho * c2o27 * (c3o1*( vx2        ));
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0M0])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S-m3+(f_N+f_S-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_N+f_S))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_0M0])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_S - f_N - c2o1 * drho * c2o27 * (c3o1*(   -vx2     ));
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0P0])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N-m3+(f_S+f_N-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_S+f_N))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_0P0])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_T - f_B - c2o1 * drho * c2o27 * (c3o1*(         vx3));
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_00M])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B-m3+(f_T+f_B-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_T+f_B))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_00M])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_B - f_T - c2o1 * drho * c2o27 * (c3o1*(        -vx3));
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_00P])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T-m3+(f_B+f_T-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_B+f_T))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_00P])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_NE - f_SW - c2o1 * drho * c1o54 * (c3o1*( vx1+vx2    ));
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MM0])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW-m3+(f_NE+f_SW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_NE+f_SW))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_MM0])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_SW - f_NE - c2o1 * drho * c1o54 * (c3o1*(-vx1-vx2    ));
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PP0])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE-m3+(f_SW+f_NE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_SW+f_NE))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_PP0])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_SE - f_NW - c2o1 * drho * c1o54 * (c3o1*( vx1-vx2    ));
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MP0])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW-m3+(f_SE+f_NW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_SE+f_NW))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_MP0])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_NW - f_SE - c2o1 * drho * c1o54 * (c3o1*(-vx1+vx2    ));
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PM0])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE-m3+(f_NW+f_SE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_NW+f_SE))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_PM0])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TE - f_BW - c2o1 * drho * c1o54 * (c3o1*( vx1    +vx3));
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_M0M])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW-m3+(f_TE+f_BW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TE+f_BW))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_M0M])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BW - f_TE - c2o1 * drho * c1o54 * (c3o1*(-vx1    -vx3));
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_P0P])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE-m3+(f_BW+f_TE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BW+f_TE))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_P0P])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BE - f_TW - c2o1 * drho * c1o54 * (c3o1*( vx1    -vx3));
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_M0P])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW-m3+(f_BE+f_TW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BE+f_TW))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_M0P])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TW - f_BE - c2o1 * drho * c1o54 * (c3o1*(-vx1    +vx3));
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_P0M])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE-m3+(f_TW+f_BE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TW+f_BE))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_P0M])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TN - f_BS - c2o1 * drho * c1o54 * (c3o1*(     vx2+vx3));
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0MM])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS-m3+(f_TN+f_BS-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TN+f_BS))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_0MM])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BS - f_TN - c2o1 * drho * c1o54 * (c3o1*(    -vx2-vx3));
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0PP])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN-m3+(f_BS+f_TN-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BS+f_TN))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_0PP])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BN - f_TS - c2o1 * drho * c1o54 * (c3o1*(     vx2-vx3));
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0MP])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS-m3+(f_BN+f_TS-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BN+f_TS))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_0MP])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TS - f_BN - c2o1 * drho * c1o54 * (c3o1*(    -vx2+vx3));
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0PM])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN-m3+(f_TS+f_BN-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TS+f_BN))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_0PM])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TNE - f_BSW - c2o1 * drho * c1o216 * (c3o1*( vx1+vx2+vx3));
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MMM])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW-m3+(f_TNE+f_BSW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TNE+f_BSW))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_MMM])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BSW - f_TNE - c2o1 * drho * c1o216 * (c3o1*(-vx1-vx2-vx3));
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PPP])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE-m3+(f_BSW+f_TNE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BSW+f_TNE))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_PPP])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BNE - f_TSW - c2o1 * drho * c1o216 * (c3o1*( vx1+vx2-vx3));
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MMP])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW-m3+(f_BNE+f_TSW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BNE+f_TSW))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_MMP])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TSW - f_BNE - c2o1 * drho * c1o216 * (c3o1*(-vx1-vx2+vx3));
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PPM])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE-m3+(f_TSW+f_BNE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TSW+f_BNE))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_PPM])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TSE - f_BNW - c2o1 * drho * c1o216 * (c3o1*( vx1-vx2+vx3));
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MPM])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW-m3+(f_TSE+f_BNW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TSE+f_BNW))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_MPM])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BNW - f_TSE - c2o1 * drho * c1o216 * (c3o1*(-vx1+vx2-vx3));
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PMP])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE-m3+(f_BNW+f_TSE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BNW+f_TSE))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_PMP])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_BSE - f_TNW - c2o1 * drho * c1o216 * (c3o1*( vx1-vx2-vx3));
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MPP])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW-m3+(f_BSE+f_TNW-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_BSE+f_TNW))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_MPP])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
		 m3 = f_TNW - f_BSE - c2o1 * drho * c1o216 * (c3o1*(-vx1+vx2+vx3));
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PMM])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE-m3+(f_TNW+f_BSE-c2o1*feq*omega)/(c1o1-omega))*c1o2+(q*(f_TNW+f_BSE))/(c1o1+q)+(m3*c1o2);
         //(D.f[DIR_PMM])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QDeviceIncompHighNu27(real* DD, 
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
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
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
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
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
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
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
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
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
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
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

      f_E   = (D.f[DIR_P00])[ke   ];
      f_W   = (D.f[DIR_M00])[kw   ];
      f_N   = (D.f[DIR_0P0])[kn   ];
      f_S   = (D.f[DIR_0M0])[ks   ];
      f_T   = (D.f[DIR_00P])[kt   ];
      f_B   = (D.f[DIR_00M])[kb   ];
      f_NE  = (D.f[DIR_PP0])[kne  ];
      f_SW  = (D.f[DIR_MM0])[ksw  ];
      f_SE  = (D.f[DIR_PM0])[kse  ];
      f_NW  = (D.f[DIR_MP0])[knw  ];
      f_TE  = (D.f[DIR_P0P])[kte  ];
      f_BW  = (D.f[DIR_M0M])[kbw  ];
      f_BE  = (D.f[DIR_P0M])[kbe  ];
      f_TW  = (D.f[DIR_M0P])[ktw  ];
      f_TN  = (D.f[DIR_0PP])[ktn  ];
      f_BS  = (D.f[DIR_0MM])[kbs  ];
      f_BN  = (D.f[DIR_0PM])[kbn  ];
      f_TS  = (D.f[DIR_0MP])[kts  ];
      f_TNE = (D.f[DIR_PPP])[ktne ];
      f_TSW = (D.f[DIR_MMP])[ktsw ];
      f_TSE = (D.f[DIR_PMP])[ktse ];
      f_TNW = (D.f[DIR_MPP])[ktnw ];
      f_BNE = (D.f[DIR_PPM])[kbne ];
      f_BSW = (D.f[DIR_MMM])[kbsw ];
      f_BSE = (D.f[DIR_PMM])[kbse ];
      f_BNW = (D.f[DIR_MPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]); 

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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
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
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
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
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
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
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
         //(D.f[DIR_000])[k]=c1o10;
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
         (D.f[DIR_M00])[kw]=((c1o1 - q) * f_E + q * ((f_E + f_W) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloX     )) / (q + c1o1) ;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_P00])[ke]=((c1o1 - q) * f_W + q * ((f_W + f_E) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloX     )) / (q + c1o1) ;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_0M0])[ks]=((c1o1 - q) * f_N + q * ((f_N + f_S) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloY     )) / (q + c1o1) ;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_0P0])[kn]=((c1o1 - q) * f_S + q * ((f_S + f_N) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloY     )) / (q + c1o1) ;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_00M])[kb]=((c1o1 - q) * f_T + q * ((f_T + f_B) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloZ     )) / (q + c1o1) ;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_00P])[kt]=((c1o1 - q) * f_B + q * ((f_B + f_T) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloZ     )) / (q + c1o1) ;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_MM0])[ksw]=((c1o1 - q) * f_NE + q * ((f_NE + f_SW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloY)) / (q + c1o1) ;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_PP0])[kne]=((c1o1 - q) * f_SW + q * ((f_SW + f_NE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloY)) / (q + c1o1) ;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_MP0])[knw]=((c1o1 - q) * f_SE + q * ((f_SE + f_NW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloY)) / (q + c1o1) ;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_PM0])[kse]=((c1o1 - q) * f_NW + q * ((f_NW + f_SE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloY)) / (q + c1o1) ;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_M0M])[kbw]=((c1o1 - q) * f_TE + q * ((f_TE + f_BW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_P0P])[kte]=((c1o1 - q) * f_BW + q * ((f_BW + f_TE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_M0P])[ktw]=((c1o1 - q) * f_BE + q * ((f_BE + f_TW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_P0M])[kbe]=((c1o1 - q) * f_TW + q * ((f_TW + f_BE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_0MM])[kbs]=((c1o1 - q) * f_TN + q * ((f_TN + f_BS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_0PP])[ktn]=((c1o1 - q) * f_BS + q * ((f_BS + f_TN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_0MP])[kts]=((c1o1 - q) * f_BN + q * ((f_BN + f_TS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_0PM])[kbn]=((c1o1 - q) * f_TS + q * ((f_TS + f_BN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_MMM])[kbsw]=((c1o1 - q) * f_TNE + q * ((f_TNE + f_BSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_PPP])[ktne]=((c1o1 - q) * f_BSW + q * ((f_BSW + f_TNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_MMP])[ktsw]=((c1o1 - q) * f_BNE + q * ((f_BNE + f_TSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_PPM])[kbne]=((c1o1 - q) * f_TSW + q * ((f_TSW + f_BNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_MPM])[kbnw]=((c1o1 - q) * f_TSE + q * ((f_TSE + f_BNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_PMP])[ktse]=((c1o1 - q) * f_BNW + q * ((f_BNW + f_TSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_MPP])[ktnw]=((c1o1 - q) * f_BSE + q * ((f_BSE + f_TNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[DIR_PMM])[kbse]=((c1o1 - q) * f_TNW + q * ((f_TNW + f_BSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY+VeloZ)) / (q + c1o1) ;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QDeviceCompHighNu27(
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
   if (isEvenTimestep==true)
   {
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
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
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
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
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
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
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
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
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
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

      f_E   = (D.f[DIR_P00])[ke   ];
      f_W   = (D.f[DIR_M00])[kw   ];
      f_N   = (D.f[DIR_0P0])[kn   ];
      f_S   = (D.f[DIR_0M0])[ks   ];
      f_T   = (D.f[DIR_00P])[kt   ];
      f_B   = (D.f[DIR_00M])[kb   ];
      f_NE  = (D.f[DIR_PP0])[kne  ];
      f_SW  = (D.f[DIR_MM0])[ksw  ];
      f_SE  = (D.f[DIR_PM0])[kse  ];
      f_NW  = (D.f[DIR_MP0])[knw  ];
      f_TE  = (D.f[DIR_P0P])[kte  ];
      f_BW  = (D.f[DIR_M0M])[kbw  ];
      f_BE  = (D.f[DIR_P0M])[kbe  ];
      f_TW  = (D.f[DIR_M0P])[ktw  ];
      f_TN  = (D.f[DIR_0PP])[ktn  ];
      f_BS  = (D.f[DIR_0MM])[kbs  ];
      f_BN  = (D.f[DIR_0PM])[kbn  ];
      f_TS  = (D.f[DIR_0MP])[kts  ];
      f_TNE = (D.f[DIR_PPP])[ktne ];
      f_TSW = (D.f[DIR_MMP])[ktsw ];
      f_TSE = (D.f[DIR_PMP])[ktse ];
      f_TNW = (D.f[DIR_MPP])[ktnw ];
      f_BNE = (D.f[DIR_PPM])[kbne ];
      f_BSW = (D.f[DIR_MMM])[kbsw ];
      f_BSE = (D.f[DIR_PMM])[kbse ];
      f_BNW = (D.f[DIR_MPM])[kbnw ];
      //f_W    = (D.f[DIR_P00])[ke   ];
      //f_E    = (D.f[DIR_M00])[kw   ];
      //f_S    = (D.f[DIR_0P0])[kn   ];
      //f_N    = (D.f[DIR_0M0])[ks   ];
      //f_B    = (D.f[DIR_00P])[kt   ];
      //f_T    = (D.f[DIR_00M])[kb   ];
      //f_SW   = (D.f[DIR_PP0])[kne  ];
      //f_NE   = (D.f[DIR_MM0])[ksw  ];
      //f_NW   = (D.f[DIR_PM0])[kse  ];
      //f_SE   = (D.f[DIR_MP0])[knw  ];
      //f_BW   = (D.f[DIR_P0P])[kte  ];
      //f_TE   = (D.f[DIR_M0M])[kbw  ];
      //f_TW   = (D.f[DIR_P0M])[kbe  ];
      //f_BE   = (D.f[DIR_M0P])[ktw  ];
      //f_BS   = (D.f[DIR_0PP])[ktn  ];
      //f_TN   = (D.f[DIR_0MM])[kbs  ];
      //f_TS   = (D.f[DIR_0PM])[kbn  ];
      //f_BN   = (D.f[DIR_0MP])[kts  ];
      //f_BSW  = (D.f[DIR_PPP])[ktne ];
      //f_BNE  = (D.f[DIR_MMP])[ktsw ];
      //f_BNW  = (D.f[DIR_PMP])[ktse ];
      //f_BSE  = (D.f[DIR_MPP])[ktnw ];
      //f_TSW  = (D.f[DIR_PPM])[kbne ];
      //f_TNE  = (D.f[DIR_MMM])[kbsw ];
      //f_TNW  = (D.f[DIR_PMM])[kbse ];
      //f_TSE  = (D.f[DIR_MPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]); 

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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
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
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
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
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
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
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
         //(D.f[DIR_000])[k]=c1o10;
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
         (D.f[DIR_M00])[kw]=((c1o1 - q) * f_E + q * ((f_E + f_W) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloX     )) / (q + c1o1) ;
         //(D.f[DIR_M00])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[DIR_M00])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_P00])[ke]=((c1o1 - q) * f_W + q * ((f_W + f_E) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloX     )) / (q + c1o1) ;
         //(D.f[DIR_P00])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[DIR_P00])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0M0])[ks]=((c1o1 - q) * f_N + q * ((f_N + f_S) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloY     )) / (q + c1o1) ;
         //(D.f[DIR_0M0])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[DIR_0M0])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0P0])[kn]=((c1o1 - q) * f_S + q * ((f_S + f_N) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloY     )) / (q + c1o1) ;
         //(D.f[DIR_0P0])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[DIR_0P0])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_00M])[kb]=((c1o1 - q) * f_T + q * ((f_T + f_B) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloZ     )) / (q + c1o1) ;
         //(D.f[DIR_00M])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[DIR_00M])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_00P])[kt]=((c1o1 - q) * f_B + q * ((f_B + f_T) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloZ     )) / (q + c1o1) ;
         //(D.f[DIR_00P])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[DIR_00P])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MM0])[ksw]=((c1o1 - q) * f_NE + q * ((f_NE + f_SW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloY)) / (q + c1o1) ;
         //(D.f[DIR_MM0])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[DIR_MM0])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PP0])[kne]=((c1o1 - q) * f_SW + q * ((f_SW + f_NE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloY)) / (q + c1o1) ;
         //(D.f[DIR_PP0])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[DIR_PP0])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MP0])[knw]=((c1o1 - q) * f_SE + q * ((f_SE + f_NW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloY)) / (q + c1o1) ;
         //(D.f[DIR_MP0])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[DIR_MP0])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PM0])[kse]=((c1o1 - q) * f_NW + q * ((f_NW + f_SE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloY)) / (q + c1o1) ;
         //(D.f[DIR_PM0])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[DIR_PM0])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_M0M])[kbw]=((c1o1 - q) * f_TE + q * ((f_TE + f_BW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_M0M])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[DIR_M0M])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_P0P])[kte]=((c1o1 - q) * f_BW + q * ((f_BW + f_TE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_P0P])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[DIR_P0P])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_M0P])[ktw]=((c1o1 - q) * f_BE + q * ((f_BE + f_TW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_M0P])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[DIR_M0P])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_P0M])[kbe]=((c1o1 - q) * f_TW + q * ((f_TW + f_BE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_P0M])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[DIR_P0M])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0MM])[kbs]=((c1o1 - q) * f_TN + q * ((f_TN + f_BS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_0MM])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[DIR_0MM])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0PP])[ktn]=((c1o1 - q) * f_BS + q * ((f_BS + f_TN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_0PP])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[DIR_0PP])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0MP])[kts]=((c1o1 - q) * f_BN + q * ((f_BN + f_TS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_0MP])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[DIR_0MP])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0PM])[kbn]=((c1o1 - q) * f_TS + q * ((f_TS + f_BN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_0PM])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[DIR_0PM])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MMM])[kbsw]=((c1o1 - q) * f_TNE + q * ((f_TNE + f_BSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_MMM])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[DIR_MMM])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PPP])[ktne]=((c1o1 - q) * f_BSW + q * ((f_BSW + f_TNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_PPP])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[DIR_PPP])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MMP])[ktsw]=((c1o1 - q) * f_BNE + q * ((f_BNE + f_TSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_MMP])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[DIR_MMP])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PPM])[kbne]=((c1o1 - q) * f_TSW + q * ((f_TSW + f_BNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_PPM])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[DIR_PPM])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MPM])[kbnw]=((c1o1 - q) * f_TSE + q * ((f_TSE + f_BNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_MPM])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[DIR_MPM])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PMP])[ktse]=((c1o1 - q) * f_BNW + q * ((f_BNW + f_TSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_PMP])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[DIR_PMP])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MPP])[ktnw]=((c1o1 - q) * f_BSE + q * ((f_BSE + f_TNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_MPP])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[DIR_MPP])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PMM])[kbse]=((c1o1 - q) * f_TNW + q * ((f_TNW + f_BSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[DIR_PMM])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[DIR_PMM])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QDeviceComp27(
										 real* distributions, 
										 int* subgridDistanceIndices, 
										 real* subgridDistances,
										 unsigned int numberOfBCnodes, 
										 real omega, 
										 unsigned int* neighborX,
										 unsigned int* neighborY,
										 unsigned int* neighborZ,
										 unsigned long long numberOfLBnodes, 
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
      real f_W    = (dist.f[DIR_P00])[ke   ];
      real f_E    = (dist.f[DIR_M00])[kw   ];
      real f_S    = (dist.f[DIR_0P0])[kn   ];
      real f_N    = (dist.f[DIR_0M0])[ks   ];
      real f_B    = (dist.f[DIR_00P])[kt   ];
      real f_T    = (dist.f[DIR_00M])[kb   ];
      real f_SW   = (dist.f[DIR_PP0])[kne  ];
      real f_NE   = (dist.f[DIR_MM0])[ksw  ];
      real f_NW   = (dist.f[DIR_PM0])[kse  ];
      real f_SE   = (dist.f[DIR_MP0])[knw  ];
      real f_BW   = (dist.f[DIR_P0P])[kte  ];
      real f_TE   = (dist.f[DIR_M0M])[kbw  ];
      real f_TW   = (dist.f[DIR_P0M])[kbe  ];
      real f_BE   = (dist.f[DIR_M0P])[ktw  ];
      real f_BS   = (dist.f[DIR_0PP])[ktn  ];
      real f_TN   = (dist.f[DIR_0MM])[kbs  ];
      real f_TS   = (dist.f[DIR_0PM])[kbn  ];
      real f_BN   = (dist.f[DIR_0MP])[kts  ];
      real f_BSW  = (dist.f[DIR_PPP])[ktne ];
      real f_BNE  = (dist.f[DIR_MMP])[ktsw ];
      real f_BNW  = (dist.f[DIR_PMP])[ktse ];
      real f_BSE  = (dist.f[DIR_MPP])[ktnw ];
      real f_TSW  = (dist.f[DIR_PPM])[kbne ];
      real f_TNE  = (dist.f[DIR_MMM])[kbsw ];
      real f_TNW  = (dist.f[DIR_PMM])[kbse ];
      real f_TSE  = (dist.f[DIR_MPM])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Calculate macroscopic quantities
      //!
      real drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                  f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[DIR_000])[kzero]); 

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
      q = (subgridD.q[DIR_P00])[k];
      if (q>=c0o1 && q<=c1o1) // only update distribution for q between zero and one
      {
         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[DIR_M00])[kw] = getInterpolatedDistributionForNoSlipBC(q, f_E, f_W, feq, omega);
      }

      q = (subgridD.q[DIR_M00])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[DIR_P00])[ke] = getInterpolatedDistributionForNoSlipBC(q, f_W, f_E, feq, omega);
      }

      q = (subgridD.q[DIR_0P0])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[DIR_0M0])[ks] = getInterpolatedDistributionForNoSlipBC(q, f_N, f_S, feq, omega);
      }

      q = (subgridD.q[DIR_0M0])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[DIR_0P0])[kn] = getInterpolatedDistributionForNoSlipBC(q, f_S, f_N, feq, omega);
      }

      q = (subgridD.q[DIR_00P])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[DIR_00M])[kb] = getInterpolatedDistributionForNoSlipBC(q, f_T, f_B, feq, omega);
      }

      q = (subgridD.q[DIR_00M])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[DIR_00P])[kt] = getInterpolatedDistributionForNoSlipBC(q, f_B, f_T, feq, omega);
      }

      q = (subgridD.q[DIR_PP0])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_MM0])[ksw] = getInterpolatedDistributionForNoSlipBC(q, f_NE, f_SW, feq, omega);
      }

      q = (subgridD.q[DIR_MM0])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_PP0])[kne] = getInterpolatedDistributionForNoSlipBC(q, f_SW, f_NE, feq, omega);
      }

      q = (subgridD.q[DIR_PM0])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_MP0])[knw] = getInterpolatedDistributionForNoSlipBC(q, f_SE, f_NW, feq, omega);
      }

      q = (subgridD.q[DIR_MP0])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_PM0])[kse] = getInterpolatedDistributionForNoSlipBC(q, f_NW, f_SE, feq, omega);
      }

      q = (subgridD.q[DIR_P0P])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_M0M])[kbw] = getInterpolatedDistributionForNoSlipBC(q, f_TE, f_BW, feq, omega);
      }

      q = (subgridD.q[DIR_M0M])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_P0P])[kte] = getInterpolatedDistributionForNoSlipBC(q, f_BW, f_TE, feq, omega);
      }

      q = (subgridD.q[DIR_P0M])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_M0P])[ktw] = getInterpolatedDistributionForNoSlipBC(q, f_BE, f_TW, feq, omega);
      }

      q = (subgridD.q[DIR_M0P])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_P0M])[kbe] = getInterpolatedDistributionForNoSlipBC(q, f_TW, f_BE, feq, omega);
      }

      q = (subgridD.q[DIR_0PP])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_0MM])[kbs] = getInterpolatedDistributionForNoSlipBC(q, f_TN, f_BS, feq, omega);
      }

      q = (subgridD.q[DIR_0MM])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_0PP])[ktn] = getInterpolatedDistributionForNoSlipBC(q, f_BS, f_TN, feq, omega);
      }

      q = (subgridD.q[DIR_0PM])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_0MP])[kts] = getInterpolatedDistributionForNoSlipBC(q, f_BN, f_TS, feq, omega);
      }

      q = (subgridD.q[DIR_0MP])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_0PM])[kbn] = getInterpolatedDistributionForNoSlipBC(q, f_TS, f_BN, feq, omega);
      }

      q = (subgridD.q[DIR_PPP])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_MMM])[kbsw] = getInterpolatedDistributionForNoSlipBC(q, f_TNE, f_BSW, feq, omega);
      }

      q = (subgridD.q[DIR_MMM])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_PPP])[ktne] = getInterpolatedDistributionForNoSlipBC(q, f_BSW, f_TNE, feq, omega);
      }

      q = (subgridD.q[DIR_PPM])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_MMP])[ktsw] = getInterpolatedDistributionForNoSlipBC(q, f_BNE, f_TSW, feq, omega);
      }

      q = (subgridD.q[DIR_MMP])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_PPM])[kbne] = getInterpolatedDistributionForNoSlipBC(q, f_TSW, f_BNE, feq, omega);
      }

      q = (subgridD.q[DIR_PMP])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_MPM])[kbnw] = getInterpolatedDistributionForNoSlipBC(q, f_TSE, f_BNW, feq, omega);
      }

      q = (subgridD.q[DIR_MPM])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_PMP])[ktse] = getInterpolatedDistributionForNoSlipBC(q, f_BNW, f_TSE, feq, omega);
      }

      q = (subgridD.q[DIR_PMM])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_MPP])[ktnw] = getInterpolatedDistributionForNoSlipBC(q, f_BSE, f_TNW, feq, omega);
      }

      q = (subgridD.q[DIR_MPP])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_PMM])[kbse] = getInterpolatedDistributionForNoSlipBC(q, f_TNW, f_BSE, feq, omega);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QDevice27(real* distributions, 
                                     int* subgridDistanceIndices, 
                                     real* subgridDistances,
                                     unsigned int numberOfBCnodes, 
                                     real omega, 
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     unsigned long long numberOfLBnodes, 
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
      real f_W    = (dist.f[DIR_P00])[ke   ];
      real f_E    = (dist.f[DIR_M00])[kw   ];
      real f_S    = (dist.f[DIR_0P0])[kn   ];
      real f_N    = (dist.f[DIR_0M0])[ks   ];
      real f_B    = (dist.f[DIR_00P])[kt   ];
      real f_T    = (dist.f[DIR_00M])[kb   ];
      real f_SW   = (dist.f[DIR_PP0])[kne  ];
      real f_NE   = (dist.f[DIR_MM0])[ksw  ];
      real f_NW   = (dist.f[DIR_PM0])[kse  ];
      real f_SE   = (dist.f[DIR_MP0])[knw  ];
      real f_BW   = (dist.f[DIR_P0P])[kte  ];
      real f_TE   = (dist.f[DIR_M0M])[kbw  ];
      real f_TW   = (dist.f[DIR_P0M])[kbe  ];
      real f_BE   = (dist.f[DIR_M0P])[ktw  ];
      real f_BS   = (dist.f[DIR_0PP])[ktn  ];
      real f_TN   = (dist.f[DIR_0MM])[kbs  ];
      real f_TS   = (dist.f[DIR_0PM])[kbn  ];
      real f_BN   = (dist.f[DIR_0MP])[kts  ];
      real f_BSW  = (dist.f[DIR_PPP])[ktne ];
      real f_BNE  = (dist.f[DIR_MMP])[ktsw ];
      real f_BNW  = (dist.f[DIR_PMP])[ktse ];
      real f_BSE  = (dist.f[DIR_MPP])[ktnw ];
      real f_TSW  = (dist.f[DIR_PPM])[kbne ];
      real f_TNE  = (dist.f[DIR_MMM])[kbsw ];
      real f_TNW  = (dist.f[DIR_PMM])[kbse ];
      real f_TSE  = (dist.f[DIR_MPM])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Calculate macroscopic quantities
      //!
      real drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                  f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[DIR_000])[kzero]); 

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
      q = (subgridD.q[DIR_P00])[k];
      if (q>=c0o1 && q<=c1o1) // only update distribution for q between zero and one
      {
         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[DIR_M00])[kw] = getInterpolatedDistributionForNoSlipBC(q, f_E, f_W, feq, omega);
      }

      q = (subgridD.q[DIR_M00])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[DIR_P00])[ke] = getInterpolatedDistributionForNoSlipBC(q, f_W, f_E, feq, omega);
      }

      q = (subgridD.q[DIR_0P0])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[DIR_0M0])[ks] = getInterpolatedDistributionForNoSlipBC(q, f_N, f_S, feq, omega);
      }

      q = (subgridD.q[DIR_0M0])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[DIR_0P0])[kn] = getInterpolatedDistributionForNoSlipBC(q, f_S, f_N, feq, omega);
      }

      q = (subgridD.q[DIR_00P])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[DIR_00M])[kb] = getInterpolatedDistributionForNoSlipBC(q, f_T, f_B, feq, omega);
      }

      q = (subgridD.q[DIR_00M])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         (dist.f[DIR_00P])[kt] = getInterpolatedDistributionForNoSlipBC(q, f_B, f_T, feq, omega);
      }

      q = (subgridD.q[DIR_PP0])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_MM0])[ksw] = getInterpolatedDistributionForNoSlipBC(q, f_NE, f_SW, feq, omega);
      }

      q = (subgridD.q[DIR_MM0])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_PP0])[kne] = getInterpolatedDistributionForNoSlipBC(q, f_SW, f_NE, feq, omega);
      }

      q = (subgridD.q[DIR_PM0])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_MP0])[knw] = getInterpolatedDistributionForNoSlipBC(q, f_SE, f_NW, feq, omega);
      }

      q = (subgridD.q[DIR_MP0])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_PM0])[kse] = getInterpolatedDistributionForNoSlipBC(q, f_NW, f_SE, feq, omega);
      }

      q = (subgridD.q[DIR_P0P])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_M0M])[kbw] = getInterpolatedDistributionForNoSlipBC(q, f_TE, f_BW, feq, omega);
      }

      q = (subgridD.q[DIR_M0M])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_P0P])[kte] = getInterpolatedDistributionForNoSlipBC(q, f_BW, f_TE, feq, omega);
      }

      q = (subgridD.q[DIR_P0M])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_M0P])[ktw] = getInterpolatedDistributionForNoSlipBC(q, f_BE, f_TW, feq, omega);
      }

      q = (subgridD.q[DIR_M0P])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_P0M])[kbe] = getInterpolatedDistributionForNoSlipBC(q, f_TW, f_BE, feq, omega);
      }

      q = (subgridD.q[DIR_0PP])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_0MM])[kbs] = getInterpolatedDistributionForNoSlipBC(q, f_TN, f_BS, feq, omega);
      }

      q = (subgridD.q[DIR_0MM])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_0PP])[ktn] = getInterpolatedDistributionForNoSlipBC(q, f_BS, f_TN, feq, omega);
      }

      q = (subgridD.q[DIR_0PM])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_0MP])[kts] = getInterpolatedDistributionForNoSlipBC(q, f_BN, f_TS, feq, omega);
      }

      q = (subgridD.q[DIR_0MP])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         (dist.f[DIR_0PM])[kbn] = getInterpolatedDistributionForNoSlipBC(q, f_TS, f_BN, feq, omega);
      }

      q = (subgridD.q[DIR_PPP])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_MMM])[kbsw] = getInterpolatedDistributionForNoSlipBC(q, f_TNE, f_BSW, feq, omega);
      }

      q = (subgridD.q[DIR_MMM])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_PPP])[ktne] = getInterpolatedDistributionForNoSlipBC(q, f_BSW, f_TNE, feq, omega);
      }

      q = (subgridD.q[DIR_PPM])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_MMP])[ktsw] = getInterpolatedDistributionForNoSlipBC(q, f_BNE, f_TSW, feq, omega);
      }

      q = (subgridD.q[DIR_MMP])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_PPM])[kbne] = getInterpolatedDistributionForNoSlipBC(q, f_TSW, f_BNE, feq, omega);
      }

      q = (subgridD.q[DIR_PMP])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_MPM])[kbnw] = getInterpolatedDistributionForNoSlipBC(q, f_TSE, f_BNW, feq, omega);
      }

      q = (subgridD.q[DIR_MPM])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_PMP])[ktse] = getInterpolatedDistributionForNoSlipBC(q, f_BNW, f_TSE, feq, omega);
      }

      q = (subgridD.q[DIR_PMM])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_MPP])[ktnw] = getInterpolatedDistributionForNoSlipBC(q, f_BSE, f_TNW, feq, omega);
      }

      q = (subgridD.q[DIR_MPP])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         (dist.f[DIR_PMM])[kbse] = getInterpolatedDistributionForNoSlipBC(q, f_TNW, f_BSE, feq, omega);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void BBDevice27(real* distributions, 
                                     int* subgridDistanceIndices, 
                                     real* subgridDistances,
                                     unsigned int numberOfBCnodes, 
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     unsigned long long numberOfLBnodes, 
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
      real f_W    = (dist.f[DIR_P00])[ke   ];
      real f_E    = (dist.f[DIR_M00])[kw   ];
      real f_S    = (dist.f[DIR_0P0])[kn   ];
      real f_N    = (dist.f[DIR_0M0])[ks   ];
      real f_B    = (dist.f[DIR_00P])[kt   ];
      real f_T    = (dist.f[DIR_00M])[kb   ];
      real f_SW   = (dist.f[DIR_PP0])[kne  ];
      real f_NE   = (dist.f[DIR_MM0])[ksw  ];
      real f_NW   = (dist.f[DIR_PM0])[kse  ];
      real f_SE   = (dist.f[DIR_MP0])[knw  ];
      real f_BW   = (dist.f[DIR_P0P])[kte  ];
      real f_TE   = (dist.f[DIR_M0M])[kbw  ];
      real f_TW   = (dist.f[DIR_P0M])[kbe  ];
      real f_BE   = (dist.f[DIR_M0P])[ktw  ];
      real f_BS   = (dist.f[DIR_0PP])[ktn  ];
      real f_TN   = (dist.f[DIR_0MM])[kbs  ];
      real f_TS   = (dist.f[DIR_0PM])[kbn  ];
      real f_BN   = (dist.f[DIR_0MP])[kts  ];
      real f_BSW  = (dist.f[DIR_PPP])[ktne ];
      real f_BNE  = (dist.f[DIR_MMP])[ktsw ];
      real f_BNW  = (dist.f[DIR_PMP])[ktse ];
      real f_BSE  = (dist.f[DIR_MPP])[ktnw ];
      real f_TSW  = (dist.f[DIR_PPM])[kbne ];
      real f_TNE  = (dist.f[DIR_MMM])[kbsw ];
      real f_TNW  = (dist.f[DIR_PMM])[kbse ];
      real f_TSE  = (dist.f[DIR_MPM])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - rewrite distributions if there is a sub-grid distance (q) in same direction
      real q;
      q = (subgridD.q[DIR_P00])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_M00])[kw  ]=f_E  ;
      q = (subgridD.q[DIR_M00])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_P00])[ke  ]=f_W  ;
      q = (subgridD.q[DIR_0P0])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_0M0])[ks  ]=f_N  ;
      q = (subgridD.q[DIR_0M0])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_0P0])[kn  ]=f_S  ;
      q = (subgridD.q[DIR_00P])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_00M])[kb  ]=f_T  ;
      q = (subgridD.q[DIR_00M])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_00P])[kt  ]=f_B  ;
      q = (subgridD.q[DIR_PP0])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_MM0])[ksw ]=f_NE ;
      q = (subgridD.q[DIR_MM0])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_PP0])[kne ]=f_SW ;
      q = (subgridD.q[DIR_PM0])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_MP0])[knw ]=f_SE ;
      q = (subgridD.q[DIR_MP0])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_PM0])[kse ]=f_NW ;
      q = (subgridD.q[DIR_P0P])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_M0M])[kbw ]=f_TE ;
      q = (subgridD.q[DIR_M0M])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_P0P])[kte ]=f_BW ;
      q = (subgridD.q[DIR_P0M])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_M0P])[ktw ]=f_BE ;
      q = (subgridD.q[DIR_M0P])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_P0M])[kbe ]=f_TW ;
      q = (subgridD.q[DIR_0PP])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_0MM])[kbs ]=f_TN ;
      q = (subgridD.q[DIR_0MM])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_0PP])[ktn ]=f_BS ;
      q = (subgridD.q[DIR_0PM])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_0MP])[kts ]=f_BN ;
      q = (subgridD.q[DIR_0MP])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_0PM])[kbn ]=f_TS ;
      q = (subgridD.q[DIR_PPP])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_MMM])[kbsw]=f_TNE;
      q = (subgridD.q[DIR_MMM])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_PPP])[ktne]=f_BSW;
      q = (subgridD.q[DIR_PPM])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_MMP])[ktsw]=f_BNE;
      q = (subgridD.q[DIR_MMP])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_PPM])[kbne]=f_TSW;
      q = (subgridD.q[DIR_PMP])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_MPM])[kbnw]=f_TSE;
      q = (subgridD.q[DIR_MPM])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_PMP])[ktse]=f_BNW;
      q = (subgridD.q[DIR_PMM])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_MPP])[ktnw]=f_BSE;
      q = (subgridD.q[DIR_MPP])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_PMM])[kbse]=f_TNW;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

