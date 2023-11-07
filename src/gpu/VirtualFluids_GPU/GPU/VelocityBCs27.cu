//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file VelocityBCs27.cu
//! \ingroup GPU
//! \author Martin Schoenherr, Anna Wellmann
//======================================================================================
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include "basics/constants/NumericConstants.h"
#include "LBM/GPUHelperFunctions/KernelUtilities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

//////////////////////////////////////////////////////////////////////////////
__global__ void QVelDeviceCompPlusSlip27(
    real* vx,
    real* vy,
    real* vz,
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
	   real VeloX = vx[k];
	   real VeloY = vy[k];
	   real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
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

	   real cu_sq = c3o2*(vx1*vx1 + vx2*vx2 + vx3*vx3) * (c1o1 + drho);

	   //////////////////////////////////////////////////////////////////////////
	   if (isEvenTimestep == false)
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
	   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	   //Test
	   //(D.f[d000])[k]=c1o10;
	   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	   //ToDo anders Klammern

	   /////To Slip Or Not To Slip?
	   // We assume slip BC if |vec(V_BC)|=1. To avoid problems we take V_BC*V_BC>0.99 (c99o100)
	   if (VeloX*VeloX + VeloY*VeloY + VeloZ*VeloZ > c99o100)
		{
		   // vt=v-(n \dot v) *n
		   // n=(VeloX,VeloY,VeloZ) a misuse of the velocity variable!
		   real normalV = VeloX*vx1 + VeloY*vx2 + VeloZ*vx3;
		   vx1 = vx1 - normalV*VeloX;
		   vx2 = vx2 - normalV*VeloY;
		   vx3 = vx3 - normalV*VeloZ;
		}
	  ////////////////

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dM00])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q);// - c2over27 * drho;
         //(D.f[dM00])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dP00])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q);// - c2over27 * drho;
         //(D.f[dP00])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0M0])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q);// - c2over27 * drho;
         //(D.f[DIR_0M0])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0P0])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q);// - c2over27 * drho;
         //(D.f[DIR_0P0])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_00M])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q);// - c2over27 * drho;
         //(D.f[DIR_00M])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_00P])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q);// - c2over27 * drho;
         //(D.f[DIR_00P])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MM0])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q);// - c1over54 * drho;
         //(D.f[DIR_MM0])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PP0])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);// - c1over54 * drho;
         //(D.f[DIR_PP0])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MP0])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);// - c1over54 * drho;
         //(D.f[DIR_MP0])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PM0])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);// - c1over54 * drho;
         //(D.f[DIR_PM0])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_M0M])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[DIR_M0M])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_P0P])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[DIR_P0P])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_M0P])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[DIR_M0P])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_P0M])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[DIR_P0M])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0MM])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[DIR_0MM])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0PP])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[DIR_0PP])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0MP])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[DIR_0MP])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0PM])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[DIR_0PM])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MMM])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[DIR_MMM])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PPP])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[DIR_PPP])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MMP])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[DIR_MMP])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PPM])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[DIR_PPM])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MPM])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[DIR_MPM])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PMP])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[DIR_PMP])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MPP])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[DIR_MPP])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PMM])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[DIR_PMM])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void QVeloDeviceEQ27(
    real* VeloX,
    real* VeloY,
    real* VeloZ,
    real* DD, 
    int* k_Q, 
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

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // based on BGK Plus Comp
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[dP00])[ke   ];
			real mfabb = (D.f[dM00])[kw   ];
			real mfbcb = (D.f[DIR_0P0])[kn   ];
			real mfbab = (D.f[DIR_0M0])[ks   ];
			real mfbbc = (D.f[DIR_00P])[kt   ];
			real mfbba = (D.f[DIR_00M])[kb   ];
			real mfccb = (D.f[DIR_PP0])[kne  ];
			real mfaab = (D.f[DIR_MM0])[ksw  ];
			real mfcab = (D.f[DIR_PM0])[kse  ];
			real mfacb = (D.f[DIR_MP0])[knw  ];
			real mfcbc = (D.f[DIR_P0P])[kte  ];
			real mfaba = (D.f[DIR_M0M])[kbw  ];
			real mfcba = (D.f[DIR_P0M])[kbe  ];
			real mfabc = (D.f[DIR_M0P])[ktw  ];
			real mfbcc = (D.f[DIR_0PP])[ktn  ];
			real mfbaa = (D.f[DIR_0MM])[kbs  ];
			real mfbca = (D.f[DIR_0PM])[kbn  ];
			real mfbac = (D.f[DIR_0MP])[kts  ];
			real mfbbb = (D.f[d000])[kzero];
			real mfccc = (D.f[DIR_PPP])[ktne ];
			real mfaac = (D.f[DIR_MMP])[ktsw ];
			real mfcac = (D.f[DIR_PMP])[ktse ];
			real mfacc = (D.f[DIR_MPP])[ktnw ];
			real mfcca = (D.f[DIR_PPM])[kbne ];
			real mfaaa = (D.f[DIR_MMM])[kbsw ];
			real mfcaa = (D.f[DIR_PMM])[kbse ];
			real mfaca = (D.f[DIR_MPM])[kbnw ];
			////////////////////////////////////////////////////////////////////////////////////
			real rho   = (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
							 mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
							 mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb + c1o1);//!!!!Achtung + one
			////////////////////////////////////////////////////////////////////////////////////
			real vvx    = ((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) / rho;
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) / rho;
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			if(VeloX[k]!=c0o1) vvx = VeloX[k];
			if(VeloY[k]!=c0o1) vvy = VeloY[k];
			if(VeloZ[k]!=c0o1) vvz = VeloZ[k];
			////////////////////////////////////////////////////////////////////////////////////
			real vx2    = vvx * vvx;
			real vy2    = vvy * vvy;
			real vz2    = vvz * vvz;
			////////////////////////////////////////////////////////////////////////////////////
            real XXb    = -c2o3 + vx2;
            real XXc    = -c1o2 * (XXb + c1o1 + vvx);
            real XXa    = XXc + vvx;
            real YYb    = -c2o3 + vy2;
            real YYc    = -c1o2 * (YYb + c1o1 + vvy);
            real YYa    = YYc + vvy;
            real ZZb    = -c2o3 + vz2;
            real ZZc    = -c1o2 * (ZZb + c1o1 + vvz);
            real ZZa    = ZZc + vvz;
			////////////////////////////////////////////////////////////////////////////////////
            mfcbb = -rho * XXc * YYb * ZZb - c2o27 ; 
			mfabb = -rho * XXa * YYb * ZZb - c2o27 ;
			mfbcb = -rho * XXb * YYc * ZZb - c2o27 ;
			mfbab = -rho * XXb * YYa * ZZb - c2o27 ;
			mfbbc = -rho * XXb * YYb * ZZc - c2o27 ;
			mfbba = -rho * XXb * YYb * ZZa - c2o27 ;
			mfccb = -rho * XXc * YYc * ZZb - c1o54 ;
			mfaab = -rho * XXa * YYa * ZZb - c1o54 ;
			mfcab = -rho * XXc * YYa * ZZb - c1o54 ;
			mfacb = -rho * XXa * YYc * ZZb - c1o54 ;
			mfcbc = -rho * XXc * YYb * ZZc - c1o54 ;
			mfaba = -rho * XXa * YYb * ZZa - c1o54 ;
			mfcba = -rho * XXc * YYb * ZZa - c1o54 ;
			mfabc = -rho * XXa * YYb * ZZc - c1o54 ;
			mfbcc = -rho * XXb * YYc * ZZc - c1o54 ;
			mfbaa = -rho * XXb * YYa * ZZa - c1o54 ;
			mfbca = -rho * XXb * YYc * ZZa - c1o54 ;
			mfbac = -rho * XXb * YYa * ZZc - c1o54 ;
			mfbbb = -rho * XXb * YYb * ZZb - c8o27 ;
			mfccc = -rho * XXc * YYc * ZZc - c1o216;
			mfaac = -rho * XXa * YYa * ZZc - c1o216;
			mfcac = -rho * XXc * YYa * ZZc - c1o216;
			mfacc = -rho * XXa * YYc * ZZc - c1o216;
			mfcca = -rho * XXc * YYc * ZZa - c1o216;
			mfaaa = -rho * XXa * YYa * ZZa - c1o216;
			mfcaa = -rho * XXc * YYa * ZZa - c1o216;
			mfaca = -rho * XXa * YYc * ZZa - c1o216;
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			(D.f[dP00])[ke   ] = mfabb;//mfcbb;
			(D.f[dM00])[kw   ] = mfcbb;//mfabb;
			(D.f[DIR_0P0])[kn   ] = mfbab;//mfbcb;
			(D.f[DIR_0M0])[ks   ] = mfbcb;//mfbab;
			(D.f[DIR_00P])[kt   ] = mfbba;//mfbbc;
			(D.f[DIR_00M])[kb   ] = mfbbc;//mfbba;
			(D.f[DIR_PP0])[kne  ] = mfaab;//mfccb;
			(D.f[DIR_MM0])[ksw  ] = mfccb;//mfaab;
			(D.f[DIR_PM0])[kse  ] = mfacb;//mfcab;
			(D.f[DIR_MP0])[knw  ] = mfcab;//mfacb;
			(D.f[DIR_P0P])[kte  ] = mfaba;//mfcbc;
			(D.f[DIR_M0M])[kbw  ] = mfcbc;//mfaba;
			(D.f[DIR_P0M])[kbe  ] = mfabc;//mfcba;
			(D.f[DIR_M0P])[ktw  ] = mfcba;//mfabc;
			(D.f[DIR_0PP])[ktn  ] = mfbaa;//mfbcc;
			(D.f[DIR_0MM])[kbs  ] = mfbcc;//mfbaa;
			(D.f[DIR_0PM])[kbn  ] = mfbac;//mfbca;
			(D.f[DIR_0MP])[kts  ] = mfbca;//mfbac;
			(D.f[d000])[kzero] = mfbbb;//mfbbb;
			(D.f[DIR_PPP])[ktne ] = mfaaa;//mfccc;
			(D.f[DIR_MMP])[ktsw ] = mfcca;//mfaac;
			(D.f[DIR_PMP])[ktse ] = mfaca;//mfcac;
			(D.f[DIR_MPP])[ktnw ] = mfcaa;//mfacc;
			(D.f[DIR_PPM])[kbne ] = mfaac;//mfcca;
			(D.f[DIR_MMM])[kbsw ] = mfccc;//mfaaa;
			(D.f[DIR_PMM])[kbse ] = mfacc;//mfcaa;
			(D.f[DIR_MPM])[kbnw ] = mfcac;//mfaca;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////













































































//////////////////////////////////////////////////////////////////////////////
__global__ void QVelDeviceIncompHighNu27(
    real* vx,
    real* vy,
    real* vz,
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
      real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
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
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_E   = (D.f[dP00])[ke   ];
      f_W   = (D.f[dM00])[kw   ];
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
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W));// / (one + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S));// / (one + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B));// / (one + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);// * (one + drho);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[d000])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //ToDo anders Klammern

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) /** (one + drho)*/-cu_sq); 
         (D.f[dM00])[kw]=((c1o1 - q) * f_E + q * ((f_E + f_W) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloX     )) / (q + c1o1) ;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) /** (one + drho)*/-cu_sq); 
         (D.f[dP00])[ke]=((c1o1 - q) * f_W + q * ((f_W + f_E) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloX     )) / (q + c1o1) ;
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
__global__ void QVelDeviceCompHighNu27(
    real* vx,
    real* vy,
    real* vz,
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
      real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
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
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_E   = (D.f[dP00])[ke   ];
      f_W   = (D.f[dM00])[kw   ];
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
      //f_W    = (D.f[dP00])[ke   ];
      //f_E    = (D.f[dM00])[kw   ];
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
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[d000])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //ToDo anders Klammern

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dM00])[kw]=((c1o1 - q) * f_E + q * ((f_E + f_W) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloX     )) / (q + c1o1) ;
         //(D.f[dM00])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[dM00])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dP00])[ke]=((c1o1 - q) * f_W + q * ((f_W + f_E) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloX     )) / (q + c1o1) ;
         //(D.f[dP00])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[dP00])[ke]=zero;
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
__global__ void QVelDeviceCompZeroPress27(
    real* velocityX,
    real* velocityY,
    real* velocityZ,
    real* distribution, 
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
   //! The velocity boundary condition is executed in the following steps
   //!
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

   //////////////////////////////////////////////////////////////////////////
   //! - Run for all indices in size of boundary condition (numberOfBCnodes)
   //!
   if(nodeIndex < numberOfBCnodes)
   {

      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist;
      getPointersToDistributions(dist, distribution, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local velocities
      //!
      real VeloX = velocityX[nodeIndex];
      real VeloY = velocityY[nodeIndex];
      real VeloZ = velocityZ[nodeIndex];


      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local subgrid distances (q's)
      //!
      SubgridDistances27 subgridD;
      getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);
     
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int KQK  = subgridDistanceIndices[nodeIndex];
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
      //! - Set local distributions
      //!
      real f_W    = (dist.f[dP00])[ke   ];
      real f_E    = (dist.f[dM00])[kw   ];
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
      real drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                     f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                     f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[d000])[kzero]); 

      real vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                      ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                      (f_E - f_W)) / (c1o1 + drho); 
         

      real vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                       ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                       (f_N - f_S)) / (c1o1 + drho); 

      real vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                       (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                       (f_T - f_B)) / (c1o1 + drho); 
    
      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distribution, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Update distributions with subgrid distance (q) between zero and one
      real feq, q, velocityLB, velocityBC;
      q = (subgridD.q[dP00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1) // only update distribution for q between zero and one
      {
         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloX;
         (dist.f[dM00])[kw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_E, f_W, feq, omega, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[dM00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloX;
         (dist.f[dP00])[ke] = getInterpolatedDistributionForVeloWithPressureBC(q, f_W, f_E, feq, omega, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[DIR_0P0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloY;
         (dist.f[DIR_0M0])[ks] = getInterpolatedDistributionForVeloWithPressureBC(q, f_N, f_S, feq, omega, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[DIR_0M0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloY;
         (dist.f[DIR_0P0])[kn] = getInterpolatedDistributionForVeloWithPressureBC(q, f_S, f_N, feq, omega, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[DIR_00P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloZ;
         (dist.f[DIR_00M])[kb] = getInterpolatedDistributionForVeloWithPressureBC(q, f_T, f_B, feq, omega, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[DIR_00M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloZ;
         (dist.f[DIR_00P])[kt] = getInterpolatedDistributionForVeloWithPressureBC(q, f_B, f_T, feq, omega, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[DIR_PP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloY;
         (dist.f[DIR_MM0])[ksw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_NE, f_SW, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_MM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloY;
         (dist.f[DIR_PP0])[kne] = getInterpolatedDistributionForVeloWithPressureBC(q, f_SW, f_NE, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_PM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloY;
         (dist.f[DIR_MP0])[knw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_SE, f_NW, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_MP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloY;
         (dist.f[DIR_PM0])[kse] = getInterpolatedDistributionForVeloWithPressureBC(q, f_NW, f_SE, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_P0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloZ;
         (dist.f[DIR_M0M])[kbw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TE, f_BW, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_M0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloZ;
         (dist.f[DIR_P0P])[kte] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BW, f_TE, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_P0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloZ;
         (dist.f[DIR_M0P])[ktw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BE, f_TW, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_M0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloZ;
         (dist.f[DIR_P0M])[kbe] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TW, f_BE, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_0PP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY + VeloZ;
         (dist.f[DIR_0MM])[kbs] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TN, f_BS, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_0MM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY - VeloZ;
         (dist.f[DIR_0PP])[ktn] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BS, f_TN, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_0PM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY - VeloZ;
         (dist.f[DIR_0MP])[kts] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BN, f_TS, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_0MP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY + VeloZ;
         (dist.f[DIR_0PM])[kbn] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TS, f_BN, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_PPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY + VeloZ;
         (dist.f[DIR_MMM])[kbsw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TNE, f_BSW, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_MMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY - VeloZ;
         (dist.f[DIR_PPP])[ktne] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BSW, f_TNE, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_PPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY - VeloZ;
         (dist.f[DIR_MMP])[ktsw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BNE, f_TSW, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_MMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY + VeloZ;
         (dist.f[DIR_PPM])[kbne] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TSW, f_BNE, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_PMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY + VeloZ;
         (dist.f[DIR_MPM])[kbnw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TSE, f_BNW, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_MPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY - VeloZ;
         (dist.f[DIR_PMP])[ktse] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BNW, f_TSE, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_PMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY - VeloZ;
         (dist.f[DIR_MPP])[ktnw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BSE, f_TNW, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_MPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY + VeloZ;
         (dist.f[DIR_PMM])[kbse] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TNW, f_BSE, feq, omega, drho, velocityBC, c1o216);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








































//////////////////////////////////////////////////////////////////////////////
__global__ void QVelDeviceCompZeroPress1h27(
    int inx,
    int iny,
    real* vx,
    real* vy,
    real* vz,
    real* DD, 
    int* k_Q, 
    real* QQ,
    unsigned int numberOfBCnodes,
    real om1, 
    real Phi,
    real angularVelocity,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    real* coordX,
    real* coordY,
    real* coordZ,
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
      //real VeloX = vx[k];
      //real VeloY = vy[k];
      //real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
		////////////////////////////////////////////////////////////////////////////////
		real VeloX = cosf(Phi)*vx[k] - sinf(Phi)*vy[k];
		real VeloY = sinf(Phi)*vx[k] + cosf(Phi)*vy[k];
		//real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
		////////////////////////////////////////////////////////////////////////////////////
		//Ship
		real coord0X = 281.125f;//7.5f;
		real coord0Y = 388.125f;//7.5f;
		real ux = - angularVelocity * (coordY[k_Q[k]] - coord0Y);
		real uy =   angularVelocity * (coordX[k_Q[k]] - coord0X);
		real VeloXpur=VeloX;
		real VeloYpur=VeloY;
		VeloX-=ux;
		VeloY-=uy;
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
      unsigned int KQK  = k_Q[k];
      //unsigned int kzero= KQK;
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
      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real vx1, vx2, vx3, drho, feq, q, cu_sq;
	  ///////// equilibrium BC
	  cu_sq=c3o2*(VeloX*VeloX +VeloY*VeloY);
	  VeloXpur*=-c1o1;
	  VeloYpur*=-c1o1;
	  vx1=VeloX;
	  vx2=VeloY;
	  vx3=c0o1;
	  drho=c0o1;

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c3o1*( VeloXpur        )+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dM00])[kw]= feq - c2o27 * drho;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c3o1*(-VeloXpur        )+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dP00])[ke]= feq - c2o27 * drho;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c3o1*(    VeloYpur     )+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0M0])[ks]= feq - c2o27 * drho;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c3o1*(   -VeloYpur     )+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0P0])[kn]= feq - c2o27 * drho;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_00M])[kb]= feq - c2o27 * drho;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_00P])[kt]= feq - c2o27 * drho;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*( VeloXpur+VeloYpur    )+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MM0])[ksw]= feq - c1o54 * drho;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(-VeloXpur-VeloYpur    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PP0])[kne]= feq - c1o54 * drho;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*( VeloXpur-VeloYpur    )+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MP0])[knw]= feq - c1o54 * drho;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(-VeloXpur+VeloYpur    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PM0])[kse]= feq - c1o54 * drho;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*( VeloXpur    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_M0M])[kbw]= feq - c1o54 * drho;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(-VeloXpur    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_P0P])[kte]= feq - c1o54 * drho;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*( VeloXpur    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_M0P])[ktw]= feq - c1o54 * drho;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(-VeloXpur    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_P0M])[kbe]= feq - c1o54 * drho;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(     VeloYpur+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0MM])[kbs]= feq - c1o54 * drho;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(    -VeloYpur-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0PP])[ktn]= feq - c1o54 * drho;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(     VeloYpur-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0MP])[kts]= feq - c1o54 * drho;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(    -VeloYpur+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_0PM])[kbn]= feq - c1o54 * drho;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*( VeloXpur+VeloYpur+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MMM])[kbsw]= feq - c1o216 * drho;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*(-VeloXpur-VeloYpur-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PPP])[ktne]= feq - c1o216 * drho;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*( VeloXpur+VeloYpur-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MMP])[ktsw]= feq - c1o216 * drho;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*(-VeloXpur-VeloYpur+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PPM])[kbne]= feq - c1o216 * drho;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*( VeloXpur-VeloYpur+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MPM])[kbnw]= feq - c1o216 * drho;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*(-VeloXpur+VeloYpur-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PMP])[ktse]= feq - c1o216 * drho;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*( VeloXpur-VeloYpur-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_MPP])[ktnw]= feq - c1o216 * drho;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*(-VeloXpur+VeloYpur+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[DIR_PMM])[kbse]= feq - c1o216 * drho;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void LB_BC_Vel_West_27(
    int nx, 
    int ny, 
    int nz, 
    int itz, 
    unsigned int* bcMatD, 
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    real* DD, 
    unsigned long long numberOfLBnodes, 
    bool isEvenTimestep, 
    real u0x, 
    unsigned int grid_nx, 
    unsigned int grid_ny, 
    real om) 
{
   //thread-index
   unsigned int ity = blockIdx.x;
   unsigned int itx = threadIdx.x;

   unsigned int  k, nxny;                   // Zugriff auf arrays im device

   unsigned int  x = itx + STARTOFFX;  // Globaler x-Index 
   unsigned int  y = ity + STARTOFFY;  // Globaler y-Index 
   unsigned int  z = itz + STARTOFFZ;  // Globaler z-Index 

   k = nx*(ny*z + y) + x;
   nxny = nx*ny;
   unsigned int k1 = k+nxny;

   if( bcMatD[k] == GEO_VELO )
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
      //index
      unsigned int kzero= k;
      unsigned int ke   = k;
      unsigned int kw   = neighborX[k];
      unsigned int kn   = k;
      unsigned int ks   = neighborY[k];
      unsigned int kt   = k;
      unsigned int kb   = neighborZ[k];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = k;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = k;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = k;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = k;
      unsigned int kbsw = neighborZ[ksw];
      //unsigned int kzero= k;
      //unsigned int ke   = k;
      //unsigned int kw   = k + 1;
      //unsigned int kn   = k;
      //unsigned int ks   = k + nx;
      //unsigned int kt   = k;
      //unsigned int kb   = k + nxny;
      //unsigned int ksw  = k + nx + 1;
      //unsigned int kne  = k;
      //unsigned int kse  = k + nx;
      //unsigned int knw  = k + 1;
      //unsigned int kbw  = k + nxny + 1;
      //unsigned int kte  = k;
      //unsigned int kbe  = k + nxny;
      //unsigned int ktw  = k + 1;
      //unsigned int kbs  = k + nxny + nx;
      //unsigned int ktn  = k;
      //unsigned int kbn  = k + nxny;
      //unsigned int kts  = k + nx;
      //unsigned int ktse = k + nx;
      //unsigned int kbnw = k + nxny + 1;
      //unsigned int ktnw = k + 1;
      //unsigned int kbse = k + nxny + nx;
      //unsigned int ktsw = k + nx + 1;
      //unsigned int kbne = k + nxny;
      //unsigned int ktne = k;
      //unsigned int kbsw = k + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
      //index1
      unsigned int k1zero= k1;
      unsigned int k1e   = k1;
      unsigned int k1w   = neighborX[k1];
      unsigned int k1n   = k1;
      unsigned int k1s   = neighborY[k1];
      unsigned int k1t   = k1;
      unsigned int k1b   = neighborZ[k1];
      unsigned int k1sw  = neighborY[k1w];
      unsigned int k1ne  = k1;
      unsigned int k1se  = k1s;
      unsigned int k1nw  = k1w;
      unsigned int k1bw  = neighborZ[k1w];
      unsigned int k1te  = k1;
      unsigned int k1be  = k1b;
      unsigned int k1tw  = k1w;
      unsigned int k1bs  = neighborZ[k1s];
      unsigned int k1tn  = k1;
      unsigned int k1bn  = k1b;
      unsigned int k1ts  = k1s;
      unsigned int k1tse = k1s;
      unsigned int k1bnw = k1bw;
      unsigned int k1tnw = k1w;
      unsigned int k1bse = k1bs;
      unsigned int k1tsw = k1sw;
      unsigned int k1bne = k1b;
      unsigned int k1tne = k1;
      unsigned int k1bsw = neighborZ[k1sw];
      //unsigned int k1zero= k1;
      //unsigned int k1e   = k1;
      //unsigned int k1w   = k1 + 1;
      //unsigned int k1n   = k1;
      //unsigned int k1s   = k1 + nx;
      //unsigned int k1t   = k1;
      //unsigned int k1b   = k1 + nxny;
      //unsigned int k1sw  = k1 + nx + 1;
      //unsigned int k1ne  = k1;
      //unsigned int k1se  = k1 + nx;
      //unsigned int k1nw  = k1 + 1;
      //unsigned int k1bw  = k1 + nxny + 1;
      //unsigned int k1te  = k1;
      //unsigned int k1be  = k1 + nxny;
      //unsigned int k1tw  = k1 + 1;
      //unsigned int k1bs  = k1 + nxny + nx;
      //unsigned int k1tn  = k1;
      //unsigned int k1bn  = k1 + nxny;
      //unsigned int k1ts  = k1 + nx;
      //unsigned int k1tse = k1 + nx;
      //unsigned int k1bnw = k1 + nxny + 1;
      //unsigned int k1tnw = k1 + 1;
      //unsigned int k1bse = k1 + nxny + nx;
      //unsigned int k1tsw = k1 + nx + 1;
      //unsigned int k1bne = k1 + nxny;
      //unsigned int k1tne = k1;
      //unsigned int k1bsw = k1 + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
         f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[dP00])[k1e   ];
      f1_E    = (D.f[dM00])[k1w   ];
      f1_S    = (D.f[DIR_0P0])[k1n   ];
      f1_N    = (D.f[DIR_0M0])[k1s   ];
      f1_B    = (D.f[DIR_00P])[k1t   ];
      f1_T    = (D.f[DIR_00M])[k1b   ];
      f1_SW   = (D.f[DIR_PP0])[k1ne  ];
      f1_NE   = (D.f[DIR_MM0])[k1sw  ];
      f1_NW   = (D.f[DIR_PM0])[k1se  ];
      f1_SE   = (D.f[DIR_MP0])[k1nw  ];
      f1_BW   = (D.f[DIR_P0P])[k1te  ];
      f1_TE   = (D.f[DIR_M0M])[k1bw  ];
      f1_TW   = (D.f[DIR_P0M])[k1be  ];
      f1_BE   = (D.f[DIR_M0P])[k1tw  ];
      f1_BS   = (D.f[DIR_0PP])[k1tn  ];
      f1_TN   = (D.f[DIR_0MM])[k1bs  ];
      f1_TS   = (D.f[DIR_0PM])[k1bn  ];
      f1_BN   = (D.f[DIR_0MP])[k1ts  ];
      f1_ZERO = (D.f[d000])[k1zero];
      f1_BSW  = (D.f[DIR_PPP])[k1tne ];
      f1_BNE  = (D.f[DIR_MMP])[k1tsw ];
      f1_BNW  = (D.f[DIR_PMP])[k1tse ];
      f1_BSE  = (D.f[DIR_MPP])[k1tnw ];
      f1_TSW  = (D.f[DIR_PPM])[k1bne ];
      f1_TNE  = (D.f[DIR_MMM])[k1bsw ];
      f1_TNW  = (D.f[DIR_PMM])[k1bse ];
      f1_TSE  = (D.f[DIR_MPM])[k1bnw ];

      real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
         f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

      __syncthreads();

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real drho = drho1;
      real  vx1 = c0o1;
      real  vx2 = c0o1;
      real  vx3 = u0x;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      (D.f[d000])[kzero] =   c8o27* (drho-cu_sq);
      (D.f[dP00])[ke   ] =   c2o27* (drho+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      (D.f[dM00])[kw   ] =   c2o27* (drho+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      (D.f[DIR_0P0])[kn   ] =   c2o27* (drho+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      (D.f[DIR_0M0])[ks   ] =   c2o27* (drho+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      (D.f[DIR_00P])[kt   ] =   c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      (D.f[DIR_00M])[kb   ] =   c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      (D.f[DIR_PP0])[kne  ] =   c1o54* (drho+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      (D.f[DIR_MM0])[ksw  ] =   c1o54* (drho+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      (D.f[DIR_PM0])[kse  ] =   c1o54* (drho+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      (D.f[DIR_MP0])[knw  ] =   c1o54* (drho+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      (D.f[DIR_P0P])[kte  ] =   c1o54* (drho+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      (D.f[DIR_M0M])[kbw  ] =   c1o54* (drho+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      (D.f[DIR_P0M])[kbe  ] =   c1o54* (drho+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      (D.f[DIR_M0P])[ktw  ] =   c1o54* (drho+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      (D.f[DIR_0PP])[ktn  ] =   c1o54* (drho+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      (D.f[DIR_0MM])[kbs  ] =   c1o54* (drho+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      (D.f[DIR_0PM])[kbn  ] =   c1o54* (drho+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      (D.f[DIR_0MP])[kts  ] =   c1o54* (drho+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      (D.f[DIR_PPP])[ktne ] =   c1o216*(drho+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      (D.f[DIR_MMM])[kbsw ] =   c1o216*(drho+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      (D.f[DIR_PPM])[kbne ] =   c1o216*(drho+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      (D.f[DIR_MMP])[ktsw ] =   c1o216*(drho+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      (D.f[DIR_PMP])[ktse ] =   c1o216*(drho+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      (D.f[DIR_MPM])[kbnw ] =   c1o216*(drho+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      (D.f[DIR_PMM])[kbse ] =   c1o216*(drho+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      (D.f[DIR_MPP])[ktnw ] =   c1o216*(drho+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
   }
   __syncthreads();
}          
//////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QVelDevPlainBB27(
    real* velocityX,
    real* velocityY,
    real* velocityZ,
    real* distributions,
    int* subgridDistanceIndices,
    real* subgridDistances,
    uint numberOfBCnodes,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
   //////////////////////////////////////////////////////////////////////////
   //! The velocity boundary condition is executed in the following steps
   //!
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

   //////////////////////////////////////////////////////////////////////////
   // run for all indices in size of boundary condition (numberOfBCnodes)
   if(nodeIndex < numberOfBCnodes)
   {
       //////////////////////////////////////////////////////////////////////////
       //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
       //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
       //!
      Distributions27 dist;
      getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local velocities
      //!
      real VeloX = velocityX[nodeIndex];
      real VeloY = velocityY[nodeIndex];
      real VeloZ = velocityZ[nodeIndex];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local subgrid distances (q's)
      //!
      SubgridDistances27 subgridD;
      getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      uint indexOfBCnode = subgridDistanceIndices[nodeIndex];
      uint ke   = indexOfBCnode;
      uint kw   = neighborX[indexOfBCnode];
      uint kn   = indexOfBCnode;
      uint ks   = neighborY[indexOfBCnode];
      uint kt   = indexOfBCnode;
      uint kb   = neighborZ[indexOfBCnode];
      uint ksw  = neighborY[kw];
      uint kne  = indexOfBCnode;
      uint kse  = ks;
      uint knw  = kw;
      uint kbw  = neighborZ[kw];
      uint kte  = indexOfBCnode;
      uint kbe  = kb;
      uint ktw  = kw;
      uint kbs  = neighborZ[ks];
      uint ktn  = indexOfBCnode;
      uint kbn  = kb;
      uint kts  = ks;
      uint ktse = ks;
      uint kbnw = kbw;
      uint ktnw = kw;
      uint kbse = kbs;
      uint ktsw = ksw;
      uint kbne = kb;
      uint ktne = indexOfBCnode;
      uint kbsw = neighborZ[ksw];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local distributions
      //!
      real f_W    = (dist.f[dP00])[ke   ];
      real f_E    = (dist.f[dM00])[kw   ];
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
      q = (subgridD.q[dP00])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dM00])[kw  ]=f_E   + c4o9  * (-VeloX);
      q = (subgridD.q[dM00])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dP00])[ke  ]=f_W   + c4o9  * ( VeloX);
      q = (subgridD.q[DIR_0P0])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_0M0])[ks  ]=f_N   + c4o9  * (-VeloY);
      q = (subgridD.q[DIR_0M0])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_0P0])[kn  ]=f_S   + c4o9  * ( VeloY);
      q = (subgridD.q[DIR_00P])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_00M])[kb  ]=f_T   + c4o9  * (-VeloZ);
      q = (subgridD.q[DIR_00M])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_00P])[kt  ]=f_B   + c4o9  * ( VeloZ);
      q = (subgridD.q[DIR_PP0])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_MM0])[ksw ]=f_NE  + c1o9  * (-VeloX - VeloY);
      q = (subgridD.q[DIR_MM0])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_PP0])[kne ]=f_SW  + c1o9  * ( VeloX + VeloY);
      q = (subgridD.q[DIR_PM0])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_MP0])[knw ]=f_SE  + c1o9  * (-VeloX + VeloY);
      q = (subgridD.q[DIR_MP0])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_PM0])[kse ]=f_NW  + c1o9  * ( VeloX - VeloY);
      q = (subgridD.q[DIR_P0P])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_M0M])[kbw ]=f_TE  + c1o9  * (-VeloX - VeloZ);
      q = (subgridD.q[DIR_M0M])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_P0P])[kte ]=f_BW  + c1o9  * ( VeloX + VeloZ);
      q = (subgridD.q[DIR_P0M])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_M0P])[ktw ]=f_BE  + c1o9  * (-VeloX + VeloZ);
      q = (subgridD.q[DIR_M0P])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_P0M])[kbe ]=f_TW  + c1o9  * ( VeloX - VeloZ);
      q = (subgridD.q[DIR_0PP])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_0MM])[kbs ]=f_TN  + c1o9  * (-VeloY - VeloZ);
      q = (subgridD.q[DIR_0MM])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_0PP])[ktn ]=f_BS  + c1o9  * ( VeloY + VeloZ);
      q = (subgridD.q[DIR_0PM])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_0MP])[kts ]=f_BN  + c1o9  * (-VeloY + VeloZ);
      q = (subgridD.q[DIR_0MP])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_0PM])[kbn ]=f_TS  + c1o9  * ( VeloY - VeloZ);
      q = (subgridD.q[DIR_PPP])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_MMM])[kbsw]=f_TNE + c1o36 * (-VeloX - VeloY - VeloZ);
      q = (subgridD.q[DIR_MMM])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_PPP])[ktne]=f_BSW + c1o36 * ( VeloX + VeloY + VeloZ);
      q = (subgridD.q[DIR_PPM])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_MMP])[ktsw]=f_BNE + c1o36 * (-VeloX - VeloY + VeloZ);
      q = (subgridD.q[DIR_MMP])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_PPM])[kbne]=f_TSW + c1o36 * ( VeloX + VeloY - VeloZ);
      q = (subgridD.q[DIR_PMP])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_MPM])[kbnw]=f_TSE + c1o36 * (-VeloX + VeloY - VeloZ);
      q = (subgridD.q[DIR_MPM])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_PMP])[ktse]=f_BNW + c1o36 * ( VeloX - VeloY + VeloZ);
      q = (subgridD.q[DIR_PMM])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_MPP])[ktnw]=f_BSE + c1o36 * (-VeloX + VeloY + VeloZ);
      q = (subgridD.q[DIR_MPP])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[DIR_PMM])[kbse]=f_TNW + c1o36 * ( VeloX - VeloY - VeloZ);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QVelDevCouette27(
    real* vx,
    real* vy,
    real* vz,
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
      real*q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
      unsigned int KQK  = k_Q[k];
      //unsigned int kzero= KQK;
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
     
      ////////////////////////////////////////////////////////////////////////////////
      real f_W    = (D.f[dP00])[ke   ];
      real f_E    = (D.f[dM00])[kw   ];
      real f_S    = (D.f[DIR_0P0])[kn   ];
      real f_N    = (D.f[DIR_0M0])[ks   ];
      real f_B    = (D.f[DIR_00P])[kt   ];
      real f_T    = (D.f[DIR_00M])[kb   ];
      real f_SW   = (D.f[DIR_PP0])[kne  ];
      real f_NE   = (D.f[DIR_MM0])[ksw  ];
      real f_NW   = (D.f[DIR_PM0])[kse  ];
      real f_SE   = (D.f[DIR_MP0])[knw  ];
      real f_BW   = (D.f[DIR_P0P])[kte  ];
      real f_TE   = (D.f[DIR_M0M])[kbw  ];
      real f_TW   = (D.f[DIR_P0M])[kbe  ];
      real f_BE   = (D.f[DIR_M0P])[ktw  ];
      real f_BS   = (D.f[DIR_0PP])[ktn  ];
      real f_TN   = (D.f[DIR_0MM])[kbs  ];
      real f_TS   = (D.f[DIR_0PM])[kbn  ];
      real f_BN   = (D.f[DIR_0MP])[kts  ];
      real f_BSW  = (D.f[DIR_PPP])[ktne ];
      real f_BNE  = (D.f[DIR_MMP])[ktsw ];
      real f_BNW  = (D.f[DIR_PMP])[ktse ];
      real f_BSE  = (D.f[DIR_MPP])[ktnw ];
      real f_TSW  = (D.f[DIR_PPM])[kbne ];
      real f_TNE  = (D.f[DIR_MMM])[kbsw ];
      real f_TNW  = (D.f[DIR_PMM])[kbse ];
      real f_TSE  = (D.f[DIR_MPM])[kbnw ];
	  ////////////////////////////////////////////////////////////////////////////////

	  ////////////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ///////               FlowDirection Y !!!!!!!!!!                                                           ///////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //calculate velocity
	  //real vx1 = ((f_TNE-f_BSW)+(f_BSE-f_TNW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)) + (((f_NE-f_SW)+(f_TE-f_BW))+((f_SE-f_NW)+(f_BE-f_TW))) + (f_E-f_W);
	  real vx2 = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_BNE-f_TSW)+(f_BNW-f_TSE)) + (((f_NE-f_SW)+(f_TN-f_BS))+((f_BN-f_TS)+(f_NW-f_SE))) + (f_N-f_S);
	  //real vx3 = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_TSW-f_BNE)+(f_TSE-f_BNW)) + (((f_TE-f_BW)+(f_TN-f_BS))+((f_TW-f_BE)+(f_TS-f_BN))) + (f_T-f_B);
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //constant
	  real on=c0o1;//c1o2;//one;
	  real ms=-c6o1;
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //2nd order moment
	  real kxxMyyFromfcNEQ = c0o1;//-c3o2 * (f_BW+f_W+f_TW-f_BS-f_S-f_TS-f_BN-f_N-f_TN+f_BE+f_E+f_TE-(vx1*vx1-vx2*vx2));		//all dP00+dM00 minus all DIR_0P0+DIR_0M0 (no combinations of xy left)

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //set distributions
      real q;
      q = q_dirE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[dM00])[kw  ]=f_E   + ms*c2o27  * VeloX;	
      q = q_dirW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[dP00])[ke  ]=f_W   - ms*c2o27  * VeloX;	
      q = q_dirN[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_0M0])[ks  ]=f_N   + ms*c2o27  * VeloY;	
      q = q_dirS[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_0P0])[kn  ]=f_S   - ms*c2o27  * VeloY;	
	  q = q_dirT[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_00M])[kb  ]=f_T   + ms*c2o27  * VeloZ - c3o2*c2o27*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on;
      q = q_dirB[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_00P])[kt  ]=f_B   - ms*c2o27  * VeloZ;
      q = q_dirNE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_MM0])[ksw ]=f_NE  + ms*c1o54  * VeloX + ms*c1o54  * VeloY;
	  q = q_dirSW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_PP0])[kne ]=f_SW  - ms*c1o54  * VeloX - ms*c1o54  * VeloY;
	  q = q_dirSE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_MP0])[knw ]=f_SE  + ms*c1o54  * VeloX - ms*c1o54  * VeloY;
	  q = q_dirNW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_PM0])[kse ]=f_NW  - ms*c1o54  * VeloX + ms*c1o54  * VeloY;
	  q = q_dirTE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_M0M])[kbw ]=f_TE  + ms*c1o54  * VeloX + ms*c1o54  * VeloZ - c3o2*c1o54*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on-c1o12*kxxMyyFromfcNEQ;
	  q = q_dirBW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_P0P])[kte ]=f_BW  - ms*c1o54  * VeloX - ms*c1o54  * VeloZ;
	  q = q_dirBE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_M0P])[ktw ]=f_BE  + ms*c1o54  * VeloX - ms*c1o54  * VeloZ;
	  q = q_dirTW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_P0M])[kbe ]=f_TW  - ms*c1o54  * VeloX + ms*c1o54  * VeloZ - c3o2*c1o54*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on-c1o12*kxxMyyFromfcNEQ;
	  q = q_dirTN[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_0MM])[kbs ]=f_TN  + ms*c1o54  * VeloY + ms*c1o54  * VeloZ + c3o1*c1o54*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on+c1o12*kxxMyyFromfcNEQ;
	  q = q_dirBS[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_0PP])[ktn ]=f_BS  - ms*c1o54  * VeloY - ms*c1o54  * VeloZ;
	  q = q_dirBN[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_0MP])[kts ]=f_BN  + ms*c1o54  * VeloY - ms*c1o54  * VeloZ;
	  q = q_dirTS[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_0PM])[kbn ]=f_TS  - ms*c1o54  * VeloY + ms*c1o54  * VeloZ + c3o1*c1o54*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on+c1o12*kxxMyyFromfcNEQ;
      q = q_dirTNE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_MMM])[kbsw]=f_TNE + ms*c1o216 * VeloX + ms*c1o216 * VeloY + ms*c1o216 * VeloZ + c3o1*c1o216*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on;
      q = q_dirBSW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_PPP])[ktne]=f_BSW - ms*c1o216 * VeloX - ms*c1o216 * VeloY - ms*c1o216 * VeloZ;
      q = q_dirBNE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_MMP])[ktsw]=f_BNE + ms*c1o216 * VeloX + ms*c1o216 * VeloY - ms*c1o216 * VeloZ;
      q = q_dirTSW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_PPM])[kbne]=f_TSW - ms*c1o216 * VeloX - ms*c1o216 * VeloY + ms*c1o216 * VeloZ + c3o1*c1o216*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on;
      q = q_dirTSE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_MPM])[kbnw]=f_TSE + ms*c1o216 * VeloX - ms*c1o216 * VeloY + ms*c1o216 * VeloZ + c3o1*c1o216*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on;
      q = q_dirBNW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_PMP])[ktse]=f_BNW - ms*c1o216 * VeloX + ms*c1o216 * VeloY - ms*c1o216 * VeloZ;
      q = q_dirBSE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_MPP])[ktnw]=f_BSE + ms*c1o216 * VeloX - ms*c1o216 * VeloY - ms*c1o216 * VeloZ;
      q = q_dirTNW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[DIR_PMM])[kbse]=f_TNW - ms*c1o216 * VeloX + ms*c1o216 * VeloY + ms*c1o216 * VeloZ + c3o1*c1o216*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on;
      //q = q_dirE[k];	if (q>=zero && q<=one)	(D.f[dM00])[kw  ]=f_E   + ms*c2over27  * VeloX;	
   //   q = q_dirW[k];	if (q>=zero && q<=one)	(D.f[dP00])[ke  ]=f_W   - ms*c2over27  * VeloX;	
   //   q = q_dirN[k];	if (q>=zero && q<=one)	(D.f[DIR_0M0])[ks  ]=f_N   + ms*c2over27  * VeloY;	
   //   q = q_dirS[k];	if (q>=zero && q<=one)	(D.f[DIR_0P0])[kn  ]=f_S   - ms*c2over27  * VeloY;	
	  //q = q_dirT[k];	if (q>=zero && q<=one)	(D.f[DIR_00M])[kb  ]=f_T   + ms*c2over27  * VeloZ - c1o9*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirB[k];	if (q>=zero && q<=one)	(D.f[DIR_00P])[kt  ]=f_B   - ms*c2over27  * VeloZ;
   //   q = q_dirNE[k];	if (q>=zero && q<=one)	(D.f[DIR_MM0])[ksw ]=f_NE  + ms*c1over54  * VeloX + ms*c1over54  * VeloY;
	  //q = q_dirSW[k];	if (q>=zero && q<=one)	(D.f[DIR_PP0])[kne ]=f_SW  - ms*c1over54  * VeloX - ms*c1over54  * VeloY;
	  //q = q_dirSE[k];	if (q>=zero && q<=one)	(D.f[DIR_MP0])[knw ]=f_SE  + ms*c1over54  * VeloX - ms*c1over54  * VeloY;
	  //q = q_dirNW[k];	if (q>=zero && q<=one)	(D.f[DIR_PM0])[kse ]=f_NW  - ms*c1over54  * VeloX + ms*c1over54  * VeloY;
	  //q = q_dirTE[k];	if (q>=zero && q<=one)	(D.f[DIR_M0M])[kbw ]=f_TE  + ms*c1over54  * VeloX + ms*c1over54  * VeloZ - c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
	  //q = q_dirBW[k];	if (q>=zero && q<=one)	(D.f[DIR_P0P])[kte ]=f_BW  - ms*c1over54  * VeloX - ms*c1over54  * VeloZ;
	  //q = q_dirBE[k];	if (q>=zero && q<=one)	(D.f[DIR_M0P])[ktw ]=f_BE  + ms*c1over54  * VeloX - ms*c1over54  * VeloZ;
	  //q = q_dirTW[k];	if (q>=zero && q<=one)	(D.f[DIR_P0M])[kbe ]=f_TW  - ms*c1over54  * VeloX + ms*c1over54  * VeloZ - c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
	  //q = q_dirTN[k];	if (q>=zero && q<=one)	(D.f[DIR_0MM])[kbs ]=f_TN  + ms*c1over54  * VeloY + ms*c1over54  * VeloZ + c1o2*c1o9*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
	  //q = q_dirBS[k];	if (q>=zero && q<=one)	(D.f[DIR_0PP])[ktn ]=f_BS  - ms*c1over54  * VeloY - ms*c1over54  * VeloZ;
	  //q = q_dirBN[k];	if (q>=zero && q<=one)	(D.f[DIR_0MP])[kts ]=f_BN  + ms*c1over54  * VeloY - ms*c1over54  * VeloZ;
	  //q = q_dirTS[k];	if (q>=zero && q<=one)	(D.f[DIR_0PM])[kbn ]=f_TS  - ms*c1over54  * VeloY + ms*c1over54  * VeloZ + c1o2*c1o9*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirTNE[k];	if (q>=zero && q<=one)	(D.f[DIR_MMM])[kbsw]=f_TNE + ms*c1over216 * VeloX + ms*c1over216 * VeloY + ms*c1over216 * VeloZ + c1o2*c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirBSW[k];	if (q>=zero && q<=one)	(D.f[DIR_PPP])[ktne]=f_BSW - ms*c1over216 * VeloX - ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
   //   q = q_dirBNE[k];	if (q>=zero && q<=one)	(D.f[DIR_MMP])[ktsw]=f_BNE + ms*c1over216 * VeloX + ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
   //   q = q_dirTSW[k];	if (q>=zero && q<=one)	(D.f[DIR_PPM])[kbne]=f_TSW - ms*c1over216 * VeloX - ms*c1over216 * VeloY + ms*c1over216 * VeloZ + c1o2*c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirTSE[k];	if (q>=zero && q<=one)	(D.f[DIR_MPM])[kbnw]=f_TSE + ms*c1over216 * VeloX - ms*c1over216 * VeloY + ms*c1over216 * VeloZ + c1o2*c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirBNW[k];	if (q>=zero && q<=one)	(D.f[DIR_PMP])[ktse]=f_BNW - ms*c1over216 * VeloX + ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
   //   q = q_dirBSE[k];	if (q>=zero && q<=one)	(D.f[DIR_MPP])[ktnw]=f_BSE + ms*c1over216 * VeloX - ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
   //   q = q_dirTNW[k];	if (q>=zero && q<=one)	(D.f[DIR_PMM])[kbse]=f_TNW - ms*c1over216 * VeloX + ms*c1over216 * VeloY + ms*c1over216 * VeloZ + c1o2*c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QVelDev1h27(
    int inx,
    int iny,
    real* vx,
    real* vy,
    real* vz,
    real* DD, 
    int* k_Q, 
    real* QQ,
    unsigned int numberOfBCnodes, 
    real om1,
    real Phi,
    real angularVelocity,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    real* coordX,
    real* coordY,
    real* coordZ,
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
		real VeloX = cosf(Phi)*vx[k] - sinf(Phi)*vy[k];
		real VeloY = sinf(Phi)*vx[k] + cosf(Phi)*vy[k];
		//real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
		////////////////////////////////////////////////////////////////////////////////////
		//Ship
		real coord0X = 281.125f;//7.5f;
		real coord0Y = 388.125f;//7.5f;
		real ux = - angularVelocity * (coordY[k_Q[k]] - coord0Y);
		real uy =   angularVelocity * (coordX[k_Q[k]] - coord0X);
		real VeloXpur=VeloX;
		real VeloYpur=VeloY;
		VeloX-=ux;
		VeloY-=uy;
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
		unsigned int KQK  = k_Q[k];
		//unsigned int kzero= KQK;
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
		//unsigned int nxny = nx*ny;
		//unsigned int kzero= KQK;
		//unsigned int ke   = KQK;
		//unsigned int kw   = KQK + 1;
		//unsigned int kn   = KQK;
		//unsigned int ks   = KQK + nx;
		//unsigned int kt   = KQK;
		//unsigned int kb   = KQK + nxny;
		//unsigned int ksw  = KQK + nx + 1;
		//unsigned int kne  = KQK;
		//unsigned int kse  = KQK + nx;
		//unsigned int knw  = KQK + 1;
		//unsigned int kbw  = KQK + nxny + 1;
		//unsigned int kte  = KQK;
		//unsigned int kbe  = KQK + nxny;
		//unsigned int ktw  = KQK + 1;
		//unsigned int kbs  = KQK + nxny + nx;
		//unsigned int ktn  = KQK;
		//unsigned int kbn  = KQK + nxny;
		//unsigned int kts  = KQK + nx;
		//unsigned int ktse = KQK + nx;
		//unsigned int kbnw = KQK + nxny + 1;
		//unsigned int ktnw = KQK + 1;
		//unsigned int kbse = KQK + nxny + nx;
		//unsigned int ktsw = KQK + nx + 1;
		//unsigned int kbne = KQK + nxny;
		//unsigned int ktne = KQK;
		//unsigned int kbsw = KQK + nxny + nx + 1;
		////////////////////////////////////////////////////////////////////////////////
		//real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
		//	f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

		//f_W    = (D.f[dP00])[ke   ];
		//f_E    = (D.f[dM00])[kw   ];
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
		real /*vx1, vx2,*/ vx3, drho, feq, q, cu_sq;
		//drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
		//	f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
		//	f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]); 

		//vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
		//	((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
		//	(f_E - f_W); 


		//vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
		//	((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
		//	(f_N - f_S); 

		//vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
		//	(-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
		//	(f_T - f_B); 

		//cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

		//////////////////////////////////////////////////////////////////////////
		if (isEvenTimestep==false)
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
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Test
		//(D.f[d000])[k]=c1o10;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//ToDo anders Klammern

		//q = q_dirE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        )-cu_sq); 
		//	(D.f[dM00])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);
		//	//(D.f[dM00])[kw]=zero;
		//}

		//q = q_dirW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
		//	(D.f[dP00])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);
		//	//(D.f[dP00])[ke]=zero;
		//}

		//q = q_dirN[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
		//	(D.f[DIR_0M0])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);
		//	//(D.f[DIR_0M0])[ks]=zero;
		//}

		//q = q_dirS[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
		//	(D.f[DIR_0P0])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);
		//	//(D.f[DIR_0P0])[kn]=zero;
		//}

		//q = q_dirT[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3)-cu_sq); 
		//	(D.f[DIR_00M])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);
		//	//(D.f[DIR_00M])[kb]=one;
		//}

		//q = q_dirB[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
		//	(D.f[DIR_00P])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);
		//	//(D.f[DIR_00P])[kt]=zero;
		//}

		//q = q_dirNE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
		//	(D.f[DIR_MM0])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);
		//	//(D.f[DIR_MM0])[ksw]=zero;
		//}

		//q = q_dirSW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
		//	(D.f[DIR_PP0])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);
		//	//(D.f[DIR_PP0])[kne]=zero;
		//}

		//q = q_dirSE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
		//	(D.f[DIR_MP0])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);
		//	//(D.f[DIR_MP0])[knw]=zero;
		//}

		//q = q_dirNW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
		//	(D.f[DIR_PM0])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);
		//	//(D.f[DIR_PM0])[kse]=zero;
		//}

		//q = q_dirTE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
		//	(D.f[DIR_M0M])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);
		//	//(D.f[DIR_M0M])[kbw]=zero;
		//}

		//q = q_dirBW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
		//	(D.f[DIR_P0P])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);
		//	//(D.f[DIR_P0P])[kte]=zero;
		//}

		//q = q_dirBE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
		//	(D.f[DIR_M0P])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);
		//	//(D.f[DIR_M0P])[ktw]=zero;
		//}

		//q = q_dirTW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
		//	(D.f[DIR_P0M])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);
		//	//(D.f[DIR_P0M])[kbe]=zero;
		//}

		//q = q_dirTN[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
		//	(D.f[DIR_0MM])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);
		//	//(D.f[DIR_0MM])[kbs]=zero;
		//}

		//q = q_dirBS[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
		//	(D.f[DIR_0PP])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);
		//	//(D.f[DIR_0PP])[ktn]=zero;
		//}

		//q = q_dirBN[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
		//	(D.f[DIR_0MP])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);
		//	//(D.f[DIR_0MP])[kts]=zero;
		//}

		//q = q_dirTS[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
		//	(D.f[DIR_0PM])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);
		//	//(D.f[DIR_0PM])[kbn]=zero;
		//}

		//q = q_dirTNE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
		//	(D.f[DIR_MMM])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);
		//	//(D.f[DIR_MMM])[kbsw]=zero;
		//}

		//q = q_dirBSW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
		//	(D.f[DIR_PPP])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);
		//	//(D.f[DIR_PPP])[ktne]=zero;
		//}

		//q = q_dirBNE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
		//	(D.f[DIR_MMP])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);
		//	//(D.f[DIR_MMP])[ktsw]=zero;
		//}

		//q = q_dirTSW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
		//	(D.f[DIR_PPM])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);
		//	//(D.f[DIR_PPM])[kbne]=zero;
		//}

		//q = q_dirTSE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
		//	(D.f[DIR_MPM])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);
		//	//(D.f[DIR_MPM])[kbnw]=zero;
		//}

		//q = q_dirBNW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
		//	(D.f[DIR_PMP])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);
		//	//(D.f[DIR_PMP])[ktse]=zero;
		//}

		//q = q_dirBSE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
		//	(D.f[DIR_MPP])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);
		//	//(D.f[DIR_MPP])[ktnw]=zero;
		//}

		//q = q_dirTNW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
		//	(D.f[DIR_PMM])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);
		//	//(D.f[DIR_PMM])[kbse]=zero;
		//}

		///////// equilibrium BC
		cu_sq=c3o2*(VeloX*VeloX +VeloY*VeloY);
		VeloXpur*=-c1o1;
		VeloYpur*=-c1o1;
		vx3=c0o1;
		drho=c0o1;
		q = q_dirE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c3o1*( VeloXpur        )+c9o2*( VeloX        )*( VeloX        )-cu_sq); 
			(D.f[dM00])[kw]=feq;
			//(D.f[dM00])[kw]=zero;
		}

		q = q_dirW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c3o1*(-VeloXpur        )+c9o2*(-VeloX        )*(-VeloX        )-cu_sq); 
			(D.f[dP00])[ke]=feq;
			//(D.f[dP00])[ke]=zero;
		}

		q = q_dirN[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c3o1*(    VeloYpur     )+c9o2*(     VeloY    )*(     VeloY    )-cu_sq); 
			(D.f[DIR_0M0])[ks]=feq;
			//(D.f[DIR_0M0])[ks]=zero;
		}

		q = q_dirS[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c3o1*(   -VeloYpur     )+c9o2*(    -VeloY    )*(    -VeloY    )-cu_sq); 
			(D.f[DIR_0P0])[kn]=feq;
			//(D.f[DIR_0P0])[kn]=zero;
		}

		q = q_dirT[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq); 
			(D.f[DIR_00M])[kb]=feq;
			//(D.f[DIR_00M])[kb]=one;
		}

		q = q_dirB[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq); 
			(D.f[DIR_00P])[kt]=feq;
			//(D.f[DIR_00P])[kt]=zero;
		}

		q = q_dirNE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*( VeloXpur+VeloYpur    )+c9o2*( VeloX+VeloY    )*( VeloX+VeloY    )-cu_sq); 
			(D.f[DIR_MM0])[ksw]=feq;
			//(D.f[DIR_MM0])[ksw]=zero;
		}

		q = q_dirSW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(-VeloXpur-VeloYpur    )+c9o2*(-VeloX-VeloY    )*(-VeloX-VeloY    )-cu_sq); 
			(D.f[DIR_PP0])[kne]=feq;
			//(D.f[DIR_PP0])[kne]=zero;
		}

		q = q_dirSE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*( VeloXpur-VeloYpur    )+c9o2*( VeloX-VeloY    )*( VeloX-VeloY    )-cu_sq); 
			(D.f[DIR_MP0])[knw]=feq;
			//(D.f[DIR_MP0])[knw]=zero;
		}

		q = q_dirNW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(-VeloXpur+VeloYpur    )+c9o2*(-VeloX+VeloY    )*(-VeloX+VeloY    )-cu_sq); 
			(D.f[DIR_PM0])[kse]=feq;
			//(D.f[DIR_PM0])[kse]=zero;
		}

		q = q_dirTE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*( VeloXpur    +vx3)+c9o2*( VeloX    +vx3)*( VeloX    +vx3)-cu_sq); 
			(D.f[DIR_M0M])[kbw]=feq;
			//(D.f[DIR_M0M])[kbw]=zero;
		}

		q = q_dirBW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(-VeloXpur    -vx3)+c9o2*(-VeloX    -vx3)*(-VeloX    -vx3)-cu_sq); 
			(D.f[DIR_P0P])[kte]=feq;
			//(D.f[DIR_P0P])[kte]=zero;
		}

		q = q_dirBE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*( VeloXpur    -vx3)+c9o2*( VeloX    -vx3)*( VeloX    -vx3)-cu_sq); 
			(D.f[DIR_M0P])[ktw]=feq;
			//(D.f[DIR_M0P])[ktw]=zero;
		}

		q = q_dirTW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(-VeloXpur    +vx3)+c9o2*(-VeloX    +vx3)*(-VeloX    +vx3)-cu_sq); 
			(D.f[DIR_P0M])[kbe]=feq;
			//(D.f[DIR_P0M])[kbe]=zero;
		}

		q = q_dirTN[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(     VeloYpur+vx3)+c9o2*(     VeloY+vx3)*(     VeloY+vx3)-cu_sq); 
			(D.f[DIR_0MM])[kbs]=feq;
			//(D.f[DIR_0MM])[kbs]=zero;
		}

		q = q_dirBS[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(    -VeloYpur-vx3)+c9o2*(    -VeloY-vx3)*(    -VeloY-vx3)-cu_sq); 
			(D.f[DIR_0PP])[ktn]=feq;
			//(D.f[DIR_0PP])[ktn]=zero;
		}

		q = q_dirBN[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(     VeloYpur-vx3)+c9o2*(     VeloY-vx3)*(     VeloY-vx3)-cu_sq); 
			(D.f[DIR_0MP])[kts]=feq;
			//(D.f[DIR_0MP])[kts]=zero;
		}

		q = q_dirTS[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(    -VeloYpur+vx3)+c9o2*(    -VeloY+vx3)*(    -VeloY+vx3)-cu_sq); 
			(D.f[DIR_0PM])[kbn]=feq;
			//(D.f[DIR_0PM])[kbn]=zero;
		}

		q = q_dirTNE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*( VeloXpur+VeloYpur+vx3)+c9o2*( VeloX+VeloY+vx3)*( VeloX+VeloY+vx3)-cu_sq); 
			(D.f[DIR_MMM])[kbsw]=feq;
			//(D.f[DIR_MMM])[kbsw]=zero;
		}

		q = q_dirBSW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*(-VeloXpur-VeloYpur-vx3)+c9o2*(-VeloX-VeloY-vx3)*(-VeloX-VeloY-vx3)-cu_sq); 
			(D.f[DIR_PPP])[ktne]=feq;
			//(D.f[DIR_PPP])[ktne]=zero;
		}

		q = q_dirBNE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*( VeloXpur+VeloYpur-vx3)+c9o2*( VeloX+VeloY-vx3)*( VeloX+VeloY-vx3)-cu_sq); 
			(D.f[DIR_MMP])[ktsw]=feq;
			//(D.f[DIR_MMP])[ktsw]=zero;
		}

		q = q_dirTSW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*(-VeloXpur-VeloYpur+vx3)+c9o2*(-VeloX-VeloY+vx3)*(-VeloX-VeloY+vx3)-cu_sq); 
			(D.f[DIR_PPM])[kbne]=feq;
			//(D.f[DIR_PPM])[kbne]=zero;
		}

		q = q_dirTSE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*( VeloXpur-VeloYpur+vx3)+c9o2*( VeloX-VeloY+vx3)*( VeloX-VeloY+vx3)-cu_sq); 
			(D.f[DIR_MPM])[kbnw]=feq;
			//(D.f[DIR_MPM])[kbnw]=zero;
		}

		q = q_dirBNW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*(-VeloXpur+VeloYpur-vx3)+c9o2*(-VeloX+VeloY-vx3)*(-VeloX+VeloY-vx3)-cu_sq); 
			(D.f[DIR_PMP])[ktse]=feq;
			//(D.f[DIR_PMP])[ktse]=zero;
		}

		q = q_dirBSE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*( VeloXpur-VeloYpur-vx3)+c9o2*( VeloX-VeloY-vx3)*( VeloX-VeloY-vx3)-cu_sq); 
			(D.f[DIR_MPP])[ktnw]=feq;
			//(D.f[DIR_MPP])[ktnw]=zero;
		}

		q = q_dirTNW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*(-VeloXpur+VeloYpur+vx3)+c9o2*(-VeloX+VeloY+vx3)*(-VeloX+VeloY+vx3)-cu_sq); 
			(D.f[DIR_PMM])[kbse]=feq;
			//(D.f[DIR_PMM])[kbse]=zero;
		}
	
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QVelDeviceComp27(
    real* velocityX,
    real* velocityY,
    real* velocityZ,
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
   //! The velocity boundary condition is executed in the following steps
   //!
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

   //////////////////////////////////////////////////////////////////////////
   //! - Run for all indices in size of boundary condition (numberOfBCnodes)
   //!
   if(nodeIndex < numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist;
      getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local velocities
      //!
      real VeloX = velocityX[nodeIndex];
      real VeloY = velocityY[nodeIndex];
      real VeloZ = velocityZ[nodeIndex];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local subgrid distances (q's)
      //!
      SubgridDistances27 subgridD;
      getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);
      
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int indexOfBCnode  = subgridDistanceIndices[nodeIndex];
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
      real f_W    = (dist.f[dP00])[ke   ];
      real f_E    = (dist.f[dM00])[kw   ];
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
                  f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[d000])[kzero]); 

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
      //!
      real feq, q, velocityLB, velocityBC;
      q = (subgridD.q[dP00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1) // only update distribution for q between zero and one
      {
         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloX;
         (dist.f[dM00])[kw] = getInterpolatedDistributionForVeloBC(q, f_E, f_W, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[dM00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloX;
         (dist.f[dP00])[ke] = getInterpolatedDistributionForVeloBC(q, f_W, f_E, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[DIR_0P0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloY;
         (dist.f[DIR_0M0])[ks] = getInterpolatedDistributionForVeloBC(q, f_N, f_S, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[DIR_0M0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloY;
         (dist.f[DIR_0P0])[kn] = getInterpolatedDistributionForVeloBC(q, f_S, f_N, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[DIR_00P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloZ;
         (dist.f[DIR_00M])[kb] = getInterpolatedDistributionForVeloBC(q, f_T, f_B, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[DIR_00M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloZ;
         (dist.f[DIR_00P])[kt] = getInterpolatedDistributionForVeloBC(q, f_B, f_T, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[DIR_PP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloY;
         (dist.f[DIR_MM0])[ksw] = getInterpolatedDistributionForVeloBC(q, f_NE, f_SW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_MM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloY;
         (dist.f[DIR_PP0])[kne] = getInterpolatedDistributionForVeloBC(q, f_SW, f_NE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_PM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloY;
         (dist.f[DIR_MP0])[knw] = getInterpolatedDistributionForVeloBC(q, f_SE, f_NW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_MP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloY;
         (dist.f[DIR_PM0])[kse] = getInterpolatedDistributionForVeloBC(q, f_NW, f_SE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_P0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloZ;
         (dist.f[DIR_M0M])[kbw] = getInterpolatedDistributionForVeloBC(q, f_TE, f_BW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_M0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloZ;
         (dist.f[DIR_P0P])[kte] = getInterpolatedDistributionForVeloBC(q, f_BW, f_TE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_P0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloZ;
         (dist.f[DIR_M0P])[ktw] = getInterpolatedDistributionForVeloBC(q, f_BE, f_TW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_M0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloZ;
         (dist.f[DIR_P0M])[kbe] = getInterpolatedDistributionForVeloBC(q, f_TW, f_BE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_0PP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY + VeloZ;
         (dist.f[DIR_0MM])[kbs] = getInterpolatedDistributionForVeloBC(q, f_TN, f_BS, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_0MM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY - VeloZ;
         (dist.f[DIR_0PP])[ktn] = getInterpolatedDistributionForVeloBC(q, f_BS, f_TN, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_0PM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY - VeloZ;
         (dist.f[DIR_0MP])[kts] = getInterpolatedDistributionForVeloBC(q, f_BN, f_TS, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_0MP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY + VeloZ;
         (dist.f[DIR_0PM])[kbn] = getInterpolatedDistributionForVeloBC(q, f_TS, f_BN, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[DIR_PPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY + VeloZ;
         (dist.f[DIR_MMM])[kbsw] = getInterpolatedDistributionForVeloBC(q, f_TNE, f_BSW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_MMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY - VeloZ;
         (dist.f[DIR_PPP])[ktne] = getInterpolatedDistributionForVeloBC(q, f_BSW, f_TNE, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_PPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY - VeloZ;
         (dist.f[DIR_MMP])[ktsw] = getInterpolatedDistributionForVeloBC(q, f_BNE, f_TSW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_MMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY + VeloZ;
         (dist.f[DIR_PPM])[kbne] = getInterpolatedDistributionForVeloBC(q, f_TSW, f_BNE, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_PMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY + VeloZ;
         (dist.f[DIR_MPM])[kbnw] = getInterpolatedDistributionForVeloBC(q, f_TSE, f_BNW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_MPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY - VeloZ;
         (dist.f[DIR_PMP])[ktse] = getInterpolatedDistributionForVeloBC(q, f_BNW, f_TSE, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_PMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY - VeloZ;
         (dist.f[DIR_MPP])[ktnw] = getInterpolatedDistributionForVeloBC(q, f_BSE, f_TNW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[DIR_MPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY + VeloZ;
         (dist.f[DIR_PMM])[kbse] = getInterpolatedDistributionForVeloBC(q, f_TNW, f_BSE, feq, omega, velocityBC, c1o216);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QVelDevice27(
    int inx,
    int iny,
    real* vx,
    real* vy,
    real* vz,
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
      real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
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
      //unsigned int nxny = nx*ny;
      //unsigned int kzero= KQK;
      //unsigned int ke   = KQK;
      //unsigned int kw   = KQK + 1;
      //unsigned int kn   = KQK;
      //unsigned int ks   = KQK + nx;
      //unsigned int kt   = KQK;
      //unsigned int kb   = KQK + nxny;
      //unsigned int ksw  = KQK + nx + 1;
      //unsigned int kne  = KQK;
      //unsigned int kse  = KQK + nx;
      //unsigned int knw  = KQK + 1;
      //unsigned int kbw  = KQK + nxny + 1;
      //unsigned int kte  = KQK;
      //unsigned int kbe  = KQK + nxny;
      //unsigned int ktw  = KQK + 1;
      //unsigned int kbs  = KQK + nxny + nx;
      //unsigned int ktn  = KQK;
      //unsigned int kbn  = KQK + nxny;
      //unsigned int kts  = KQK + nx;
      //unsigned int ktse = KQK + nx;
      //unsigned int kbnw = KQK + nxny + 1;
      //unsigned int ktnw = KQK + 1;
      //unsigned int kbse = KQK + nxny + nx;
      //unsigned int ktsw = KQK + nx + 1;
      //unsigned int kbne = KQK + nxny;
      //unsigned int ktne = KQK;
      //unsigned int kbsw = KQK + nxny + nx + 1;
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

      vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W); 
         

      vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S); 

      vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[d000])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //ToDo anders Klammern

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        )-cu_sq); 
         (D.f[dM00])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q);
         //(D.f[dM00])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        )-cu_sq); 
         (D.f[dP00])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q);
         //(D.f[dP00])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    )-cu_sq); 
         (D.f[DIR_0M0])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q);
         //(D.f[DIR_0M0])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    )-cu_sq); 
         (D.f[DIR_0P0])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q);
         //(D.f[DIR_0P0])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3)-cu_sq); 
         (D.f[DIR_00M])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q);
         //(D.f[DIR_00M])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3)-cu_sq); 
         (D.f[DIR_00P])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q);
         //(D.f[DIR_00P])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         (D.f[DIR_MM0])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q);
         //(D.f[DIR_MM0])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         (D.f[DIR_PP0])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
         //(D.f[DIR_PP0])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         (D.f[DIR_MP0])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
         //(D.f[DIR_MP0])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         (D.f[DIR_PM0])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
         //(D.f[DIR_PM0])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         (D.f[DIR_M0M])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
         //(D.f[DIR_M0M])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         (D.f[DIR_P0P])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
         //(D.f[DIR_P0P])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         (D.f[DIR_M0P])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
         //(D.f[DIR_M0P])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         (D.f[DIR_P0M])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
         //(D.f[DIR_P0M])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         (D.f[DIR_0MM])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
         //(D.f[DIR_0MM])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         (D.f[DIR_0PP])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q);
         //(D.f[DIR_0PP])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         (D.f[DIR_0MP])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
         //(D.f[DIR_0MP])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         (D.f[DIR_0PM])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q);
         //(D.f[DIR_0PM])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         (D.f[DIR_MMM])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
         //(D.f[DIR_MMM])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         (D.f[DIR_PPP])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
         //(D.f[DIR_PPP])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         (D.f[DIR_MMP])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
         //(D.f[DIR_MMP])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         (D.f[DIR_PPM])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
         //(D.f[DIR_PPM])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         (D.f[DIR_MPM])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
         //(D.f[DIR_MPM])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         (D.f[DIR_PMP])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
         //(D.f[DIR_PMP])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         (D.f[DIR_MPP])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
         //(D.f[DIR_MPP])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         (D.f[DIR_PMM])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);
         //(D.f[DIR_PMM])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
__global__ void PropellerBC(
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    real* rho,
    real* ux,
    real* uy,
    real* uz,
    int* k_Q, 
    unsigned int size_Prop,
    unsigned long long numberOfLBnodes,
    unsigned int* bcMatD,
    real* DD,
    bool EvenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Prop)
   {
    ////////////////////////////////////////////////////////////////////////////////
        Distributions27 D;
        if (EvenOrOdd==true)
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
			D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
			D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
			D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
			D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
			D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
			D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
			D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
			D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
        }
        //////////////////////////////////////////////////////////////////////////
		unsigned int KQK = k_Q[k];
		unsigned int BC  = bcMatD[KQK];
		if( (BC != GEO_SOLID) && (BC != GEO_VOID))
		{		
		//////////////////////////////////////////////////////////////////////////
        real  vx1 = ux[k];
        real  vx2 = uy[k];
        real  vx3 = uz[k];
        //real  vx1 = -c1o100;
        //real  vx2 = zero;
        //real  vx3 = zero;
        //////////////////////////////////////////////////////////////////////////
        //index
        //////////////////////////////////////////////////////////////////////////
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
        //////////////////////////////////////////////////////////////////////////
		real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
		f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW, f_ZERO;

		f_ZERO= (D.f[d000])[kzero];
		f_E   = (D.f[dP00])[ke   ];
		f_W   = (D.f[dM00])[kw   ];
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
		f_BSW = (D.f[DIR_MMM])[kbsw ];
		f_BNE = (D.f[DIR_PPM])[kbne ];
		f_TSW = (D.f[DIR_MMP])[ktsw ];
		f_TSE = (D.f[DIR_PMP])[ktse ];
		f_BNW = (D.f[DIR_MPM])[kbnw ];
		f_BSE = (D.f[DIR_PMM])[kbse ];
		f_TNW = (D.f[DIR_MPP])[ktnw ];
		//f_W    = (D.f[dP00])[ke   ];
		//f_E    = (D.f[dM00])[kw   ];
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
		//f_TNE  = (D.f[DIR_MMM])[kbsw ];
		//f_TSW  = (D.f[DIR_PPM])[kbne ];
		//f_BNE  = (D.f[DIR_MMP])[ktsw ];
		//f_BNW  = (D.f[DIR_PMP])[ktse ];
		//f_TSE  = (D.f[DIR_MPM])[kbnw ];
		//f_TNW  = (D.f[DIR_PMM])[kbse ];
		//f_BSE  = (D.f[DIR_MPP])[ktnw ];
		//////////////////////////////////////////////////////////////////////////////////
		real vxo1, vxo2, vxo3, drho;
		drho   =  /*zero;*/f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				  f_T + f_B + f_N + f_S + f_E + f_W + f_ZERO; 

		vxo1   =   (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
					((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
					(f_E - f_W) )/ (c1o1 + drho); 
        

		vxo2   =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
					((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
					(f_N - f_S) )/ (c1o1 + drho); 

		vxo3   =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
		 			(-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
					(f_T - f_B) )/ (c1o1 + drho); 

		real cusq=c3o2*(vxo1*vxo1+vxo2*vxo2+vxo3*vxo3);
		//vx1 = vx1 * two - vxo1;
		//vx2 = vx2 * two - vxo2;
		//vx3 = vx3 * two - vxo3;
		real cusq2=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         //f_ZERO = ((one+drho) * (   c8over27 *(one+(-cusq2)))) - c8over27;
         //f_E    = ((one+drho) * (   c2over27 *(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq2))) - c2over27 ;
         //f_W    = ((one+drho) * (   c2over27 *(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq2))) - c2over27 ;
         //f_N    = ((one+drho) * (   c2over27 *(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq2))) - c2over27 ;
         //f_S    = ((one+drho) * (   c2over27 *(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq2))) - c2over27 ;
         //f_T    = ((one+drho) * (   c2over27 *(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq2))) - c2over27 ;
         //f_B    = ((one+drho) * (   c2over27 *(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq2))) - c2over27 ;
         //f_NE   = ((one+drho) * (   c1over54 *(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq2))) - c1over54 ;
         //f_SW   = ((one+drho) * (   c1over54 *(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq2))) - c1over54 ;
         //f_SE   = ((one+drho) * (   c1over54 *(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq2))) - c1over54 ;
         //f_NW   = ((one+drho) * (   c1over54 *(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq2))) - c1over54 ;
         //f_TE   = ((one+drho) * (   c1over54 *(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq2))) - c1over54 ;
         //f_BW   = ((one+drho) * (   c1over54 *(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq2))) - c1over54 ;
         //f_BE   = ((one+drho) * (   c1over54 *(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq2))) - c1over54 ;
         //f_TW   = ((one+drho) * (   c1over54 *(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq2))) - c1over54 ;
         //f_TN   = ((one+drho) * (   c1over54 *(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq2))) - c1over54 ;
         //f_BS   = ((one+drho) * (   c1over54 *(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq2))) - c1over54 ;
         //f_BN   = ((one+drho) * (   c1over54 *(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq2))) - c1over54 ;
         //f_TS   = ((one+drho) * (   c1over54 *(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq2))) - c1over54 ;
         //f_TNE  = ((one+drho) * (   c1over216*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq2))) - c1over216;
         //f_BSW  = ((one+drho) * (   c1over216*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq2))) - c1over216;
         //f_BNE  = ((one+drho) * (   c1over216*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq2))) - c1over216;
         //f_TSW  = ((one+drho) * (   c1over216*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq2))) - c1over216;
         //f_TSE  = ((one+drho) * (   c1over216*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq2))) - c1over216;
         //f_BNW  = ((one+drho) * (   c1over216*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq2))) - c1over216;
         //f_BSE  = ((one+drho) * (   c1over216*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq2))) - c1over216;
         //f_TNW  = ((one+drho) * (   c1over216*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq2))) - c1over216;
         f_ZERO = f_ZERO + ((c1o1+drho) * (-  c8o27* (-cusq)																   +   c8o27* (-cusq2)));
         f_E    = f_E    + ((c1o1+drho) * (-  c2o27* (c3o1*( vxo1          )+c9o2*( vxo1          )*( vxo1          )-cusq) +   c2o27* (c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq2)));
         f_W    = f_W    + ((c1o1+drho) * (-  c2o27* (c3o1*(-vxo1          )+c9o2*(-vxo1          )*(-vxo1          )-cusq) +   c2o27* (c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq2)));
         f_N    = f_N    + ((c1o1+drho) * (-  c2o27* (c3o1*(      vxo2     )+c9o2*(      vxo2     )*(      vxo2     )-cusq) +   c2o27* (c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq2)));
         f_S    = f_S    + ((c1o1+drho) * (-  c2o27* (c3o1*(     -vxo2     )+c9o2*(     -vxo2     )*(     -vxo2     )-cusq) +   c2o27* (c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq2)));
         f_T    = f_T    + ((c1o1+drho) * (-  c2o27* (c3o1*(           vxo3)+c9o2*(           vxo3)*(           vxo3)-cusq) +   c2o27* (c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq2)));
         f_B    = f_B    + ((c1o1+drho) * (-  c2o27* (c3o1*(          -vxo3)+c9o2*(          -vxo3)*(          -vxo3)-cusq) +   c2o27* (c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq2)));
         f_NE   = f_NE   + ((c1o1+drho) * (-  c1o54* (c3o1*( vxo1+vxo2     )+c9o2*( vxo1+vxo2     )*( vxo1+vxo2     )-cusq) +   c1o54* (c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq2)));
         f_SW   = f_SW   + ((c1o1+drho) * (-  c1o54* (c3o1*(-vxo1-vxo2     )+c9o2*(-vxo1-vxo2     )*(-vxo1-vxo2     )-cusq) +   c1o54* (c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq2)));
         f_SE   = f_SE   + ((c1o1+drho) * (-  c1o54* (c3o1*( vxo1-vxo2     )+c9o2*( vxo1-vxo2     )*( vxo1-vxo2     )-cusq) +   c1o54* (c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq2)));
         f_NW   = f_NW   + ((c1o1+drho) * (-  c1o54* (c3o1*(-vxo1+vxo2     )+c9o2*(-vxo1+vxo2     )*(-vxo1+vxo2     )-cusq) +   c1o54* (c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq2)));
         f_TE   = f_TE   + ((c1o1+drho) * (-  c1o54* (c3o1*( vxo1     +vxo3)+c9o2*( vxo1     +vxo3)*( vxo1     +vxo3)-cusq) +   c1o54* (c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq2)));
         f_BW   = f_BW   + ((c1o1+drho) * (-  c1o54* (c3o1*(-vxo1     -vxo3)+c9o2*(-vxo1     -vxo3)*(-vxo1     -vxo3)-cusq) +   c1o54* (c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq2)));
         f_BE   = f_BE   + ((c1o1+drho) * (-  c1o54* (c3o1*( vxo1     -vxo3)+c9o2*( vxo1     -vxo3)*( vxo1     -vxo3)-cusq) +   c1o54* (c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq2)));
         f_TW   = f_TW   + ((c1o1+drho) * (-  c1o54* (c3o1*(-vxo1     +vxo3)+c9o2*(-vxo1     +vxo3)*(-vxo1     +vxo3)-cusq) +   c1o54* (c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq2)));
         f_TN   = f_TN   + ((c1o1+drho) * (-  c1o54* (c3o1*(      vxo2+vxo3)+c9o2*(      vxo2+vxo3)*(      vxo2+vxo3)-cusq) +   c1o54* (c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq2)));
         f_BS   = f_BS   + ((c1o1+drho) * (-  c1o54* (c3o1*(     -vxo2-vxo3)+c9o2*(     -vxo2-vxo3)*(     -vxo2-vxo3)-cusq) +   c1o54* (c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq2)));
         f_BN   = f_BN   + ((c1o1+drho) * (-  c1o54* (c3o1*(      vxo2-vxo3)+c9o2*(      vxo2-vxo3)*(      vxo2-vxo3)-cusq) +   c1o54* (c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq2)));
         f_TS   = f_TS   + ((c1o1+drho) * (-  c1o54* (c3o1*(     -vxo2+vxo3)+c9o2*(     -vxo2+vxo3)*(     -vxo2+vxo3)-cusq) +   c1o54* (c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq2)));
         f_TNE  = f_TNE  + ((c1o1+drho) * (-  c1o216*(c3o1*( vxo1+vxo2+vxo3)+c9o2*( vxo1+vxo2+vxo3)*( vxo1+vxo2+vxo3)-cusq) +   c1o216*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq2)));
         f_BSW  = f_BSW  + ((c1o1+drho) * (-  c1o216*(c3o1*(-vxo1-vxo2-vxo3)+c9o2*(-vxo1-vxo2-vxo3)*(-vxo1-vxo2-vxo3)-cusq) +   c1o216*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq2)));
         f_BNE  = f_BNE  + ((c1o1+drho) * (-  c1o216*(c3o1*( vxo1+vxo2-vxo3)+c9o2*( vxo1+vxo2-vxo3)*( vxo1+vxo2-vxo3)-cusq) +   c1o216*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq2)));
         f_TSW  = f_TSW  + ((c1o1+drho) * (-  c1o216*(c3o1*(-vxo1-vxo2+vxo3)+c9o2*(-vxo1-vxo2+vxo3)*(-vxo1-vxo2+vxo3)-cusq) +   c1o216*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq2)));
         f_TSE  = f_TSE  + ((c1o1+drho) * (-  c1o216*(c3o1*( vxo1-vxo2+vxo3)+c9o2*( vxo1-vxo2+vxo3)*( vxo1-vxo2+vxo3)-cusq) +   c1o216*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq2)));
         f_BNW  = f_BNW  + ((c1o1+drho) * (-  c1o216*(c3o1*(-vxo1+vxo2-vxo3)+c9o2*(-vxo1+vxo2-vxo3)*(-vxo1+vxo2-vxo3)-cusq) +   c1o216*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq2)));
         f_BSE  = f_BSE  + ((c1o1+drho) * (-  c1o216*(c3o1*( vxo1-vxo2-vxo3)+c9o2*( vxo1-vxo2-vxo3)*( vxo1-vxo2-vxo3)-cusq) +   c1o216*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq2)));
         f_TNW  = f_TNW  + ((c1o1+drho) * (-  c1o216*(c3o1*(-vxo1+vxo2+vxo3)+c9o2*(-vxo1+vxo2+vxo3)*(-vxo1+vxo2+vxo3)-cusq) +   c1o216*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq2)));

		(D.f[d000])[kzero] =  f_ZERO;
        (D.f[dP00])[ke   ] =  f_E   ;	// f_W   ;//    	
        (D.f[dM00])[kw   ] =  f_W   ;	// f_E   ;//    	
        (D.f[DIR_0P0])[kn   ] =  f_N   ;	// f_S   ;//    	
        (D.f[DIR_0M0])[ks   ] =  f_S   ;	// f_N   ;//    	
        (D.f[DIR_00P])[kt   ] =  f_T   ;	// f_B   ;//    	
        (D.f[DIR_00M])[kb   ] =  f_B   ;	// f_T   ;//    	
        (D.f[DIR_PP0])[kne  ] =  f_NE  ;	// f_SW  ;//    	
        (D.f[DIR_MM0])[ksw  ] =  f_SW  ;	// f_NE  ;//    	
        (D.f[DIR_PM0])[kse  ] =  f_SE  ;	// f_NW  ;//    	
        (D.f[DIR_MP0])[knw  ] =  f_NW  ;	// f_SE  ;//    	
        (D.f[DIR_P0P])[kte  ] =  f_TE  ;	// f_BW  ;//    	
        (D.f[DIR_M0M])[kbw  ] =  f_BW  ;	// f_TE  ;//    	
        (D.f[DIR_P0M])[kbe  ] =  f_BE  ;	// f_TW  ;//    	
        (D.f[DIR_M0P])[ktw  ] =  f_TW  ;	// f_BE  ;//    	
        (D.f[DIR_0PP])[ktn  ] =  f_TN  ;	// f_BS  ;//    	
        (D.f[DIR_0MM])[kbs  ] =  f_BS  ;	// f_TN  ;//    	
        (D.f[DIR_0PM])[kbn  ] =  f_BN  ;	// f_TS  ;//    	
        (D.f[DIR_0MP])[kts  ] =  f_TS  ;	// f_BN  ;//    	
        (D.f[DIR_PPP])[ktne ] =  f_TNE ;	// f_BSW ;//    	
        (D.f[DIR_MMM])[kbsw ] =  f_BSW ;	// f_BNE ;//    	
        (D.f[DIR_PPM])[kbne ] =  f_BNE ;	// f_BNW ;//    	
        (D.f[DIR_MMP])[ktsw ] =  f_TSW ;	// f_BSE ;//    	
        (D.f[DIR_PMP])[ktse ] =  f_TSE ;	// f_TSW ;//    	
        (D.f[DIR_MPM])[kbnw ] =  f_BNW ;	// f_TNE ;//    	
        (D.f[DIR_PMM])[kbse ] =  f_BSE ;	// f_TNW ;//    	
        (D.f[DIR_MPP])[ktnw ] =  f_TNW ;	// f_TSE ;//    	

		//////////////////////////////////////////////////////////////////////////
        ////(D.f[d000])[kzero] =   c8over27* (drho-cu_sq);
        //(D.f[dP00])[ke   ] =   three*c2over27* ( vx1        );		//six
        //(D.f[dM00])[kw   ] =   three*c2over27* (-vx1        );		//six
        //(D.f[DIR_0P0])[kn   ] =   three*c2over27* (     vx2    );		//six
        //(D.f[DIR_0M0])[ks   ] =   three*c2over27* (    -vx2    );		//six
        //(D.f[DIR_00P])[kt   ] =   three*c2over27* (         vx3);		//six
        //(D.f[DIR_00M])[kb   ] =   three*c2over27* (        -vx3);		//six
        //(D.f[DIR_PP0])[kne  ] =   three*c1over54* ( vx1+vx2    );		//six
        //(D.f[DIR_MM0])[ksw  ] =   three*c1over54* (-vx1-vx2    );		//six
        //(D.f[DIR_PM0])[kse  ] =   three*c1over54* ( vx1-vx2    );		//six
        //(D.f[DIR_MP0])[knw  ] =   three*c1over54* (-vx1+vx2    );		//six
        //(D.f[DIR_P0P])[kte  ] =   three*c1over54* ( vx1    +vx3);		//six
        //(D.f[DIR_M0M])[kbw  ] =   three*c1over54* (-vx1    -vx3);		//six
        //(D.f[DIR_P0M])[kbe  ] =   three*c1over54* ( vx1    -vx3);		//six
        //(D.f[DIR_M0P])[ktw  ] =   three*c1over54* (-vx1    +vx3);		//six
        //(D.f[DIR_0PP])[ktn  ] =   three*c1over54* (     vx2+vx3);		//six
        //(D.f[DIR_0MM])[kbs  ] =   three*c1over54* (    -vx2-vx3);		//six
        //(D.f[DIR_0PM])[kbn  ] =   three*c1over54* (     vx2-vx3);		//six
        //(D.f[DIR_0MP])[kts  ] =   three*c1over54* (    -vx2+vx3);		//six
        //(D.f[DIR_PPP])[ktne ] =   three*c1over216*( vx1+vx2+vx3);		//six
        //(D.f[DIR_MMM])[kbsw ] =   three*c1over216*(-vx1-vx2-vx3);		//six
        //(D.f[DIR_PPM])[kbne ] =   three*c1over216*( vx1+vx2-vx3);		//six
        //(D.f[DIR_MMP])[ktsw ] =   three*c1over216*(-vx1-vx2+vx3);		//six
        //(D.f[DIR_PMP])[ktse ] =   three*c1over216*( vx1-vx2+vx3);		//six
        //(D.f[DIR_MPM])[kbnw ] =   three*c1over216*(-vx1+vx2-vx3);		//six
        //(D.f[DIR_PMM])[kbse ] =   three*c1over216*( vx1-vx2-vx3);		//six
        //(D.f[DIR_MPP])[ktnw ] =   three*c1over216*(-vx1+vx2+vx3);		//six
        //(D.f[d000])[kzero] =   c8over27* (drho-cu_sq);
        //(D.f[dP00])[ke   ] =   c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
        //(D.f[dM00])[kw   ] =   c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
        //(D.f[DIR_0P0])[kn   ] =   c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
        //(D.f[DIR_0M0])[ks   ] =   c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
        //(D.f[DIR_00P])[kt   ] =   c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
        //(D.f[DIR_00M])[kb   ] =   c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
        //(D.f[DIR_PP0])[kne  ] =   c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
        //(D.f[DIR_MM0])[ksw  ] =   c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
        //(D.f[DIR_PM0])[kse  ] =   c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
        //(D.f[DIR_MP0])[knw  ] =   c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
        //(D.f[DIR_P0P])[kte  ] =   c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
        //(D.f[DIR_M0M])[kbw  ] =   c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
        //(D.f[DIR_P0M])[kbe  ] =   c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
        //(D.f[DIR_M0P])[ktw  ] =   c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
        //(D.f[DIR_0PP])[ktn  ] =   c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
        //(D.f[DIR_0MM])[kbs  ] =   c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
        //(D.f[DIR_0PM])[kbn  ] =   c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
        //(D.f[DIR_0MP])[kts  ] =   c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
        //(D.f[DIR_PPP])[ktne ] =   c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
        //(D.f[DIR_MMM])[kbsw ] =   c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
        //(D.f[DIR_PPM])[kbne ] =   c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
        //(D.f[DIR_MMP])[ktsw ] =   c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
        //(D.f[DIR_PMP])[ktse ] =   c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
        //(D.f[DIR_MPM])[kbnw ] =   c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
        //(D.f[DIR_PMM])[kbse ] =   c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
        //(D.f[DIR_MPP])[ktnw ] =   c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
		}
    }
}
//////////////////////////////////////////////////////////////////////////


