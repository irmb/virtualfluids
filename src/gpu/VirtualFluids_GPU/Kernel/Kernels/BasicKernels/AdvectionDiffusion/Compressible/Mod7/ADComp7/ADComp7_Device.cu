#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

#include "math.h"

__global__ void LB_Kernel_AD_Comp_7(real diffusivity,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	real* DD7,
	int size_Mat,
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

	if (k<size_Mat)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if ((BC != GEO_SOLID) && (BC != GEO_VOID))
		{
			Distributions27 D;
			if (EvenOrOdd == true)
			{
				D.f[DIR_P00] = &DDStart[DIR_P00 * size_Mat];
				D.f[DIR_M00] = &DDStart[DIR_M00 * size_Mat];
				D.f[DIR_0P0] = &DDStart[DIR_0P0 * size_Mat];
				D.f[DIR_0M0] = &DDStart[DIR_0M0 * size_Mat];
				D.f[DIR_00P] = &DDStart[DIR_00P * size_Mat];
				D.f[DIR_00M] = &DDStart[DIR_00M * size_Mat];
				D.f[DIR_PP0] = &DDStart[DIR_PP0 * size_Mat];
				D.f[DIR_MM0] = &DDStart[DIR_MM0 * size_Mat];
				D.f[DIR_PM0] = &DDStart[DIR_PM0 * size_Mat];
				D.f[DIR_MP0] = &DDStart[DIR_MP0 * size_Mat];
				D.f[DIR_P0P] = &DDStart[DIR_P0P * size_Mat];
				D.f[DIR_M0M] = &DDStart[DIR_M0M * size_Mat];
				D.f[DIR_P0M] = &DDStart[DIR_P0M * size_Mat];
				D.f[DIR_M0P] = &DDStart[DIR_M0P * size_Mat];
				D.f[DIR_0PP] = &DDStart[DIR_0PP * size_Mat];
				D.f[DIR_0MM] = &DDStart[DIR_0MM * size_Mat];
				D.f[DIR_0PM] = &DDStart[DIR_0PM * size_Mat];
				D.f[DIR_0MP] = &DDStart[DIR_0MP * size_Mat];
				D.f[DIR_000] = &DDStart[DIR_000 * size_Mat];
				D.f[DIR_PPP] = &DDStart[DIR_PPP * size_Mat];
				D.f[DIR_MMP] = &DDStart[DIR_MMP * size_Mat];
				D.f[DIR_PMP] = &DDStart[DIR_PMP * size_Mat];
				D.f[DIR_MPP] = &DDStart[DIR_MPP * size_Mat];
				D.f[DIR_PPM] = &DDStart[DIR_PPM * size_Mat];
				D.f[DIR_MMM] = &DDStart[DIR_MMM * size_Mat];
				D.f[DIR_PMM] = &DDStart[DIR_PMM * size_Mat];
				D.f[DIR_MPM] = &DDStart[DIR_MPM * size_Mat];
			}
			else
			{
				D.f[DIR_M00] = &DDStart[DIR_P00 * size_Mat];
				D.f[DIR_P00] = &DDStart[DIR_M00 * size_Mat];
				D.f[DIR_0M0] = &DDStart[DIR_0P0 * size_Mat];
				D.f[DIR_0P0] = &DDStart[DIR_0M0 * size_Mat];
				D.f[DIR_00M] = &DDStart[DIR_00P * size_Mat];
				D.f[DIR_00P] = &DDStart[DIR_00M * size_Mat];
				D.f[DIR_MM0] = &DDStart[DIR_PP0 * size_Mat];
				D.f[DIR_PP0] = &DDStart[DIR_MM0 * size_Mat];
				D.f[DIR_MP0] = &DDStart[DIR_PM0 * size_Mat];
				D.f[DIR_PM0] = &DDStart[DIR_MP0 * size_Mat];
				D.f[DIR_M0M] = &DDStart[DIR_P0P * size_Mat];
				D.f[DIR_P0P] = &DDStart[DIR_M0M * size_Mat];
				D.f[DIR_M0P] = &DDStart[DIR_P0M * size_Mat];
				D.f[DIR_P0M] = &DDStart[DIR_M0P * size_Mat];
				D.f[DIR_0MM] = &DDStart[DIR_0PP * size_Mat];
				D.f[DIR_0PP] = &DDStart[DIR_0MM * size_Mat];
				D.f[DIR_0MP] = &DDStart[DIR_0PM * size_Mat];
				D.f[DIR_0PM] = &DDStart[DIR_0MP * size_Mat];
				D.f[DIR_000] = &DDStart[DIR_000 * size_Mat];
				D.f[DIR_MMM] = &DDStart[DIR_PPP * size_Mat];
				D.f[DIR_PPM] = &DDStart[DIR_MMP * size_Mat];
				D.f[DIR_MPM] = &DDStart[DIR_PMP * size_Mat];
				D.f[DIR_PMM] = &DDStart[DIR_MPP * size_Mat];
				D.f[DIR_MMP] = &DDStart[DIR_PPM * size_Mat];
				D.f[DIR_PPP] = &DDStart[DIR_MMM * size_Mat];
				D.f[DIR_MPP] = &DDStart[DIR_PMM * size_Mat];
				D.f[DIR_PMP] = &DDStart[DIR_MPM * size_Mat];
			}

			Distributions7 D7;
			if (EvenOrOdd == true)
			{
				D7.f[0] = &DD7[0 * size_Mat];
				D7.f[1] = &DD7[1 * size_Mat];
				D7.f[2] = &DD7[2 * size_Mat];
				D7.f[3] = &DD7[3 * size_Mat];
				D7.f[4] = &DD7[4 * size_Mat];
				D7.f[5] = &DD7[5 * size_Mat];
				D7.f[6] = &DD7[6 * size_Mat];
			}
			else
			{
				D7.f[0] = &DD7[0 * size_Mat];
				D7.f[2] = &DD7[1 * size_Mat];
				D7.f[1] = &DD7[2 * size_Mat];
				D7.f[4] = &DD7[3 * size_Mat];
				D7.f[3] = &DD7[4 * size_Mat];
				D7.f[6] = &DD7[5 * size_Mat];
				D7.f[5] = &DD7[6 * size_Mat];
			}

			////////////////////////////////////////////////////////////////////////////////
			//index
			unsigned int kw = neighborX[k];
			unsigned int ks = neighborY[k];
			unsigned int kb = neighborZ[k];
			unsigned int ksw = neighborY[kw];
			unsigned int kbw = neighborZ[kw];
			unsigned int kbs = neighborZ[ks];
			unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real fW = (D.f[DIR_P00])[k];//ke
			real fE = (D.f[DIR_M00])[kw];
			real fS = (D.f[DIR_0P0])[k];//kn
			real fN = (D.f[DIR_0M0])[ks];
			real fB = (D.f[DIR_00P])[k];//kt
			real fT = (D.f[DIR_00M])[kb];
			real fSW = (D.f[DIR_PP0])[k];//kne
			real fNE = (D.f[DIR_MM0])[ksw];
			real fNW = (D.f[DIR_PM0])[ks];//kse
			real fSE = (D.f[DIR_MP0])[kw];//knw
			real fBW = (D.f[DIR_P0P])[k];//kte
			real fTE = (D.f[DIR_M0M])[kbw];
			real fTW = (D.f[DIR_P0M])[kb];//kbe
			real fBE = (D.f[DIR_M0P])[kw];//ktw
			real fBS = (D.f[DIR_0PP])[k];//ktn
			real fTN = (D.f[DIR_0MM])[kbs];
			real fTS = (D.f[DIR_0PM])[kb];//kbn
			real fBN = (D.f[DIR_0MP])[ks];//kts
			real fZERO = (D.f[DIR_000])[k];//kzero
			real fBSW = (D.f[DIR_PPP])[k];//ktne
			real fBNE = (D.f[DIR_MMP])[ksw];//ktsw
			real fBNW = (D.f[DIR_PMP])[ks];//ktse
			real fBSE = (D.f[DIR_MPP])[kw];//ktnw
			real fTSW = (D.f[DIR_PPM])[kb];//kbne
			real fTNE = (D.f[DIR_MMM])[kbsw];
			real fTNW = (D.f[DIR_PMM])[kbs];//kbse
			real fTSE = (D.f[DIR_MPM])[kbw];//kbnw
										   //real fE    =  (D.f[DIR_P00])[k  ];//ke
										   //real fW    =  (D.f[DIR_M00])[kw ];
										   //real fN    =  (D.f[DIR_0P0])[k  ];//kn
										   //real fS    =  (D.f[DIR_0M0])[ks ];
										   //real fT    =  (D.f[DIR_00P])[k  ];//kt
										   //real fB    =  (D.f[DIR_00M])[kb ];
										   //real fNE   =  (D.f[DIR_PP0])[k  ];//kne
										   //real fSW   =  (D.f[DIR_MM0])[ksw];
										   //real fSE   =  (D.f[DIR_PM0])[ks ];//kse
										   //real fNW   =  (D.f[DIR_MP0])[kw ];//knw
										   //real fTE   =  (D.f[DIR_P0P])[k  ];//kte
										   //real fBW   =  (D.f[DIR_M0M])[kbw];
										   //real fBE   =  (D.f[DIR_P0M])[kb ];//kbe
										   //real fTW   =  (D.f[DIR_M0P])[kw ];//ktw
										   //real fTN   =  (D.f[DIR_0PP])[k  ];//ktn
										   //real fBS   =  (D.f[DIR_0MM])[kbs];
										   //real fBN   =  (D.f[DIR_0PM])[kb ];//kbn
										   //real fTS   =  (D.f[DIR_0MP])[ks ];//kts
										   //real fZERO =  (D.f[DIR_000])[k  ];//kzero
										   //real fTNE   = (D.f[DIR_PPP])[k  ];//ktne
										   //real fTSW   = (D.f[DIR_MMP])[ksw];//ktsw
										   //real fTSE   = (D.f[DIR_PMP])[ks ];//ktse
										   //real fTNW   = (D.f[DIR_MPP])[kw ];//ktnw
										   //real fBNE   = (D.f[DIR_PPM])[kb ];//kbne
										   //real fBSW   = (D.f[DIR_MMM])[kbsw];
										   //real fBSE   = (D.f[DIR_PMM])[kbs];//kbse
										   //real fBNW   = (D.f[DIR_MPM])[kbw];//kbnw
										   ////////////////////////////////////////////////////////////////////////////////
			real f7ZERO = (D7.f[0])[k];
			real f7E = (D7.f[1])[k];
			real f7W = (D7.f[2])[kw];
			real f7N = (D7.f[3])[k];
			real f7S = (D7.f[4])[ks];
			real f7T = (D7.f[5])[k];
			real f7B = (D7.f[6])[kb];
			////////////////////////////////////////////////////////////////////////////////
			real rho0 = (fTNE + fBSW) + (fTSW + fBNE) + (fTSE + fBNW) + (fTNW + fBSE) + (fNE + fSW) + (fNW + fSE) + (fTE + fBW) + (fBE + fTW) + (fTN + fBS) + (fBN + fTS) + (fE + fW) + (fN + fS) + (fT + fB) + fZERO;
			real rho = rho0 + c1o1;
			real OORho = c1o1 / rho;
			real vx = OORho*((fTNE - fBSW) + (fBNE - fTSW) + (fTSE - fBNW) + (fBSE - fTNW) + (fNE - fSW) + (fSE - fNW) + (fTE - fBW) + (fBE - fTW) + (fE - fW));
			real vy = OORho*((fTNE - fBSW) + (fBNE - fTSW) + (fBNW - fTSE) + (fTNW - fBSE) + (fNE - fSW) + (fNW - fSE) + (fTN - fBS) + (fBN - fTS) + (fN - fS));
			real vz = OORho*((fTNE - fBSW) + (fTSW - fBNE) + (fTSE - fBNW) + (fTNW - fBSE) + (fTE - fBW) + (fTW - fBE) + (fTN - fBS) + (fTS - fBN) + (fT - fB));
			////////////////////////////////////////////////////////////////////////////////
			real omegaD = -c3o1 + sqrt(c3o1);
			real Lam = -(c1o2 + c1o1 / omegaD);
			real nue_d = Lam / c3o1;
			real ae = diffusivity / nue_d - c1o1;
			real ux_sq = vx * vx;
			real uy_sq = vy * vy;
			real uz_sq = vz * vz;

			real ConcD = f7ZERO + f7E + f7W + f7N + f7S + f7T + f7B;

			(D7.f[0])[k] = f7ZERO*(c1o1 + omegaD) - omegaD*ConcD*(c1o3*(ae*(-c3o1)) - (ux_sq + uy_sq + uz_sq));
			(D7.f[2])[kw] = f7E   *(c1o1 + omegaD) - omegaD*ConcD*(c1o6*(ae + c1o1) + c1o2*(ux_sq)+vx*c1o2);
			(D7.f[1])[k] = f7W   *(c1o1 + omegaD) - omegaD*ConcD*(c1o6*(ae + c1o1) + c1o2*(ux_sq)-vx*c1o2);
			(D7.f[4])[ks] = f7N   *(c1o1 + omegaD) - omegaD*ConcD*(c1o6*(ae + c1o1) + c1o2*(uy_sq)+vy*c1o2);
			(D7.f[3])[k] = f7S   *(c1o1 + omegaD) - omegaD*ConcD*(c1o6*(ae + c1o1) + c1o2*(uy_sq)-vy*c1o2);
			(D7.f[6])[kb] = f7T   *(c1o1 + omegaD) - omegaD*ConcD*(c1o6*(ae + c1o1) + c1o2*(uz_sq)+vz*c1o2);
			(D7.f[5])[k] = f7B   *(c1o1 + omegaD) - omegaD*ConcD*(c1o6*(ae + c1o1) + c1o2*(uz_sq)-vz*c1o2);
		}
	}
}