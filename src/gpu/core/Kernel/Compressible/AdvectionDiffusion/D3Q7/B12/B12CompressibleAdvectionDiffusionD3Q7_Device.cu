#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

#include "math.h"

__global__ void B12CompressibleAdvectionDiffusionD3Q7_Device(real diffusivity,
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
				D.f[dP00] = &DDStart[dP00 * size_Mat];
				D.f[dM00] = &DDStart[dM00 * size_Mat];
				D.f[d0P0] = &DDStart[d0P0 * size_Mat];
				D.f[d0M0] = &DDStart[d0M0 * size_Mat];
				D.f[d00P] = &DDStart[d00P * size_Mat];
				D.f[d00M] = &DDStart[d00M * size_Mat];
				D.f[dPP0] = &DDStart[dPP0 * size_Mat];
				D.f[dMM0] = &DDStart[dMM0 * size_Mat];
				D.f[dPM0] = &DDStart[dPM0 * size_Mat];
				D.f[dMP0] = &DDStart[dMP0 * size_Mat];
				D.f[dP0P] = &DDStart[dP0P * size_Mat];
				D.f[dM0M] = &DDStart[dM0M * size_Mat];
				D.f[dP0M] = &DDStart[dP0M * size_Mat];
				D.f[dM0P] = &DDStart[dM0P * size_Mat];
				D.f[d0PP] = &DDStart[d0PP * size_Mat];
				D.f[d0MM] = &DDStart[d0MM * size_Mat];
				D.f[d0PM] = &DDStart[d0PM * size_Mat];
				D.f[d0MP] = &DDStart[d0MP * size_Mat];
				D.f[d000] = &DDStart[d000 * size_Mat];
				D.f[dPPP] = &DDStart[dPPP * size_Mat];
				D.f[dMMP] = &DDStart[dMMP * size_Mat];
				D.f[dPMP] = &DDStart[dPMP * size_Mat];
				D.f[dMPP] = &DDStart[dMPP * size_Mat];
				D.f[dPPM] = &DDStart[dPPM * size_Mat];
				D.f[dMMM] = &DDStart[dMMM * size_Mat];
				D.f[dPMM] = &DDStart[dPMM * size_Mat];
				D.f[dMPM] = &DDStart[dMPM * size_Mat];
			}
			else
			{
				D.f[dM00] = &DDStart[dP00 * size_Mat];
				D.f[dP00] = &DDStart[dM00 * size_Mat];
				D.f[d0M0] = &DDStart[d0P0 * size_Mat];
				D.f[d0P0] = &DDStart[d0M0 * size_Mat];
				D.f[d00M] = &DDStart[d00P * size_Mat];
				D.f[d00P] = &DDStart[d00M * size_Mat];
				D.f[dMM0] = &DDStart[dPP0 * size_Mat];
				D.f[dPP0] = &DDStart[dMM0 * size_Mat];
				D.f[dMP0] = &DDStart[dPM0 * size_Mat];
				D.f[dPM0] = &DDStart[dMP0 * size_Mat];
				D.f[dM0M] = &DDStart[dP0P * size_Mat];
				D.f[dP0P] = &DDStart[dM0M * size_Mat];
				D.f[dM0P] = &DDStart[dP0M * size_Mat];
				D.f[dP0M] = &DDStart[dM0P * size_Mat];
				D.f[d0MM] = &DDStart[d0PP * size_Mat];
				D.f[d0PP] = &DDStart[d0MM * size_Mat];
				D.f[d0MP] = &DDStart[d0PM * size_Mat];
				D.f[d0PM] = &DDStart[d0MP * size_Mat];
				D.f[d000] = &DDStart[d000 * size_Mat];
				D.f[dMMM] = &DDStart[dPPP * size_Mat];
				D.f[dPPM] = &DDStart[dMMP * size_Mat];
				D.f[dMPM] = &DDStart[dPMP * size_Mat];
				D.f[dPMM] = &DDStart[dMPP * size_Mat];
				D.f[dMMP] = &DDStart[dPPM * size_Mat];
				D.f[dPPP] = &DDStart[dMMM * size_Mat];
				D.f[dMPP] = &DDStart[dPMM * size_Mat];
				D.f[dPMP] = &DDStart[dMPM * size_Mat];
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
			real fW = (D.f[dP00])[k];//ke
			real fE = (D.f[dM00])[kw];
			real fS = (D.f[d0P0])[k];//kn
			real fN = (D.f[d0M0])[ks];
			real fB = (D.f[d00P])[k];//kt
			real fT = (D.f[d00M])[kb];
			real fSW = (D.f[dPP0])[k];//kne
			real fNE = (D.f[dMM0])[ksw];
			real fNW = (D.f[dPM0])[ks];//kse
			real fSE = (D.f[dMP0])[kw];//knw
			real fBW = (D.f[dP0P])[k];//kte
			real fTE = (D.f[dM0M])[kbw];
			real fTW = (D.f[dP0M])[kb];//kbe
			real fBE = (D.f[dM0P])[kw];//ktw
			real fBS = (D.f[d0PP])[k];//ktn
			real fTN = (D.f[d0MM])[kbs];
			real fTS = (D.f[d0PM])[kb];//kbn
			real fBN = (D.f[d0MP])[ks];//kts
			real fZERO = (D.f[d000])[k];//kzero
			real fBSW = (D.f[dPPP])[k];//ktne
			real fBNE = (D.f[dMMP])[ksw];//ktsw
			real fBNW = (D.f[dPMP])[ks];//ktse
			real fBSE = (D.f[dMPP])[kw];//ktnw
			real fTSW = (D.f[dPPM])[kb];//kbne
			real fTNE = (D.f[dMMM])[kbsw];
			real fTNW = (D.f[dPMM])[kbs];//kbse
			real fTSE = (D.f[dMPM])[kbw];//kbnw
										   //real fE    =  (D.f[dP00])[k  ];//ke
										   //real fW    =  (D.f[dM00])[kw ];
										   //real fN    =  (D.f[d0P0])[k  ];//kn
										   //real fS    =  (D.f[d0M0])[ks ];
										   //real fT    =  (D.f[d00P])[k  ];//kt
										   //real fB    =  (D.f[d00M])[kb ];
										   //real fNE   =  (D.f[dPP0])[k  ];//kne
										   //real fSW   =  (D.f[dMM0])[ksw];
										   //real fSE   =  (D.f[dPM0])[ks ];//kse
										   //real fNW   =  (D.f[dMP0])[kw ];//knw
										   //real fTE   =  (D.f[dP0P])[k  ];//kte
										   //real fBW   =  (D.f[dM0M])[kbw];
										   //real fBE   =  (D.f[dP0M])[kb ];//kbe
										   //real fTW   =  (D.f[dM0P])[kw ];//ktw
										   //real fTN   =  (D.f[d0PP])[k  ];//ktn
										   //real fBS   =  (D.f[d0MM])[kbs];
										   //real fBN   =  (D.f[d0PM])[kb ];//kbn
										   //real fTS   =  (D.f[d0MP])[ks ];//kts
										   //real fZERO =  (D.f[d000])[k  ];//kzero
										   //real fTNE   = (D.f[dPPP])[k  ];//ktne
										   //real fTSW   = (D.f[dMMP])[ksw];//ktsw
										   //real fTSE   = (D.f[dPMP])[ks ];//ktse
										   //real fTNW   = (D.f[dMPP])[kw ];//ktnw
										   //real fBNE   = (D.f[dPPM])[kb ];//kbne
										   //real fBSW   = (D.f[dMMM])[kbsw];
										   //real fBSE   = (D.f[dPMM])[kbs];//kbse
										   //real fBNW   = (D.f[dMPM])[kbw];//kbnw
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