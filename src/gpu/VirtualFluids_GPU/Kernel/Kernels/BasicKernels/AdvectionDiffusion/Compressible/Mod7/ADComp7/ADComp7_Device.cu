#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
#include "math.h"

extern "C" __global__ void LB_Kernel_AD_Comp_7(real diffusivity,
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
				D.f[dirE] = &DDStart[dirE   *size_Mat];
				D.f[dirW] = &DDStart[dirW   *size_Mat];
				D.f[dirN] = &DDStart[dirN   *size_Mat];
				D.f[dirS] = &DDStart[dirS   *size_Mat];
				D.f[dirT] = &DDStart[dirT   *size_Mat];
				D.f[dirB] = &DDStart[dirB   *size_Mat];
				D.f[dirNE] = &DDStart[dirNE  *size_Mat];
				D.f[dirSW] = &DDStart[dirSW  *size_Mat];
				D.f[dirSE] = &DDStart[dirSE  *size_Mat];
				D.f[dirNW] = &DDStart[dirNW  *size_Mat];
				D.f[dirTE] = &DDStart[dirTE  *size_Mat];
				D.f[dirBW] = &DDStart[dirBW  *size_Mat];
				D.f[dirBE] = &DDStart[dirBE  *size_Mat];
				D.f[dirTW] = &DDStart[dirTW  *size_Mat];
				D.f[dirTN] = &DDStart[dirTN  *size_Mat];
				D.f[dirBS] = &DDStart[dirBS  *size_Mat];
				D.f[dirBN] = &DDStart[dirBN  *size_Mat];
				D.f[dirTS] = &DDStart[dirTS  *size_Mat];
				D.f[dirREST] = &DDStart[dirREST*size_Mat];
				D.f[dirTNE] = &DDStart[dirTNE *size_Mat];
				D.f[dirTSW] = &DDStart[dirTSW *size_Mat];
				D.f[dirTSE] = &DDStart[dirTSE *size_Mat];
				D.f[dirTNW] = &DDStart[dirTNW *size_Mat];
				D.f[dirBNE] = &DDStart[dirBNE *size_Mat];
				D.f[dirBSW] = &DDStart[dirBSW *size_Mat];
				D.f[dirBSE] = &DDStart[dirBSE *size_Mat];
				D.f[dirBNW] = &DDStart[dirBNW *size_Mat];
			}
			else
			{
				D.f[dirW] = &DDStart[dirE   *size_Mat];
				D.f[dirE] = &DDStart[dirW   *size_Mat];
				D.f[dirS] = &DDStart[dirN   *size_Mat];
				D.f[dirN] = &DDStart[dirS   *size_Mat];
				D.f[dirB] = &DDStart[dirT   *size_Mat];
				D.f[dirT] = &DDStart[dirB   *size_Mat];
				D.f[dirSW] = &DDStart[dirNE  *size_Mat];
				D.f[dirNE] = &DDStart[dirSW  *size_Mat];
				D.f[dirNW] = &DDStart[dirSE  *size_Mat];
				D.f[dirSE] = &DDStart[dirNW  *size_Mat];
				D.f[dirBW] = &DDStart[dirTE  *size_Mat];
				D.f[dirTE] = &DDStart[dirBW  *size_Mat];
				D.f[dirTW] = &DDStart[dirBE  *size_Mat];
				D.f[dirBE] = &DDStart[dirTW  *size_Mat];
				D.f[dirBS] = &DDStart[dirTN  *size_Mat];
				D.f[dirTN] = &DDStart[dirBS  *size_Mat];
				D.f[dirTS] = &DDStart[dirBN  *size_Mat];
				D.f[dirBN] = &DDStart[dirTS  *size_Mat];
				D.f[dirREST] = &DDStart[dirREST*size_Mat];
				D.f[dirBSW] = &DDStart[dirTNE *size_Mat];
				D.f[dirBNE] = &DDStart[dirTSW *size_Mat];
				D.f[dirBNW] = &DDStart[dirTSE *size_Mat];
				D.f[dirBSE] = &DDStart[dirTNW *size_Mat];
				D.f[dirTSW] = &DDStart[dirBNE *size_Mat];
				D.f[dirTNE] = &DDStart[dirBSW *size_Mat];
				D.f[dirTNW] = &DDStart[dirBSE *size_Mat];
				D.f[dirTSE] = &DDStart[dirBNW *size_Mat];
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
			real fW = (D.f[dirE])[k];//ke
			real fE = (D.f[dirW])[kw];
			real fS = (D.f[dirN])[k];//kn
			real fN = (D.f[dirS])[ks];
			real fB = (D.f[dirT])[k];//kt
			real fT = (D.f[dirB])[kb];
			real fSW = (D.f[dirNE])[k];//kne
			real fNE = (D.f[dirSW])[ksw];
			real fNW = (D.f[dirSE])[ks];//kse
			real fSE = (D.f[dirNW])[kw];//knw
			real fBW = (D.f[dirTE])[k];//kte
			real fTE = (D.f[dirBW])[kbw];
			real fTW = (D.f[dirBE])[kb];//kbe
			real fBE = (D.f[dirTW])[kw];//ktw
			real fBS = (D.f[dirTN])[k];//ktn
			real fTN = (D.f[dirBS])[kbs];
			real fTS = (D.f[dirBN])[kb];//kbn
			real fBN = (D.f[dirTS])[ks];//kts
			real fZERO = (D.f[dirREST])[k];//kzero
			real fBSW = (D.f[dirTNE])[k];//ktne
			real fBNE = (D.f[dirTSW])[ksw];//ktsw
			real fBNW = (D.f[dirTSE])[ks];//ktse
			real fBSE = (D.f[dirTNW])[kw];//ktnw
			real fTSW = (D.f[dirBNE])[kb];//kbne
			real fTNE = (D.f[dirBSW])[kbsw];
			real fTNW = (D.f[dirBSE])[kbs];//kbse
			real fTSE = (D.f[dirBNW])[kbw];//kbnw
										   //real fE    =  (D.f[dirE   ])[k  ];//ke
										   //real fW    =  (D.f[dirW   ])[kw ];
										   //real fN    =  (D.f[dirN   ])[k  ];//kn
										   //real fS    =  (D.f[dirS   ])[ks ];
										   //real fT    =  (D.f[dirT   ])[k  ];//kt
										   //real fB    =  (D.f[dirB   ])[kb ];
										   //real fNE   =  (D.f[dirNE  ])[k  ];//kne
										   //real fSW   =  (D.f[dirSW  ])[ksw];
										   //real fSE   =  (D.f[dirSE  ])[ks ];//kse
										   //real fNW   =  (D.f[dirNW  ])[kw ];//knw
										   //real fTE   =  (D.f[dirTE  ])[k  ];//kte
										   //real fBW   =  (D.f[dirBW  ])[kbw];
										   //real fBE   =  (D.f[dirBE  ])[kb ];//kbe
										   //real fTW   =  (D.f[dirTW  ])[kw ];//ktw
										   //real fTN   =  (D.f[dirTN  ])[k  ];//ktn
										   //real fBS   =  (D.f[dirBS  ])[kbs];
										   //real fBN   =  (D.f[dirBN  ])[kb ];//kbn
										   //real fTS   =  (D.f[dirTS  ])[ks ];//kts
										   //real fZERO =  (D.f[dirREST])[k  ];//kzero
										   //real fTNE   = (D.f[dirTNE ])[k  ];//ktne
										   //real fTSW   = (D.f[dirTSW ])[ksw];//ktsw
										   //real fTSE   = (D.f[dirTSE ])[ks ];//ktse
										   //real fTNW   = (D.f[dirTNW ])[kw ];//ktnw
										   //real fBNE   = (D.f[dirBNE ])[kb ];//kbne
										   //real fBSW   = (D.f[dirBSW ])[kbsw];
										   //real fBSE   = (D.f[dirBSE ])[kbs];//kbse
										   //real fBNW   = (D.f[dirBNW ])[kbw];//kbnw
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