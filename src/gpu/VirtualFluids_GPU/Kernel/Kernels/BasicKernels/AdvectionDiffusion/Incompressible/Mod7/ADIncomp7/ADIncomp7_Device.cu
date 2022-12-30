#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void LB_Kernel_AD_Incomp_7(real diffusivity,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	real* DD7,
	unsigned long long numberOfLBnodes,
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

	if (k<numberOfLBnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if ((BC != GEO_SOLID) && (BC != GEO_VOID))
		{
			Distributions27 D;
			if (EvenOrOdd == true)
			{
				D.f[DIR_P00] = &DDStart[DIR_P00   *numberOfLBnodes];
				D.f[DIR_M00] = &DDStart[DIR_M00   *numberOfLBnodes];
				D.f[DIR_0P0] = &DDStart[DIR_0P0   *numberOfLBnodes];
				D.f[DIR_0M0] = &DDStart[DIR_0M0   *numberOfLBnodes];
				D.f[DIR_00P] = &DDStart[DIR_00P   *numberOfLBnodes];
				D.f[DIR_00M] = &DDStart[DIR_00M   *numberOfLBnodes];
				D.f[DIR_PP0] = &DDStart[DIR_PP0  *numberOfLBnodes];
				D.f[DIR_MM0] = &DDStart[DIR_MM0  *numberOfLBnodes];
				D.f[DIR_PM0] = &DDStart[DIR_PM0  *numberOfLBnodes];
				D.f[DIR_MP0] = &DDStart[DIR_MP0  *numberOfLBnodes];
				D.f[DIR_P0P] = &DDStart[DIR_P0P  *numberOfLBnodes];
				D.f[DIR_M0M] = &DDStart[DIR_M0M  *numberOfLBnodes];
				D.f[DIR_P0M] = &DDStart[DIR_P0M  *numberOfLBnodes];
				D.f[DIR_M0P] = &DDStart[DIR_M0P  *numberOfLBnodes];
				D.f[DIR_0PP] = &DDStart[DIR_0PP  *numberOfLBnodes];
				D.f[DIR_0MM] = &DDStart[DIR_0MM  *numberOfLBnodes];
				D.f[DIR_0PM] = &DDStart[DIR_0PM  *numberOfLBnodes];
				D.f[DIR_0MP] = &DDStart[DIR_0MP  *numberOfLBnodes];
				D.f[DIR_000] = &DDStart[DIR_000*numberOfLBnodes];
				D.f[DIR_PPP] = &DDStart[DIR_PPP *numberOfLBnodes];
				D.f[DIR_MMP] = &DDStart[DIR_MMP *numberOfLBnodes];
				D.f[DIR_PMP] = &DDStart[DIR_PMP *numberOfLBnodes];
				D.f[DIR_MPP] = &DDStart[DIR_MPP *numberOfLBnodes];
				D.f[DIR_PPM] = &DDStart[DIR_PPM *numberOfLBnodes];
				D.f[DIR_MMM] = &DDStart[DIR_MMM *numberOfLBnodes];
				D.f[DIR_PMM] = &DDStart[DIR_PMM *numberOfLBnodes];
				D.f[DIR_MPM] = &DDStart[DIR_MPM *numberOfLBnodes];
			}
			else
			{
				D.f[DIR_M00] = &DDStart[DIR_P00   *numberOfLBnodes];
				D.f[DIR_P00] = &DDStart[DIR_M00   *numberOfLBnodes];
				D.f[DIR_0M0] = &DDStart[DIR_0P0   *numberOfLBnodes];
				D.f[DIR_0P0] = &DDStart[DIR_0M0   *numberOfLBnodes];
				D.f[DIR_00M] = &DDStart[DIR_00P   *numberOfLBnodes];
				D.f[DIR_00P] = &DDStart[DIR_00M   *numberOfLBnodes];
				D.f[DIR_MM0] = &DDStart[DIR_PP0  *numberOfLBnodes];
				D.f[DIR_PP0] = &DDStart[DIR_MM0  *numberOfLBnodes];
				D.f[DIR_MP0] = &DDStart[DIR_PM0  *numberOfLBnodes];
				D.f[DIR_PM0] = &DDStart[DIR_MP0  *numberOfLBnodes];
				D.f[DIR_M0M] = &DDStart[DIR_P0P  *numberOfLBnodes];
				D.f[DIR_P0P] = &DDStart[DIR_M0M  *numberOfLBnodes];
				D.f[DIR_M0P] = &DDStart[DIR_P0M  *numberOfLBnodes];
				D.f[DIR_P0M] = &DDStart[DIR_M0P  *numberOfLBnodes];
				D.f[DIR_0MM] = &DDStart[DIR_0PP  *numberOfLBnodes];
				D.f[DIR_0PP] = &DDStart[DIR_0MM  *numberOfLBnodes];
				D.f[DIR_0MP] = &DDStart[DIR_0PM  *numberOfLBnodes];
				D.f[DIR_0PM] = &DDStart[DIR_0MP  *numberOfLBnodes];
				D.f[DIR_000] = &DDStart[DIR_000*numberOfLBnodes];
				D.f[DIR_MMM] = &DDStart[DIR_PPP *numberOfLBnodes];
				D.f[DIR_PPM] = &DDStart[DIR_MMP *numberOfLBnodes];
				D.f[DIR_MPM] = &DDStart[DIR_PMP *numberOfLBnodes];
				D.f[DIR_PMM] = &DDStart[DIR_MPP *numberOfLBnodes];
				D.f[DIR_MMP] = &DDStart[DIR_PPM *numberOfLBnodes];
				D.f[DIR_PPP] = &DDStart[DIR_MMM *numberOfLBnodes];
				D.f[DIR_MPP] = &DDStart[DIR_PMM *numberOfLBnodes];
				D.f[DIR_PMP] = &DDStart[DIR_MPM *numberOfLBnodes];
			}

			Distributions7 D7;
			if (EvenOrOdd == true)
			{
				D7.f[0] = &DD7[0 * numberOfLBnodes];
				D7.f[1] = &DD7[1 * numberOfLBnodes];
				D7.f[2] = &DD7[2 * numberOfLBnodes];
				D7.f[3] = &DD7[3 * numberOfLBnodes];
				D7.f[4] = &DD7[4 * numberOfLBnodes];
				D7.f[5] = &DD7[5 * numberOfLBnodes];
				D7.f[6] = &DD7[6 * numberOfLBnodes];
			}
			else
			{
				D7.f[0] = &DD7[0 * numberOfLBnodes];
				D7.f[2] = &DD7[1 * numberOfLBnodes];
				D7.f[1] = &DD7[2 * numberOfLBnodes];
				D7.f[4] = &DD7[3 * numberOfLBnodes];
				D7.f[3] = &DD7[4 * numberOfLBnodes];
				D7.f[6] = &DD7[5 * numberOfLBnodes];
				D7.f[5] = &DD7[6 * numberOfLBnodes];
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
			//real fZERO = (D.f[DIR_000])[k];//kzero
			real fBSW = (D.f[DIR_PPP])[k];//ktne
			real fBNE = (D.f[DIR_MMP])[ksw];//ktsw
			real fBNW = (D.f[DIR_PMP])[ks];//ktse
			real fBSE = (D.f[DIR_MPP])[kw];//ktnw
			real fTSW = (D.f[DIR_PPM])[kb];//kbne
			real fTNE = (D.f[DIR_MMM])[kbsw];
			real fTNW = (D.f[DIR_PMM])[kbs];//kbse
			real fTSE = (D.f[DIR_MPM])[kbw];//kbnw
										   //real fE    =  (D.f[DIR_P00   ])[k  ];//ke
										   //real fW    =  (D.f[DIR_M00   ])[kw ];
										   //real fN    =  (D.f[DIR_0P0   ])[k  ];//kn
										   //real fS    =  (D.f[DIR_0M0   ])[ks ];
										   //real fT    =  (D.f[DIR_00P   ])[k  ];//kt
										   //real fB    =  (D.f[DIR_00M   ])[kb ];
										   //real fNE   =  (D.f[DIR_PP0  ])[k  ];//kne
										   //real fSW   =  (D.f[DIR_MM0  ])[ksw];
										   //real fSE   =  (D.f[DIR_PM0  ])[ks ];//kse
										   //real fNW   =  (D.f[DIR_MP0  ])[kw ];//knw
										   //real fTE   =  (D.f[DIR_P0P  ])[k  ];//kte
										   //real fBW   =  (D.f[DIR_M0M  ])[kbw];
										   //real fBE   =  (D.f[DIR_P0M  ])[kb ];//kbe
										   //real fTW   =  (D.f[DIR_M0P  ])[kw ];//ktw
										   //real fTN   =  (D.f[DIR_0PP  ])[k  ];//ktn
										   //real fBS   =  (D.f[DIR_0MM  ])[kbs];
										   //real fBN   =  (D.f[DIR_0PM  ])[kb ];//kbn
										   //real fTS   =  (D.f[DIR_0MP  ])[ks ];//kts
										   //real fZERO =  (D.f[DIR_000])[k  ];//kzero
										   //real fTNE   = (D.f[DIR_PPP ])[k  ];//ktne
										   //real fTSW   = (D.f[DIR_MMP ])[ksw];//ktsw
										   //real fTSE   = (D.f[DIR_PMP ])[ks ];//ktse
										   //real fTNW   = (D.f[DIR_MPP ])[kw ];//ktnw
										   //real fBNE   = (D.f[DIR_PPM ])[kb ];//kbne
										   //real fBSW   = (D.f[DIR_MMM ])[kbsw];
										   //real fBSE   = (D.f[DIR_PMM ])[kbs];//kbse
										   //real fBNW   = (D.f[DIR_MPM ])[kbw];//kbnw
										   ////////////////////////////////////////////////////////////////////////////////
			real f7ZERO = (D7.f[0])[k];
			real f7E = (D7.f[1])[k];
			real f7W = (D7.f[2])[kw];
			real f7N = (D7.f[3])[k];
			real f7S = (D7.f[4])[ks];
			real f7T = (D7.f[5])[k];
			real f7B = (D7.f[6])[kb];
			////////////////////////////////////////////////////////////////////////////////
			real vx = ((fTNE - fBSW) + (fBNE - fTSW) + (fTSE - fBNW) + (fBSE - fTNW) + (fNE - fSW) + (fSE - fNW) + (fTE - fBW) + (fBE - fTW) + (fE - fW));
			real vy = ((fTNE - fBSW) + (fBNE - fTSW) + (fBNW - fTSE) + (fTNW - fBSE) + (fNE - fSW) + (fNW - fSE) + (fTN - fBS) + (fBN - fTS) + (fN - fS));
			real vz = ((fTNE - fBSW) + (fTSW - fBNE) + (fTSE - fBNW) + (fTNW - fBSE) + (fTE - fBW) + (fTW - fBE) + (fTN - fBS) + (fTS - fBN) + (fT - fB));
			////d�rrrrrty !!!!!!!!!!!!!
			//      real vx     =  ten * ((fTNE-fBSW)+(fBNE-fTSW)+(fTSE-fBNW)+(fBSE-fTNW) +(fNE-fSW)+(fSE-fNW)+(fTE-fBW)+(fBE-fTW)+(fE-fW));
			//      real vy     =  ten * ((fTNE-fBSW)+(fBNE-fTSW)+(fBNW-fTSE)+(fTNW-fBSE) +(fNE-fSW)+(fNW-fSE)+(fTN-fBS)+(fBN-fTS)+(fN-fS));
			//      real vz     =  ten * ((fTNE-fBSW)+(fTSW-fBNE)+(fTSE-fBNW)+(fTNW-fBSE) +(fTE-fBW)+(fTW-fBE)+(fTN-fBS)+(fTS-fBN)+(fT-fB));
			////////////////////////////////////////////////////////////////////////////////
			//real ux_sq = vx * vx;
			//real uy_sq = vy * vy;
			//real uz_sq = vz * vz;
			////////////////////////////////////////////////////////////////////////////////
			//BGK
			//real omegaD     = -three + sqrt(three); !!!!!!!!!!!!!!Achtung!!!!!!!!!!!!!!!!!! anderes Vorzeichen als in den Randbedingungen
			//real Lam         = -(c1o2+one/omegaD);
			//real nue_d       = Lam/three;
			//real ae          = diffusivity/nue_d - one;

			//real ConcD       = f7ZERO+f7E+f7W+f7N+f7S+f7T+f7B;

			//(D7.f[0])[k  ] = f7ZERO*(one+omegaD)-omegaD*ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
			//(D7.f[2])[kw ] = f7E   *(one+omegaD)-omegaD*ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx*c1o2);
			//(D7.f[1])[k  ] = f7W   *(one+omegaD)-omegaD*ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx*c1o2);
			//(D7.f[4])[ks ] = f7N   *(one+omegaD)-omegaD*ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vy*c1o2);
			//(D7.f[3])[k  ] = f7S   *(one+omegaD)-omegaD*ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vy*c1o2);
			//(D7.f[6])[kb ] = f7T   *(one+omegaD)-omegaD*ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vz*c1o2);
			//(D7.f[5])[k  ] = f7B   *(one+omegaD)-omegaD*ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vz*c1o2);

			////////////////////////////////////////////////////////////////////////////////
			//TRT  Yoshida Kernel - based on Ying
			//real cs2 = c1o4;
			real Lam = diffusivity*c4o1;//diffusivity/(one)/cs2;
			real omegaD = -c1o1 / (Lam + c1o2);
			//real ae = c0o1;
			////////////////////////////////////////////////////////////////////////////////
			real ConcD = f7ZERO + f7E + f7W + f7N + f7S + f7T + f7B;

			real Mom000 = f7ZERO + f7W + f7E + f7N + f7S + f7T + f7B; //1
			real Mom100 = f7E - f7W;
			real Mom010 = f7N - f7S;
			real Mom001 = f7T - f7B;
			real Mom222 = c6o1*f7ZERO - f7W - f7E - f7N - f7S - f7T - f7B;
			real Mom200 = c2o1*f7W + c2o1*f7E - f7N - f7S - f7T - f7B;
			real Mom022 = f7N + f7S - f7T - f7B;

			real Meq000 = ConcD;
			real Meq100 = ConcD*vx;
			real Meq010 = ConcD*vy;
			real Meq001 = ConcD*vz;
			real Meq222 = c3o4*ConcD;
			real Meq200 = c0o1;
			real Meq022 = c0o1;

			// relaxation TRT Yoshida

			// odd 
			Mom100 = omegaD * (Mom100 - Meq100);
			Mom010 = omegaD * (Mom010 - Meq010);
			Mom001 = omegaD * (Mom001 - Meq001);

			// even
			Mom000 = -c1o1*(Mom000 - Meq000);
			Mom222 = -c1o1*(Mom222 - Meq222);
			Mom200 = -c1o1*(Mom200 - Meq200);
			Mom022 = -c1o1*(Mom022 - Meq022);

			//Back transformation to distributions
			f7ZERO = f7ZERO + c1o7*Mom000 + c1o7*Mom222;                                                  //1
			f7E = f7E + c1o7*Mom000 + c1o2*Mom100 - c1o6*c1o7*Mom222 + c1o6*Mom200;                 //2
			f7W = f7W + c1o7*Mom000 - c1o2*Mom100 - c1o6*c1o7*Mom222 + c1o6*Mom200;                 //3
			f7N = f7N + c1o7*Mom000 + c1o2*Mom010 - c1o6*c1o7*Mom222 - c1o12*Mom200 + c1o4 *Mom022; //4
			f7S = f7S + c1o7*Mom000 - c1o2*Mom010 - c1o6*c1o7*Mom222 - c1o12*Mom200 + c1o4 *Mom022; //5
			f7T = f7T + c1o7*Mom000 + c1o2*Mom001 - c1o6*c1o7*Mom222 - c1o12*Mom200 - c1o4 *Mom022; //6
			f7B = f7B + c1o7*Mom000 - c1o2*Mom001 - c1o6*c1o7*Mom222 - c1o12*Mom200 - c1o4 *Mom022; //7

			(D7.f[0])[k] = f7ZERO;
			(D7.f[2])[kw] = f7E;
			(D7.f[1])[k] = f7W;
			(D7.f[4])[ks] = f7N;
			(D7.f[3])[k] = f7S;
			(D7.f[6])[kb] = f7T;
			(D7.f[5])[k] = f7B;
		}
	}
}