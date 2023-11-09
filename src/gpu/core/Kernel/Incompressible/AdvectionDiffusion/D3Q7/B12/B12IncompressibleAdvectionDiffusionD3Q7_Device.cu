#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void B12IncompressibleAdvectionDiffusionD3Q7_Device(real diffusivity,
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
				D.f[dP00] = &DDStart[dP00 * numberOfLBnodes];
				D.f[dM00] = &DDStart[dM00 * numberOfLBnodes];
				D.f[d0P0] = &DDStart[d0P0 * numberOfLBnodes];
				D.f[d0M0] = &DDStart[d0M0 * numberOfLBnodes];
				D.f[d00P] = &DDStart[d00P * numberOfLBnodes];
				D.f[d00M] = &DDStart[d00M * numberOfLBnodes];
				D.f[dPP0] = &DDStart[dPP0 * numberOfLBnodes];
				D.f[dMM0] = &DDStart[dMM0 * numberOfLBnodes];
				D.f[dPM0] = &DDStart[dPM0 * numberOfLBnodes];
				D.f[dMP0] = &DDStart[dMP0 * numberOfLBnodes];
				D.f[dP0P] = &DDStart[dP0P * numberOfLBnodes];
				D.f[dM0M] = &DDStart[dM0M * numberOfLBnodes];
				D.f[dP0M] = &DDStart[dP0M * numberOfLBnodes];
				D.f[dM0P] = &DDStart[dM0P * numberOfLBnodes];
				D.f[d0PP] = &DDStart[d0PP * numberOfLBnodes];
				D.f[d0MM] = &DDStart[d0MM * numberOfLBnodes];
				D.f[d0PM] = &DDStart[d0PM * numberOfLBnodes];
				D.f[d0MP] = &DDStart[d0MP * numberOfLBnodes];
				D.f[d000] = &DDStart[d000 * numberOfLBnodes];
				D.f[dPPP] = &DDStart[dPPP * numberOfLBnodes];
				D.f[dMMP] = &DDStart[dMMP * numberOfLBnodes];
				D.f[dPMP] = &DDStart[dPMP * numberOfLBnodes];
				D.f[dMPP] = &DDStart[dMPP * numberOfLBnodes];
				D.f[dPPM] = &DDStart[dPPM * numberOfLBnodes];
				D.f[dMMM] = &DDStart[dMMM * numberOfLBnodes];
				D.f[dPMM] = &DDStart[dPMM * numberOfLBnodes];
				D.f[dMPM] = &DDStart[dMPM * numberOfLBnodes];
			}
			else
			{
				D.f[dM00] = &DDStart[dP00 * numberOfLBnodes];
				D.f[dP00] = &DDStart[dM00 * numberOfLBnodes];
				D.f[d0M0] = &DDStart[d0P0 * numberOfLBnodes];
				D.f[d0P0] = &DDStart[d0M0 * numberOfLBnodes];
				D.f[d00M] = &DDStart[d00P * numberOfLBnodes];
				D.f[d00P] = &DDStart[d00M * numberOfLBnodes];
				D.f[dMM0] = &DDStart[dPP0 * numberOfLBnodes];
				D.f[dPP0] = &DDStart[dMM0 * numberOfLBnodes];
				D.f[dMP0] = &DDStart[dPM0 * numberOfLBnodes];
				D.f[dPM0] = &DDStart[dMP0 * numberOfLBnodes];
				D.f[dM0M] = &DDStart[dP0P * numberOfLBnodes];
				D.f[dP0P] = &DDStart[dM0M * numberOfLBnodes];
				D.f[dM0P] = &DDStart[dP0M * numberOfLBnodes];
				D.f[dP0M] = &DDStart[dM0P * numberOfLBnodes];
				D.f[d0MM] = &DDStart[d0PP * numberOfLBnodes];
				D.f[d0PP] = &DDStart[d0MM * numberOfLBnodes];
				D.f[d0MP] = &DDStart[d0PM * numberOfLBnodes];
				D.f[d0PM] = &DDStart[d0MP * numberOfLBnodes];
				D.f[d000] = &DDStart[d000 * numberOfLBnodes];
				D.f[dMMM] = &DDStart[dPPP * numberOfLBnodes];
				D.f[dPPM] = &DDStart[dMMP * numberOfLBnodes];
				D.f[dMPM] = &DDStart[dPMP * numberOfLBnodes];
				D.f[dPMM] = &DDStart[dMPP * numberOfLBnodes];
				D.f[dMMP] = &DDStart[dPPM * numberOfLBnodes];
				D.f[dPPP] = &DDStart[dMMM * numberOfLBnodes];
				D.f[dMPP] = &DDStart[dPMM * numberOfLBnodes];
				D.f[dPMP] = &DDStart[dMPM * numberOfLBnodes];
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
			//real fZERO = (D.f[d000])[k];//kzero
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