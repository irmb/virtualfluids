#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
#include "math.h"

extern "C" __global__ void LB_Kernel_AD_Incomp_7(real diffusivity,
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
				D.f[E] = &DDStart[E   *size_Mat];
				D.f[W] = &DDStart[W   *size_Mat];
				D.f[N] = &DDStart[N   *size_Mat];
				D.f[S] = &DDStart[S   *size_Mat];
				D.f[T] = &DDStart[T   *size_Mat];
				D.f[B] = &DDStart[B   *size_Mat];
				D.f[NE] = &DDStart[NE  *size_Mat];
				D.f[SW] = &DDStart[SW  *size_Mat];
				D.f[SE] = &DDStart[SE  *size_Mat];
				D.f[NW] = &DDStart[NW  *size_Mat];
				D.f[TE] = &DDStart[TE  *size_Mat];
				D.f[BW] = &DDStart[BW  *size_Mat];
				D.f[BE] = &DDStart[BE  *size_Mat];
				D.f[TW] = &DDStart[TW  *size_Mat];
				D.f[TN] = &DDStart[TN  *size_Mat];
				D.f[BS] = &DDStart[BS  *size_Mat];
				D.f[BN] = &DDStart[BN  *size_Mat];
				D.f[TS] = &DDStart[TS  *size_Mat];
				D.f[dirREST] = &DDStart[dirREST*size_Mat];
				D.f[TNE] = &DDStart[TNE *size_Mat];
				D.f[TSW] = &DDStart[TSW *size_Mat];
				D.f[TSE] = &DDStart[TSE *size_Mat];
				D.f[TNW] = &DDStart[TNW *size_Mat];
				D.f[BNE] = &DDStart[BNE *size_Mat];
				D.f[BSW] = &DDStart[BSW *size_Mat];
				D.f[BSE] = &DDStart[BSE *size_Mat];
				D.f[BNW] = &DDStart[BNW *size_Mat];
			}
			else
			{
				D.f[W] = &DDStart[E   *size_Mat];
				D.f[E] = &DDStart[W   *size_Mat];
				D.f[S] = &DDStart[N   *size_Mat];
				D.f[N] = &DDStart[S   *size_Mat];
				D.f[B] = &DDStart[T   *size_Mat];
				D.f[T] = &DDStart[B   *size_Mat];
				D.f[SW] = &DDStart[NE  *size_Mat];
				D.f[NE] = &DDStart[SW  *size_Mat];
				D.f[NW] = &DDStart[SE  *size_Mat];
				D.f[SE] = &DDStart[NW  *size_Mat];
				D.f[BW] = &DDStart[TE  *size_Mat];
				D.f[TE] = &DDStart[BW  *size_Mat];
				D.f[TW] = &DDStart[BE  *size_Mat];
				D.f[BE] = &DDStart[TW  *size_Mat];
				D.f[BS] = &DDStart[TN  *size_Mat];
				D.f[TN] = &DDStart[BS  *size_Mat];
				D.f[TS] = &DDStart[BN  *size_Mat];
				D.f[BN] = &DDStart[TS  *size_Mat];
				D.f[dirREST] = &DDStart[dirREST*size_Mat];
				D.f[BSW] = &DDStart[TNE *size_Mat];
				D.f[BNE] = &DDStart[TSW *size_Mat];
				D.f[BNW] = &DDStart[TSE *size_Mat];
				D.f[BSE] = &DDStart[TNW *size_Mat];
				D.f[TSW] = &DDStart[BNE *size_Mat];
				D.f[TNE] = &DDStart[BSW *size_Mat];
				D.f[TNW] = &DDStart[BSE *size_Mat];
				D.f[TSE] = &DDStart[BNW *size_Mat];
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
			real fW = (D.f[E])[k];//ke
			real fE = (D.f[W])[kw];
			real fS = (D.f[N])[k];//kn
			real fN = (D.f[S])[ks];
			real fB = (D.f[T])[k];//kt
			real fT = (D.f[B])[kb];
			real fSW = (D.f[NE])[k];//kne
			real fNE = (D.f[SW])[ksw];
			real fNW = (D.f[SE])[ks];//kse
			real fSE = (D.f[NW])[kw];//knw
			real fBW = (D.f[TE])[k];//kte
			real fTE = (D.f[BW])[kbw];
			real fTW = (D.f[BE])[kb];//kbe
			real fBE = (D.f[TW])[kw];//ktw
			real fBS = (D.f[TN])[k];//ktn
			real fTN = (D.f[BS])[kbs];
			real fTS = (D.f[BN])[kb];//kbn
			real fBN = (D.f[TS])[ks];//kts
			//real fZERO = (D.f[dirREST])[k];//kzero
			real fBSW = (D.f[TNE])[k];//ktne
			real fBNE = (D.f[TSW])[ksw];//ktsw
			real fBNW = (D.f[TSE])[ks];//ktse
			real fBSE = (D.f[TNW])[kw];//ktnw
			real fTSW = (D.f[BNE])[kb];//kbne
			real fTNE = (D.f[BSW])[kbsw];
			real fTNW = (D.f[BSE])[kbs];//kbse
			real fTSE = (D.f[BNW])[kbw];//kbnw
										   //real fE    =  (D.f[E   ])[k  ];//ke
										   //real fW    =  (D.f[W   ])[kw ];
										   //real fN    =  (D.f[N   ])[k  ];//kn
										   //real fS    =  (D.f[S   ])[ks ];
										   //real fT    =  (D.f[T   ])[k  ];//kt
										   //real fB    =  (D.f[B   ])[kb ];
										   //real fNE   =  (D.f[NE  ])[k  ];//kne
										   //real fSW   =  (D.f[SW  ])[ksw];
										   //real fSE   =  (D.f[SE  ])[ks ];//kse
										   //real fNW   =  (D.f[NW  ])[kw ];//knw
										   //real fTE   =  (D.f[TE  ])[k  ];//kte
										   //real fBW   =  (D.f[BW  ])[kbw];
										   //real fBE   =  (D.f[BE  ])[kb ];//kbe
										   //real fTW   =  (D.f[TW  ])[kw ];//ktw
										   //real fTN   =  (D.f[TN  ])[k  ];//ktn
										   //real fBS   =  (D.f[BS  ])[kbs];
										   //real fBN   =  (D.f[BN  ])[kb ];//kbn
										   //real fTS   =  (D.f[TS  ])[ks ];//kts
										   //real fZERO =  (D.f[dirREST])[k  ];//kzero
										   //real fTNE   = (D.f[TNE ])[k  ];//ktne
										   //real fTSW   = (D.f[TSW ])[ksw];//ktsw
										   //real fTSE   = (D.f[TSE ])[ks ];//ktse
										   //real fTNW   = (D.f[TNW ])[kw ];//ktnw
										   //real fBNE   = (D.f[BNE ])[kb ];//kbne
										   //real fBSW   = (D.f[BSW ])[kbsw];
										   //real fBSE   = (D.f[BSE ])[kbs];//kbse
										   //real fBNW   = (D.f[BNW ])[kbw];//kbnw
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