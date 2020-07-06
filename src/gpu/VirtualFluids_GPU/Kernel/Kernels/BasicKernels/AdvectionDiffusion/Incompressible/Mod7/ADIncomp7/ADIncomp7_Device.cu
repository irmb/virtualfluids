#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include "Core/RealConstants.h"
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
				D.f[dirZERO] = &DDStart[dirZERO*size_Mat];
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
				D.f[dirZERO] = &DDStart[dirZERO*size_Mat];
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
			real fZERO = (D.f[dirZERO])[k];//kzero
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
										   //real fZERO =  (D.f[dirZERO])[k  ];//kzero
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
			real vx = ((fTNE - fBSW) + (fBNE - fTSW) + (fTSE - fBNW) + (fBSE - fTNW) + (fNE - fSW) + (fSE - fNW) + (fTE - fBW) + (fBE - fTW) + (fE - fW));
			real vy = ((fTNE - fBSW) + (fBNE - fTSW) + (fBNW - fTSE) + (fTNW - fBSE) + (fNE - fSW) + (fNW - fSE) + (fTN - fBS) + (fBN - fTS) + (fN - fS));
			real vz = ((fTNE - fBSW) + (fTSW - fBNE) + (fTSE - fBNW) + (fTNW - fBSE) + (fTE - fBW) + (fTW - fBE) + (fTN - fBS) + (fTS - fBN) + (fT - fB));
			////dörrrrrty !!!!!!!!!!!!!
			//      real vx     =  ten * ((fTNE-fBSW)+(fBNE-fTSW)+(fTSE-fBNW)+(fBSE-fTNW) +(fNE-fSW)+(fSE-fNW)+(fTE-fBW)+(fBE-fTW)+(fE-fW));
			//      real vy     =  ten * ((fTNE-fBSW)+(fBNE-fTSW)+(fBNW-fTSE)+(fTNW-fBSE) +(fNE-fSW)+(fNW-fSE)+(fTN-fBS)+(fBN-fTS)+(fN-fS));
			//      real vz     =  ten * ((fTNE-fBSW)+(fTSW-fBNE)+(fTSE-fBNW)+(fTNW-fBSE) +(fTE-fBW)+(fTW-fBE)+(fTN-fBS)+(fTS-fBN)+(fT-fB));
			////////////////////////////////////////////////////////////////////////////////
			real ux_sq = vx * vx;
			real uy_sq = vy * vy;
			real uz_sq = vz * vz;
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
			real cs2 = c1o4;
			real Lam = diffusivity*c4o1;//diffusivity/(one)/cs2;
			real omegaD = -c1o1 / (Lam + c1o2);
			real ae = c0o1;
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