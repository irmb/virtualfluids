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
#include "Core/RealConstants.h"
#include "math.h"


////////////////////////////////////////////////////////////////////////////////
inline __device__ void forwardInverseChimeraWithK(real &mfa, real &mfb, real &mfc, real vv, real v2, real Kinverse, real K) {
	real m2 = mfa + mfc;
	real m1 = mfc - mfa;
	real m0 = m2 + mfb;
	mfa = m0;
	m0 *= Kinverse;
	m0 += c1o1;
	mfb = (m1*Kinverse - m0 * vv) * K;
	mfc = ((m2 - c2o1*	m1 * vv)*Kinverse + v2 * m0) * K;
}

inline __device__ void backwardInverseChimeraWithK(real &mfa, real &mfb, real &mfc, real vv, real v2, real Kinverse, real K) {
	real m0 = (((mfc - mfb) * c1o2 + mfb *  vv)*Kinverse + (mfa*Kinverse + c1o1) * (v2 - vv) * c1o2) * K;
	real m1 = (((mfa - mfc) -  c2o1 * mfb *  vv)*Kinverse + (mfa*Kinverse + c1o1) * (           -v2)) * K;
	mfc     = (((mfc + mfb) * c1o2 + mfb *  vv)*Kinverse + (mfa*Kinverse + c1o1) * (v2 + vv) * c1o2) * K;
	mfa = m0;
	mfb = m1;
}
////////////////////////////////////////////////////////////////////////////////





inline __device__ void forwardChimeraWithK(real &mfa, real &mfb, real &mfc, real vv, real v2, real K) {

	real m2 = mfa + mfc;
	real m1 = mfc - mfa;
	real m0 = m2 + mfb;
	mfa = m0;
	//m0     += K;
	mfb = (m1 - K*vv) - m0 * vv;
	mfc = ((m2 - c2o1*	m1 * vv) + v2*K) + v2 * m0;
	//m0 += K;
	//mfb = m1 - m0 * vv;
	//mfc = m2 - two*	m1 * vv + v2 * m0;
}

inline __device__ void forwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2) {
	real m1 = (mfa + mfc) + mfb;
	real m2 = mfc - mfa;
	mfc = (mfc + mfa) + (v2*m1 - c2o1*vv*m2);
	mfb = m2 - vv*m1;
	mfa = m1;
}


inline __device__ void backwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2) {
	real ma = (mfc + mfa*(v2 - vv))*c1o2 + mfb*(vv - c1o2);
	real mb = ((mfa - mfc) - mfa*v2) - c2o1*mfb*vv;
	mfc = (mfc + mfa*(v2 + vv))*c1o2 + mfb*(vv + c1o2);
	mfb = mb;
	mfa = ma;
}


inline __device__ void backwardChimeraWithK(real &mfa, real &mfb, real &mfc, real vv, real v2, real K) {
	real  m0 = (mfc - mfb)* c1o2 + mfb * (vv)+(mfa + K) * (v2 - vv) * c1o2;
	real m1 = (mfa - mfc) - c2o1* mfb * vv + (mfa + K) * (-v2);
	mfc = (mfc + mfb)* c1o2 + mfb * (vv)+(mfa + K) * (v2 + vv) * c1o2;
	mfa = m0;
	mfb = m1;

}






////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void Cumulant_One_preconditioned_errorDiffusion_chim_Comp_SP_27(
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	int level,
	real* forces,
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

		if (BC >= GEO_FLUID/*(BC != GEO_SOLID) && (BC != GEO_VOID)*/)
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

			////////////////////////////////////////////////////////////////////////////////
			//index
			//unsigned int kzero= k;
			//unsigned int ke   = k;
			unsigned int kw = neighborX[k];
			//unsigned int kn   = k;
			unsigned int ks = neighborY[k];
			//unsigned int kt   = k;
			unsigned int kb = neighborZ[k];
			unsigned int ksw = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			unsigned int kbw = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			unsigned int kbs = neighborZ[ks];
			//unsigned int ktn  = k;
			//unsigned int kbn  = kb;
			//unsigned int kts  = ks;
			//unsigned int ktse = ks;
			//unsigned int kbnw = kbw;
			//unsigned int ktnw = kw;
			//unsigned int kbse = kbs;
			//unsigned int ktsw = ksw;
			//unsigned int kbne = kb;
			//unsigned int ktne = k;
			unsigned int kbsw = neighborZ[ksw];



			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[dirE])[k];//[ke   ];// +  c2over27 ;(D.f[dirE   ])[k  ];//ke
			real mfabb = (D.f[dirW])[kw];//[kw   ];// +  c2over27 ;(D.f[dirW   ])[kw ];
			real mfbcb = (D.f[dirN])[k];//[kn   ];// +  c2over27 ;(D.f[dirN   ])[k  ];//kn
			real mfbab = (D.f[dirS])[ks];//[ks   ];// +  c2over27 ;(D.f[dirS   ])[ks ];
			real mfbbc = (D.f[dirT])[k];//[kt   ];// +  c2over27 ;(D.f[dirT   ])[k  ];//kt
			real mfbba = (D.f[dirB])[kb];//[kb   ];// +  c2over27 ;(D.f[dirB   ])[kb ];
			real mfccb = (D.f[dirNE])[k];//[kne  ];// +  c1over54 ;(D.f[dirNE  ])[k  ];//kne
			real mfaab = (D.f[dirSW])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dirSW  ])[ksw];
			real mfcab = (D.f[dirSE])[ks];//[kse  ];// +  c1over54 ;(D.f[dirSE  ])[ks ];//kse
			real mfacb = (D.f[dirNW])[kw];//[knw  ];// +  c1over54 ;(D.f[dirNW  ])[kw ];//knw
			real mfcbc = (D.f[dirTE])[k];//[kte  ];// +  c1over54 ;(D.f[dirTE  ])[k  ];//kte
			real mfaba = (D.f[dirBW])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dirBW  ])[kbw];
			real mfcba = (D.f[dirBE])[kb];//[kbe  ];// +  c1over54 ;(D.f[dirBE  ])[kb ];//kbe
			real mfabc = (D.f[dirTW])[kw];//[ktw  ];// +  c1over54 ;(D.f[dirTW  ])[kw ];//ktw
			real mfbcc = (D.f[dirTN])[k];//[ktn  ];// +  c1over54 ;(D.f[dirTN  ])[k  ];//ktn
			real mfbaa = (D.f[dirBS])[kbs];//[kbs  ];// +  c1over54 ;(D.f[dirBS  ])[kbs];
			real mfbca = (D.f[dirBN])[kb];//[kbn  ];// +  c1over54 ;(D.f[dirBN  ])[kb ];//kbn
			real mfbac = (D.f[dirTS])[ks];//[kts  ];// +  c1over54 ;(D.f[dirTS  ])[ks ];//kts
			real mfbbb = (D.f[dirZERO])[k];//[kzero];// +  c8over27 ;(D.f[dirZERO])[k  ];//kzero
			real mfccc = (D.f[dirTNE])[k];//[ktne ];// +  c1over216;(D.f[dirTNE ])[k  ];//ktne
			real mfaac = (D.f[dirTSW])[ksw];//[ktsw ];// +  c1over216;(D.f[dirTSW ])[ksw];//ktsw
			real mfcac = (D.f[dirTSE])[ks];//[ktse ];// +  c1over216;(D.f[dirTSE ])[ks ];//ktse
			real mfacc = (D.f[dirTNW])[kw];//[ktnw ];// +  c1over216;(D.f[dirTNW ])[kw ];//ktnw
			real mfcca = (D.f[dirBNE])[kb];//[kbne ];// +  c1over216;(D.f[dirBNE ])[kb ];//kbne
			real mfaaa = (D.f[dirBSW])[kbsw];//[kbsw ];// +  c1over216;(D.f[dirBSW ])[kbsw];
			real mfcaa = (D.f[dirBSE])[kbs];//[kbse ];// +  c1over216;(D.f[dirBSE ])[kbs];//kbse
			real mfaca = (D.f[dirBNW])[kbw];//[kbnw ];// +  c1over216;(D.f[dirBNW ])[kbw];//kbnw
											   ////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
				((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

			real rho = c1o1 + drho;
			////////////////////////////////////////////////////////////////////////////////////
			real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
				(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
				(mfcbb - mfabb)) / rho;
			real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
				(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
				(mfbcb - mfbab)) / rho;
			real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
				(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
				(mfbbc - mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			//the force be with you
			real fx = forces[0] / (pow((double)c2o1, (double)level)); //zero;//0.0032653/(pow(two,level)); //0.000000005;//(two/1600000.0) / 120.0; //
			real fy = forces[1] / (pow((double)c2o1, (double)level)); //zero;
			real fz = forces[2] / (pow((double)c2o1, (double)level)); //zero;
			vvx += fx*c1o2;
			vvy += fy*c1o2;
			vvz += fz*c1o2;
			////////////////////////////////////////////////////////////////////////////////////
			//real omega = omega_in;
			////////////////////////////////////////////////////////////////////////////////////
			//fast
			real oMdrho = c1o1; // comp special
			real m0, m1, m2;
			real vx2;
			real vy2;
			real vz2;
			vx2 = vvx*vvx;
			vy2 = vvy*vvy;
			vz2 = vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real wadjust;
			real qudricLimitP = c1o100;// * 0.0001f;
			real qudricLimitM = c1o100;// * 0.0001f;
			real qudricLimitD = c1o100;// * 0.001f;
			//real s9 = minusomega;
			//test
			//s9 = 0.;


			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real EQcbb = c0o1;
			real EQabb = c0o1;
			real EQbcb = c0o1;
			real EQbab = c0o1;
			real EQbbc = c0o1;
			real EQbba = c0o1;
			real EQccb = c0o1;
			real EQaab = c0o1;
			real EQcab = c0o1;
			real EQacb = c0o1;
			real EQcbc = c0o1;
			real EQaba = c0o1;
			real EQcba = c0o1;
			real EQabc = c0o1;
			real EQbcc = c0o1;
			real EQbaa = c0o1;
			real EQbca = c0o1;
			real EQbac = c0o1;
			real EQbbb = c0o1;
			real EQccc = drho * c1o27;
			real EQaac = drho * c1o3;
			real EQcac = drho * c1o9;
			real EQacc = drho * c1o9;
			real EQcca = drho * c1o9;
			real EQaaa = drho;
			real EQcaa = drho * c1o3;
			real EQaca = drho * c1o3;
			////////////////////////////////////////////////////////////////////////////////////
			backwardChimeraWithK(EQaaa, EQaab, EQaac, vvz, vz2, c1o1);
			backwardChimeraWithK(EQaca, EQacb, EQacc, vvz, vz2, c1o3);
			///////////////////////////////////////////////////////////
			EQcaa = EQaca; EQcab = EQacb; EQcac = EQacc;
			///////////////////////////////////////////////////////////
			backwardChimeraWithK(EQcca, EQccb, EQccc, vvz, vz2, c1o9);

			backwardChimeraWithK(EQaaa, EQaba, EQaca, vvy, vy2, c1o6);
			backwardChimeraWithK(EQaab, EQabb, EQacb, vvy, vy2, c2o3);
			backwardChimeraWithK(EQaac, EQabc, EQacc, vvy, vy2, c1o6);
			backwardChimeraWithK(EQcaa, EQcba, EQcca, vvy, vy2, c1o18);
			backwardChimeraWithK(EQcab, EQcbb, EQccb, vvy, vy2, c2o9);
			backwardChimeraWithK(EQcac, EQcbc, EQccc, vvy, vy2, c1o18);

			backwardChimeraWithK(EQaaa, EQbaa, EQcaa, vvx, vx2, c1o36);
			backwardChimeraWithK(EQaab, EQbab, EQcab, vvx, vx2, c1o9);
			backwardChimeraWithK(EQaac, EQbac, EQcac, vvx, vx2, c1o36);
			backwardChimeraWithK(EQaba, EQbba, EQcba, vvx, vx2, c1o9);
			backwardChimeraWithK(EQabb, EQbbb, EQcbb, vvx, vx2, c4o9);
			backwardChimeraWithK(EQabc, EQbbc, EQcbc, vvx, vx2, c1o9);
			backwardChimeraWithK(EQaca, EQbca, EQcca, vvx, vx2, c1o36);
			backwardChimeraWithK(EQacb, EQbcb, EQccb, vvx, vx2, c1o9);
			backwardChimeraWithK(EQacc, EQbcc, EQccc, vvx, vx2, c1o36);

			////////////////////////////////////////////////////////////////////////////////////
			//Pre-condition
			mfcbb -= EQcbb;
			mfabb -= EQabb;
			mfbcb -= EQbcb;
			mfbab -= EQbab;
			mfbbc -= EQbbc;
			mfbba -= EQbba;
			mfccb -= EQccb;
			mfaab -= EQaab;
			mfcab -= EQcab;
			mfacb -= EQacb;
			mfcbc -= EQcbc;
			mfaba -= EQaba;
			mfcba -= EQcba;
			mfabc -= EQabc;
			mfbcc -= EQbcc;
			mfbaa -= EQbaa;
			mfbca -= EQbca;
			mfbac -= EQbac;
			mfbbb -= EQbbb;
			mfccc -= EQccc;
			mfaac -= EQaac;
			mfcac -= EQcac;
			mfacc -= EQacc;
			mfcca -= EQcca;
			mfaaa -= EQaaa;
			mfcaa -= EQcaa;
			mfaca -= EQaca;

			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			forwardChimera(mfaaa, mfaab, mfaac, vvz, vz2);
			forwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
			forwardChimera(mfaca, mfacb, mfacc, vvz, vz2);
			forwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
			forwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
			forwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
			forwardChimera(mfcaa, mfcab, mfcac, vvz, vz2);
			forwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
			forwardChimera(mfcca, mfccb, mfccc, vvz, vz2);

			forwardChimera(mfaaa, mfaba, mfaca, vvy, vy2);
			forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
			forwardChimera(mfaac, mfabc, mfacc, vvy, vy2);
			forwardChimera(mfbaa, mfbba, mfbca, vvy, vy2);
			forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
			forwardChimera(mfbac, mfbbc, mfbcc, vvy, vy2);
			forwardChimera(mfcaa, mfcba, mfcca, vvy, vy2);
			forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
			forwardChimera(mfcac, mfcbc, mfccc, vvy, vy2);

			forwardChimera(mfaaa, mfbaa, mfcaa, vvx, vx2);
			forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
			forwardChimera(mfaac, mfbac, mfcac, vvx, vx2);
			forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
			forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
			forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
			forwardChimera(mfaca, mfbca, mfcca, vvx, vx2);
			forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
			forwardChimera(mfacc, mfbcc, mfccc, vvx, vx2);

			//////////////////////////////////////////////////////////////////////////////////////
			////Hin
			//////////////////////////////////////////////////////////////////////////////////////
			//// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			//////////////////////////////////////////////////////////////////////////////////////
			//// Z - Dir
			//forwardChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2, c4o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			//////////////////////////////////////////////////////////////////////////////////////
			//// Y - Dir
			//forwardChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, c1o6);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c1o18);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2, c2o3);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2, c2o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, c1o6);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c1o18);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			//////////////////////////////////////////////////////////////////////////////////////
			//// X - Dir
			//forwardChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, one);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c1o3);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c1o3);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = c1o1; //omega; // one;	//set the bulk viscosity one is high / two is very low and zero is (too) high

			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			real OxyyPxzz = c1o1;
			real OxyyMxzz = c1o1;
			real Oxyz = c1o1;
			////////////////////////////////////////////////////////////
			//4.
			//////////////////////////////
			real O4 = c1o1;
			////////////////////////////////////////////////////////////
			//5.
			//////////////////////////////
			real O5 = c1o1;
			////////////////////////////////////////////////////////////
			//6.
			//////////////////////////////
			real O6 = c1o1;
			////////////////////////////////////////////////////////////


			//central moments to cumulants
			//4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho;
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho;

			real CUMcca = mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
			real CUMcac = mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
			real CUMacc = mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

			//5.
			real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)) / rho;
			real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)) / rho;
			real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)) / rho;

			//6.

			real CUMccc = mfccc + ((-c4o1 *  mfbbb * mfbbb
				- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				- c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				- c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
				+ (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
					+ c2o1 * (mfcaa * mfaca * mfaac)
					+ c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
				- c1o3 * (mfacc + mfcac + mfcca) / rho
				- c1o9 * (mfcaa + mfaca + mfaac) / rho
				+ (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
					+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3
				+ c1o27*((drho * drho - drho) / (rho*rho)));




			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy = mfcaa - mfaca;
			real mxxMzz = mfcaa - mfaac;

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
			{
				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
				real dyuy = dxux + omega * c3o2 * mxxMyy;
				real dzuz = dxux + omega * c3o2 * mxxMzz;

				//relax
				mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
				mxxMyy += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
				mxxMzz += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

			}
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////no correction
			//mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz);//-magicBulk*OxxPyyPzz;
			//mxxMyy += -(-omega) * (-mxxMyy);
			//mxxMzz += -(-omega) * (-mxxMzz);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb += omega * (-mfabb);
			mfbab += omega * (-mfbab);
			mfbba += omega * (-mfbba);

			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-c2o1*  mxxMyy + mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (mxxMyy - c2o1* mxxMzz + mxxPyyPzz);

			//3.
			// linear combinations

			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

			//relax
			//////////////////////////////////////////////////////////////////////////
			mfbbb += OxyyMxzz * (-mfbbb);
			mxxyPyzz += OxyyPxzz * (-mxxyPyzz);
			mxxyMyzz += OxyyMxzz * (-mxxyMyzz);
			mxxzPyyz += OxyyPxzz * (-mxxzPyyz);
			mxxzMyyz += OxyyMxzz * (-mxxzMyyz);
			mxyyPxzz += OxyyPxzz * (-mxyyPxzz);
			mxyyMxzz += OxyyMxzz * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////

			mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			//////////////////////////////////////////////////////////////////////////
			CUMacc += O4 * (-CUMacc);
			CUMcac += O4 * (-CUMcac);
			CUMcca += O4 * (-CUMcca);

			CUMbbc += O4 * (-CUMbbc);
			CUMbcb += O4 * (-CUMbcb);
			CUMcbb += O4 * (-CUMcbb);
			//////////////////////////////////////////////////////////////////////////


			//5.
			CUMbcc += O5 * (-CUMbcc);
			CUMcbc += O5 * (-CUMcbc);
			CUMccb += O5 * (-CUMccb);

			//6.
			CUMccc += O6 * (-CUMccc);



			//back cumulants to central moments
			//4.
			mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;
			mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho;
			mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho;

			mfcca = CUMcca + (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
			mfcac = CUMcac + (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
			mfacc = CUMacc + (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

			//5.
			mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)) / rho;
			mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)) / rho;
			mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)) / rho;

			//6.
			mfccc = CUMccc - ((-c4o1 *  mfbbb * mfbbb
				- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				- c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				- c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
				+ (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
					+ c2o1 * (mfcaa * mfaca * mfaac)
					+ c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
				- c1o3 * (mfacc + mfcac + mfcca) / rho
				- c1o9 * (mfcaa + mfaca + mfaac) / rho
				+ (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
					+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3
				+ c1o27*((drho * drho - drho) / (rho*rho)));

			////////////////////////////////////////////////////////////////////////////////////
			//the force be with you
			mfbaa = -mfbaa;
			mfaba = -mfaba;
			mfaab = -mfaab;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			//back
			////////////////////////////////////////////////////////////////////////////////////
			backwardChimera(mfaaa, mfaab, mfaac, vvz, vz2);
			backwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
			backwardChimera(mfaca, mfacb, mfacc, vvz, vz2);
			backwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
			backwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
			backwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
			backwardChimera(mfcaa, mfcab, mfcac, vvz, vz2);
			backwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
			backwardChimera(mfcca, mfccb, mfccc, vvz, vz2);

			backwardChimera(mfaaa, mfaba, mfaca, vvy, vy2);
			backwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
			backwardChimera(mfaac, mfabc, mfacc, vvy, vy2);
			backwardChimera(mfbaa, mfbba, mfbca, vvy, vy2);
			backwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
			backwardChimera(mfbac, mfbbc, mfbcc, vvy, vy2);
			backwardChimera(mfcaa, mfcba, mfcca, vvy, vy2);
			backwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
			backwardChimera(mfcac, mfcbc, mfccc, vvy, vy2);

			backwardChimera(mfaaa, mfbaa, mfcaa, vvx, vx2);
			backwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
			backwardChimera(mfaac, mfbac, mfcac, vvx, vx2);
			backwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
			backwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
			backwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
			backwardChimera(mfaca, mfbca, mfcca, vvx, vx2);
			backwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
			backwardChimera(mfacc, mfbcc, mfccc, vvx, vx2);

			////////////////////////////////////////////////////////////////////////////////////
			//mfcbb += EQcbb;
			//mfabb += EQabb;
			//mfbcb += EQbcb;
			//mfbab += EQbab;
			//mfbbc += EQbbc;
			//mfbba += EQbba;
			//mfccb += EQccb;
			//mfaab += EQaab;
			//mfcab += EQcab;
			//mfacb += EQacb;
			//mfcbc += EQcbc;
			//mfaba += EQaba;
			//mfcba += EQcba;
			//mfabc += EQabc;
			//mfbcc += EQbcc;
			//mfbaa += EQbaa;
			//mfbca += EQbca;
			//mfbac += EQbac;
			//mfbbb += EQbbb;
			//mfccc += EQccc;
			//mfaac += EQaac;
			//mfcac += EQcac;
			//mfacc += EQacc;
			//mfcca += EQcca;
			//mfaaa += EQaaa;
			//mfcaa += EQcaa;
			//mfaca += EQaca;
			////////////////////////////////////////////////////////////////////////////////////
			////Error diffusion
			real fTEMP = mfbbb + EQbbb;
			real delta0 = mfbbb - (fTEMP - EQbbb);
			delta0 *= c1o4;
			mfbbb = fTEMP;


			fTEMP = mfcbb + EQcbb;
			real deltacbb = mfcbb - (fTEMP - EQcbb);
			mfcbb = fTEMP;
			//mfcbb+=EQcbb;

			fTEMP = mfabb + EQabb;
			real deltaabb = mfabb - (fTEMP - EQabb);
			mfabb = fTEMP;
			//mfabb+=EQabb;

			fTEMP = mfbcb + EQbcb;
			real deltabcb = mfbcb - (fTEMP - EQbcb);
			mfbcb = fTEMP;
			//mfbcb+=EQbcb;

			fTEMP = mfbab + EQbab;
			real deltabab = mfbab - (fTEMP - EQbab);
			mfbab = fTEMP;
			//mfbab+=EQbab;

			fTEMP = mfbbc + EQbbc;
			real deltabbc = mfbbc - (fTEMP - EQbbc);
			mfbbc = fTEMP;
			//mfbbc+=EQbbc;

			fTEMP = mfbba + EQbba;
			real deltabba = mfbba - (fTEMP - EQbba);
			mfbba = fTEMP;
			//mfbba+=EQbba;

			EQccb += (delta0 + c1o2*(deltacbb + deltabcb));
			fTEMP = mfccb + EQccb;
			real deltaccb = mfccb - (fTEMP - EQccb);
			mfccb = fTEMP;
			//mfccb+=EQccb+(delta0+c1o2*(deltacbb+deltabcb));

			EQaab += (delta0 + c1o2*(deltaabb + deltabab));
			fTEMP = mfaab + EQaab;
			real deltaaab = mfaab - (fTEMP - EQaab);
			mfaab = fTEMP;
			//mfaab+=EQaab+(delta0+c1o2*(deltaabb+deltabab));

			EQcab += (delta0 + c1o2*(deltacbb + deltabab));
			fTEMP = mfcab + EQcab;
			real deltacab = mfcab - (fTEMP - EQcab);
			mfcab = fTEMP;
			//mfcab+=EQcab+(delta0+c1o2*(deltacbb+deltabab));

			EQacb += (delta0 + c1o2*(deltaabb + deltabcb));
			fTEMP = mfacb + EQacb;
			real deltaacb = mfacb - (fTEMP - EQacb);
			mfacb = fTEMP;
			//mfacb+=EQacb+(delta0+c1o2*(deltaabb+deltabcb));

			EQcbc += (delta0 + c1o2*(deltacbb + deltabbc));
			fTEMP = mfcbc + EQcbc;
			real deltacbc = mfcbc - (fTEMP - EQcbc);
			mfcbc = fTEMP;
			//mfcbc+=EQcbc+(delta0+c1o2*(deltacbb+deltabbc));

			EQaba += (delta0 + c1o2*(deltaabb + deltabba));
			fTEMP = mfaba + EQaba;
			real deltaaba = mfaba - (fTEMP - EQaba);
			mfaba = fTEMP;
			//mfaba+=EQaba+(delta0+c1o2*(deltaabb+deltabba));

			EQcba += (delta0 + c1o2*(deltacbb + deltabba));
			fTEMP = mfcba + EQcba;
			real deltacba = mfcba - (fTEMP - EQcba);
			mfcba = fTEMP;
			//mfcba+=EQcba+(delta0+c1o2*(deltacbb+deltabba));

			EQabc += (delta0 + c1o2*(deltaabb + deltabbc));
			fTEMP = mfabc + EQabc;
			real deltaabc = mfabc - (fTEMP - EQabc);
			mfabc = fTEMP;
			//mfabc+=EQabc+(delta0+c1o2*(deltaabb+deltabbc));

			EQbcc += (delta0 + c1o2*(deltabcb + deltabbc));
			fTEMP = mfbcc + EQbcc;
			real deltabcc = mfbcc - (fTEMP - EQbcc);
			mfbcc = fTEMP;
			//mfbcc+=EQbcc+(delta0+c1o2*(deltabcb+deltabbc));

			EQbaa += (delta0 + c1o2*(deltabab + deltabba));
			fTEMP = mfbaa + EQbaa;
			real deltabaa = mfbaa - (fTEMP - EQbaa);
			mfbaa = fTEMP;
			//mfbaa+=EQbaa+(delta0+c1o2*(deltabab+deltabba));

			EQbca += (delta0 + c1o2*(deltabcb + deltabba));
			fTEMP = mfbca + EQbca;
			real deltabca = mfbca - (fTEMP - EQbca);
			mfbca = fTEMP;
			//mfbca+=EQbca+(delta0+c1o2*(deltabcb+deltabba));

			EQbac += (delta0 + c1o2*(deltabab + deltabbc));
			fTEMP = mfbac + EQbac;
			real deltabac = mfbac - (fTEMP - EQbac);
			mfbac = fTEMP;
			//mfbac+=EQbac+(delta0+c1o2*(deltabab+deltabbc));

			mfccc += EQccc - (delta0 + c1o4*(deltacbb + deltabcb + deltabbc) - c1o2*(deltabcc + deltacbc + deltaccb));
			mfaac += EQaac - (delta0 + c1o4*(deltaabb + deltabab + deltabbc) - c1o2*(deltabac + deltaabc + deltaaab));
			mfcac += EQcac - (delta0 + c1o4*(deltacbb + deltabab + deltabbc) - c1o2*(deltabac + deltacbc + deltacab));
			mfacc += EQacc - (delta0 + c1o4*(deltaabb + deltabcb + deltabbc) - c1o2*(deltabcc + deltaabc + deltaacb));
			mfcca += EQcca - (delta0 + c1o4*(deltacbb + deltabcb + deltabba) - c1o2*(deltabca + deltacba + deltaccb));
			mfaaa += EQaaa - (delta0 + c1o4*(deltaabb + deltabab + deltabba) - c1o2*(deltabaa + deltaaba + deltaaab));
			mfcaa += EQcaa - (delta0 + c1o4*(deltacbb + deltabab + deltabba) - c1o2*(deltabaa + deltacba + deltacab));
			mfaca += EQaca - (delta0 + c1o4*(deltaabb + deltabcb + deltabba) - c1o2*(deltabca + deltaaba + deltaacb));



			//////////////////////////////////////////////////////////////////////////////////////
			////back
			//////////////////////////////////////////////////////////////////////////////////////
			////mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
			//////////////////////////////////////////////////////////////////////////////////////
			//// Z - Dir
			//backwardChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, one);
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c1o3);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
			///////////b//////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
			///////////b//////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c1o3);
			///////////c//////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
			///////////c//////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			////mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			//////////////////////////////////////////////////////////////////////////////////////
			//// Y - Dir
			//backwardChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, c1o6);
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaab, mfabb, mfacb, vvy, vy2, c2o3);
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c1o6);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfbaa, mfbba, mfbca, vvz, vz2);
			///////////b//////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfbab, mfbbb, mfbcb, vvz, vz2);
			///////////b//////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfbac, mfbbc, mfbcc, vvz, vz2);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, c1o18);
			///////////c//////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfcab, mfcbb, mfccb, vvy, vy2, c2o9);
			///////////c//////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c1o18);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			////mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			//////////////////////////////////////////////////////////////////////////////////////
			//// X - Dir
			//backwardChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaba, mfbba, mfcba, vvx, vx2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaab, mfbab, mfcab, vvx, vx2, c1o9);
			/////////////b////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfabb, mfbbb, mfcbb, vvx, vx2, c4o9);
			/////////////b////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfacb, mfbcb, mfccb, vvx, vx2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c1o36);
			/////////////c////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfabc, mfbbc, mfcbc, vvx, vx2, c1o9);
			/////////////c////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////////
			//real drhoPost =
			//	((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
			//	(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
			//		((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
			//mfbbb += drho - drhoPost;
			////////////////////////////////////////////////////////////////////////////////////
			(D.f[dirE])[k] = mfabb;//(D.f[ dirE   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dirE   ])[k   ]                                                                     
			(D.f[dirW])[kw] = mfcbb;//(D.f[ dirW   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dirW   ])[kw  ]                                                                   
			(D.f[dirN])[k] = mfbab;//(D.f[ dirN   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ dirN   ])[k   ]
			(D.f[dirS])[ks] = mfbcb;//(D.f[ dirS   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ dirS   ])[ks  ]
			(D.f[dirT])[k] = mfbba;//(D.f[ dirT   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ dirT   ])[k   ]
			(D.f[dirB])[kb] = mfbbc;//(D.f[ dirB   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ dirB   ])[kb  ]
			(D.f[dirNE])[k] = mfaab;//(D.f[ dirNE  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dirNE  ])[k   ]
			(D.f[dirSW])[ksw] = mfccb;//(D.f[ dirSW  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dirSW  ])[ksw ]
			(D.f[dirSE])[ks] = mfacb;//(D.f[ dirSE  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dirSE  ])[ks  ]
			(D.f[dirNW])[kw] = mfcab;//(D.f[ dirNW  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dirNW  ])[kw  ]
			(D.f[dirTE])[k] = mfaba;//(D.f[ dirTE  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dirTE  ])[k   ]
			(D.f[dirBW])[kbw] = mfcbc;//(D.f[ dirBW  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dirBW  ])[kbw ]
			(D.f[dirBE])[kb] = mfabc;//(D.f[ dirBE  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dirBE  ])[kb  ]
			(D.f[dirTW])[kw] = mfcba;//(D.f[ dirTW  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dirTW  ])[kw  ]
			(D.f[dirTN])[k] = mfbaa;//(D.f[ dirTN  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ dirTN  ])[k   ]
			(D.f[dirBS])[kbs] = mfbcc;//(D.f[ dirBS  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ dirBS  ])[kbs ]
			(D.f[dirBN])[kb] = mfbac;//(D.f[ dirBN  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ dirBN  ])[kb  ]
			(D.f[dirTS])[ks] = mfbca;//(D.f[ dirTS  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ dirTS  ])[ks  ]
			(D.f[dirZERO])[k] = mfbbb;//(D.f[ dirZERO])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ dirZERO])[k   ]
			(D.f[dirTNE])[k] = mfaaa;//(D.f[ dirTNE ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dirTNE ])[k   ]
			(D.f[dirTSE])[ks] = mfaca;//(D.f[ dirTSE ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dirTSE ])[ks  ]
			(D.f[dirBNE])[kb] = mfaac;//(D.f[ dirBNE ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dirBNE ])[kb  ]
			(D.f[dirBSE])[kbs] = mfacc;//(D.f[ dirBSE ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dirBSE ])[kbs ]
			(D.f[dirTNW])[kw] = mfcaa;//(D.f[ dirTNW ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dirTNW ])[kw  ]
			(D.f[dirTSW])[ksw] = mfcca;//(D.f[ dirTSW ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dirTSW ])[ksw ]
			(D.f[dirBNW])[kbw] = mfcac;//(D.f[ dirBNW ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dirBNW ])[kbw ]
			(D.f[dirBSW])[kbsw] = mfccc;//(D.f[ dirBSW ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dirBSW ])[kbsw]
										////////////////////////////////////////////////////////////////////////////////////
		}
	}
}
////////////////////////////////////////////////////////////////////////////////








































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void Cumulant_One_preconditioned_chim_Comp_SP_27(
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	int level,
	real* forces,
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

		if (BC >= GEO_FLUID/*(BC != GEO_SOLID) && (BC != GEO_VOID)*/)
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

			////////////////////////////////////////////////////////////////////////////////
			//index
			//unsigned int kzero= k;
			//unsigned int ke   = k;
			unsigned int kw = neighborX[k];
			//unsigned int kn   = k;
			unsigned int ks = neighborY[k];
			//unsigned int kt   = k;
			unsigned int kb = neighborZ[k];
			unsigned int ksw = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			unsigned int kbw = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			unsigned int kbs = neighborZ[ks];
			//unsigned int ktn  = k;
			//unsigned int kbn  = kb;
			//unsigned int kts  = ks;
			//unsigned int ktse = ks;
			//unsigned int kbnw = kbw;
			//unsigned int ktnw = kw;
			//unsigned int kbse = kbs;
			//unsigned int ktsw = ksw;
			//unsigned int kbne = kb;
			//unsigned int ktne = k;
			unsigned int kbsw = neighborZ[ksw];



			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[dirE])[k];//[ke   ];// +  c2over27 ;(D.f[dirE   ])[k  ];//ke
			real mfabb = (D.f[dirW])[kw];//[kw   ];// +  c2over27 ;(D.f[dirW   ])[kw ];
			real mfbcb = (D.f[dirN])[k];//[kn   ];// +  c2over27 ;(D.f[dirN   ])[k  ];//kn
			real mfbab = (D.f[dirS])[ks];//[ks   ];// +  c2over27 ;(D.f[dirS   ])[ks ];
			real mfbbc = (D.f[dirT])[k];//[kt   ];// +  c2over27 ;(D.f[dirT   ])[k  ];//kt
			real mfbba = (D.f[dirB])[kb];//[kb   ];// +  c2over27 ;(D.f[dirB   ])[kb ];
			real mfccb = (D.f[dirNE])[k];//[kne  ];// +  c1over54 ;(D.f[dirNE  ])[k  ];//kne
			real mfaab = (D.f[dirSW])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dirSW  ])[ksw];
			real mfcab = (D.f[dirSE])[ks];//[kse  ];// +  c1over54 ;(D.f[dirSE  ])[ks ];//kse
			real mfacb = (D.f[dirNW])[kw];//[knw  ];// +  c1over54 ;(D.f[dirNW  ])[kw ];//knw
			real mfcbc = (D.f[dirTE])[k];//[kte  ];// +  c1over54 ;(D.f[dirTE  ])[k  ];//kte
			real mfaba = (D.f[dirBW])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dirBW  ])[kbw];
			real mfcba = (D.f[dirBE])[kb];//[kbe  ];// +  c1over54 ;(D.f[dirBE  ])[kb ];//kbe
			real mfabc = (D.f[dirTW])[kw];//[ktw  ];// +  c1over54 ;(D.f[dirTW  ])[kw ];//ktw
			real mfbcc = (D.f[dirTN])[k];//[ktn  ];// +  c1over54 ;(D.f[dirTN  ])[k  ];//ktn
			real mfbaa = (D.f[dirBS])[kbs];//[kbs  ];// +  c1over54 ;(D.f[dirBS  ])[kbs];
			real mfbca = (D.f[dirBN])[kb];//[kbn  ];// +  c1over54 ;(D.f[dirBN  ])[kb ];//kbn
			real mfbac = (D.f[dirTS])[ks];//[kts  ];// +  c1over54 ;(D.f[dirTS  ])[ks ];//kts
			real mfbbb = (D.f[dirZERO])[k];//[kzero];// +  c8over27 ;(D.f[dirZERO])[k  ];//kzero
			real mfccc = (D.f[dirTNE])[k];//[ktne ];// +  c1over216;(D.f[dirTNE ])[k  ];//ktne
			real mfaac = (D.f[dirTSW])[ksw];//[ktsw ];// +  c1over216;(D.f[dirTSW ])[ksw];//ktsw
			real mfcac = (D.f[dirTSE])[ks];//[ktse ];// +  c1over216;(D.f[dirTSE ])[ks ];//ktse
			real mfacc = (D.f[dirTNW])[kw];//[ktnw ];// +  c1over216;(D.f[dirTNW ])[kw ];//ktnw
			real mfcca = (D.f[dirBNE])[kb];//[kbne ];// +  c1over216;(D.f[dirBNE ])[kb ];//kbne
			real mfaaa = (D.f[dirBSW])[kbsw];//[kbsw ];// +  c1over216;(D.f[dirBSW ])[kbsw];
			real mfcaa = (D.f[dirBSE])[kbs];//[kbse ];// +  c1over216;(D.f[dirBSE ])[kbs];//kbse
			real mfaca = (D.f[dirBNW])[kbw];//[kbnw ];// +  c1over216;(D.f[dirBNW ])[kbw];//kbnw
											   ////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
				((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

			real rho = c1o1 + drho;
			////////////////////////////////////////////////////////////////////////////////////
			real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
				(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
				(mfcbb - mfabb)) / rho;
			real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
				(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
				(mfbcb - mfbab)) / rho;
			real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
				(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
				(mfbbc - mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			//the force be with you
			real fx = forces[0] / (pow((double)c2o1, (double)level)); //zero;//0.0032653/(pow(two,level)); //0.000000005;//(two/1600000.0) / 120.0; //
			real fy = forces[1] / (pow((double)c2o1, (double)level)); //zero;
			real fz = forces[2] / (pow((double)c2o1, (double)level)); //zero;
			vvx += fx*c1o2;
			vvy += fy*c1o2;
			vvz += fz*c1o2;
			////////////////////////////////////////////////////////////////////////////////////
			//real omega = omega_in;
			////////////////////////////////////////////////////////////////////////////////////
			//fast
			real oMdrho = c1o1; // comp special
			real m0, m1, m2;
			real vx2;
			real vy2;
			real vz2;
			vx2 = vvx*vvx;
			vy2 = vvy*vvy;
			vz2 = vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real wadjust;
			real qudricLimitP = c1o100;// * 0.0001f;
			real qudricLimitM = c1o100;// * 0.0001f;
			real qudricLimitD = c1o100;// * 0.001f;
			//real s9 = minusomega;
			//test
			//s9 = 0.;


			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real EQcbb = c0o1;
			real EQabb = c0o1;
			real EQbcb = c0o1;
			real EQbab = c0o1;
			real EQbbc = c0o1;
			real EQbba = c0o1;
			real EQccb = c0o1;
			real EQaab = c0o1;
			real EQcab = c0o1;
			real EQacb = c0o1;
			real EQcbc = c0o1;
			real EQaba = c0o1;
			real EQcba = c0o1;
			real EQabc = c0o1;
			real EQbcc = c0o1;
			real EQbaa = c0o1;
			real EQbca = c0o1;
			real EQbac = c0o1;
			real EQbbb = c0o1;
			real EQccc = drho * c1o27;
			real EQaac = drho * c1o3;
			real EQcac = drho * c1o9;
			real EQacc = drho * c1o9;
			real EQcca = drho * c1o9;
			real EQaaa = drho;
			real EQcaa = drho * c1o3;
			real EQaca = drho * c1o3;
			////////////////////////////////////////////////////////////////////////////////////
			backwardChimeraWithK(EQaaa, EQaab, EQaac, vvz, vz2, c1o1);
			backwardChimeraWithK(EQaca, EQacb, EQacc, vvz, vz2, c1o3);
			///////////////////////////////////////////////////////////
			EQcaa = EQaca; EQcab = EQacb; EQcac = EQacc;
			///////////////////////////////////////////////////////////
			backwardChimeraWithK(EQcca, EQccb, EQccc, vvz, vz2, c1o9);

			backwardChimeraWithK(EQaaa, EQaba, EQaca, vvy, vy2, c1o6);
			backwardChimeraWithK(EQaab, EQabb, EQacb, vvy, vy2, c2o3);
			backwardChimeraWithK(EQaac, EQabc, EQacc, vvy, vy2, c1o6);
			backwardChimeraWithK(EQcaa, EQcba, EQcca, vvy, vy2, c1o18);
			backwardChimeraWithK(EQcab, EQcbb, EQccb, vvy, vy2, c2o9);
			backwardChimeraWithK(EQcac, EQcbc, EQccc, vvy, vy2, c1o18);

			backwardChimeraWithK(EQaaa, EQbaa, EQcaa, vvx, vx2, c1o36);
			backwardChimeraWithK(EQaab, EQbab, EQcab, vvx, vx2, c1o9);
			backwardChimeraWithK(EQaac, EQbac, EQcac, vvx, vx2, c1o36);
			backwardChimeraWithK(EQaba, EQbba, EQcba, vvx, vx2, c1o9);
			backwardChimeraWithK(EQabb, EQbbb, EQcbb, vvx, vx2, c4o9);
			backwardChimeraWithK(EQabc, EQbbc, EQcbc, vvx, vx2, c1o9);
			backwardChimeraWithK(EQaca, EQbca, EQcca, vvx, vx2, c1o36);
			backwardChimeraWithK(EQacb, EQbcb, EQccb, vvx, vx2, c1o9);
			backwardChimeraWithK(EQacc, EQbcc, EQccc, vvx, vx2, c1o36);

			////////////////////////////////////////////////////////////////////////////////////
			//Pre-condition
			mfcbb -= EQcbb;
			mfabb -= EQabb;
			mfbcb -= EQbcb;
			mfbab -= EQbab;
			mfbbc -= EQbbc;
			mfbba -= EQbba;
			mfccb -= EQccb;
			mfaab -= EQaab;
			mfcab -= EQcab;
			mfacb -= EQacb;
			mfcbc -= EQcbc;
			mfaba -= EQaba;
			mfcba -= EQcba;
			mfabc -= EQabc;
			mfbcc -= EQbcc;
			mfbaa -= EQbaa;
			mfbca -= EQbca;
			mfbac -= EQbac;
			mfbbb -= EQbbb;
			mfccc -= EQccc;
			mfaac -= EQaac;
			mfcac -= EQcac;
			mfacc -= EQacc;
			mfcca -= EQcca;
			mfaaa -= EQaaa;
			mfcaa -= EQcaa;
			mfaca -= EQaca;

			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			forwardChimera(mfaaa, mfaab, mfaac, vvz, vz2);
			forwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
			forwardChimera(mfaca, mfacb, mfacc, vvz, vz2);
			forwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
			forwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
			forwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
			forwardChimera(mfcaa, mfcab, mfcac, vvz, vz2);
			forwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
			forwardChimera(mfcca, mfccb, mfccc, vvz, vz2);

			forwardChimera(mfaaa, mfaba, mfaca, vvy, vy2);
			forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
			forwardChimera(mfaac, mfabc, mfacc, vvy, vy2);
			forwardChimera(mfbaa, mfbba, mfbca, vvy, vy2);
			forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
			forwardChimera(mfbac, mfbbc, mfbcc, vvy, vy2);
			forwardChimera(mfcaa, mfcba, mfcca, vvy, vy2);
			forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
			forwardChimera(mfcac, mfcbc, mfccc, vvy, vy2);

			forwardChimera(mfaaa, mfbaa, mfcaa, vvx, vx2);
			forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
			forwardChimera(mfaac, mfbac, mfcac, vvx, vx2);
			forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
			forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
			forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
			forwardChimera(mfaca, mfbca, mfcca, vvx, vx2);
			forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
			forwardChimera(mfacc, mfbcc, mfccc, vvx, vx2);

			//////////////////////////////////////////////////////////////////////////////////////
			////Hin
			//////////////////////////////////////////////////////////////////////////////////////
			//// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			//////////////////////////////////////////////////////////////////////////////////////
			//// Z - Dir
			//forwardChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2, c4o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			//////////////////////////////////////////////////////////////////////////////////////
			//// Y - Dir
			//forwardChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, c1o6);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c1o18);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2, c2o3);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2, c2o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, c1o6);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c1o18);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			//////////////////////////////////////////////////////////////////////////////////////
			//// X - Dir
			//forwardChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, one);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c1o3);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c1o3);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
			//////////////////////////////////////////////////////////////////////////////////////
			//forwardChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = c1o1; //omega; // one;	//set the bulk viscosity one is high / two is very low and zero is (too) high

			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			real OxyyPxzz = c1o1;
			real OxyyMxzz = c1o1;
			real Oxyz = c1o1;
			////////////////////////////////////////////////////////////
			//4.
			//////////////////////////////
			real O4 = c1o1;
			////////////////////////////////////////////////////////////
			//5.
			//////////////////////////////
			real O5 = c1o1;
			////////////////////////////////////////////////////////////
			//6.
			//////////////////////////////
			real O6 = c1o1;
			////////////////////////////////////////////////////////////


			//central moments to cumulants
			//4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho;
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho;

			real CUMcca = mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
			real CUMcac = mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
			real CUMacc = mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

			//5.
			real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)) / rho;
			real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)) / rho;
			real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)) / rho;

			//6.

			real CUMccc = mfccc + ((-c4o1 *  mfbbb * mfbbb
				- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				- c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				- c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
				+ (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
					+ c2o1 * (mfcaa * mfaca * mfaac)
					+ c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
				- c1o3 * (mfacc + mfcac + mfcca) / rho
				- c1o9 * (mfcaa + mfaca + mfaac) / rho
				+ (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
					+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3
				+ c1o27*((drho * drho - drho) / (rho*rho)));




			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy = mfcaa - mfaca;
			real mxxMzz = mfcaa - mfaac;

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
			{
				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
				real dyuy = dxux + omega * c3o2 * mxxMyy;
				real dzuz = dxux + omega * c3o2 * mxxMzz;

				//relax
				mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
				mxxMyy += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
				mxxMzz += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

			}
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////no correction
			//mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz);//-magicBulk*OxxPyyPzz;
			//mxxMyy += -(-omega) * (-mxxMyy);
			//mxxMzz += -(-omega) * (-mxxMzz);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb += omega * (-mfabb);
			mfbab += omega * (-mfbab);
			mfbba += omega * (-mfbba);

			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-c2o1*  mxxMyy + mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (mxxMyy - c2o1* mxxMzz + mxxPyyPzz);

			//3.
			// linear combinations

			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

			//relax
			//////////////////////////////////////////////////////////////////////////
			mfbbb += OxyyMxzz * (-mfbbb);
			mxxyPyzz += OxyyPxzz * (-mxxyPyzz);
			mxxyMyzz += OxyyMxzz * (-mxxyMyzz);
			mxxzPyyz += OxyyPxzz * (-mxxzPyyz);
			mxxzMyyz += OxyyMxzz * (-mxxzMyyz);
			mxyyPxzz += OxyyPxzz * (-mxyyPxzz);
			mxyyMxzz += OxyyMxzz * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////

			mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			//////////////////////////////////////////////////////////////////////////
			CUMacc += O4 * (-CUMacc);
			CUMcac += O4 * (-CUMcac);
			CUMcca += O4 * (-CUMcca);

			CUMbbc += O4 * (-CUMbbc);
			CUMbcb += O4 * (-CUMbcb);
			CUMcbb += O4 * (-CUMcbb);
			//////////////////////////////////////////////////////////////////////////


			//5.
			CUMbcc += O5 * (-CUMbcc);
			CUMcbc += O5 * (-CUMcbc);
			CUMccb += O5 * (-CUMccb);

			//6.
			CUMccc += O6 * (-CUMccc);



			//back cumulants to central moments
			//4.
			mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;
			mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho;
			mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho;

			mfcca = CUMcca + (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
			mfcac = CUMcac + (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
			mfacc = CUMacc + (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

			//5.
			mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)) / rho;
			mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)) / rho;
			mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)) / rho;

			//6.
			mfccc = CUMccc - ((-c4o1 *  mfbbb * mfbbb
				- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				- c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				- c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
				+ (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
					+ c2o1 * (mfcaa * mfaca * mfaac)
					+ c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
				- c1o3 * (mfacc + mfcac + mfcca) / rho
				- c1o9 * (mfcaa + mfaca + mfaac) / rho
				+ (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
					+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3
				+ c1o27*((drho * drho - drho) / (rho*rho)));

			////////////////////////////////////////////////////////////////////////////////////
			//the force be with you
			mfbaa = -mfbaa;
			mfaba = -mfaba;
			mfaab = -mfaab;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			//back
			////////////////////////////////////////////////////////////////////////////////////
			backwardChimera(mfaaa, mfaab, mfaac, vvz, vz2);
			backwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
			backwardChimera(mfaca, mfacb, mfacc, vvz, vz2);
			backwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
			backwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
			backwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
			backwardChimera(mfcaa, mfcab, mfcac, vvz, vz2);
			backwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
			backwardChimera(mfcca, mfccb, mfccc, vvz, vz2);

			backwardChimera(mfaaa, mfaba, mfaca, vvy, vy2);
			backwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
			backwardChimera(mfaac, mfabc, mfacc, vvy, vy2);
			backwardChimera(mfbaa, mfbba, mfbca, vvy, vy2);
			backwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
			backwardChimera(mfbac, mfbbc, mfbcc, vvy, vy2);
			backwardChimera(mfcaa, mfcba, mfcca, vvy, vy2);
			backwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
			backwardChimera(mfcac, mfcbc, mfccc, vvy, vy2);

			backwardChimera(mfaaa, mfbaa, mfcaa, vvx, vx2);
			backwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
			backwardChimera(mfaac, mfbac, mfcac, vvx, vx2);
			backwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
			backwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
			backwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
			backwardChimera(mfaca, mfbca, mfcca, vvx, vx2);
			backwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
			backwardChimera(mfacc, mfbcc, mfccc, vvx, vx2);

			////////////////////////////////////////////////////////////////////////////////////
			mfcbb+=EQcbb;
			mfabb+=EQabb;
			mfbcb+=EQbcb;
			mfbab+=EQbab;
			mfbbc+=EQbbc;
			mfbba+=EQbba;
			mfccb+=EQccb;
			mfaab+=EQaab;
			mfcab+=EQcab;
			mfacb+=EQacb;
			mfcbc+=EQcbc;
			mfaba+=EQaba;
			mfcba+=EQcba;
			mfabc+=EQabc;
			mfbcc+=EQbcc;
			mfbaa+=EQbaa;
			mfbca+=EQbca;
			mfbac+=EQbac;
			mfbbb+=EQbbb;
			mfccc+=EQccc;
			mfaac+=EQaac;
			mfcac+=EQcac;
			mfacc+=EQacc;
			mfcca+=EQcca;
			mfaaa+=EQaaa;
			mfcaa+=EQcaa;
			mfaca+=EQaca;


			//////////////////////////////////////////////////////////////////////////////////////
			////back
			//////////////////////////////////////////////////////////////////////////////////////
			////mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
			//////////////////////////////////////////////////////////////////////////////////////
			//// Z - Dir
			//backwardChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, one);
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c1o3);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
			///////////b//////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
			///////////b//////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c1o3);
			///////////c//////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
			///////////c//////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			////mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			//////////////////////////////////////////////////////////////////////////////////////
			//// Y - Dir
			//backwardChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, c1o6);
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaab, mfabb, mfacb, vvy, vy2, c2o3);
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c1o6);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfbaa, mfbba, mfbca, vvz, vz2);
			///////////b//////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfbab, mfbbb, mfbcb, vvz, vz2);
			///////////b//////////////////////////////////////////////////////////////////////////
			//backwardChimera(mfbac, mfbbc, mfbcc, vvz, vz2);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, c1o18);
			///////////c//////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfcab, mfcbb, mfccb, vvy, vy2, c2o9);
			///////////c//////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c1o18);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			////mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			//////////////////////////////////////////////////////////////////////////////////////
			//// X - Dir
			//backwardChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaba, mfbba, mfcba, vvx, vx2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaab, mfbab, mfcab, vvx, vx2, c1o9);
			/////////////b////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfabb, mfbbb, mfcbb, vvx, vx2, c4o9);
			/////////////b////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfacb, mfbcb, mfccb, vvx, vx2, c1o9);
			//////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c1o36);
			/////////////c////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfabc, mfbbc, mfcbc, vvx, vx2, c1o9);
			/////////////c////////////////////////////////////////////////////////////////////////
			//backwardChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c1o36);
			//////////////////////////////////////////////////////////////////////////////////////

			//////////////////////////////////////////////////////////////////////////////////////
			real drhoPost =
				((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
					((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
			mfbbb += drho - drhoPost;
			////////////////////////////////////////////////////////////////////////////////////
			(D.f[dirE])[k] = mfabb;//(D.f[ dirE   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dirE   ])[k   ]                                                                     
			(D.f[dirW])[kw] = mfcbb;//(D.f[ dirW   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dirW   ])[kw  ]                                                                   
			(D.f[dirN])[k] = mfbab;//(D.f[ dirN   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ dirN   ])[k   ]
			(D.f[dirS])[ks] = mfbcb;//(D.f[ dirS   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ dirS   ])[ks  ]
			(D.f[dirT])[k] = mfbba;//(D.f[ dirT   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ dirT   ])[k   ]
			(D.f[dirB])[kb] = mfbbc;//(D.f[ dirB   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ dirB   ])[kb  ]
			(D.f[dirNE])[k] = mfaab;//(D.f[ dirNE  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dirNE  ])[k   ]
			(D.f[dirSW])[ksw] = mfccb;//(D.f[ dirSW  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dirSW  ])[ksw ]
			(D.f[dirSE])[ks] = mfacb;//(D.f[ dirSE  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dirSE  ])[ks  ]
			(D.f[dirNW])[kw] = mfcab;//(D.f[ dirNW  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dirNW  ])[kw  ]
			(D.f[dirTE])[k] = mfaba;//(D.f[ dirTE  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dirTE  ])[k   ]
			(D.f[dirBW])[kbw] = mfcbc;//(D.f[ dirBW  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dirBW  ])[kbw ]
			(D.f[dirBE])[kb] = mfabc;//(D.f[ dirBE  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dirBE  ])[kb  ]
			(D.f[dirTW])[kw] = mfcba;//(D.f[ dirTW  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dirTW  ])[kw  ]
			(D.f[dirTN])[k] = mfbaa;//(D.f[ dirTN  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ dirTN  ])[k   ]
			(D.f[dirBS])[kbs] = mfbcc;//(D.f[ dirBS  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ dirBS  ])[kbs ]
			(D.f[dirBN])[kb] = mfbac;//(D.f[ dirBN  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ dirBN  ])[kb  ]
			(D.f[dirTS])[ks] = mfbca;//(D.f[ dirTS  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ dirTS  ])[ks  ]
			(D.f[dirZERO])[k] = mfbbb;//(D.f[ dirZERO])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ dirZERO])[k   ]
			(D.f[dirTNE])[k] = mfaaa;//(D.f[ dirTNE ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dirTNE ])[k   ]
			(D.f[dirTSE])[ks] = mfaca;//(D.f[ dirTSE ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dirTSE ])[ks  ]
			(D.f[dirBNE])[kb] = mfaac;//(D.f[ dirBNE ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dirBNE ])[kb  ]
			(D.f[dirBSE])[kbs] = mfacc;//(D.f[ dirBSE ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dirBSE ])[kbs ]
			(D.f[dirTNW])[kw] = mfcaa;//(D.f[ dirTNW ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dirTNW ])[kw  ]
			(D.f[dirTSW])[ksw] = mfcca;//(D.f[ dirTSW ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dirTSW ])[ksw ]
			(D.f[dirBNW])[kbw] = mfcac;//(D.f[ dirBNW ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dirBNW ])[kbw ]
			(D.f[dirBSW])[kbsw] = mfccc;//(D.f[ dirBSW ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dirBSW ])[kbsw]
			////////////////////////////////////////////////////////////////////////////////////
		}
	}
}
////////////////////////////////////////////////////////////////////////////////








































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void Cumulant_One_chim_Comp_SP_27(
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	int level,
	real* forces,
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

		if (BC >= GEO_FLUID/*(BC != GEO_SOLID) && (BC != GEO_VOID)*/)
		{
			Distributions27 D;
			if (EvenOrOdd == true)
			{
				D.f[dirE   ] = &DDStart[dirE   *size_Mat];
				D.f[dirW   ] = &DDStart[dirW   *size_Mat];
				D.f[dirN   ] = &DDStart[dirN   *size_Mat];
				D.f[dirS   ] = &DDStart[dirS   *size_Mat];
				D.f[dirT   ] = &DDStart[dirT   *size_Mat];
				D.f[dirB   ] = &DDStart[dirB   *size_Mat];
				D.f[dirNE  ] = &DDStart[dirNE  *size_Mat];
				D.f[dirSW  ] = &DDStart[dirSW  *size_Mat];
				D.f[dirSE  ] = &DDStart[dirSE  *size_Mat];
				D.f[dirNW  ] = &DDStart[dirNW  *size_Mat];
				D.f[dirTE  ] = &DDStart[dirTE  *size_Mat];
				D.f[dirBW  ] = &DDStart[dirBW  *size_Mat];
				D.f[dirBE  ] = &DDStart[dirBE  *size_Mat];
				D.f[dirTW  ] = &DDStart[dirTW  *size_Mat];
				D.f[dirTN  ] = &DDStart[dirTN  *size_Mat];
				D.f[dirBS  ] = &DDStart[dirBS  *size_Mat];
				D.f[dirBN  ] = &DDStart[dirBN  *size_Mat];
				D.f[dirTS  ] = &DDStart[dirTS  *size_Mat];
				D.f[dirZERO] = &DDStart[dirZERO*size_Mat];
				D.f[dirTNE ] = &DDStart[dirTNE *size_Mat];
				D.f[dirTSW ] = &DDStart[dirTSW *size_Mat];
				D.f[dirTSE ] = &DDStart[dirTSE *size_Mat];
				D.f[dirTNW ] = &DDStart[dirTNW *size_Mat];
				D.f[dirBNE ] = &DDStart[dirBNE *size_Mat];
				D.f[dirBSW ] = &DDStart[dirBSW *size_Mat];
				D.f[dirBSE ] = &DDStart[dirBSE *size_Mat];
				D.f[dirBNW ] = &DDStart[dirBNW *size_Mat];
			}
			else
			{
				D.f[dirW   ] = &DDStart[dirE   *size_Mat];
				D.f[dirE   ] = &DDStart[dirW   *size_Mat];
				D.f[dirS   ] = &DDStart[dirN   *size_Mat];
				D.f[dirN   ] = &DDStart[dirS   *size_Mat];
				D.f[dirB   ] = &DDStart[dirT   *size_Mat];
				D.f[dirT   ] = &DDStart[dirB   *size_Mat];
				D.f[dirSW  ] = &DDStart[dirNE  *size_Mat];
				D.f[dirNE  ] = &DDStart[dirSW  *size_Mat];
				D.f[dirNW  ] = &DDStart[dirSE  *size_Mat];
				D.f[dirSE  ] = &DDStart[dirNW  *size_Mat];
				D.f[dirBW  ] = &DDStart[dirTE  *size_Mat];
				D.f[dirTE  ] = &DDStart[dirBW  *size_Mat];
				D.f[dirTW  ] = &DDStart[dirBE  *size_Mat];
				D.f[dirBE  ] = &DDStart[dirTW  *size_Mat];
				D.f[dirBS  ] = &DDStart[dirTN  *size_Mat];
				D.f[dirTN  ] = &DDStart[dirBS  *size_Mat];
				D.f[dirTS  ] = &DDStart[dirBN  *size_Mat];
				D.f[dirBN  ] = &DDStart[dirTS  *size_Mat];
				D.f[dirZERO] = &DDStart[dirZERO*size_Mat];
				D.f[dirBSW ] = &DDStart[dirTNE *size_Mat];
				D.f[dirBNE ] = &DDStart[dirTSW *size_Mat];
				D.f[dirBNW ] = &DDStart[dirTSE *size_Mat];
				D.f[dirBSE ] = &DDStart[dirTNW *size_Mat];
				D.f[dirTSW ] = &DDStart[dirBNE *size_Mat];
				D.f[dirTNE ] = &DDStart[dirBSW *size_Mat];
				D.f[dirTNW ] = &DDStart[dirBSE *size_Mat];
				D.f[dirTSE ] = &DDStart[dirBNW *size_Mat];
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
			////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[dirE   ])[k   ];
			real mfabb = (D.f[dirW   ])[kw  ];
			real mfbcb = (D.f[dirN   ])[k   ];
			real mfbab = (D.f[dirS   ])[ks  ];
			real mfbbc = (D.f[dirT   ])[k   ];
			real mfbba = (D.f[dirB   ])[kb  ];
			real mfccb = (D.f[dirNE  ])[k   ];
			real mfaab = (D.f[dirSW  ])[ksw ];
			real mfcab = (D.f[dirSE  ])[ks  ];
			real mfacb = (D.f[dirNW  ])[kw  ];
			real mfcbc = (D.f[dirTE  ])[k   ];
			real mfaba = (D.f[dirBW  ])[kbw ];
			real mfcba = (D.f[dirBE  ])[kb  ];
			real mfabc = (D.f[dirTW  ])[kw  ];
			real mfbcc = (D.f[dirTN  ])[k   ];
			real mfbaa = (D.f[dirBS  ])[kbs ];
			real mfbca = (D.f[dirBN  ])[kb  ];
			real mfbac = (D.f[dirTS  ])[ks  ];
			real mfbbb = (D.f[dirZERO])[k   ];
			real mfccc = (D.f[dirTNE ])[k   ];
			real mfaac = (D.f[dirTSW ])[ksw ];
			real mfcac = (D.f[dirTSE ])[ks  ];
			real mfacc = (D.f[dirTNW ])[kw  ];
			real mfcca = (D.f[dirBNE ])[kb  ];
			real mfaaa = (D.f[dirBSW ])[kbsw];
			real mfcaa = (D.f[dirBSE ])[kbs ];
			real mfaca = (D.f[dirBNW ])[kbw ];
			////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
				((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

			real rho = c1o1 + drho;
			real OOrho = c1o1 / rho;
			////////////////////////////////////////////////////////////////////////////////////
			real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
				(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
				(mfcbb - mfabb)) * OOrho;
			real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
				(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
				(mfbcb - mfbab)) * OOrho;
			real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
				(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
				(mfbbc - mfbba)) * OOrho;
			////////////////////////////////////////////////////////////////////////////////////
			//the force be with you
			real fx = forces[0] / (pow((double)c2o1, (double)level)); //zero;//0.0032653/(pow(two,level)); //0.000000005;//(two/1600000.0) / 120.0; //
			real fy = forces[1] / (pow((double)c2o1, (double)level)); //zero;
			real fz = forces[2] / (pow((double)c2o1, (double)level)); //zero;
			vvx += fx*c1o2;
			vvy += fy*c1o2;
			vvz += fz*c1o2;
			////////////////////////////////////////////////////////////////////////////////////
			//real omega = omega_in;
			////////////////////////////////////////////////////////////////////////////////////
			//fast
			real oMdrho = c1o1; // comp special
			real m0, m1, m2;
			real vx2;
			real vy2;
			real vz2;
			vx2 = vvx*vvx;
			vy2 = vvy*vvy;
			vz2 = vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real wadjust;
			real qudricLimitP = c1o100;// * 0.0001f;
			real qudricLimitM = c1o100;// * 0.0001f;
			real qudricLimitD = c1o100;// * 0.001f;
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			forwardInverseChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, 36.0f, c1o36);
			forwardInverseChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2, 9.0f , c1o9 );
			forwardInverseChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, 36.0f, c1o36);
			forwardInverseChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2, 9.0f , c1o9 );
			forwardInverseChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2, 2.25f, c4o9 );
			forwardInverseChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2, 9.0f , c1o9 );
			forwardInverseChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, 36.0f, c1o36);
			forwardInverseChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2, 9.0f , c1o9 );
			forwardInverseChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, 36.0f, c1o36);

			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			forwardInverseChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, 6.0f , c1o6 );
			forwardChimera(     mfaab, mfabb, mfacb, vvy, vy2);
			forwardInverseChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, 18.0f, c1o18);
			forwardInverseChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2, 1.5f , c2o3 );
			forwardChimera(     mfbab, mfbbb, mfbcb, vvy, vy2);
			forwardInverseChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2, 4.5f , c2o9 );
			forwardInverseChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, 6.0f , c1o6 );
			forwardChimera(     mfcab, mfcbb, mfccb, vvy, vy2);
			forwardInverseChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, 18.0f, c1o18);

			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			forwardInverseChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o1, c1o1);
			forwardChimera(     mfaba, mfbba, mfcba, vvx, vx2);
			forwardInverseChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, 3.0f, c1o3);
			forwardChimera(     mfaab, mfbab, mfcab, vvx, vx2);
			forwardChimera(     mfabb, mfbbb, mfcbb, vvx, vx2);
			forwardChimera(     mfacb, mfbcb, mfccb, vvx, vx2);
			forwardInverseChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, 3.0f, c1o3);
			forwardChimera(     mfabc, mfbbc, mfcbc, vvx, vx2);
			forwardInverseChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, 9.0f, c1o9);

			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = c1o1;
			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			real OxyyPxzz = c1o1;
			real OxyyMxzz = c1o1;
			real Oxyz = c1o1;
			////////////////////////////////////////////////////////////
			//4.
			//////////////////////////////
			real O4 = c1o1;
			////////////////////////////////////////////////////////////
			//5.
			//////////////////////////////
			real O5 = c1o1;
			////////////////////////////////////////////////////////////
			//6.
			//////////////////////////////
			real O6 = c1o1;
			////////////////////////////////////////////////////////////


			//central moments to cumulants
			//4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) * OOrho;
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) * OOrho;
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) * OOrho;

			real CUMcca = mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) * OOrho - c1o9*(drho * OOrho));
			real CUMcac = mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) * OOrho - c1o9*(drho * OOrho));
			real CUMacc = mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) * OOrho - c1o9*(drho * OOrho));

			//5.
			real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)) * OOrho;
			real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)) * OOrho;
			real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)) * OOrho;

			//6.
			real CUMccc = mfccc + ((-c4o1 *  mfbbb * mfbbb
				- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				- c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				- c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) * OOrho
				+ (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
				+ c2o1 * (mfcaa * mfaca * mfaac)
				+ c16o1 *  mfbba * mfbab * mfabb) * OOrho * OOrho
				- c1o3 * (mfacc + mfcac + mfcca) * OOrho
				- c1o9 * (mfcaa + mfaca + mfaac) * OOrho
				+ (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
				+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) * OOrho * OOrho  * c2o3
				+ c1o27*((drho * drho - drho) * OOrho * OOrho ));


			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy = mfcaa - mfaca;
			real mxxMzz = mfcaa - mfaac;

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
			{
				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
				real dyuy = dxux + omega * c3o2 * mxxMyy;
				real dzuz = dxux + omega * c3o2 * mxxMzz;

				//relax
				mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
				mxxMyy    += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
				mxxMzz    += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

			}
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////no correction
			//mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz);//-magicBulk*OxxPyyPzz;
			//mxxMyy += -(-omega) * (-mxxMyy);
			//mxxMzz += -(-omega) * (-mxxMzz);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb += omega * (-mfabb);
			mfbab += omega * (-mfbab);
			mfbba += omega * (-mfbba);

			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-c2o1*  mxxMyy + mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (mxxMyy - c2o1* mxxMzz + mxxPyyPzz);

			//3.
			// linear combinations

			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

			//relax
			//////////////////////////////////////////////////////////////////////////
			mfbbb     += OxyyMxzz * (-mfbbb);
			mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
			mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
			mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
			mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
			mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
			mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////

			mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			//////////////////////////////////////////////////////////////////////////
			CUMacc += O4 * (-CUMacc);
			CUMcac += O4 * (-CUMcac);
			CUMcca += O4 * (-CUMcca);

			CUMbbc += O4 * (-CUMbbc);
			CUMbcb += O4 * (-CUMbcb);
			CUMcbb += O4 * (-CUMcbb);
			//////////////////////////////////////////////////////////////////////////


			//5.
			CUMbcc += O5 * (-CUMbcc);
			CUMcbc += O5 * (-CUMcbc);
			CUMccb += O5 * (-CUMccb);

			//6.
			CUMccc += O6 * (-CUMccc);



			//back cumulants to central moments
			//4.
			mfcbb = CUMcbb + c1o3*((c3o1*mfcaa + c1o1) * mfabb + c6o1 * mfbba * mfbab) * OOrho; 
			mfbcb = CUMbcb + c1o3*((c3o1*mfaca + c1o1) * mfbab + c6o1 * mfbba * mfabb) * OOrho;
			mfbbc = CUMbbc + c1o3*((c3o1*mfaac + c1o1) * mfbba + c6o1 * mfbab * mfabb) * OOrho;

			mfcca = CUMcca + (((mfcaa * mfaca + c2o1 * mfbba * mfbba)*c9o1 + c3o1 * (mfcaa + mfaca)) * OOrho - (drho * OOrho))*c1o9;
			mfcac = CUMcac + (((mfcaa * mfaac + c2o1 * mfbab * mfbab)*c9o1 + c3o1 * (mfcaa + mfaac)) * OOrho - (drho * OOrho))*c1o9;
			mfacc = CUMacc + (((mfaac * mfaca + c2o1 * mfabb * mfabb)*c9o1 + c3o1 * (mfaac + mfaca)) * OOrho - (drho * OOrho))*c1o9;

			//5.
			mfbcc = CUMbcc + c1o3 *(c3o1*(mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + (mfbca + mfbac)) * OOrho;
			mfcbc = CUMcbc + c1o3 *(c3o1*(mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + (mfcba + mfabc)) * OOrho;
			mfccb = CUMccb + c1o3 *(c3o1*(mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) +  (mfacb + mfcab)) * OOrho;

			//6.
			mfccc = 
				CUMccc - ((-c4o1 *  mfbbb * mfbbb
				- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				- c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				- c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) * OOrho
				+ (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
				+ c2o1 * (mfcaa * mfaca * mfaac)
				+ c16o1 *  mfbba * mfbab * mfabb) * OOrho * OOrho
				- c1o3 * (mfacc + mfcac + mfcca) * OOrho
				- c1o9 * (mfcaa + mfaca + mfaac) * OOrho
				+ (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
				+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) * OOrho * OOrho * c2o3
				+ c1o27*((drho * drho - drho) * OOrho * OOrho ));

			////////////////////////////////////////////////////////////////////////////////////
			//the force be with you
			mfbaa = -mfbaa;
			mfaba = -mfaba;
			mfaab = -mfaab;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			//back
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			backwardInverseChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o1, c1o1);
			backwardChimera(			mfaba, mfbba, mfcba, vvx, vx2);
			backwardInverseChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, 3.0f, c1o3);
			backwardChimera(			mfaab, mfbab, mfcab, vvx, vx2);
			backwardChimera(			mfabb, mfbbb, mfcbb, vvx, vx2);
			backwardChimera(			mfacb, mfbcb, mfccb, vvx, vx2);
			backwardInverseChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, 3.0f, c1o3);
			backwardChimera(			mfabc, mfbbc, mfcbc, vvx, vx2);
			backwardInverseChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, 9.0f, c1o9);

			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			backwardInverseChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, 6.0f , c1o6 );
			backwardChimera(			mfaab, mfabb, mfacb, vvy, vy2);
			backwardInverseChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, 18.0f, c1o18);
			backwardInverseChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2, 1.5f , c2o3 );
			backwardChimera(			mfbab, mfbbb, mfbcb, vvy, vy2);
			backwardInverseChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2, 4.5f , c2o9 );
			backwardInverseChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, 6.0f , c1o6 );
			backwardChimera(			mfcab, mfcbb, mfccb, vvy, vy2);
			backwardInverseChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, 18.0f, c1o18);

			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			backwardInverseChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, 36.0f, c1o36);
			backwardInverseChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2, 9.0f , c1o9 );
			backwardInverseChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, 36.0f, c1o36);
			backwardInverseChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2, 9.0f , c1o9 );
			backwardInverseChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2, 2.25f, c4o9 );
			backwardInverseChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2, 9.0f , c1o9 );
			backwardInverseChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, 36.0f, c1o36);
			backwardInverseChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2, 9.0f , c1o9 );
			backwardInverseChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, 36.0f, c1o36);

			//////////////////////////////////////////////////////////////////////////////////////
			real drhoPost =
				((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
					((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
			mfbbb += drho - drhoPost;
			////////////////////////////////////////////////////////////////////////////////////
			(D.f[dirE   ])[k   ] = mfabb;                                                                   
			(D.f[dirW   ])[kw  ] = mfcbb;                                                                 
			(D.f[dirN   ])[k   ] = mfbab;
			(D.f[dirS   ])[ks  ] = mfbcb;
			(D.f[dirT   ])[k   ] = mfbba;
			(D.f[dirB   ])[kb  ] = mfbbc;
			(D.f[dirNE  ])[k   ] = mfaab;
			(D.f[dirSW  ])[ksw ] = mfccb;
			(D.f[dirSE  ])[ks  ] = mfacb;
			(D.f[dirNW  ])[kw  ] = mfcab;
			(D.f[dirTE  ])[k   ] = mfaba;
			(D.f[dirBW  ])[kbw ] = mfcbc;
			(D.f[dirBE  ])[kb  ] = mfabc;
			(D.f[dirTW  ])[kw  ] = mfcba;
			(D.f[dirTN  ])[k   ] = mfbaa;
			(D.f[dirBS  ])[kbs ] = mfbcc;
			(D.f[dirBN  ])[kb  ] = mfbac;
			(D.f[dirTS  ])[ks  ] = mfbca;
			(D.f[dirZERO])[k   ] = mfbbb;
			(D.f[dirTNE ])[k   ] = mfaaa;
			(D.f[dirTSE ])[ks  ] = mfaca;
			(D.f[dirBNE ])[kb  ] = mfaac;
			(D.f[dirBSE ])[kbs ] = mfacc;
			(D.f[dirTNW ])[kw  ] = mfcaa;
			(D.f[dirTSW ])[ksw ] = mfcca;
			(D.f[dirBNW ])[kbw ] = mfcac;
			(D.f[dirBSW ])[kbsw] = mfccc;
		}
	}
}
////////////////////////////////////////////////////////////////////////////////








































