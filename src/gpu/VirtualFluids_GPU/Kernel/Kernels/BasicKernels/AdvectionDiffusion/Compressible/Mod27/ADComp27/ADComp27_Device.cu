#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void LB_KERNEL_AD_COMP_27(real diffusivity,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	real* DD27,
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

			Distributions27 D27;
			if (EvenOrOdd == true)
			{
				D27.f[dP00] = &DD27[dP00 * size_Mat];
				D27.f[dM00] = &DD27[dM00 * size_Mat];
				D27.f[d0P0] = &DD27[d0P0 * size_Mat];
				D27.f[d0M0] = &DD27[d0M0 * size_Mat];
				D27.f[d00P] = &DD27[d00P * size_Mat];
				D27.f[d00M] = &DD27[d00M * size_Mat];
				D27.f[dPP0] = &DD27[dPP0 * size_Mat];
				D27.f[dMM0] = &DD27[dMM0 * size_Mat];
				D27.f[dPM0] = &DD27[dPM0 * size_Mat];
				D27.f[dMP0] = &DD27[dMP0 * size_Mat];
				D27.f[dP0P] = &DD27[dP0P * size_Mat];
				D27.f[dM0M] = &DD27[dM0M * size_Mat];
				D27.f[dP0M] = &DD27[dP0M * size_Mat];
				D27.f[dM0P] = &DD27[dM0P * size_Mat];
				D27.f[d0PP] = &DD27[d0PP * size_Mat];
				D27.f[d0MM] = &DD27[d0MM * size_Mat];
				D27.f[d0PM] = &DD27[d0PM * size_Mat];
				D27.f[d0MP] = &DD27[d0MP * size_Mat];
				D27.f[d000] = &DD27[d000 * size_Mat];
				D27.f[dPPP] = &DD27[dPPP * size_Mat];
				D27.f[dMMP] = &DD27[dMMP * size_Mat];
				D27.f[dPMP] = &DD27[dPMP * size_Mat];
				D27.f[dMPP] = &DD27[dMPP * size_Mat];
				D27.f[dPPM] = &DD27[dPPM * size_Mat];
				D27.f[dMMM] = &DD27[dMMM * size_Mat];
				D27.f[dPMM] = &DD27[dPMM * size_Mat];
				D27.f[dMPM] = &DD27[dMPM * size_Mat];
			}
			else
			{
				D27.f[dM00] = &DD27[dP00 * size_Mat];
				D27.f[dP00] = &DD27[dM00 * size_Mat];
				D27.f[d0M0] = &DD27[d0P0 * size_Mat];
				D27.f[d0P0] = &DD27[d0M0 * size_Mat];
				D27.f[d00M] = &DD27[d00P * size_Mat];
				D27.f[d00P] = &DD27[d00M * size_Mat];
				D27.f[dMM0] = &DD27[dPP0 * size_Mat];
				D27.f[dPP0] = &DD27[dMM0 * size_Mat];
				D27.f[dMP0] = &DD27[dPM0 * size_Mat];
				D27.f[dPM0] = &DD27[dMP0 * size_Mat];
				D27.f[dM0M] = &DD27[dP0P * size_Mat];
				D27.f[dP0P] = &DD27[dM0M * size_Mat];
				D27.f[dM0P] = &DD27[dP0M * size_Mat];
				D27.f[dP0M] = &DD27[dM0P * size_Mat];
				D27.f[d0MM] = &DD27[d0PP * size_Mat];
				D27.f[d0PP] = &DD27[d0MM * size_Mat];
				D27.f[d0MP] = &DD27[d0PM * size_Mat];
				D27.f[d0PM] = &DD27[d0MP * size_Mat];
				D27.f[d000] = &DD27[d000 * size_Mat];
				D27.f[dMMM] = &DD27[dPPP * size_Mat];
				D27.f[dPPM] = &DD27[dMMP * size_Mat];
				D27.f[dMPM] = &DD27[dPMP * size_Mat];
				D27.f[dPMM] = &DD27[dMPP * size_Mat];
				D27.f[dMMP] = &DD27[dPPM * size_Mat];
				D27.f[dPPP] = &DD27[dMMM * size_Mat];
				D27.f[dMPP] = &DD27[dPMM * size_Mat];
				D27.f[dPMP] = &DD27[dMPM * size_Mat];
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
										   ////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D27.f[dP00])[k];
			real mfabb = (D27.f[dM00])[kw];
			real mfbcb = (D27.f[d0P0])[k];
			real mfbab = (D27.f[d0M0])[ks];
			real mfbbc = (D27.f[d00P])[k];
			real mfbba = (D27.f[d00M])[kb];
			real mfccb = (D27.f[dPP0])[k];
			real mfaab = (D27.f[dMM0])[ksw];
			real mfcab = (D27.f[dPM0])[ks];
			real mfacb = (D27.f[dMP0])[kw];
			real mfcbc = (D27.f[dP0P])[k];
			real mfaba = (D27.f[dM0M])[kbw];
			real mfcba = (D27.f[dP0M])[kb];
			real mfabc = (D27.f[dM0P])[kw];
			real mfbcc = (D27.f[d0PP])[k];
			real mfbaa = (D27.f[d0MM])[kbs];
			real mfbca = (D27.f[d0PM])[kb];
			real mfbac = (D27.f[d0MP])[ks];
			real mfbbb = (D27.f[d000])[k];
			real mfccc = (D27.f[dPPP])[k];
			real mfaac = (D27.f[dMMP])[ksw];
			real mfcac = (D27.f[dPMP])[ks];
			real mfacc = (D27.f[dMPP])[kw];
			real mfcca = (D27.f[dPPM])[kb];
			real mfaaa = (D27.f[dMMM])[kbsw];
			real mfcaa = (D27.f[dPMM])[kbs];
			real mfaca = (D27.f[dMPM])[kbw];
			////////////////////////////////////////////////////////////////////////////////////
			//Conc
			real drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
				((mfabb + mfcbb) + (mfbab + mfbcb)) + (mfbba + mfbbc)) + mfbbb;

			//real rho = c1o1 + drho;
			////////////////////////////////////////////////////////////////////////////////
			real rho0fluid = (fTNE + fBSW) + (fTSW + fBNE) + (fTSE + fBNW) + (fTNW + fBSE) + (fNE + fSW) + (fNW + fSE) + (fTE + fBW) + (fBE + fTW) + (fTN + fBS) + (fBN + fTS) + (fE + fW) + (fN + fS) + (fT + fB) + fZERO;
			real rhofluid = rho0fluid + c1o1;
			real OORhofluid = c1o1 / rhofluid;
			real vvx = OORhofluid*((fTNE - fBSW) + (fBNE - fTSW) + (fTSE - fBNW) + (fBSE - fTNW) + (fNE - fSW) + (fSE - fNW) + (fTE - fBW) + (fBE - fTW) + (fE - fW));
			real vvy = OORhofluid*((fTNE - fBSW) + (fBNE - fTSW) + (fBNW - fTSE) + (fTNW - fBSE) + (fNE - fSW) + (fNW - fSE) + (fTN - fBS) + (fBN - fTS) + (fN - fS));
			real vvz = OORhofluid*((fTNE - fBSW) + (fTSW - fBNE) + (fTSE - fBNW) + (fTNW - fBSE) + (fTE - fBW) + (fTW - fBE) + (fTN - fBS) + (fTS - fBN) + (fT - fB));
			////////////////////////////////////////////////////////////////////////////////
			//real omegaD     = zero;
			real omegaD = c2o1 / (c6o1 * diffusivity + c1o1);

			real oMdrho = c0o1;//one; // comp special
			real m0, m1, m2;
			real vx2 = vvx*vvx;
			real vy2 = vvy*vvy;
			real vz2 = vvz*vvz;

			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2 = mfaaa + mfaac;
			m1 = mfaac - mfaaa;
			m0 = m2 + mfaab;
			mfaaa = m0;
			m0 += c1o36 * oMdrho;
			mfaab = m1 - m0 * vvz;
			mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaba + mfabc;
			m1 = mfabc - mfaba;
			m0 = m2 + mfabb;
			mfaba = m0;
			m0 += c1o9 * oMdrho;
			mfabb = m1 - m0 * vvz;
			mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaca + mfacc;
			m1 = mfacc - mfaca;
			m0 = m2 + mfacb;
			mfaca = m0;
			m0 += c1o36 * oMdrho;
			mfacb = m1 - m0 * vvz;
			mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbaa + mfbac;
			m1 = mfbac - mfbaa;
			m0 = m2 + mfbab;
			mfbaa = m0;
			m0 += c1o9 * oMdrho;
			mfbab = m1 - m0 * vvz;
			mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbba + mfbbc;
			m1 = mfbbc - mfbba;
			m0 = m2 + mfbbb;
			mfbba = m0;
			m0 += c4o9 * oMdrho;
			mfbbb = m1 - m0 * vvz;
			mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbca + mfbcc;
			m1 = mfbcc - mfbca;
			m0 = m2 + mfbcb;
			mfbca = m0;
			m0 += c1o9 * oMdrho;
			mfbcb = m1 - m0 * vvz;
			mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcaa + mfcac;
			m1 = mfcac - mfcaa;
			m0 = m2 + mfcab;
			mfcaa = m0;
			m0 += c1o36 * oMdrho;
			mfcab = m1 - m0 * vvz;
			mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcba + mfcbc;
			m1 = mfcbc - mfcba;
			m0 = m2 + mfcbb;
			mfcba = m0;
			m0 += c1o9 * oMdrho;
			mfcbb = m1 - m0 * vvz;
			mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcca + mfccc;
			m1 = mfccc - mfcca;
			m0 = m2 + mfccb;
			mfcca = m0;
			m0 += c1o36 * oMdrho;
			mfccb = m1 - m0 * vvz;
			mfccc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2 = mfaaa + mfaca;
			m1 = mfaca - mfaaa;
			m0 = m2 + mfaba;
			mfaaa = m0;
			m0 += c1o6 * oMdrho;
			mfaba = m1 - m0 * vvy;
			mfaca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaab + mfacb;
			m1 = mfacb - mfaab;
			m0 = m2 + mfabb;
			mfaab = m0;
			mfabb = m1 - m0 * vvy;
			mfacb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaac + mfacc;
			m1 = mfacc - mfaac;
			m0 = m2 + mfabc;
			mfaac = m0;
			m0 += c1o18 * oMdrho;
			mfabc = m1 - m0 * vvy;
			mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbaa + mfbca;
			m1 = mfbca - mfbaa;
			m0 = m2 + mfbba;
			mfbaa = m0;
			m0 += c2o3 * oMdrho;
			mfbba = m1 - m0 * vvy;
			mfbca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbab + mfbcb;
			m1 = mfbcb - mfbab;
			m0 = m2 + mfbbb;
			mfbab = m0;
			mfbbb = m1 - m0 * vvy;
			mfbcb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbac + mfbcc;
			m1 = mfbcc - mfbac;
			m0 = m2 + mfbbc;
			mfbac = m0;
			m0 += c2o9 * oMdrho;
			mfbbc = m1 - m0 * vvy;
			mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcaa + mfcca;
			m1 = mfcca - mfcaa;
			m0 = m2 + mfcba;
			mfcaa = m0;
			m0 += c1o6 * oMdrho;
			mfcba = m1 - m0 * vvy;
			mfcca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcab + mfccb;
			m1 = mfccb - mfcab;
			m0 = m2 + mfcbb;
			mfcab = m0;
			mfcbb = m1 - m0 * vvy;
			mfccb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcac + mfccc;
			m1 = mfccc - mfcac;
			m0 = m2 + mfcbc;
			mfcac = m0;
			m0 += c1o18 * oMdrho;
			mfcbc = m1 - m0 * vvy;
			mfccc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2 = mfaaa + mfcaa;
			m1 = mfcaa - mfaaa;
			m0 = m2 + mfbaa;
			mfaaa = m0;
			m0 += c1o1* oMdrho;
			mfbaa = m1 - m0 * vvx;
			mfcaa = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaba + mfcba;
			m1 = mfcba - mfaba;
			m0 = m2 + mfbba;
			mfaba = m0;
			mfbba = m1 - m0 * vvx;
			mfcba = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaca + mfcca;
			m1 = mfcca - mfaca;
			m0 = m2 + mfbca;
			mfaca = m0;
			m0 += c1o3 * oMdrho;
			mfbca = m1 - m0 * vvx;
			mfcca = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaab + mfcab;
			m1 = mfcab - mfaab;
			m0 = m2 + mfbab;
			mfaab = m0;
			mfbab = m1 - m0 * vvx;
			mfcab = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfabb + mfcbb;
			m1 = mfcbb - mfabb;
			m0 = m2 + mfbbb;
			mfabb = m0;
			mfbbb = m1 - m0 * vvx;
			mfcbb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfacb + mfccb;
			m1 = mfccb - mfacb;
			m0 = m2 + mfbcb;
			mfacb = m0;
			mfbcb = m1 - m0 * vvx;
			mfccb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaac + mfcac;
			m1 = mfcac - mfaac;
			m0 = m2 + mfbac;
			mfaac = m0;
			m0 += c1o3 * oMdrho;
			mfbac = m1 - m0 * vvx;
			mfcac = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfabc + mfcbc;
			m1 = mfcbc - mfabc;
			m0 = m2 + mfbbc;
			mfabc = m0;
			mfbbc = m1 - m0 * vvx;
			mfcbc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfacc + mfccc;
			m1 = mfccc - mfacc;
			m0 = m2 + mfbcc;
			mfacc = m0;
			m0 += c1o9 * oMdrho;
			mfbcc = m1 - m0 * vvx;
			mfccc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////

			//if(mfaaa < zero) omegaD = one;
			real limit = c9o1*omegaD*omegaD*(mfbaa*mfbaa + mfaba*mfaba + mfaab*mfaab);
			//real CC=c1o2;
			//if ((two*mfaaa*mfaaa<limit)) omegaD=two / (six * (diffusivity+((limit/(1.0e-10f+two*mfaaa*mfaaa)-one)*(c1o6-diffusivity))*c1o2) + one);
			if ((c2o1*mfaaa*mfaaa<limit)) omegaD = c1o1;
			//omegaD = two / (six * (diffusivity+CC*limit) + one);

			//mfaaa = c1o2;
			//trans 3.
			real Mabc = mfabc - mfaba*c1o3;
			real Mbca = mfbca - mfbaa*c1o3;
			real Macb = mfacb - mfaab*c1o3;
			real Mcba = mfcba - mfaba*c1o3;
			real Mcab = mfcab - mfaab*c1o3;
			real Mbac = mfbac - mfbaa*c1o3;
			//trans 5.
			real Mcbc = mfcbc - mfaba*c1o9;
			real Mbcc = mfbcc - mfbaa*c1o9;
			real Mccb = mfccb - mfaab*c1o9;

			//1.
			mfbaa *= c1o1 - omegaD;
			mfaba *= c1o1 - omegaD;
			mfaab *= c1o1 - omegaD;

			//3.
			//mfbca *= one - omegaD;
			//mfbac *= one - omegaD;
			//mfcba *= one - omegaD;
			//mfabc *= one - omegaD;
			//mfcab *= one - omegaD;
			//mfacb *= one - omegaD;

			//mfbbb *= one - omegaD; 
			Mabc = c0o1;
			Mbca = c0o1;
			Macb = c0o1;
			Mcba = c0o1;
			Mcab = c0o1;
			Mbac = c0o1;
			mfbbb = c0o1;

			//5.
			//mfbcc *= one - omegaD;
			//mfcbc *= one - omegaD;
			//mfccb *= one - omegaD;
			Mcbc = c0o1;
			Mbcc = c0o1;
			Mccb = c0o1;

			//2.
			mfbba = c0o1;
			mfbab = c0o1;
			mfabb = c0o1;

			mfcaa = c1o3 * drho;
			mfaca = c1o3 * drho;
			mfaac = c1o3 * drho;

			//4.
			mfacc = c1o9 * drho;
			mfcac = c1o9 * drho;
			mfcca = c1o9 * drho;

			mfcbb = c0o1;
			mfbcb = c0o1;
			mfbbc = c0o1;

			//6.
			mfccc = c1o27 * drho;

			//3.
			mfabc = Mabc + mfaba*c1o3;
			mfbca = Mbca + mfbaa*c1o3;
			mfacb = Macb + mfaab*c1o3;
			mfcba = Mcba + mfaba*c1o3;
			mfcab = Mcab + mfaab*c1o3;
			mfbac = Mbac + mfbaa*c1o3;
			//5.	  
			mfcbc = Mcbc + mfaba*c1o9;
			mfbcc = Mbcc + mfbaa*c1o9;
			mfccb = Mccb + mfaab*c1o9;

			////////////////////////////////////////////////////////////////////////////////////
			//back
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m0 = mfaac * c1o2 + mfaab * (vvz - c1o2) + (mfaaa + c1o1* oMdrho) * (vz2 - vvz) * c1o2;
			m1 = -mfaac - c2o1* mfaab *  vvz + mfaaa                * (c1o1 - vz2) - c1o1* oMdrho * vz2;
			m2 = mfaac * c1o2 + mfaab * (vvz + c1o2) + (mfaaa + c1o1* oMdrho) * (vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfabc * c1o2 + mfabb * (vvz - c1o2) + mfaba * (vz2 - vvz) * c1o2;
			m1 = -mfabc - c2o1* mfabb *  vvz + mfaba * (c1o1 - vz2);
			m2 = mfabc * c1o2 + mfabb * (vvz + c1o2) + mfaba * (vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfacc * c1o2 + mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
			m1 = -mfacc - c2o1* mfacb *  vvz + mfaca                  * (c1o1 - vz2) - c1o3 * oMdrho * vz2;
			m2 = mfacc * c1o2 + mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfbac * c1o2 + mfbab * (vvz - c1o2) + mfbaa * (vz2 - vvz) * c1o2;
			m1 = -mfbac - c2o1* mfbab *  vvz + mfbaa * (c1o1 - vz2);
			m2 = mfbac * c1o2 + mfbab * (vvz + c1o2) + mfbaa * (vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 = mfbbc * c1o2 + mfbbb * (vvz - c1o2) + mfbba * (vz2 - vvz) * c1o2;
			m1 = -mfbbc - c2o1* mfbbb *  vvz + mfbba * (c1o1 - vz2);
			m2 = mfbbc * c1o2 + mfbbb * (vvz + c1o2) + mfbba * (vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 = mfbcc * c1o2 + mfbcb * (vvz - c1o2) + mfbca * (vz2 - vvz) * c1o2;
			m1 = -mfbcc - c2o1* mfbcb *  vvz + mfbca * (c1o1 - vz2);
			m2 = mfbcc * c1o2 + mfbcb * (vvz + c1o2) + mfbca * (vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcac * c1o2 + mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
			m1 = -mfcac - c2o1* mfcab *  vvz + mfcaa                  * (c1o1 - vz2) - c1o3 * oMdrho * vz2;
			m2 = mfcac * c1o2 + mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 = mfcbc * c1o2 + mfcbb * (vvz - c1o2) + mfcba * (vz2 - vvz) * c1o2;
			m1 = -mfcbc - c2o1* mfcbb *  vvz + mfcba * (c1o1 - vz2);
			m2 = mfcbc * c1o2 + mfcbb * (vvz + c1o2) + mfcba * (vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 = mfccc * c1o2 + mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 - vvz) * c1o2;
			m1 = -mfccc - c2o1* mfccb *  vvz + mfcca                  * (c1o1 - vz2) - c1o9 * oMdrho * vz2;
			m2 = mfccc * c1o2 + mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 + vvz) * c1o2;
			mfcca = m0;
			mfccb = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m0 = mfaca * c1o2 + mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
			m1 = -mfaca - c2o1* mfaba *  vvy + mfaaa                  * (c1o1 - vy2) - c1o6 * oMdrho * vy2;
			m2 = mfaca * c1o2 + mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfacb * c1o2 + mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 - vvy) * c1o2;
			m1 = -mfacb - c2o1* mfabb *  vvy + mfaab                  * (c1o1 - vy2) - c2o3 * oMdrho * vy2;
			m2 = mfacb * c1o2 + mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfacc * c1o2 + mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
			m1 = -mfacc - c2o1* mfabc *  vvy + mfaac                  * (c1o1 - vy2) - c1o6 * oMdrho * vy2;
			m2 = mfacc * c1o2 + mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfbca * c1o2 + mfbba * (vvy - c1o2) + mfbaa * (vy2 - vvy) * c1o2;
			m1 = -mfbca - c2o1* mfbba *  vvy + mfbaa * (c1o1 - vy2);
			m2 = mfbca * c1o2 + mfbba * (vvy + c1o2) + mfbaa * (vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 = mfbcb * c1o2 + mfbbb * (vvy - c1o2) + mfbab * (vy2 - vvy) * c1o2;
			m1 = -mfbcb - c2o1* mfbbb *  vvy + mfbab * (c1o1 - vy2);
			m2 = mfbcb * c1o2 + mfbbb * (vvy + c1o2) + mfbab * (vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 = mfbcc * c1o2 + mfbbc * (vvy - c1o2) + mfbac * (vy2 - vvy) * c1o2;
			m1 = -mfbcc - c2o1* mfbbc *  vvy + mfbac * (c1o1 - vy2);
			m2 = mfbcc * c1o2 + mfbbc * (vvy + c1o2) + mfbac * (vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcca * c1o2 + mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
			m1 = -mfcca - c2o1* mfcba *  vvy + mfcaa                   * (c1o1 - vy2) - c1o18 * oMdrho * vy2;
			m2 = mfcca * c1o2 + mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 = mfccb * c1o2 + mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 - vvy) * c1o2;
			m1 = -mfccb - c2o1* mfcbb *  vvy + mfcab                  * (c1o1 - vy2) - c2o9 * oMdrho * vy2;
			m2 = mfccb * c1o2 + mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 = mfccc * c1o2 + mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
			m1 = -mfccc - c2o1* mfcbc *  vvy + mfcac                   * (c1o1 - vy2) - c1o18 * oMdrho * vy2;
			m2 = mfccc * c1o2 + mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
			mfcac = m0;
			mfcbc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m0 = mfcaa * c1o2 + mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcaa - c2o1* mfbaa *  vvx + mfaaa                   * (c1o1 - vx2) - c1o36 * oMdrho * vx2;
			m2 = mfcaa * c1o2 + mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcba * c1o2 + mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcba - c2o1* mfbba *  vvx + mfaba                  * (c1o1 - vx2) - c1o9 * oMdrho * vx2;
			m2 = mfcba * c1o2 + mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcca * c1o2 + mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcca - c2o1* mfbca *  vvx + mfaca                   * (c1o1 - vx2) - c1o36 * oMdrho * vx2;
			m2 = mfcca * c1o2 + mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcab * c1o2 + mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcab - c2o1* mfbab *  vvx + mfaab                  * (c1o1 - vx2) - c1o9 * oMdrho * vx2;
			m2 = mfcab * c1o2 + mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 = mfcbb * c1o2 + mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcbb - c2o1* mfbbb *  vvx + mfabb                  * (c1o1 - vx2) - c4o9 * oMdrho * vx2;
			m2 = mfcbb * c1o2 + mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 = mfccb * c1o2 + mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfccb - c2o1* mfbcb *  vvx + mfacb                  * (c1o1 - vx2) - c1o9 * oMdrho * vx2;
			m2 = mfccb * c1o2 + mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcac * c1o2 + mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcac - c2o1* mfbac *  vvx + mfaac                   * (c1o1 - vx2) - c1o36 * oMdrho * vx2;
			m2 = mfcac * c1o2 + mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 = mfcbc * c1o2 + mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcbc - c2o1* mfbbc *  vvx + mfabc                  * (c1o1 - vx2) - c1o9 * oMdrho * vx2;
			m2 = mfcbc * c1o2 + mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 = mfccc * c1o2 + mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfccc - c2o1* mfbcc *  vvx + mfacc                   * (c1o1 - vx2) - c1o36 * oMdrho * vx2;
			m2 = mfccc * c1o2 + mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			(D27.f[dP00])[k] = mfabb;
			(D27.f[dM00])[kw] = mfcbb;
			(D27.f[d0P0])[k] = mfbab;
			(D27.f[d0M0])[ks] = mfbcb;
			(D27.f[d00P])[k] = mfbba;
			(D27.f[d00M])[kb] = mfbbc;
			(D27.f[dPP0])[k] = mfaab;
			(D27.f[dMM0])[ksw] = mfccb;
			(D27.f[dPM0])[ks] = mfacb;
			(D27.f[dMP0])[kw] = mfcab;
			(D27.f[dP0P])[k] = mfaba;
			(D27.f[dM0M])[kbw] = mfcbc;
			(D27.f[dP0M])[kb] = mfabc;
			(D27.f[dM0P])[kw] = mfcba;
			(D27.f[d0PP])[k] = mfbaa;
			(D27.f[d0MM])[kbs] = mfbcc;
			(D27.f[d0PM])[kb] = mfbac;
			(D27.f[d0MP])[ks] = mfbca;
			(D27.f[d000])[k] = mfbbb;
			(D27.f[dPPP])[k] = mfaaa;
			(D27.f[dPMP])[ks] = mfaca;
			(D27.f[dPPM])[kb] = mfaac;
			(D27.f[dPMM])[kbs] = mfacc;
			(D27.f[dMPP])[kw] = mfcaa;
			(D27.f[dMMP])[ksw] = mfcca;
			(D27.f[dMPM])[kbw] = mfcac;
			(D27.f[dMMM])[kbsw] = mfccc;
			////////////////////////////////////////////////////////////////////////////////////
		}
	}
}