#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
#include "math.h"

extern "C" __global__ void LB_Kernel_Cum_Comp_SP_27(real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
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
				D.f[REST] = &DDStart[REST*size_Mat];
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
				D.f[REST] = &DDStart[REST*size_Mat];
				D.f[BSW] = &DDStart[TNE *size_Mat];
				D.f[BNE] = &DDStart[TSW *size_Mat];
				D.f[BNW] = &DDStart[TSE *size_Mat];
				D.f[BSE] = &DDStart[TNW *size_Mat];
				D.f[TSW] = &DDStart[BNE *size_Mat];
				D.f[TNE] = &DDStart[BSW *size_Mat];
				D.f[TNW] = &DDStart[BSE *size_Mat];
				D.f[TSE] = &DDStart[BNW *size_Mat];
			}

			////////////////////////////////////////////////////////////////////////////////
			//index
			unsigned int kzero = k;
			unsigned int ke = k;
			unsigned int kw = neighborX[k];
			unsigned int kn = k;
			unsigned int ks = neighborY[k];
			unsigned int kt = k;
			unsigned int kb = neighborZ[k];
			unsigned int ksw = neighborY[kw];
			unsigned int kne = k;
			unsigned int kse = ks;
			unsigned int knw = kw;
			unsigned int kbw = neighborZ[kw];
			unsigned int kte = k;
			unsigned int kbe = kb;
			unsigned int ktw = kw;
			unsigned int kbs = neighborZ[ks];
			unsigned int ktn = k;
			unsigned int kbn = kb;
			unsigned int kts = ks;
			unsigned int ktse = ks;
			unsigned int kbnw = kbw;
			unsigned int ktnw = kw;
			unsigned int kbse = kbs;
			unsigned int ktsw = ksw;
			unsigned int kbne = kb;
			unsigned int ktne = k;
			unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real f_E = (D.f[E])[ke];// +  c2over27 ;
			real f_W = (D.f[W])[kw];// +  c2over27 ;
			real f_N = (D.f[N])[kn];// +  c2over27 ;
			real f_S = (D.f[S])[ks];// +  c2over27 ;
			real f_T = (D.f[T])[kt];// +  c2over27 ;
			real f_B = (D.f[B])[kb];// +  c2over27 ;
			real f_NE = (D.f[NE])[kne];// +  c1over54 ;
			real f_SS = (D.f[SW])[ksw];// +  c1over54 ;
			real f_SE = (D.f[SE])[kse];// +  c1over54 ;
			real f_NW = (D.f[NW])[knw];// +  c1over54 ;
			real f_TE = (D.f[TE])[kte];// +  c1over54 ;
			real f_BW = (D.f[BW])[kbw];// +  c1over54 ;
			real f_EB = (D.f[BE])[kbe];// +  c1over54 ;
			real f_TW = (D.f[TW])[ktw];// +  c1over54 ;
			real f_TN = (D.f[TN])[ktn];// +  c1over54 ;
			real f_BS = (D.f[BS])[kbs];// +  c1over54 ;
			real f_BN = (D.f[BN])[kbn];// +  c1over54 ;
			real f_TS = (D.f[TS])[kts];// +  c1over54 ;
			real f_R = (D.f[REST])[kzero];// +  c8over27 ;
			real f_TNE = (D.f[TNE])[ktne];// +  c1over216;
			real f_TSW = (D.f[TSW])[ktsw];// +  c1over216;
			real f_TSE = (D.f[TSE])[ktse];// +  c1over216;
			real f_TNW = (D.f[TNW])[ktnw];// +  c1over216;
			real f_BNE = (D.f[BNE])[kbne];// +  c1over216;
			real f_BSW = (D.f[BSW])[kbsw];// +  c1over216;
			real f_BSE = (D.f[BSE])[kbse];// +  c1over216;
			real f_BNW = (D.f[BNW])[kbnw];// +  c1over216;
										   ////////////////////////////////////////////////////////////////////////////////////
			real fx = c0o1;
			real fy = c0o1;
			real fz = c0o1;
			////////////////////////////////////////////////////////////////////////////////////
			real rho = f_NW + f_W + f_SS + f_S + f_SE + f_E + f_NE + f_N + f_R + f_TN + f_BN + f_TS + f_BS + f_TE + f_EB + f_TW + f_BW + f_TNW + f_BNW + f_TNE + f_BNE + f_TSW + f_BSW + f_TSE + f_BSE + f_T + f_B + c1o1;// ACHTUNG ne EINS !!!!!!!!
			real pix = (f_NE + f_E + f_SE + f_TE + f_EB - f_NW - f_W - f_SS - f_TW - f_BW + f_TNE + f_BNE + f_TSE + f_BSE - f_TNW - f_BNW - f_TSW - f_BSW);
			real piy = (f_NE + f_N + f_NW + f_TN + f_BN - f_SE - f_S - f_SS - f_TS - f_BS + f_TNE + f_BNE + f_TNW + f_BNW - f_TSE - f_BSE - f_TSW - f_BSW);
			real piz = (f_TN + f_TS + f_TW + f_TE + f_T - f_BN - f_BS - f_BW - f_EB - f_B + f_TNE + f_TNW + f_TSE + f_TSW - f_BNE - f_BNW - f_BSE - f_BSW);
			real vvx = pix / rho + fx;
			real vvy = piy / rho + fy;
			real vvz = piz / rho + fz;
			real vx2 = vvx*vvx;
			real vy2 = vvy*vvy;
			real vz2 = vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real mfaaa = f_BSW;
			real mfaab = f_SS;
			real mfaac = f_TSW;
			real mfaba = f_BW;
			real mfabb = f_W;
			real mfabc = f_TW;
			real mfbaa = f_BS;
			real mfbab = f_S;
			real mfbac = f_TS;
			real mfbba = f_B;
			real mfbbb = f_R;
			real mfbbc = f_T;
			real mfaca = f_BNW;
			real mfacb = f_NW;
			real mfacc = f_TNW;
			real mfcaa = f_BSE;
			real mfcab = f_SE;
			real mfcac = f_TSE;
			real mfcca = f_BNE;
			real mfccb = f_NE;
			real mfccc = f_TNE;
			real mfbca = f_BN;
			real mfbcb = f_N;
			real mfbcc = f_TN;
			real mfcba = f_EB;
			real mfcbb = f_E;
			real mfcbc = f_TE;
			real m0, m1, m2;
			real wadjust;
			real qudricLimit = c1o100;
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
			m0 += c1o36;
			mfaab = m1 - m0 * vvz;
			mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaba + mfabc;
			m1 = mfabc - mfaba;
			m0 = m2 + mfabb;
			mfaba = m0;
			m0 += c1o9;
			mfabb = m1 - m0 * vvz;
			mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaca + mfacc;
			m1 = mfacc - mfaca;
			m0 = m2 + mfacb;
			mfaca = m0;
			m0 += c1o36;
			mfacb = m1 - m0 * vvz;
			mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbaa + mfbac;
			m1 = mfbac - mfbaa;
			m0 = m2 + mfbab;
			mfbaa = m0;
			m0 += c1o9;
			mfbab = m1 - m0 * vvz;
			mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbba + mfbbc;
			m1 = mfbbc - mfbba;
			m0 = m2 + mfbbb;
			mfbba = m0;
			m0 += c4o9;
			mfbbb = m1 - m0 * vvz;
			mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbca + mfbcc;
			m1 = mfbcc - mfbca;
			m0 = m2 + mfbcb;
			mfbca = m0;
			m0 += c1o9;
			mfbcb = m1 - m0 * vvz;
			mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcaa + mfcac;
			m1 = mfcac - mfcaa;
			m0 = m2 + mfcab;
			mfcaa = m0;
			m0 += c1o36;
			mfcab = m1 - m0 * vvz;
			mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcba + mfcbc;
			m1 = mfcbc - mfcba;
			m0 = m2 + mfcbb;
			mfcba = m0;
			m0 += c1o9;
			mfcbb = m1 - m0 * vvz;
			mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcca + mfccc;
			m1 = mfccc - mfcca;
			m0 = m2 + mfccb;
			mfcca = m0;
			m0 += c1o36;
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
			m0 += c1o6;
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
			m0 += c1o18;
			mfabc = m1 - m0 * vvy;
			mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbaa + mfbca;
			m1 = mfbca - mfbaa;
			m0 = m2 + mfbba;
			mfbaa = m0;
			m0 += c2o3;
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
			m0 += c2o9;
			mfbbc = m1 - m0 * vvy;
			mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcaa + mfcca;
			m1 = mfcca - mfcaa;
			m0 = m2 + mfcba;
			mfcaa = m0;
			m0 += c1o6;
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
			m0 += c1o18;
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
			m0 += c1o1;
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
			m0 += c1o3;
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
			m0 += c1o3;
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
			m0 += c1o9;
			mfbcc = m1 - m0 * vvx;
			mfccc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////


			//////////////////////////////////////////////////////////////////////////////////////
			//// BGK
			//////////////////////////////////////////////////////////////////////////////////////
			////2.
			//mfabb += -(-omega) * (-mfabb);
			//mfbab += -(-omega) * (-mfbab);
			//mfbba += -(-omega) * (-mfbba);
			//
			//mfcaa += -(-omega) * (c1o3 * mfaaa - mfcaa);
			//mfaca += -(-omega) * (c1o3 * mfaaa - mfaca);
			//mfaac += -(-omega) * (c1o3 * mfaaa - mfaac);
			//
			////3.
			//mfabc += -(-omega) * (-mfabc);
			//mfbac += -(-omega) * (-mfbac);
			//
			//mfacb += -(-omega) * (-mfacb);
			//mfbca += -(-omega) * (-mfbca);

			//mfcab += -(-omega) * (-mfcab);
			//mfcba += -(-omega) * (-mfcba);

			//mfbbb += -(-omega) * (-mfbbb);

			////4.
			//mfacc += -(-omega) * (c1o9 * mfaaa - mfacc);
			//mfcac += -(-omega) * (c1o9 * mfaaa - mfcac);
			//mfcca += -(-omega) * (c1o9 * mfaaa - mfcca);

			//mfbbc += -(-omega) * (-mfbbc);
			//mfbcb += -(-omega) * (-mfbcb);
			//mfcbb += -(-omega) * (-mfcbb);

			////5.
			//mfbcc += -(-omega) * (-mfbcc);
			//mfcbc += -(-omega) * (-mfcbc);
			//mfccb += -(-omega) * (-mfccb);

			////6.
			//mfccc += -(-omega) * (c1o27 * mfaaa - mfccc);
			//////////////////////////////////////////////////////////////////////////////////////



			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = c1o1;
			real OxyyPxzz = c1o1;//two+(-omega);//one;
			real OxyyMxzz = c1o1;//two+(-omega);//one;
			real O4 = c1o1;
			real O5 = c1o1;
			real O6 = c1o1;

			//Cum 4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3 * rho) * mfabb + c2o1* mfbba * mfbab) / rho;
			real CUMbcb = mfbcb - ((mfaca + c1o3 * rho) * mfbab + c2o1* mfbba * mfabb) / rho;
			real CUMbbc = mfbbc - ((mfaac + c1o3 * rho) * mfbba + c2o1* mfbab * mfabb) / rho;

			real CUMcca = mfcca - (mfcaa * mfaca + c2o1* mfbba * mfbba) / rho - c1o3 * (mfcaa + mfaca);
			real CUMcac = mfcac - (mfcaa * mfaac + c2o1* mfbab * mfbab) / rho - c1o3 * (mfcaa + mfaac);
			real CUMacc = mfacc - (mfaac * mfaca + c2o1* mfabb * mfabb) / rho - c1o3 * (mfaac + mfaca);

			//Cum 5.
			real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + c4o1* mfabb * mfbbb + c2o1* (mfbab * mfacb + mfbba * mfabc)) / rho - c1o3 * (mfbca + mfbac);
			real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + c4o1* mfbab * mfbbb + c2o1* (mfabb * mfcab + mfbba * mfbac)) / rho - c1o3 * (mfcba + mfabc);
			real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + c4o1* mfbba * mfbbb + c2o1* (mfbab * mfbca + mfabb * mfcba)) / rho - c1o3 * (mfacb + mfcab);

			//Cum 6.
			real CUMccc = mfccc + (-c4o1*  mfbbb * mfbbb
				- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				- c4o1* (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				- c2o1* (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
				+ (c4o1* (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
					+ c2o1* (mfcaa * mfaca * mfaac)
					+ c16o1*  mfbba * mfbab * mfabb) / (rho * rho)
				- c1o3* (mfacc + mfcac + mfcca)
				+ c1o9* (mfcaa + mfaca + mfaac)
				+ (c2o1* (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
					+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3 / rho;


			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy = mfcaa - mfaca;
			real mxxMzz = mfcaa - mfaac;

			//relax
			//hat noch nicht so gut funktioniert...Optimierungsbedarf
			//{
			//	real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) +  OxxPyyPzz * (mfaaa - mxxPyyPzz);
			//	real dyuy = dxux + omega * c3o2 * mxxMyy;
			//	real dzuz = dxux + omega * c3o2 * mxxMzz;

			//	//relax
			//	mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
			//	mxxMyy    += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vy2 * dyuy);
			//	mxxMzz    += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vz2 * dzuz);

			//	//////////////////////////////////////////////////////////////////////////
			//	//limiter-Scheise Teil 2
			//	//mxxMyy    += oxxyy * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vy2 * dyuy);
			//	//mxxMzz    += oxxzz * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vz2 * dzuz);
			//	//////////////////////////////////////////////////////////////////////////

			//}

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
			{
				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
				real dyuy = dxux + omega * c3o2 * mxxMyy;
				real dzuz = dxux + omega * c3o2 * mxxMzz;

				//relax
				mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
				mxxMyy += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
				mxxMzz += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

				//////////////////////////////////////////////////////////////////////////
				//limiter-Scheise Teil 2
				//mxxMyy    += oxxyy * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vy2 * dyuy);
				//mxxMzz    += oxxzz * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vz2 * dzuz);
				//////////////////////////////////////////////////////////////////////////

			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////no correction
			//mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);
			//mxxMyy    += -(-omega) * (-mxxMyy);
			//mxxMzz    += -(-omega) * (-mxxMzz);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb += -(-omega) * (-mfabb);
			mfbab += -(-omega) * (-mfbab);
			mfbba += -(-omega) * (-mfbba);

			// linear combinations back
			mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-c2o1 * mxxMyy + mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (mxxMyy - c2o1 * mxxMzz + mxxPyyPzz);

			//3.
			// linear combinations
			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

			//relax
			wadjust = OxyyMxzz + (c1o1 - OxyyMxzz)*abs(mfbbb) / (abs(mfbbb) + qudricLimit);
			mfbbb += wadjust * (-mfbbb);
			wadjust = OxyyPxzz + (c1o1 - OxyyPxzz)*abs(mxxyPyzz) / (abs(mxxyPyzz) + qudricLimit);
			mxxyPyzz += wadjust * (-mxxyPyzz);
			wadjust = OxyyMxzz + (c1o1 - OxyyMxzz)*abs(mxxyMyzz) / (abs(mxxyMyzz) + qudricLimit);
			mxxyMyzz += wadjust * (-mxxyMyzz);
			wadjust = OxyyPxzz + (c1o1 - OxyyPxzz)*abs(mxxzPyyz) / (abs(mxxzPyyz) + qudricLimit);
			mxxzPyyz += wadjust * (-mxxzPyyz);
			wadjust = OxyyMxzz + (c1o1 - OxyyMxzz)*abs(mxxzMyyz) / (abs(mxxzMyyz) + qudricLimit);
			mxxzMyyz += wadjust * (-mxxzMyyz);
			wadjust = OxyyPxzz + (c1o1 - OxyyPxzz)*abs(mxyyPxzz) / (abs(mxyyPxzz) + qudricLimit);
			mxyyPxzz += wadjust * (-mxyyPxzz);
			wadjust = OxyyMxzz + (c1o1 - OxyyMxzz)*abs(mxyyMxzz) / (abs(mxyyMxzz) + qudricLimit);
			mxyyMxzz += wadjust * (-mxyyMxzz);

			// linear combinations back
			mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			CUMacc += O4 * (-CUMacc);
			CUMcac += O4 * (-CUMcac);
			CUMcca += O4 * (-CUMcca);

			CUMbbc += O4 * (-CUMbbc);
			CUMbcb += O4 * (-CUMbcb);
			CUMcbb += O4 * (-CUMcbb);

			//5.
			CUMbcc += O5 * (-CUMbcc);
			CUMcbc += O5 * (-CUMcbc);
			CUMccb += O5 * (-CUMccb);

			//6.
			CUMccc += O6 * (-CUMccc);

			//back cumulants to central moments
			//4.
			mfcbb = CUMcbb + ((mfcaa + c1o3 * rho) * mfabb + c2o1* mfbba * mfbab) / rho;
			mfbcb = CUMbcb + ((mfaca + c1o3 * rho) * mfbab + c2o1* mfbba * mfabb) / rho;
			mfbbc = CUMbbc + ((mfaac + c1o3 * rho) * mfbba + c2o1* mfbab * mfabb) / rho;

			mfcca = CUMcca + (mfcaa * mfaca + c2o1* mfbba * mfbba) / rho + c1o3 * (mfcaa + mfaca);
			mfcac = CUMcac + (mfcaa * mfaac + c2o1* mfbab * mfbab) / rho + c1o3 * (mfcaa + mfaac);
			mfacc = CUMacc + (mfaac * mfaca + c2o1* mfabb * mfabb) / rho + c1o3 * (mfaac + mfaca);

			//5.
			mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + c4o1* mfabb * mfbbb + c2o1* (mfbab * mfacb + mfbba * mfabc)) / rho + c1o3 * (mfbca + mfbac);
			mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + c4o1* mfbab * mfbbb + c2o1* (mfabb * mfcab + mfbba * mfbac)) / rho + c1o3 * (mfcba + mfabc);
			mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + c4o1* mfbba * mfbbb + c2o1* (mfbab * mfbca + mfabb * mfcba)) / rho + c1o3 * (mfacb + mfcab);

			//6.
			mfccc = CUMccc - ((-c4o1*  mfbbb * mfbbb
				- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				- c4o1* (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				- c2o1* (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
				+ (c4o1* (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
					+ c2o1* (mfcaa * mfaca * mfaac)
					+ c16o1*  mfbba * mfbab * mfabb) / (rho * rho)
				- c1o3* (mfacc + mfcac + mfcca)
				+ c1o9* (mfcaa + mfaca + mfaac)
				+ (c2o1* (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
					+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3 / rho);
			////////////////////////////////////////////////////////////////////////////////////

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
			// Z - Dir
			m0 = mfaac * c1o2 + mfaab * (vvz - c1o2) + (mfaaa + 1.) * (vz2 - vvz) * c1o2;
			m1 = -mfaac - c2o1* mfaab *  vvz + mfaaa       * (c1o1 - vz2) - c1o1* vz2;
			m2 = mfaac * c1o2 + mfaab * (vvz + c1o2) + (mfaaa + 1.) * (vz2 + vvz) * c1o2;
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
			m0 = mfacc * c1o2 + mfacb * (vvz - c1o2) + (mfaca + c1o3) * (vz2 - vvz) * c1o2;
			m1 = -mfacc - c2o1* mfacb *  vvz + mfaca         * (c1o1 - vz2) - c1o3 * vz2;
			m2 = mfacc * c1o2 + mfacb * (vvz + c1o2) + (mfaca + c1o3) * (vz2 + vvz) * c1o2;
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
			m0 = mfcac * c1o2 + mfcab * (vvz - c1o2) + (mfcaa + c1o3) * (vz2 - vvz) * c1o2;
			m1 = -mfcac - c2o1* mfcab *  vvz + mfcaa         * (c1o1 - vz2) - c1o3 * vz2;
			m2 = mfcac * c1o2 + mfcab * (vvz + c1o2) + (mfcaa + c1o3) * (vz2 + vvz) * c1o2;
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
			m0 = mfccc * c1o2 + mfccb * (vvz - c1o2) + (mfcca + c1o9) * (vz2 - vvz) * c1o2;
			m1 = -mfccc - c2o1* mfccb *  vvz + mfcca         * (c1o1 - vz2) - c1o9 * vz2;
			m2 = mfccc * c1o2 + mfccb * (vvz + c1o2) + (mfcca + c1o9) * (vz2 + vvz) * c1o2;
			mfcca = m0;
			mfccb = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m0 = mfaca * c1o2 + mfaba * (vvy - c1o2) + (mfaaa + c1o6) * (vy2 - vvy) * c1o2;
			m1 = -mfaca - c2o1* mfaba *  vvy + mfaaa         * (c1o1 - vy2) - c1o6 * vy2;
			m2 = mfaca * c1o2 + mfaba * (vvy + c1o2) + (mfaaa + c1o6) * (vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfacb * c1o2 + mfabb * (vvy - c1o2) + (mfaab + c2o3) * (vy2 - vvy) * c1o2;
			m1 = -mfacb - c2o1* mfabb *  vvy + mfaab         * (c1o1 - vy2) - c2o3 * vy2;
			m2 = mfacb * c1o2 + mfabb * (vvy + c1o2) + (mfaab + c2o3) * (vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfacc * c1o2 + mfabc * (vvy - c1o2) + (mfaac + c1o6) * (vy2 - vvy) * c1o2;
			m1 = -mfacc - c2o1* mfabc *  vvy + mfaac         * (c1o1 - vy2) - c1o6 * vy2;
			m2 = mfacc * c1o2 + mfabc * (vvy + c1o2) + (mfaac + c1o6) * (vy2 + vvy) * c1o2;
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
			m0 = mfcca * c1o2 + mfcba * (vvy - c1o2) + (mfcaa + c1o18) * (vy2 - vvy) * c1o2;
			m1 = -mfcca - c2o1* mfcba *  vvy + mfcaa          * (c1o1 - vy2) - c1o18 * vy2;
			m2 = mfcca * c1o2 + mfcba * (vvy + c1o2) + (mfcaa + c1o18) * (vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 = mfccb * c1o2 + mfcbb * (vvy - c1o2) + (mfcab + c2o9) * (vy2 - vvy) * c1o2;
			m1 = -mfccb - c2o1* mfcbb *  vvy + mfcab         * (c1o1 - vy2) - c2o9 * vy2;
			m2 = mfccb * c1o2 + mfcbb * (vvy + c1o2) + (mfcab + c2o9) * (vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 = mfccc * c1o2 + mfcbc * (vvy - c1o2) + (mfcac + c1o18) * (vy2 - vvy) * c1o2;
			m1 = -mfccc - c2o1* mfcbc *  vvy + mfcac          * (c1o1 - vy2) - c1o18 * vy2;
			m2 = mfccc * c1o2 + mfcbc * (vvy + c1o2) + (mfcac + c1o18) * (vy2 + vvy) * c1o2;
			mfcac = m0;
			mfcbc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m0 = mfcaa * c1o2 + mfbaa * (vvx - c1o2) + (mfaaa + c1o36) * (vx2 - vvx) * c1o2;
			m1 = -mfcaa - c2o1* mfbaa *  vvx + mfaaa          * (c1o1 - vx2) - c1o36 * vx2;
			m2 = mfcaa * c1o2 + mfbaa * (vvx + c1o2) + (mfaaa + c1o36) * (vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcba * c1o2 + mfbba * (vvx - c1o2) + (mfaba + c1o9) * (vx2 - vvx) * c1o2;
			m1 = -mfcba - c2o1* mfbba *  vvx + mfaba         * (c1o1 - vx2) - c1o9 * vx2;
			m2 = mfcba * c1o2 + mfbba * (vvx + c1o2) + (mfaba + c1o9) * (vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcca * c1o2 + mfbca * (vvx - c1o2) + (mfaca + c1o36) * (vx2 - vvx) * c1o2;
			m1 = -mfcca - c2o1* mfbca *  vvx + mfaca          * (c1o1 - vx2) - c1o36 * vx2;
			m2 = mfcca * c1o2 + mfbca * (vvx + c1o2) + (mfaca + c1o36) * (vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcab * c1o2 + mfbab * (vvx - c1o2) + (mfaab + c1o9) * (vx2 - vvx) * c1o2;
			m1 = -mfcab - c2o1* mfbab *  vvx + mfaab         * (c1o1 - vx2) - c1o9 * vx2;
			m2 = mfcab * c1o2 + mfbab * (vvx + c1o2) + (mfaab + c1o9) * (vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 = mfcbb * c1o2 + mfbbb * (vvx - c1o2) + (mfabb + c4o9) * (vx2 - vvx) * c1o2;
			m1 = -mfcbb - c2o1* mfbbb *  vvx + mfabb         * (c1o1 - vx2) - c4o9 * vx2;
			m2 = mfcbb * c1o2 + mfbbb * (vvx + c1o2) + (mfabb + c4o9) * (vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 = mfccb * c1o2 + mfbcb * (vvx - c1o2) + (mfacb + c1o9) * (vx2 - vvx) * c1o2;
			m1 = -mfccb - c2o1* mfbcb *  vvx + mfacb         * (c1o1 - vx2) - c1o9 * vx2;
			m2 = mfccb * c1o2 + mfbcb * (vvx + c1o2) + (mfacb + c1o9) * (vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcac * c1o2 + mfbac * (vvx - c1o2) + (mfaac + c1o36) * (vx2 - vvx) * c1o2;
			m1 = -mfcac - c2o1* mfbac *  vvx + mfaac          * (c1o1 - vx2) - c1o36 * vx2;
			m2 = mfcac * c1o2 + mfbac * (vvx + c1o2) + (mfaac + c1o36) * (vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 = mfcbc * c1o2 + mfbbc * (vvx - c1o2) + (mfabc + c1o9) * (vx2 - vvx) * c1o2;
			m1 = -mfcbc - c2o1* mfbbc *  vvx + mfabc         * (c1o1 - vx2) - c1o9 * vx2;
			m2 = mfcbc * c1o2 + mfbbc * (vvx + c1o2) + (mfabc + c1o9) * (vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 = mfccc * c1o2 + mfbcc * (vvx - c1o2) + (mfacc + c1o36) * (vx2 - vvx) * c1o2;
			m1 = -mfccc - c2o1* mfbcc *  vvx + mfacc          * (c1o1 - vx2) - c1o36 * vx2;
			m2 = mfccc * c1o2 + mfbcc * (vvx + c1o2) + (mfacc + c1o36) * (vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			(D.f[E])[ke] = mfabb;// -  c2over27 ;//                                                                     
			(D.f[W])[kw] = mfcbb;// -  c2over27 ;                                                                     
			(D.f[N])[kn] = mfbab;// -  c2over27 ;
			(D.f[S])[ks] = mfbcb;// -  c2over27 ;
			(D.f[T])[kt] = mfbba;// -  c2over27 ;
			(D.f[B])[kb] = mfbbc;// -  c2over27 ;
			(D.f[NE])[kne] = mfaab;// -  c1over54 ;
			(D.f[SW])[ksw] = mfccb;// -  c1over54 ;
			(D.f[SE])[kse] = mfacb;// -  c1over54 ;
			(D.f[NW])[knw] = mfcab;// -  c1over54 ;
			(D.f[TE])[kte] = mfaba;// -  c1over54 ;
			(D.f[BW])[kbw] = mfcbc;// -  c1over54 ;
			(D.f[BE])[kbe] = mfabc;// -  c1over54 ;
			(D.f[TW])[ktw] = mfcba;// -  c1over54 ;
			(D.f[TN])[ktn] = mfbaa;// -  c1over54 ;
			(D.f[BS])[kbs] = mfbcc;// -  c1over54 ;
			(D.f[BN])[kbn] = mfbac;// -  c1over54 ;
			(D.f[TS])[kts] = mfbca;// -  c1over54 ;
			(D.f[REST])[kzero] = mfbbb;// -  c8over27 ;
			(D.f[TNE])[ktne] = mfaaa;// -  c1over216;
			(D.f[TSE])[ktse] = mfaca;// -  c1over216;
			(D.f[BNE])[kbne] = mfaac;// -  c1over216;
			(D.f[BSE])[kbse] = mfacc;// -  c1over216;
			(D.f[TNW])[ktnw] = mfcaa;// -  c1over216;
			(D.f[TSW])[ktsw] = mfcca;// -  c1over216;
			(D.f[BNW])[kbnw] = mfcac;// -  c1over216;
			(D.f[BSW])[kbsw] = mfccc;// -  c1over216;
										////////////////////////////////////////////////////////////////////////////////////
		}
	}
}