//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
/* Device code */
#include "LBM/D3Q27.h"
#include "math.h"
#include "GPU/constant.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Cumulant_D3Q27All4(real omega,
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

	if (k < size_Mat)
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

			//unsigned int kzero= k;
			//unsigned int ke   = k;
			//unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			//unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			//unsigned int kb   = neighborZ[k];
			//unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			//unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			//unsigned int kbs  = neighborZ[ks];
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
			//unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

			real rho = one + drho;
			////////////////////////////////////////////////////////////////////////////////////
			//slow
			//real oMdrho = one - ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
			//					   (((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
			//						((mfabb+mfcbb) + (mfbab+mfbcb)  +  (mfbba+mfbbc)));//fehlt mfbbb
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
			real fx = forces[0] / (pow(two, level)); //zero;//0.0032653/(pow(two,level)); //0.000000005;//(two/1600000.0) / 120.0; //
			real fy = forces[1] / (pow(two, level)); //zero;
			real fz = forces[2] / (pow(two, level)); //zero;
			vvx += fx;
			vvy += fy;
			vvz += fz;
			////////////////////////////////////////////////////////////////////////////////////
			//real omega = omega_in;
			////////////////////////////////////////////////////////////////////////////////////
			//fast
			real oMdrho = one; // comp special
								  //real oMdrho = one - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
								  //					   mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
								  //					   mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb);//fehlt mfbbb nicht mehr
								  //real vvx    =mfccc-mfaaa + mfcac-mfaca + mfcaa-mfacc + mfcca-mfaac + 
								  //				mfcba-mfabc + mfcbc-mfaba + mfcab-mfacb + mfccb-mfaab +
								  //				mfcbb-mfabb;
								  //real vvy    =mfccc-mfaaa + mfaca-mfcac + mfacc-mfcaa + mfcca-mfaac + 
								  //				mfbca-mfbac + mfbcc-mfbaa + mfacb-mfcab + mfccb-mfaab +
								  //				mfbcb-mfbab;
								  //real vvz    =mfccc-mfaaa + mfcac-mfaca + mfacc-mfcaa + mfaac-mfcca + 
								  //				mfbac-mfbca + mfbcc-mfbaa + mfabc-mfcba + mfcbc-mfaba +
								  //				mfbbc-mfbba;
								  ////////////////////////////////////////////////////////////////////////////////////
								  // oMdrho assembler style -------> faaaaaastaaaa
								  // or much sloooowaaaa ... it depändssssss on sadaku
			real m0, m1, m2;
			//real oMdrho;
			//{
			//	oMdrho=mfccc+mfaaa;
			//	m0=mfaca+mfcac;
			//	m1=mfacc+mfcaa;
			//	m2=mfaac+mfcca;
			//	oMdrho+=m0;
			//	m1+=m2;
			//	oMdrho+=m1;
			//	m0=mfbac+mfbca;
			//	m1=mfbaa+mfbcc;
			//	m0+=m1;
			//	m1=mfabc+mfcba;
			//	m2=mfaba+mfcbc;
			//	m1+=m2;
			//	m0+=m1;
			//	m1=mfacb+mfcab;
			//	m2=mfaab+mfccb;
			//	m1+=m2;
			//	m0+=m1;
			//	oMdrho+=m0;
			//	m0=mfabb+mfcbb;
			//	m1=mfbab+mfbcb;
			//	m2=mfbba+mfbbc;
			//	m0+=m1+m2;
			//	m0+=mfbbb; //hat gefehlt
			//	oMdrho = one - (oMdrho + m0);
			//}
			//real vvx;
			real vx2;
			//{
			//	vvx = mfccc-mfaaa;
			//	m0  = mfcac-mfaca;
			//	m1  = mfcaa-mfacc;
			//	m2  = mfcca-mfaac;
			//	vvx+= m0;
			//	m1 += m2;
			//	vvx+= m1;
			//	vx2 = mfcba-mfabc;
			//	m0  = mfcbc-mfaba;
			//	m1  = mfcab-mfacb;
			//	m2  = mfccb-mfaab;
			//	vx2+= m0;
			//	m1 += m2;
			//	vx2+= m1;
			//	vvx+= vx2;
			//	vx2 = mfcbb-mfabb;
			//	vvx+= vx2;
			//}
			//real vvy;
			real vy2;
			//{
			//	vvy = mfccc-mfaaa;
			//	m0  = mfaca-mfcac;
			//	m1  = mfacc-mfcaa;
			//	m2  = mfcca-mfaac;
			//	vvy+= m0;
			//	m1 += m2;
			//	vvy+= m1;
			//	vy2 = mfbca-mfbac;
			//	m0  = mfbcc-mfbaa;
			//	m1  = mfacb-mfcab;
			//	m2  = mfccb-mfaab;
			//	vy2+= m0;
			//	m1 += m2;
			//	vy2+= m1;
			//	vvy+= vy2;
			//	vy2 = mfbcb-mfbab;
			//	vvy+= vy2;
			//}
			//real vvz;
			real vz2;
			//{
			//	vvz = mfccc-mfaaa;
			//	m0  = mfcac-mfaca;
			//	m1  = mfacc-mfcaa;
			//	m2  = mfaac-mfcca;
			//	vvz+= m0;
			//	m1 += m2;
			//	vvz+= m1;
			//	vz2 = mfbac-mfbca;
			//	m0  = mfbcc-mfbaa;
			//	m1  = mfabc-mfcba;
			//	m2  = mfcbc-mfaba;
			//	vz2+= m0;
			//	m1 += m2;
			//	vz2+= m1;
			//	vvz+= vz2;
			//	vz2 = mfbbc-mfbba;
			//	vvz+= vz2;
			//}
			vx2 = vvx*vvx;
			vy2 = vvy*vvy;
			vz2 = vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real wadjust;
			real qudricLimitP = 1.0e10;// 0.01f;// * 0.0001f;
			real qudricLimitM = 1.0e10;//0.01f;// * 0.0001f;
			real qudricLimitD = 1.0e10;//0.01f;// * 0.001f;
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
			mfaac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaba + mfabc;
			m1 = mfabc - mfaba;
			m0 = m2 + mfabb;
			mfaba = m0;
			m0 += c1o9 * oMdrho;
			mfabb = m1 - m0 * vvz;
			mfabc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaca + mfacc;
			m1 = mfacc - mfaca;
			m0 = m2 + mfacb;
			mfaca = m0;
			m0 += c1o36 * oMdrho;
			mfacb = m1 - m0 * vvz;
			mfacc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbaa + mfbac;
			m1 = mfbac - mfbaa;
			m0 = m2 + mfbab;
			mfbaa = m0;
			m0 += c1o9 * oMdrho;
			mfbab = m1 - m0 * vvz;
			mfbac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbba + mfbbc;
			m1 = mfbbc - mfbba;
			m0 = m2 + mfbbb;
			mfbba = m0;
			m0 += c4o9 * oMdrho;
			mfbbb = m1 - m0 * vvz;
			mfbbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbca + mfbcc;
			m1 = mfbcc - mfbca;
			m0 = m2 + mfbcb;
			mfbca = m0;
			m0 += c1o9 * oMdrho;
			mfbcb = m1 - m0 * vvz;
			mfbcc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcaa + mfcac;
			m1 = mfcac - mfcaa;
			m0 = m2 + mfcab;
			mfcaa = m0;
			m0 += c1o36 * oMdrho;
			mfcab = m1 - m0 * vvz;
			mfcac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcba + mfcbc;
			m1 = mfcbc - mfcba;
			m0 = m2 + mfcbb;
			mfcba = m0;
			m0 += c1o9 * oMdrho;
			mfcbb = m1 - m0 * vvz;
			mfcbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcca + mfccc;
			m1 = mfccc - mfcca;
			m0 = m2 + mfccb;
			mfcca = m0;
			m0 += c1o36 * oMdrho;
			mfccb = m1 - m0 * vvz;
			mfccc = m2 - two*	m1 * vvz + vz2 * m0;
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
			mfaca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaab + mfacb;
			m1 = mfacb - mfaab;
			m0 = m2 + mfabb;
			mfaab = m0;
			mfabb = m1 - m0 * vvy;
			mfacb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaac + mfacc;
			m1 = mfacc - mfaac;
			m0 = m2 + mfabc;
			mfaac = m0;
			m0 += c1o18 * oMdrho;
			mfabc = m1 - m0 * vvy;
			mfacc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbaa + mfbca;
			m1 = mfbca - mfbaa;
			m0 = m2 + mfbba;
			mfbaa = m0;
			m0 += c2o3 * oMdrho;
			mfbba = m1 - m0 * vvy;
			mfbca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbab + mfbcb;
			m1 = mfbcb - mfbab;
			m0 = m2 + mfbbb;
			mfbab = m0;
			mfbbb = m1 - m0 * vvy;
			mfbcb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbac + mfbcc;
			m1 = mfbcc - mfbac;
			m0 = m2 + mfbbc;
			mfbac = m0;
			m0 += c2o9 * oMdrho;
			mfbbc = m1 - m0 * vvy;
			mfbcc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcaa + mfcca;
			m1 = mfcca - mfcaa;
			m0 = m2 + mfcba;
			mfcaa = m0;
			m0 += c1o6 * oMdrho;
			mfcba = m1 - m0 * vvy;
			mfcca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcab + mfccb;
			m1 = mfccb - mfcab;
			m0 = m2 + mfcbb;
			mfcab = m0;
			mfcbb = m1 - m0 * vvy;
			mfccb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcac + mfccc;
			m1 = mfccc - mfcac;
			m0 = m2 + mfcbc;
			mfcac = m0;
			m0 += c1o18 * oMdrho;
			mfcbc = m1 - m0 * vvy;
			mfccc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2 = mfaaa + mfcaa;
			m1 = mfcaa - mfaaa;
			m0 = m2 + mfbaa;
			mfaaa = m0;
			m0 += one* oMdrho;
			mfbaa = m1 - m0 * vvx;
			mfcaa = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaba + mfcba;
			m1 = mfcba - mfaba;
			m0 = m2 + mfbba;
			mfaba = m0;
			mfbba = m1 - m0 * vvx;
			mfcba = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaca + mfcca;
			m1 = mfcca - mfaca;
			m0 = m2 + mfbca;
			mfaca = m0;
			m0 += c1o3 * oMdrho;
			mfbca = m1 - m0 * vvx;
			mfcca = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaab + mfcab;
			m1 = mfcab - mfaab;
			m0 = m2 + mfbab;
			mfaab = m0;
			mfbab = m1 - m0 * vvx;
			mfcab = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfabb + mfcbb;
			m1 = mfcbb - mfabb;
			m0 = m2 + mfbbb;
			mfabb = m0;
			mfbbb = m1 - m0 * vvx;
			mfcbb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfacb + mfccb;
			m1 = mfccb - mfacb;
			m0 = m2 + mfbcb;
			mfacb = m0;
			mfbcb = m1 - m0 * vvx;
			mfccb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaac + mfcac;
			m1 = mfcac - mfaac;
			m0 = m2 + mfbac;
			mfaac = m0;
			m0 += c1o3 * oMdrho;
			mfbac = m1 - m0 * vvx;
			mfcac = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfabc + mfcbc;
			m1 = mfcbc - mfabc;
			m0 = m2 + mfbbc;
			mfabc = m0;
			mfbbc = m1 - m0 * vvx;
			mfcbc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfacc + mfccc;
			m1 = mfccc - mfacc;
			m0 = m2 + mfbcc;
			mfacc = m0;
			m0 += c1o9 * oMdrho;
			mfbcc = m1 - m0 * vvx;
			mfccc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////
			//limiter advect
			real lambdaAdvect = 1000;//100000;
			real limAdvect = lambdaAdvect * (one / omega - c1o2) * (one / omega - c1o2);
			////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = one;	//set the bulk viscosity one is high / two is very low and zero is (too) high ... (also called omega 2)

			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			real OxyyPxzz = eight*(-two + omega)*(one + two*omega) / (-eight - fourteen*omega + seven*omega*omega);//one;
			real OxyyMxzz = eight*(-two + omega)*(-seven + four*omega) / (fiftysix - fifty*omega + nine*omega*omega);//one;
			real Oxyz = twentyfour*(-two + omega)*(-two - seven*omega + three*omega*omega) / (fourtyeight + c152*omega - c130*omega*omega + twentynine*omega*omega*omega);//one;
			////////////////////////////////////////////////////////////
			//4.
			//////////////////////////////
			//real O4 = one;
			//////////////////////////////
			real O4        = one/(100*(one/ omega-c1o2)+c1o2);//TRT
			////////////////////////////////////////////////////////////
			//5.
			//////////////////////////////
			real O5 =  OxyyPxzz;// two - omega; // one;
			////////////////////////////////////////////////////////////
			//6.
			//////////////////////////////
			real O6 = one;
			////////////////////////////////////////////////////////////

			real A = (four * omega * omega + two * omega * OxxPyyPzz * (omega - six) + OxxPyyPzz * OxxPyyPzz * (omega * (ten - three * omega) - four)) /
				((omega - OxxPyyPzz) * (OxxPyyPzz * (two + three * omega) - eight * omega));
			real B = (four * omega * OxxPyyPzz * (nine * omega - sixteen) - four * omega * omega - two * OxxPyyPzz * OxxPyyPzz * (two + nine * omega * (omega - two))) /
				(three * (omega - OxxPyyPzz) * (OxxPyyPzz * (two + three * omega) - eight * omega));


			//central moments to cumulants
			//4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho;
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho;
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho;

			real CUMcca = mfcca - (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
			real CUMcac = mfcac - (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
			real CUMacc = mfacc - (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

			//5.
			real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)*(one + rho*six*B / (two + three * B))) / rho;
			real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)*(one + rho*six*B / (two + three * B))) / rho;
			real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)*(one + rho*six*B / (two + three * B))) / rho;

			real scaleD2 = 
				(-twentyseven * omega*(omega - OxxPyyPzz)*(two * OxxPyyPzz + omega*(-eight + three * OxxPyyPzz))*(-four * OxxPyyPzz + omega*(-two + three * OxxPyyPzz))*O5) /
				((-two + omega)*(six * OxxPyyPzz*(thirtytwo + (c100 - sixtynine * OxxPyyPzz)*OxxPyyPzz)*(omega*omega) +
				(fourty + three * OxxPyyPzz*(-fiftytwo + three * OxxPyyPzz*(-eight + nine * OxxPyyPzz)))*(omega*omega*omega) +
					c264 * omega*(-three + two * OxxPyyPzz)*(OxxPyyPzz*OxxPyyPzz) - eightyeight * (OxxPyyPzz*OxxPyyPzz*OxxPyyPzz)));

			real dxxux = CUMbcc*scaleD2;
			real dyyuy = CUMcbc*scaleD2;
			real dzzuz = CUMccb*scaleD2;

			//dxxux *= limAdvect / (limAdvect + sqrt(abs(vvx * dxxux)));
			//dyyuy *= limAdvect / (limAdvect + sqrt(abs(vvy * dyyuy)));
			//dzzuz *= limAdvect / (limAdvect + sqrt(abs(vvz * dzzuz)));

			//real limD2 = 1.0e10;
			//dxxux /= one + limD2*abs(dxxux);
			//dyyuy /= one + limD2*abs(dyyuy);
			//dzzuz /= one + limD2*abs(dzzuz);

			//6.

			real CUMccc = mfccc + ((-four *  mfbbb * mfbbb
				- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				- four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				- two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
				+ (four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
					+ two * (mfcaa * mfaca * mfaac)
					+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
				- c1o3 * (mfacc + mfcac + mfcca) / rho
				- c1o9 * (mfcaa + mfaca + mfaac) / rho
				+ (two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
					+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3
				+ c1o27*((drho * drho - drho) / (rho*rho)));

			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy = mfcaa - mfaca;
			real mxxMzz = mfcaa - mfaac;


			////////////////////////////////////////////////////////////////////////////
			real Dxy = -three*omega*mfbba;
			real Dxz = -three*omega*mfbab;
			real Dyz = -three*omega*mfabb;

			//3.
			// linear combinations

			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)

			real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
			real dyuy = dxux + omega * c3o2 * mxxMyy;
			real dzuz = dxux + omega * c3o2 * mxxMzz;

			//////////////////////////////////////////////////////////////////////////////
			real divTest = (mxxPyyPzz - (mfaaa + (-three * (one - OxxPyyPzz*c1o2)*(vx2*dxux + vy2*dyuy + vz2*dzuz) + (six - three * (omega + OxxPyyPzz) + omega*OxxPyyPzz) / 
				(three * omega)*(dxux*dxux / rho + dyuy*dyuy / rho + dzuz*dzuz / rho + vvx*dxxux + vvy*dyyuy + vvz*dzzuz)) / OxxPyyPzz)) / OxxPyyPzz;

			//dxxux *= limAdvect / (limAdvect + sqrt(abs(divTest)));
			//dyyuy *= limAdvect / (limAdvect + sqrt(abs(divTest)));
			//dzzuz *= limAdvect / (limAdvect + sqrt(abs(divTest)));

			//dxxux *= exp(-sqrt(fabs(divTest))/limAdvect);
			//dyyuy *= exp(-sqrt(fabs(divTest)) / limAdvect);
			//dzzuz *= exp(-sqrt(fabs(divTest)) / limAdvect);

			////////////////////////////////////////////////////////////////////////////
			//relax
			//mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz) - three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
			//mxxMyy += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
			//mxxMzz += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
			mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz) - three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz)
					  + (six - three * (omega + OxxPyyPzz) + omega * OxxPyyPzz) / (three * omega) * ((dxux * dxux + dyuy * dyuy + dzuz * dzuz) / rho + vvx * dxxux + vvy * dyyuy + vvz * dzzuz);
			mxxMyy += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy)
				      +  omega * (two*(one / omega - c1o2) * (one / omega - c1o2) - c1o6) * ((dxux * dxux - dyuy * dyuy) / rho + vvx * dxxux - vvy * dyyuy);
			mxxMzz += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz)
					  + omega * (two*(one / omega - c1o2) * (one / omega - c1o2) - c1o6) * ((dxux * dxux - dzuz * dzuz) / rho + vvx * dxxux - vvz * dzzuz);

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////no correction
			//mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);//-magicBulk*OxxPyyPzz;
			//mxxMyy    += -(-omega) * (-mxxMyy);
			//mxxMzz    += -(-omega) * (-mxxMzz);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb += omega * (-mfabb);
			mfbab += omega * (-mfbab);
			mfbba += omega * (-mfbba);
			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-two*  mxxMyy + mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (mxxMyy - two* mxxMzz + mxxPyyPzz);


			//relax
			//////////////////////////////////////////////////////////////////////////
			//das ist der limiter
			wadjust = Oxyz + (one - Oxyz)*abs(mfbbb) / (abs(mfbbb) + qudricLimitD);
			mfbbb += wadjust * (-mfbbb);
			wadjust = OxyyPxzz + (one - OxyyPxzz)*abs(mxxyPyzz) / (abs(mxxyPyzz) + qudricLimitP);
			mxxyPyzz += wadjust * (-mxxyPyzz);
			wadjust = OxyyMxzz + (one - OxyyMxzz)*abs(mxxyMyzz) / (abs(mxxyMyzz) + qudricLimitM);
			mxxyMyzz += wadjust * (-mxxyMyzz);
			wadjust = OxyyPxzz + (one - OxyyPxzz)*abs(mxxzPyyz) / (abs(mxxzPyyz) + qudricLimitP);
			mxxzPyyz += wadjust * (-mxxzPyyz);
			wadjust = OxyyMxzz + (one - OxyyMxzz)*abs(mxxzMyyz) / (abs(mxxzMyyz) + qudricLimitM);
			mxxzMyyz += wadjust * (-mxxzMyyz);
			wadjust = OxyyPxzz + (one - OxyyPxzz)*abs(mxyyPxzz) / (abs(mxyyPxzz) + qudricLimitP);
			mxyyPxzz += wadjust * (-mxyyPxzz);
			wadjust = OxyyMxzz + (one - OxyyMxzz)*abs(mxyyMxzz) / (abs(mxyyMxzz) + qudricLimitM);
			mxyyMxzz += wadjust * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
			//mfbbb     += OxyyMxzz * (-mfbbb);
			//mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
			//mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
			//mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
			//mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
			//mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
			//mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			//////////////////////////////////////////////////////////////////////////
			//mit limiter
			//	wadjust    = O4+(one-O4)*abs(CUMacc)/(abs(CUMacc)+qudricLimit);
			//CUMacc    += wadjust * (-CUMacc);
			//	wadjust    = O4+(one-O4)*abs(CUMcac)/(abs(CUMcac)+qudricLimit);
			//CUMcac    += wadjust * (-CUMcac); 
			//	wadjust    = O4+(one-O4)*abs(CUMcca)/(abs(CUMcca)+qudricLimit);
			//CUMcca    += wadjust * (-CUMcca); 

			//	wadjust    = O4+(one-O4)*abs(CUMbbc)/(abs(CUMbbc)+qudricLimit);
			//CUMbbc    += wadjust * (-CUMbbc); 
			//	wadjust    = O4+(one-O4)*abs(CUMbcb)/(abs(CUMbcb)+qudricLimit);
			//CUMbcb    += wadjust * (-CUMbcb); 
			//	wadjust    = O4+(one-O4)*abs(CUMcbb)/(abs(CUMcbb)+qudricLimit);
			//CUMcbb    += wadjust * (-CUMcbb); 
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
			//CUMacc += O4 * (-CUMacc);
			//CUMcac += O4 * (-CUMcac);
			//CUMcca += O4 * (-CUMcca);
			//CUMbbc += O4 * (-CUMbbc);
			//CUMbcb += O4 * (-CUMbcb);
			//CUMcbb += O4 * (-CUMcbb);
			//CUMacc = -O4*(one / omega - c1o2)*(dyuy + dzuz)*c2o3 *(four + two*omega - three*omega*omega) / (two - seven*omega + five*omega*omega) + (one - O4) * (CUMacc);
			//CUMcac = -O4*(one / omega - c1o2)*(dxux + dzuz)*c2o3 *(four + two*omega - three*omega*omega) / (two - seven*omega + five*omega*omega) + (one - O4) * (CUMcac);
			//CUMcca = -O4*(one / omega - c1o2)*(dyuy + dxux)*c2o3 *(four + two*omega - three*omega*omega) / (two - seven*omega + five*omega*omega) + (one - O4) * (CUMcca);
			//CUMbbc = -O4*(one / omega - c1o2)*Dxy*c1o3 *(four + twentyeight*omega - fourteen*omega*omega) / (six - twentyone*omega + fiveteen*omega*omega) + (one - O4) * (CUMbbc);
			//CUMbcb = -O4*(one / omega - c1o2)*Dxz*c1o3 *(four + twentyeight*omega - fourteen*omega*omega) / (six - twentyone*omega + fiveteen*omega*omega) + (one - O4) * (CUMbcb);
			//CUMcbb = -O4*(one / omega - c1o2)*Dyz*c1o3 *(four + twentyeight*omega - fourteen*omega*omega) / (six - twentyone*omega + fiveteen*omega*omega) + (one - O4) * (CUMcbb);
			CUMacc = -O4 * (one / omega - c1o2) * (dyuy + dzuz) * c2o3 * A + (one - O4) * (CUMacc);
			CUMcac = -O4 * (one / omega - c1o2) * (dxux + dzuz) * c2o3 * A + (one - O4) * (CUMcac);
			CUMcca = -O4 * (one / omega - c1o2) * (dyuy + dxux) * c2o3 * A + (one - O4) * (CUMcca);
			CUMbbc = -O4 * (one / omega - c1o2) * Dxy           * c1o3 * B + (one - O4) * (CUMbbc);
			CUMbcb = -O4 * (one / omega - c1o2) * Dxz           * c1o3 * B + (one - O4) * (CUMbcb);
			CUMcbb = -O4 * (one / omega - c1o2) * Dyz           * c1o3 * B + (one - O4) * (CUMcbb);
			//////////////////////////////////////////////////////////////////////////


			//5.
			CUMbcc += O5 * (-CUMbcc);
			CUMcbc += O5 * (-CUMcbc);
			CUMccb += O5 * (-CUMccb);

			//6.
			CUMccc += O6 * (-CUMccc);



			//back cumulants to central moments
			//4.
			mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho;
			mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho;
			mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho;

			mfcca = CUMcca + (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
			mfcac = CUMcac + (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
			mfacc = CUMacc + (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

			//5.
			mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)*(one + rho*six*B / (two + three * B))) / rho;
			mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)*(one + rho*six*B / (two + three * B))) / rho;
			mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)*(one + rho*six*B / (two + three * B))) / rho;

			//6.

			mfccc = CUMccc - ((-four *  mfbbb * mfbbb
				- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				- four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				- two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
				+ (four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
					+ two * (mfcaa * mfaca * mfaac)
					+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
				- c1o3 * (mfacc + mfcac + mfcca) / rho
				- c1o9 * (mfcaa + mfaca + mfaac) / rho
				+ (two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
					+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3
				+ c1o27*((drho * drho - drho) / (rho*rho)));
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
			m0 = mfaac * c1o2 + mfaab * (vvz - c1o2) + (mfaaa + one* oMdrho) * (vz2 - vvz) * c1o2;
			m1 = -mfaac - two* mfaab *  vvz + mfaaa                * (one - vz2) - one* oMdrho * vz2;
			m2 = mfaac * c1o2 + mfaab * (vvz + c1o2) + (mfaaa + one* oMdrho) * (vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfabc * c1o2 + mfabb * (vvz - c1o2) + mfaba * (vz2 - vvz) * c1o2;
			m1 = -mfabc - two* mfabb *  vvz + mfaba * (one - vz2);
			m2 = mfabc * c1o2 + mfabb * (vvz + c1o2) + mfaba * (vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfacc * c1o2 + mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
			m1 = -mfacc - two* mfacb *  vvz + mfaca                  * (one - vz2) - c1o3 * oMdrho * vz2;
			m2 = mfacc * c1o2 + mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfbac * c1o2 + mfbab * (vvz - c1o2) + mfbaa * (vz2 - vvz) * c1o2;
			m1 = -mfbac - two* mfbab *  vvz + mfbaa * (one - vz2);
			m2 = mfbac * c1o2 + mfbab * (vvz + c1o2) + mfbaa * (vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 = mfbbc * c1o2 + mfbbb * (vvz - c1o2) + mfbba * (vz2 - vvz) * c1o2;
			m1 = -mfbbc - two* mfbbb *  vvz + mfbba * (one - vz2);
			m2 = mfbbc * c1o2 + mfbbb * (vvz + c1o2) + mfbba * (vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 = mfbcc * c1o2 + mfbcb * (vvz - c1o2) + mfbca * (vz2 - vvz) * c1o2;
			m1 = -mfbcc - two* mfbcb *  vvz + mfbca * (one - vz2);
			m2 = mfbcc * c1o2 + mfbcb * (vvz + c1o2) + mfbca * (vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcac * c1o2 + mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
			m1 = -mfcac - two* mfcab *  vvz + mfcaa                  * (one - vz2) - c1o3 * oMdrho * vz2;
			m2 = mfcac * c1o2 + mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 = mfcbc * c1o2 + mfcbb * (vvz - c1o2) + mfcba * (vz2 - vvz) * c1o2;
			m1 = -mfcbc - two* mfcbb *  vvz + mfcba * (one - vz2);
			m2 = mfcbc * c1o2 + mfcbb * (vvz + c1o2) + mfcba * (vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 = mfccc * c1o2 + mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 - vvz) * c1o2;
			m1 = -mfccc - two* mfccb *  vvz + mfcca                  * (one - vz2) - c1o9 * oMdrho * vz2;
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
			m1 = -mfaca - two* mfaba *  vvy + mfaaa                  * (one - vy2) - c1o6 * oMdrho * vy2;
			m2 = mfaca * c1o2 + mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfacb * c1o2 + mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 - vvy) * c1o2;
			m1 = -mfacb - two* mfabb *  vvy + mfaab                  * (one - vy2) - c2o3 * oMdrho * vy2;
			m2 = mfacb * c1o2 + mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfacc * c1o2 + mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
			m1 = -mfacc - two* mfabc *  vvy + mfaac                  * (one - vy2) - c1o6 * oMdrho * vy2;
			m2 = mfacc * c1o2 + mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfbca * c1o2 + mfbba * (vvy - c1o2) + mfbaa * (vy2 - vvy) * c1o2;
			m1 = -mfbca - two* mfbba *  vvy + mfbaa * (one - vy2);
			m2 = mfbca * c1o2 + mfbba * (vvy + c1o2) + mfbaa * (vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 = mfbcb * c1o2 + mfbbb * (vvy - c1o2) + mfbab * (vy2 - vvy) * c1o2;
			m1 = -mfbcb - two* mfbbb *  vvy + mfbab * (one - vy2);
			m2 = mfbcb * c1o2 + mfbbb * (vvy + c1o2) + mfbab * (vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 = mfbcc * c1o2 + mfbbc * (vvy - c1o2) + mfbac * (vy2 - vvy) * c1o2;
			m1 = -mfbcc - two* mfbbc *  vvy + mfbac * (one - vy2);
			m2 = mfbcc * c1o2 + mfbbc * (vvy + c1o2) + mfbac * (vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcca * c1o2 + mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
			m1 = -mfcca - two* mfcba *  vvy + mfcaa                   * (one - vy2) - c1o18 * oMdrho * vy2;
			m2 = mfcca * c1o2 + mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 = mfccb * c1o2 + mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 - vvy) * c1o2;
			m1 = -mfccb - two* mfcbb *  vvy + mfcab                  * (one - vy2) - c2o9 * oMdrho * vy2;
			m2 = mfccb * c1o2 + mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 = mfccc * c1o2 + mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
			m1 = -mfccc - two* mfcbc *  vvy + mfcac                   * (one - vy2) - c1o18 * oMdrho * vy2;
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
			m1 = -mfcaa - two* mfbaa *  vvx + mfaaa                   * (one - vx2) - c1o36 * oMdrho * vx2;
			m2 = mfcaa * c1o2 + mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcba * c1o2 + mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcba - two* mfbba *  vvx + mfaba                  * (one - vx2) - c1o9 * oMdrho * vx2;
			m2 = mfcba * c1o2 + mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcca * c1o2 + mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcca - two* mfbca *  vvx + mfaca                   * (one - vx2) - c1o36 * oMdrho * vx2;
			m2 = mfcca * c1o2 + mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcab * c1o2 + mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcab - two* mfbab *  vvx + mfaab                  * (one - vx2) - c1o9 * oMdrho * vx2;
			m2 = mfcab * c1o2 + mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 = mfcbb * c1o2 + mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcbb - two* mfbbb *  vvx + mfabb                  * (one - vx2) - c4o9 * oMdrho * vx2;
			m2 = mfcbb * c1o2 + mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 = mfccb * c1o2 + mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfccb - two* mfbcb *  vvx + mfacb                  * (one - vx2) - c1o9 * oMdrho * vx2;
			m2 = mfccb * c1o2 + mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 = mfcac * c1o2 + mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcac - two* mfbac *  vvx + mfaac                   * (one - vx2) - c1o36 * oMdrho * vx2;
			m2 = mfcac * c1o2 + mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 = mfcbc * c1o2 + mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfcbc - two* mfbbc *  vvx + mfabc                  * (one - vx2) - c1o9 * oMdrho * vx2;
			m2 = mfcbc * c1o2 + mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 = mfccc * c1o2 + mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
			m1 = -mfccc - two* mfbcc *  vvx + mfacc                   * (one - vx2) - c1o36 * oMdrho * vx2;
			m2 = mfccc * c1o2 + mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////

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
extern "C" __global__ void LB_Kernel_Kum_AA2016_Comp_Bulk_SP_27(real omega,
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

	if(k<size_Mat)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			//unsigned int kzero= k;
			//unsigned int ke   = k;
			unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			unsigned int kbs  = neighborZ[ks];
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
			real mfcbb = (D.f[dirE   ])[k  ];//[ke   ];// +  c2over27 ;(D.f[dirE   ])[k  ];//ke
			real mfabb = (D.f[dirW   ])[kw ];//[kw   ];// +  c2over27 ;(D.f[dirW   ])[kw ];
			real mfbcb = (D.f[dirN   ])[k  ];//[kn   ];// +  c2over27 ;(D.f[dirN   ])[k  ];//kn
			real mfbab = (D.f[dirS   ])[ks ];//[ks   ];// +  c2over27 ;(D.f[dirS   ])[ks ];
			real mfbbc = (D.f[dirT   ])[k  ];//[kt   ];// +  c2over27 ;(D.f[dirT   ])[k  ];//kt
			real mfbba = (D.f[dirB   ])[kb ];//[kb   ];// +  c2over27 ;(D.f[dirB   ])[kb ];
			real mfccb = (D.f[dirNE  ])[k  ];//[kne  ];// +  c1over54 ;(D.f[dirNE  ])[k  ];//kne
			real mfaab = (D.f[dirSW  ])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dirSW  ])[ksw];
			real mfcab = (D.f[dirSE  ])[ks ];//[kse  ];// +  c1over54 ;(D.f[dirSE  ])[ks ];//kse
			real mfacb = (D.f[dirNW  ])[kw ];//[knw  ];// +  c1over54 ;(D.f[dirNW  ])[kw ];//knw
			real mfcbc = (D.f[dirTE  ])[k  ];//[kte  ];// +  c1over54 ;(D.f[dirTE  ])[k  ];//kte
			real mfaba = (D.f[dirBW  ])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dirBW  ])[kbw];
			real mfcba = (D.f[dirBE  ])[kb ];//[kbe  ];// +  c1over54 ;(D.f[dirBE  ])[kb ];//kbe
			real mfabc = (D.f[dirTW  ])[kw ];//[ktw  ];// +  c1over54 ;(D.f[dirTW  ])[kw ];//ktw
			real mfbcc = (D.f[dirTN  ])[k  ];//[ktn  ];// +  c1over54 ;(D.f[dirTN  ])[k  ];//ktn
			real mfbaa = (D.f[dirBS  ])[kbs];//[kbs  ];// +  c1over54 ;(D.f[dirBS  ])[kbs];
			real mfbca = (D.f[dirBN  ])[kb ];//[kbn  ];// +  c1over54 ;(D.f[dirBN  ])[kb ];//kbn
			real mfbac = (D.f[dirTS  ])[ks ];//[kts  ];// +  c1over54 ;(D.f[dirTS  ])[ks ];//kts
			real mfbbb = (D.f[dirZERO])[k  ];//[kzero];// +  c8over27 ;(D.f[dirZERO])[k  ];//kzero
			real mfccc = (D.f[dirTNE ])[k  ];//[ktne ];// +  c1over216;(D.f[dirTNE ])[k  ];//ktne
			real mfaac = (D.f[dirTSW ])[ksw];//[ktsw ];// +  c1over216;(D.f[dirTSW ])[ksw];//ktsw
			real mfcac = (D.f[dirTSE ])[ks ];//[ktse ];// +  c1over216;(D.f[dirTSE ])[ks ];//ktse
			real mfacc = (D.f[dirTNW ])[kw ];//[ktnw ];// +  c1over216;(D.f[dirTNW ])[kw ];//ktnw
			real mfcca = (D.f[dirBNE ])[kb ];//[kbne ];// +  c1over216;(D.f[dirBNE ])[kb ];//kbne
			real mfaaa = (D.f[dirBSW ])[kbsw];//[kbsw ];// +  c1over216;(D.f[dirBSW ])[kbsw];
			real mfcaa = (D.f[dirBSE ])[kbs];//[kbse ];// +  c1over216;(D.f[dirBSE ])[kbs];//kbse
			real mfaca = (D.f[dirBNW ])[kbw];//[kbnw ];// +  c1over216;(D.f[dirBNW ])[kbw];//kbnw
			////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
							(((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
							((mfabb+mfcbb) + (mfbab+mfbcb) + (mfbba+mfbbc))) + mfbbb;

			real rho = one+drho;
			////////////////////////////////////////////////////////////////////////////////////
			//slow
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) / rho;
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) / rho;
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			//the force be with you
			real fx = forces[0]/(pow(two,level)); //zero;//0.0032653/(pow(two,level)); //0.000000005;//(two/1600000.0) / 120.0; //
			real fy = forces[1]/(pow(two,level)); //zero;
			real fz = forces[2]/(pow(two,level)); //zero;
			vvx += fx;
			vvy += fy;
			vvz += fz;
			////////////////////////////////////////////////////////////////////////////////////
			//real omega = omega_in;
			////////////////////////////////////////////////////////////////////////////////////
			real oMdrho = one; // comp special
			////////////////////////////////////////////////////////////////////////////////////
			real m0, m1, m2;	
			real vx2;
			real vy2;
			real vz2;
			vx2=vvx*vvx;
			vy2=vvy*vvy;
			vz2=vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real wadjust;
			real qudricLimitP = 0.01f;// * 0.0001f;
			real qudricLimitM = 0.01f;// * 0.0001f;
			real qudricLimitD = 0.01f;// * 0.001f;
			real qudricLimitOmega2 = 0.01f;
			const real gamma = 3400.f; //air
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += one* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = four * omega / (four + three * gamma * (two - omega));//one;	//set the bulk viscosity one is high / two is very low and zero is (too) high ... (also called omega 2)

			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			//no bulk viscosity
			//real OxyyPxzz  = eight*(-two+omega)*(one+two*omega)/(-eight-fourteen*omega+seven*omega*omega);//one;
			//real OxyyMxzz  = eight*(-two+omega)*(-seven+four*omega)/(fiftysix-fifty*omega+nine*omega*omega);//one;
			//real Oxyz      = twentyfour*(-two+omega)*(-two-seven*omega+three*omega*omega)/(fourtyeight+c152*omega-c130*omega*omega+twentynine*omega*omega*omega);//one;
			//with bulk viscosity
			real OxyyPxzz  = eight*(-two+omega)*(OxxPyyPzz * (-one+three*omega) - five * omega)/(eight * (five - two * omega) * omega + OxxPyyPzz * ( eight + omega * ( nine * omega - twentysix)));//one;
			real OxyyMxzz  = eight*(-two+omega)*(omega + OxxPyyPzz*(three * omega - seven))/(OxxPyyPzz * (fiftysix - fortytwo * omega + nine * omega * omega) - eight * omega);//one;
			real Oxyz      = twentyfour*(-two+omega)*(four*omega*omega + omega * OxxPyyPzz * (eighteen - thirteen * omega) + OxxPyyPzz * OxxPyyPzz * (two + omega * ( six * omega - eleven)))/
								(sixteen * omega * omega * ( omega - six ) - two * omega * OxxPyyPzz *( c216 + five  * omega * ( nine * omega - fourtysix )) + 
								OxxPyyPzz * OxxPyyPzz * ( omega * ( three * omega - ten ) * ( fiveteen * omega - twentyeight ) - fourtyeight ));//one;
			////////////////////////////////////////////////////////////
			//4.
			//////////////////////////////
			real O4        = one;
			//////////////////////////////
			//real O4        = omega;//TRT
			////////////////////////////////////////////////////////////
			//5.
			//////////////////////////////
			real O5        = one;
			////////////////////////////////////////////////////////////
			//6.
			//////////////////////////////
			real O6        = one;
			////////////////////////////////////////////////////////////


			//central moments to cumulants
			//4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho;
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho;
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho;
				  	 		
			real CUMcca = mfcca - (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			real CUMcac = mfcac - (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			real CUMacc = mfacc - (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			real CUMccc = mfccc + ((-four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));




			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;
			
			////////////////////////////////////////////////////////////////////////////
            real Dxy =-three*omega*mfbba;
            real Dxz =-three*omega*mfbab;
            real Dyz =-three*omega*mfabb;

			//3.
			// linear combinations

			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			//incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
 			
 				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
 				real dyuy = dxux + omega * c3o2 * mxxMyy;
 				real dzuz = dxux + omega * c3o2 * mxxMzz;
 
 				//relax
				//with limiter (bulk)
				wadjust    = OxxPyyPzz+(one-OxxPyyPzz)*abs(mfaaa  - mxxPyyPzz- three * (one/OxxPyyPzz - c1o2 ) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz))/
					         (abs(mfaaa  - mxxPyyPzz- three * (one/OxxPyyPzz - c1o2 ) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz))+qudricLimitOmega2);
 				mxxPyyPzz += wadjust*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 				//without limiter (no bulk)
				//mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 				mxxMyy    += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
 				mxxMzz    += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
 
 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			////no correction
 			//mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);//-magicBulk*OxxPyyPzz;
 			//mxxMyy    += -(-omega) * (-mxxMyy);
 			//mxxMzz    += -(-omega) * (-mxxMzz);
 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb     += omega * (-mfabb);
			mfbab     += omega * (-mfbab);
			mfbba     += omega * (-mfbba);

			//////////////////////////////////////////////////////////////////////////
			// linear combinations back
			mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-two*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - two* mxxMzz + mxxPyyPzz);


			//relax
			//////////////////////////////////////////////////////////////////////////
			//das ist der limiter
 			wadjust    = Oxyz+(one-Oxyz)*abs(mfbbb)/(abs(mfbbb)+qudricLimitD);
 			mfbbb     += wadjust * (-mfbbb);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimitP);
 			mxxyPyzz  += wadjust * (-mxxyPyzz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimitM);
 			mxxyMyzz  += wadjust * (-mxxyMyzz);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimitP);
 			mxxzPyyz  += wadjust * (-mxxzPyyz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimitM);
 			mxxzMyyz  += wadjust * (-mxxzMyyz);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimitP);
 			mxyyPxzz  += wadjust * (-mxyyPxzz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimitM);
 			mxyyMxzz  += wadjust * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
			//mfbbb     += OxyyMxzz * (-mfbbb);
			//mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
			//mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
			//mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
			//mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
			//mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
			//mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			//////////////////////////////////////////////////////////////////////////
			//mit limiter
			real A = ( four * omega * omega + two * omega * OxxPyyPzz * ( omega - six ) + OxxPyyPzz * OxxPyyPzz * ( omega * ( ten - three * omega ) - four )) /
				        ( ( omega - OxxPyyPzz ) * ( OxxPyyPzz * ( two + three * omega ) - eight * omega ) );
			real B = ( four * omega * OxxPyyPzz * ( nine * omega - sixteen ) - four * omega * omega - two * OxxPyyPzz * OxxPyyPzz * ( two + nine * omega * ( omega - two ))) /
				        ( three * ( omega - OxxPyyPzz ) * ( OxxPyyPzz * ( two + three * omega ) - eight * omega ) );
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
			//CUMacc += O4 * (-CUMacc); 
			//CUMcac += O4 * (-CUMcac); 
			//CUMcca += O4 * (-CUMcca); 
			//CUMbbc += O4 * (-CUMbbc); 
			//CUMbcb += O4 * (-CUMbcb); 
			//CUMcbb += O4 * (-CUMcbb); 
			//no bulk viscosity 
			//CUMacc = -O4*(one/omega-c1o2)*(dyuy+dzuz)*c2o3 *(four+two*omega-three*omega*omega)/(two-seven*omega+five*omega*omega)+(one-O4) * (CUMacc);
			//CUMcac = -O4*(one/omega-c1o2)*(dxux+dzuz)*c2o3 *(four+two*omega-three*omega*omega)/(two-seven*omega+five*omega*omega)+(one-O4) * (CUMcac);
			//CUMcca = -O4*(one/omega-c1o2)*(dyuy+dxux)*c2o3 *(four+two*omega-three*omega*omega)/(two-seven*omega+five*omega*omega)+(one-O4) * (CUMcca);
			//CUMbbc = -O4*(one/omega-c1o2)*Dxy*c1o3 *(four+twentyeight*omega-fourteen*omega*omega)/(six-twentyone*omega+fiveteen*omega*omega)+(one-O4) * (CUMbbc);
			//CUMbcb = -O4*(one/omega-c1o2)*Dxz*c1o3 *(four+twentyeight*omega-fourteen*omega*omega)/(six-twentyone*omega+fiveteen*omega*omega)+(one-O4) * (CUMbcb);
			//CUMcbb = -O4*(one/omega-c1o2)*Dyz*c1o3 *(four+twentyeight*omega-fourteen*omega*omega)/(six-twentyone*omega+fiveteen*omega*omega)+(one-O4) * (CUMcbb);
			//with bulk viscosity 
			CUMacc = -O4 * (one / omega - c1o2) * (dyuy + dzuz) * c2o3 * A + (one - O4) * (CUMacc);
			CUMcac = -O4 * (one / omega - c1o2) * (dxux + dzuz) * c2o3 * A + (one - O4) * (CUMcac);
			CUMcca = -O4 * (one / omega - c1o2) * (dyuy + dxux) * c2o3 * A + (one - O4) * (CUMcca);
			CUMbbc = -O4 * (one / omega - c1o2) * Dxy           * c1o3 * B + (one - O4) * (CUMbbc);
			CUMbcb = -O4 * (one / omega - c1o2) * Dxz           * c1o3 * B + (one - O4) * (CUMbcb);
			CUMcbb = -O4 * (one / omega - c1o2) * Dyz           * c1o3 * B + (one - O4) * (CUMcbb);


			//////////////////////////////////////////////////////////////////////////
			
					
			//5.
			CUMbcc += O5 * (-CUMbcc);
			CUMcbc += O5 * (-CUMcbc);
			CUMccb += O5 * (-CUMccb);

			//6.
			CUMccc += O6 * (-CUMccc);
			


			//back cumulants to central moments
			//4.
			mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho;
			mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho;
			mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho;
						   
			mfcca = CUMcca + (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			mfcac = CUMcac + (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			mfacc = CUMacc + (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			mfccc = CUMccc - ((-four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));

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
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + one* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - two* mfaab *  vvz         +  mfaaa                * (one- vz2)              - one* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + one* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - two* mfabb *  vvz         + mfaba * (one- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - two* mfacb *  vvz         +  mfaca                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - two* mfbab *  vvz         + mfbaa * (one- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - two* mfbbb *  vvz         + mfbba * (one- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - two* mfbcb *  vvz         + mfbca * (one- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - two* mfcab *  vvz         +  mfcaa                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - two* mfcbb *  vvz         + mfcba * (one- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - two* mfccb *  vvz         +  mfcca                  * (one- vz2)              - c1o9 * oMdrho * vz2; 
			m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcca = m0;
			mfccb = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfaca        - two* mfaba *  vvy         +  mfaaa                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - two* mfabb *  vvy         +  mfaab                  * (one- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - two* mfabc *  vvy         +  mfaac                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - two* mfbba *  vvy         + mfbaa * (one- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - two* mfbbb *  vvy         + mfbab * (one- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - two* mfbbc *  vvy         + mfbac * (one- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - two* mfcba *  vvy         +  mfcaa                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - two* mfcbb *  vvy         +  mfcab                  * (one- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - two* mfcbc *  vvy         +  mfcac                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcac = m0;
			mfcbc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcaa        - two* mfbaa *  vvx         +  mfaaa                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - two* mfbba *  vvx         +  mfaba                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - two* mfbca *  vvx         +  mfaca                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - two* mfbab *  vvx         +  mfaab                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - two* mfbbb *  vvx         +  mfabb                  * (one- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - two* mfbcb *  vvx         +  mfacb                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - two* mfbac *  vvx         +  mfaac                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - two* mfbbc *  vvx         +  mfabc                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - two* mfbcc *  vvx         +  mfacc                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dirE   ])[k   ] = mfabb;//(D.f[ dirE   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dirE   ])[k   ]                                                                     
			(D.f[ dirW   ])[kw  ] = mfcbb;//(D.f[ dirW   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dirW   ])[kw  ]                                                                   
			(D.f[ dirN   ])[k   ] = mfbab;//(D.f[ dirN   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ dirN   ])[k   ]
			(D.f[ dirS   ])[ks  ] = mfbcb;//(D.f[ dirS   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ dirS   ])[ks  ]
			(D.f[ dirT   ])[k   ] = mfbba;//(D.f[ dirT   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ dirT   ])[k   ]
			(D.f[ dirB   ])[kb  ] = mfbbc;//(D.f[ dirB   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ dirB   ])[kb  ]
			(D.f[ dirNE  ])[k   ] = mfaab;//(D.f[ dirNE  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dirNE  ])[k   ]
			(D.f[ dirSW  ])[ksw ] = mfccb;//(D.f[ dirSW  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dirSW  ])[ksw ]
			(D.f[ dirSE  ])[ks  ] = mfacb;//(D.f[ dirSE  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dirSE  ])[ks  ]
			(D.f[ dirNW  ])[kw  ] = mfcab;//(D.f[ dirNW  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dirNW  ])[kw  ]
			(D.f[ dirTE  ])[k   ] = mfaba;//(D.f[ dirTE  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dirTE  ])[k   ]
			(D.f[ dirBW  ])[kbw ] = mfcbc;//(D.f[ dirBW  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dirBW  ])[kbw ]
			(D.f[ dirBE  ])[kb  ] = mfabc;//(D.f[ dirBE  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dirBE  ])[kb  ]
			(D.f[ dirTW  ])[kw  ] = mfcba;//(D.f[ dirTW  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dirTW  ])[kw  ]
			(D.f[ dirTN  ])[k   ] = mfbaa;//(D.f[ dirTN  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ dirTN  ])[k   ]
			(D.f[ dirBS  ])[kbs ] = mfbcc;//(D.f[ dirBS  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ dirBS  ])[kbs ]
			(D.f[ dirBN  ])[kb  ] = mfbac;//(D.f[ dirBN  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ dirBN  ])[kb  ]
			(D.f[ dirTS  ])[ks  ] = mfbca;//(D.f[ dirTS  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ dirTS  ])[ks  ]
			(D.f[ dirZERO])[k   ] = mfbbb;//(D.f[ dirZERO])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ dirZERO])[k   ]
			(D.f[ dirTNE ])[k   ] = mfaaa;//(D.f[ dirTNE ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dirTNE ])[k   ]
			(D.f[ dirTSE ])[ks  ] = mfaca;//(D.f[ dirTSE ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dirTSE ])[ks  ]
			(D.f[ dirBNE ])[kb  ] = mfaac;//(D.f[ dirBNE ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dirBNE ])[kb  ]
			(D.f[ dirBSE ])[kbs ] = mfacc;//(D.f[ dirBSE ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dirBSE ])[kbs ]
			(D.f[ dirTNW ])[kw  ] = mfcaa;//(D.f[ dirTNW ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dirTNW ])[kw  ]
			(D.f[ dirTSW ])[ksw ] = mfcca;//(D.f[ dirTSW ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dirTSW ])[ksw ]
			(D.f[ dirBNW ])[kbw ] = mfcac;//(D.f[ dirBNW ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dirBNW ])[kbw ]
			(D.f[ dirBSW ])[kbsw] = mfccc;//(D.f[ dirBSW ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dirBSW ])[kbsw]
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Kum_AA2016_Comp_SP_27( real omega,
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

	if(k<size_Mat)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			//unsigned int kzero= k;
			//unsigned int ke   = k;
			unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			unsigned int kbs  = neighborZ[ks];
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

			//unsigned int kzero= k;
			//unsigned int ke   = k;
			//unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			//unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			//unsigned int kb   = neighborZ[k];
			//unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			//unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			//unsigned int kbs  = neighborZ[ks];
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
			//unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[dirE   ])[k  ];//[ke   ];// +  c2over27 ;(D.f[dirE   ])[k  ];//ke
			real mfabb = (D.f[dirW   ])[kw ];//[kw   ];// +  c2over27 ;(D.f[dirW   ])[kw ];
			real mfbcb = (D.f[dirN   ])[k  ];//[kn   ];// +  c2over27 ;(D.f[dirN   ])[k  ];//kn
			real mfbab = (D.f[dirS   ])[ks ];//[ks   ];// +  c2over27 ;(D.f[dirS   ])[ks ];
			real mfbbc = (D.f[dirT   ])[k  ];//[kt   ];// +  c2over27 ;(D.f[dirT   ])[k  ];//kt
			real mfbba = (D.f[dirB   ])[kb ];//[kb   ];// +  c2over27 ;(D.f[dirB   ])[kb ];
			real mfccb = (D.f[dirNE  ])[k  ];//[kne  ];// +  c1over54 ;(D.f[dirNE  ])[k  ];//kne
			real mfaab = (D.f[dirSW  ])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dirSW  ])[ksw];
			real mfcab = (D.f[dirSE  ])[ks ];//[kse  ];// +  c1over54 ;(D.f[dirSE  ])[ks ];//kse
			real mfacb = (D.f[dirNW  ])[kw ];//[knw  ];// +  c1over54 ;(D.f[dirNW  ])[kw ];//knw
			real mfcbc = (D.f[dirTE  ])[k  ];//[kte  ];// +  c1over54 ;(D.f[dirTE  ])[k  ];//kte
			real mfaba = (D.f[dirBW  ])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dirBW  ])[kbw];
			real mfcba = (D.f[dirBE  ])[kb ];//[kbe  ];// +  c1over54 ;(D.f[dirBE  ])[kb ];//kbe
			real mfabc = (D.f[dirTW  ])[kw ];//[ktw  ];// +  c1over54 ;(D.f[dirTW  ])[kw ];//ktw
			real mfbcc = (D.f[dirTN  ])[k  ];//[ktn  ];// +  c1over54 ;(D.f[dirTN  ])[k  ];//ktn
			real mfbaa = (D.f[dirBS  ])[kbs];//[kbs  ];// +  c1over54 ;(D.f[dirBS  ])[kbs];
			real mfbca = (D.f[dirBN  ])[kb ];//[kbn  ];// +  c1over54 ;(D.f[dirBN  ])[kb ];//kbn
			real mfbac = (D.f[dirTS  ])[ks ];//[kts  ];// +  c1over54 ;(D.f[dirTS  ])[ks ];//kts
			real mfbbb = (D.f[dirZERO])[k  ];//[kzero];// +  c8over27 ;(D.f[dirZERO])[k  ];//kzero
			real mfccc = (D.f[dirTNE ])[k  ];//[ktne ];// +  c1over216;(D.f[dirTNE ])[k  ];//ktne
			real mfaac = (D.f[dirTSW ])[ksw];//[ktsw ];// +  c1over216;(D.f[dirTSW ])[ksw];//ktsw
			real mfcac = (D.f[dirTSE ])[ks ];//[ktse ];// +  c1over216;(D.f[dirTSE ])[ks ];//ktse
			real mfacc = (D.f[dirTNW ])[kw ];//[ktnw ];// +  c1over216;(D.f[dirTNW ])[kw ];//ktnw
			real mfcca = (D.f[dirBNE ])[kb ];//[kbne ];// +  c1over216;(D.f[dirBNE ])[kb ];//kbne
			real mfaaa = (D.f[dirBSW ])[kbsw];//[kbsw ];// +  c1over216;(D.f[dirBSW ])[kbsw];
			real mfcaa = (D.f[dirBSE ])[kbs];//[kbse ];// +  c1over216;(D.f[dirBSE ])[kbs];//kbse
			real mfaca = (D.f[dirBNW ])[kbw];//[kbnw ];// +  c1over216;(D.f[dirBNW ])[kbw];//kbnw
			////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
							(((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
							((mfabb+mfcbb) + (mfbab+mfbcb) + (mfbba+mfbbc))) + mfbbb;

			real rho = one+drho;
			////////////////////////////////////////////////////////////////////////////////////
			//slow
			//real oMdrho = one - ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
			//					   (((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
			//						((mfabb+mfcbb) + (mfbab+mfbcb)  +  (mfbba+mfbbc)));//fehlt mfbbb
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) / rho;
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) / rho;
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			//the force be with you
			real fx = forces[0]/(pow(two,level)); //zero;//0.0032653/(pow(two,level)); //0.000000005;//(two/1600000.0) / 120.0; //
			real fy = forces[1]/(pow(two,level)); //zero;
			real fz = forces[2]/(pow(two,level)); //zero;
			vvx += fx*c1o2;
			vvy += fy*c1o2;
			vvz += fz*c1o2;
			////////////////////////////////////////////////////////////////////////////////////
			//real omega = omega_in;
			////////////////////////////////////////////////////////////////////////////////////
			//fast
			real oMdrho = one; // comp special
			//real oMdrho = one - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
			//					   mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
			//					   mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb);//fehlt mfbbb nicht mehr
			//real vvx    =mfccc-mfaaa + mfcac-mfaca + mfcaa-mfacc + mfcca-mfaac + 
			//				mfcba-mfabc + mfcbc-mfaba + mfcab-mfacb + mfccb-mfaab +
			//				mfcbb-mfabb;
			//real vvy    =mfccc-mfaaa + mfaca-mfcac + mfacc-mfcaa + mfcca-mfaac + 
			//				mfbca-mfbac + mfbcc-mfbaa + mfacb-mfcab + mfccb-mfaab +
			//				mfbcb-mfbab;
			//real vvz    =mfccc-mfaaa + mfcac-mfaca + mfacc-mfcaa + mfaac-mfcca + 
			//				mfbac-mfbca + mfbcc-mfbaa + mfabc-mfcba + mfcbc-mfaba +
			//				mfbbc-mfbba;
			////////////////////////////////////////////////////////////////////////////////////
			// oMdrho assembler style -------> faaaaaastaaaa
			// or much sloooowaaaa ... it depändssssss on sadaku
			real m0, m1, m2;	
			//real oMdrho;
			//{
			//	oMdrho=mfccc+mfaaa;
			//	m0=mfaca+mfcac;
			//	m1=mfacc+mfcaa;
			//	m2=mfaac+mfcca;
			//	oMdrho+=m0;
			//	m1+=m2;
			//	oMdrho+=m1;
			//	m0=mfbac+mfbca;
			//	m1=mfbaa+mfbcc;
			//	m0+=m1;
			//	m1=mfabc+mfcba;
			//	m2=mfaba+mfcbc;
			//	m1+=m2;
			//	m0+=m1;
			//	m1=mfacb+mfcab;
			//	m2=mfaab+mfccb;
			//	m1+=m2;
			//	m0+=m1;
			//	oMdrho+=m0;
			//	m0=mfabb+mfcbb;
			//	m1=mfbab+mfbcb;
			//	m2=mfbba+mfbbc;
			//	m0+=m1+m2;
			//	m0+=mfbbb; //hat gefehlt
			//	oMdrho = one - (oMdrho + m0);
			//}
			//real vvx;
			real vx2;
			//{
			//	vvx = mfccc-mfaaa;
			//	m0  = mfcac-mfaca;
			//	m1  = mfcaa-mfacc;
			//	m2  = mfcca-mfaac;
			//	vvx+= m0;
			//	m1 += m2;
			//	vvx+= m1;
			//	vx2 = mfcba-mfabc;
			//	m0  = mfcbc-mfaba;
			//	m1  = mfcab-mfacb;
			//	m2  = mfccb-mfaab;
			//	vx2+= m0;
			//	m1 += m2;
			//	vx2+= m1;
			//	vvx+= vx2;
			//	vx2 = mfcbb-mfabb;
			//	vvx+= vx2;
			//}
			//real vvy;
			real vy2;
			//{
			//	vvy = mfccc-mfaaa;
			//	m0  = mfaca-mfcac;
			//	m1  = mfacc-mfcaa;
			//	m2  = mfcca-mfaac;
			//	vvy+= m0;
			//	m1 += m2;
			//	vvy+= m1;
			//	vy2 = mfbca-mfbac;
			//	m0  = mfbcc-mfbaa;
			//	m1  = mfacb-mfcab;
			//	m2  = mfccb-mfaab;
			//	vy2+= m0;
			//	m1 += m2;
			//	vy2+= m1;
			//	vvy+= vy2;
			//	vy2 = mfbcb-mfbab;
			//	vvy+= vy2;
			//}
			//real vvz;
			real vz2;
			//{
			//	vvz = mfccc-mfaaa;
			//	m0  = mfcac-mfaca;
			//	m1  = mfacc-mfcaa;
			//	m2  = mfaac-mfcca;
			//	vvz+= m0;
			//	m1 += m2;
			//	vvz+= m1;
			//	vz2 = mfbac-mfbca;
			//	m0  = mfbcc-mfbaa;
			//	m1  = mfabc-mfcba;
			//	m2  = mfcbc-mfaba;
			//	vz2+= m0;
			//	m1 += m2;
			//	vz2+= m1;
			//	vvz+= vz2;
			//	vz2 = mfbbc-mfbba;
			//	vvz+= vz2;
			//}
			vx2=vvx*vvx;
			vy2=vvy*vvy;
			vz2=vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real wadjust;
			real qudricLimitP = 0.01f;// * 0.0001f; // 1000000.0f; // 1000000.0f; //
			real qudricLimitM = 0.01f;// * 0.0001f; // 1000000.0f; // 1000000.0f; //
			real qudricLimitD = 0.01f;// * 0.001f;  // 1000000.0f; // 1000000.0f; //
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += one* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = one;	//set the bulk viscosity one is high / two is very low and zero is (too) high ... (also called omega 2)

			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			real OxyyPxzz  = eight*(-two+omega)*(one+two*omega)/(-eight-fourteen*omega+seven*omega*omega);//one;
			real OxyyMxzz  = eight*(-two+omega)*(-seven+four*omega)/(fiftysix-fifty*omega+nine*omega*omega);//one;
			real Oxyz      = twentyfour*(-two+omega)*(-two-seven*omega+three*omega*omega)/(fourtyeight+c152*omega-c130*omega*omega+twentynine*omega*omega*omega);//one;
			////////////////////////////////////////////////////////////
			//4.
			//////////////////////////////
			real O4        = one;
			//////////////////////////////
			//real O4        = omega;//TRT
			////////////////////////////////////////////////////////////
			//5.
			//////////////////////////////
			real O5        = one;
			////////////////////////////////////////////////////////////
			//6.
			//////////////////////////////
			real O6        = one;
			////////////////////////////////////////////////////////////


			//central moments to cumulants
			//4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho;
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho;
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho;
				  	 		
			real CUMcca = mfcca - (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			real CUMcac = mfcac - (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			real CUMacc = mfacc - (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			real CUMccc = mfccc + ((-four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));

			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;
			
			////////////////////////////////////////////////////////////////////////////
            real Dxy =-three*omega*mfbba;
            real Dxz =-three*omega*mfbab;
            real Dyz =-three*omega*mfabb;

			//3.
			// linear combinations

			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			//incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
 			
 				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
 				real dyuy = dxux + omega * c3o2 * mxxMyy;
 				real dzuz = dxux + omega * c3o2 * mxxMzz;
 
 				//relax
				mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 				mxxMyy    += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
 				mxxMzz    += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
 
 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			////no correction
 			//mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);//-magicBulk*OxxPyyPzz;
 			//mxxMyy    += -(-omega) * (-mxxMyy);
 			//mxxMzz    += -(-omega) * (-mxxMzz);
 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb     += omega * (-mfabb);
			mfbab     += omega * (-mfbab);
			mfbba     += omega * (-mfbba);
			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-two*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - two* mxxMzz + mxxPyyPzz);


			//relax
			//////////////////////////////////////////////////////////////////////////
			//das ist der limiter
 			wadjust    = Oxyz+(one-Oxyz)*abs(mfbbb)/(abs(mfbbb)+qudricLimitD);
 			mfbbb     += wadjust * (-mfbbb);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimitP);
 			mxxyPyzz  += wadjust * (-mxxyPyzz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimitM);
 			mxxyMyzz  += wadjust * (-mxxyMyzz);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimitP);
 			mxxzPyyz  += wadjust * (-mxxzPyyz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimitM);
 			mxxzMyyz  += wadjust * (-mxxzMyyz);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimitP);
 			mxyyPxzz  += wadjust * (-mxyyPxzz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimitM);
 			mxyyMxzz  += wadjust * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
			//mfbbb     += OxyyMxzz * (-mfbbb);
			//mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
			//mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
			//mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
			//mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
			//mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
			//mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			//////////////////////////////////////////////////////////////////////////
			//mit limiter
 		//	wadjust    = O4+(one-O4)*abs(CUMacc)/(abs(CUMacc)+qudricLimit);
			//CUMacc    += wadjust * (-CUMacc);
 		//	wadjust    = O4+(one-O4)*abs(CUMcac)/(abs(CUMcac)+qudricLimit);
			//CUMcac    += wadjust * (-CUMcac); 
 		//	wadjust    = O4+(one-O4)*abs(CUMcca)/(abs(CUMcca)+qudricLimit);
			//CUMcca    += wadjust * (-CUMcca); 

 		//	wadjust    = O4+(one-O4)*abs(CUMbbc)/(abs(CUMbbc)+qudricLimit);
			//CUMbbc    += wadjust * (-CUMbbc); 
 		//	wadjust    = O4+(one-O4)*abs(CUMbcb)/(abs(CUMbcb)+qudricLimit);
			//CUMbcb    += wadjust * (-CUMbcb); 
 		//	wadjust    = O4+(one-O4)*abs(CUMcbb)/(abs(CUMcbb)+qudricLimit);
			//CUMcbb    += wadjust * (-CUMcbb); 
			//////////////////////////////////////////////////////////////////////////
			real A = (four + two*omega - three*omega*omega) / (two - seven*omega + five*omega*omega);
			real B = (four + twentyeight*omega - fourteen*omega*omega) / (six - twentyone*omega + fiveteen*omega*omega);
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
			//CUMacc += O4 * (-CUMacc); 
			//CUMcac += O4 * (-CUMcac); 
			//CUMcca += O4 * (-CUMcca); 
			//CUMbbc += O4 * (-CUMbbc); 
			//CUMbcb += O4 * (-CUMbcb); 
			//CUMcbb += O4 * (-CUMcbb); 
			CUMacc = -O4*(one / omega - c1o2) * (dyuy + dzuz) * c2o3 * A + (one - O4) * (CUMacc);
			CUMcac = -O4*(one / omega - c1o2) * (dxux + dzuz) * c2o3 * A + (one - O4) * (CUMcac);
			CUMcca = -O4*(one / omega - c1o2) * (dyuy + dxux) * c2o3 * A + (one - O4) * (CUMcca);
			CUMbbc = -O4*(one / omega - c1o2) * Dxy           * c1o3 * B + (one - O4) * (CUMbbc);
			CUMbcb = -O4*(one / omega - c1o2) * Dxz           * c1o3 * B + (one - O4) * (CUMbcb);
			CUMcbb = -O4*(one / omega - c1o2) * Dyz           * c1o3 * B + (one - O4) * (CUMcbb);
			//////////////////////////////////////////////////////////////////////////
			
					
			//5.
			CUMbcc += O5 * (-CUMbcc);
			CUMcbc += O5 * (-CUMcbc);
			CUMccb += O5 * (-CUMccb);

			//6.
			CUMccc += O6 * (-CUMccc);
			


			//back cumulants to central moments
			//4.
			mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho;
			mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho;
			mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho;
						   
			mfcca = CUMcca + (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			mfcac = CUMcac + (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			mfacc = CUMacc + (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			mfccc = CUMccc - ((-four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));
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
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + one* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - two* mfaab *  vvz         +  mfaaa                * (one- vz2)              - one* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + one* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - two* mfabb *  vvz         + mfaba * (one- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - two* mfacb *  vvz         +  mfaca                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - two* mfbab *  vvz         + mfbaa * (one- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - two* mfbbb *  vvz         + mfbba * (one- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - two* mfbcb *  vvz         + mfbca * (one- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - two* mfcab *  vvz         +  mfcaa                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - two* mfcbb *  vvz         + mfcba * (one- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - two* mfccb *  vvz         +  mfcca                  * (one- vz2)              - c1o9 * oMdrho * vz2; 
			m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcca = m0;
			mfccb = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfaca        - two* mfaba *  vvy         +  mfaaa                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - two* mfabb *  vvy         +  mfaab                  * (one- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - two* mfabc *  vvy         +  mfaac                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - two* mfbba *  vvy         + mfbaa * (one- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - two* mfbbb *  vvy         + mfbab * (one- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - two* mfbbc *  vvy         + mfbac * (one- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - two* mfcba *  vvy         +  mfcaa                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - two* mfcbb *  vvy         +  mfcab                  * (one- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - two* mfcbc *  vvy         +  mfcac                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcac = m0;
			mfcbc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcaa        - two* mfbaa *  vvx         +  mfaaa                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - two* mfbba *  vvx         +  mfaba                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - two* mfbca *  vvx         +  mfaca                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - two* mfbab *  vvx         +  mfaab                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - two* mfbbb *  vvx         +  mfabb                  * (one- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - two* mfbcb *  vvx         +  mfacb                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - two* mfbac *  vvx         +  mfaac                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - two* mfbbc *  vvx         +  mfabc                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - two* mfbcc *  vvx         +  mfacc                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dirE   ])[k   ] = mfabb;//(D.f[ dirE   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dirE   ])[k   ]                                                                     
			(D.f[ dirW   ])[kw  ] = mfcbb;//(D.f[ dirW   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dirW   ])[kw  ]                                                                   
			(D.f[ dirN   ])[k   ] = mfbab;//(D.f[ dirN   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ dirN   ])[k   ]
			(D.f[ dirS   ])[ks  ] = mfbcb;//(D.f[ dirS   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ dirS   ])[ks  ]
			(D.f[ dirT   ])[k   ] = mfbba;//(D.f[ dirT   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ dirT   ])[k   ]
			(D.f[ dirB   ])[kb  ] = mfbbc;//(D.f[ dirB   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ dirB   ])[kb  ]
			(D.f[ dirNE  ])[k   ] = mfaab;//(D.f[ dirNE  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dirNE  ])[k   ]
			(D.f[ dirSW  ])[ksw ] = mfccb;//(D.f[ dirSW  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dirSW  ])[ksw ]
			(D.f[ dirSE  ])[ks  ] = mfacb;//(D.f[ dirSE  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dirSE  ])[ks  ]
			(D.f[ dirNW  ])[kw  ] = mfcab;//(D.f[ dirNW  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dirNW  ])[kw  ]
			(D.f[ dirTE  ])[k   ] = mfaba;//(D.f[ dirTE  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dirTE  ])[k   ]
			(D.f[ dirBW  ])[kbw ] = mfcbc;//(D.f[ dirBW  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dirBW  ])[kbw ]
			(D.f[ dirBE  ])[kb  ] = mfabc;//(D.f[ dirBE  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dirBE  ])[kb  ]
			(D.f[ dirTW  ])[kw  ] = mfcba;//(D.f[ dirTW  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dirTW  ])[kw  ]
			(D.f[ dirTN  ])[k   ] = mfbaa;//(D.f[ dirTN  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ dirTN  ])[k   ]
			(D.f[ dirBS  ])[kbs ] = mfbcc;//(D.f[ dirBS  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ dirBS  ])[kbs ]
			(D.f[ dirBN  ])[kb  ] = mfbac;//(D.f[ dirBN  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ dirBN  ])[kb  ]
			(D.f[ dirTS  ])[ks  ] = mfbca;//(D.f[ dirTS  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ dirTS  ])[ks  ]
			(D.f[ dirZERO])[k   ] = mfbbb;//(D.f[ dirZERO])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ dirZERO])[k   ]
			(D.f[ dirTNE ])[k   ] = mfaaa;//(D.f[ dirTNE ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dirTNE ])[k   ]
			(D.f[ dirTSE ])[ks  ] = mfaca;//(D.f[ dirTSE ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dirTSE ])[ks  ]
			(D.f[ dirBNE ])[kb  ] = mfaac;//(D.f[ dirBNE ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dirBNE ])[kb  ]
			(D.f[ dirBSE ])[kbs ] = mfacc;//(D.f[ dirBSE ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dirBSE ])[kbs ]
			(D.f[ dirTNW ])[kw  ] = mfcaa;//(D.f[ dirTNW ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dirTNW ])[kw  ]
			(D.f[ dirTSW ])[ksw ] = mfcca;//(D.f[ dirTSW ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dirTSW ])[ksw ]
			(D.f[ dirBNW ])[kbw ] = mfcac;//(D.f[ dirBNW ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dirBNW ])[kbw ]
			(D.f[ dirBSW ])[kbsw] = mfccc;//(D.f[ dirBSW ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dirBSW ])[kbsw]
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Kum_New_Comp_Bulk_SP_27(real omega,
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

	if(k<size_Mat)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			//unsigned int kzero= k;
			//unsigned int ke   = k;
			unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			unsigned int kbs  = neighborZ[ks];
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

			//unsigned int kzero= k;
			//unsigned int ke   = k;
			//unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			//unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			//unsigned int kb   = neighborZ[k];
			//unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			//unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			//unsigned int kbs  = neighborZ[ks];
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
			//unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[dirE   ])[k  ];//[ke   ];// +  c2over27 ;(D.f[dirE   ])[k  ];//ke
			real mfabb = (D.f[dirW   ])[kw ];//[kw   ];// +  c2over27 ;(D.f[dirW   ])[kw ];
			real mfbcb = (D.f[dirN   ])[k  ];//[kn   ];// +  c2over27 ;(D.f[dirN   ])[k  ];//kn
			real mfbab = (D.f[dirS   ])[ks ];//[ks   ];// +  c2over27 ;(D.f[dirS   ])[ks ];
			real mfbbc = (D.f[dirT   ])[k  ];//[kt   ];// +  c2over27 ;(D.f[dirT   ])[k  ];//kt
			real mfbba = (D.f[dirB   ])[kb ];//[kb   ];// +  c2over27 ;(D.f[dirB   ])[kb ];
			real mfccb = (D.f[dirNE  ])[k  ];//[kne  ];// +  c1over54 ;(D.f[dirNE  ])[k  ];//kne
			real mfaab = (D.f[dirSW  ])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dirSW  ])[ksw];
			real mfcab = (D.f[dirSE  ])[ks ];//[kse  ];// +  c1over54 ;(D.f[dirSE  ])[ks ];//kse
			real mfacb = (D.f[dirNW  ])[kw ];//[knw  ];// +  c1over54 ;(D.f[dirNW  ])[kw ];//knw
			real mfcbc = (D.f[dirTE  ])[k  ];//[kte  ];// +  c1over54 ;(D.f[dirTE  ])[k  ];//kte
			real mfaba = (D.f[dirBW  ])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dirBW  ])[kbw];
			real mfcba = (D.f[dirBE  ])[kb ];//[kbe  ];// +  c1over54 ;(D.f[dirBE  ])[kb ];//kbe
			real mfabc = (D.f[dirTW  ])[kw ];//[ktw  ];// +  c1over54 ;(D.f[dirTW  ])[kw ];//ktw
			real mfbcc = (D.f[dirTN  ])[k  ];//[ktn  ];// +  c1over54 ;(D.f[dirTN  ])[k  ];//ktn
			real mfbaa = (D.f[dirBS  ])[kbs];//[kbs  ];// +  c1over54 ;(D.f[dirBS  ])[kbs];
			real mfbca = (D.f[dirBN  ])[kb ];//[kbn  ];// +  c1over54 ;(D.f[dirBN  ])[kb ];//kbn
			real mfbac = (D.f[dirTS  ])[ks ];//[kts  ];// +  c1over54 ;(D.f[dirTS  ])[ks ];//kts
			real mfbbb = (D.f[dirZERO])[k  ];//[kzero];// +  c8over27 ;(D.f[dirZERO])[k  ];//kzero
			real mfccc = (D.f[dirTNE ])[k  ];//[ktne ];// +  c1over216;(D.f[dirTNE ])[k  ];//ktne
			real mfaac = (D.f[dirTSW ])[ksw];//[ktsw ];// +  c1over216;(D.f[dirTSW ])[ksw];//ktsw
			real mfcac = (D.f[dirTSE ])[ks ];//[ktse ];// +  c1over216;(D.f[dirTSE ])[ks ];//ktse
			real mfacc = (D.f[dirTNW ])[kw ];//[ktnw ];// +  c1over216;(D.f[dirTNW ])[kw ];//ktnw
			real mfcca = (D.f[dirBNE ])[kb ];//[kbne ];// +  c1over216;(D.f[dirBNE ])[kb ];//kbne
			real mfaaa = (D.f[dirBSW ])[kbsw];//[kbsw ];// +  c1over216;(D.f[dirBSW ])[kbsw];
			real mfcaa = (D.f[dirBSE ])[kbs];//[kbse ];// +  c1over216;(D.f[dirBSE ])[kbs];//kbse
			real mfaca = (D.f[dirBNW ])[kbw];//[kbnw ];// +  c1over216;(D.f[dirBNW ])[kbw];//kbnw
			////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
							(((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
							((mfabb+mfcbb) + (mfbab+mfbcb) + (mfbba+mfbbc))) + mfbbb;

			real rho = one+drho;
			////////////////////////////////////////////////////////////////////////////////////
			//slow
			//real oMdrho = one - ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
			//					   (((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
			//						((mfabb+mfcbb) + (mfbab+mfbcb)  +  (mfbba+mfbbc)));//fehlt mfbbb
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) / rho;
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) / rho;
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			//the force be with you
			real fx = forces[0]/(pow(two,level)); //zero;//0.0032653/(pow(two,level)); //0.000000005;//(two/1600000.0) / 120.0; //
			real fy = forces[1]/(pow(two,level)); //zero;
			real fz = forces[2]/(pow(two,level)); //zero;
			vvx += fx;
			vvy += fy;
			vvz += fz;
			////////////////////////////////////////////////////////////////////////////////////
			//real omega = omega_in;
			////////////////////////////////////////////////////////////////////////////////////
			//fast
			real oMdrho = one; // comp special
			//real oMdrho = one - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
			//					   mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
			//					   mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb);//fehlt mfbbb nicht mehr
			//real vvx    =mfccc-mfaaa + mfcac-mfaca + mfcaa-mfacc + mfcca-mfaac + 
			//				mfcba-mfabc + mfcbc-mfaba + mfcab-mfacb + mfccb-mfaab +
			//				mfcbb-mfabb;
			//real vvy    =mfccc-mfaaa + mfaca-mfcac + mfacc-mfcaa + mfcca-mfaac + 
			//				mfbca-mfbac + mfbcc-mfbaa + mfacb-mfcab + mfccb-mfaab +
			//				mfbcb-mfbab;
			//real vvz    =mfccc-mfaaa + mfcac-mfaca + mfacc-mfcaa + mfaac-mfcca + 
			//				mfbac-mfbca + mfbcc-mfbaa + mfabc-mfcba + mfcbc-mfaba +
			//				mfbbc-mfbba;
			////////////////////////////////////////////////////////////////////////////////////
			// oMdrho assembler style -------> faaaaaastaaaa
			// or much sloooowaaaa ... it depändssssss on sadaku
			real m0, m1, m2;	
			real vx2;
			real vy2;
			real vz2;
			vx2=vvx*vvx;
			vy2=vvy*vvy;
			vz2=vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real wadjust;
			real qudricLimitP = 0.01f;// * 0.0001f;
			real qudricLimitM = 0.01f;// * 0.0001f;
			real qudricLimitD = 0.01f;// * 0.001f;
			real qudricLimitOmega2 = 0.01f;
			const real gamma = 3400.f; //air
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += one* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = four * omega / (four + three * gamma * (two - omega));//one;	//set the bulk viscosity one is high / two is very low and zero is (too) high ... (also called omega 2)

			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			real OxyyPxzz  = one;
			real OxyyMxzz  = one;
			real Oxyz      = one;
			////////////////////////////////////////////////////////////
			//4.
			//////////////////////////////
			real O4        = one;
			//////////////////////////////
			//real O4        = omega;//TRT
			////////////////////////////////////////////////////////////
			//5.
			//////////////////////////////
			real O5        = one;
			////////////////////////////////////////////////////////////
			//6.
			//////////////////////////////
			real O6        = one;
			////////////////////////////////////////////////////////////


			//central moments to cumulants
			//4.
			//real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab) / rho;  //bis 15.05.2015 verwendet
			//real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb) / rho;  //bis 15.05.2015 verwendet
			//real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb) / rho;  //bis 15.05.2015 verwendet
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho;	//ab 15.05.2015 verwendet
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho; //ab 15.05.2015 verwendet
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho; //ab 15.05.2015 verwendet
				  	 		
			real CUMcca = mfcca - (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			real CUMcac = mfcac - (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			real CUMacc = mfacc - (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			real CUMccc = mfccc + ((-four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));

			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;
			
 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			//incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
 			{
 				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
 				real dyuy = dxux + omega * c3o2 * mxxMyy;
 				real dzuz = dxux + omega * c3o2 * mxxMzz;
 
 				//relax
				//with limiter (bulk)
				wadjust    = OxxPyyPzz+(one-OxxPyyPzz)*abs(mfaaa  - mxxPyyPzz- three * (one/OxxPyyPzz - c1o2 ) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz))/
					         (abs(mfaaa  - mxxPyyPzz- three * (one/OxxPyyPzz - c1o2 ) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz))+qudricLimitOmega2);
 				mxxPyyPzz += wadjust*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 				//without limiter (no bulk)
				//mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 				mxxMyy    += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
 				mxxMzz    += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
  			}
 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			////no correction
 			//mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);//-magicBulk*OxxPyyPzz;
 			//mxxMyy    += -(-omega) * (-mxxMyy);
 			//mxxMzz    += -(-omega) * (-mxxMzz);
 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb     += omega * (-mfabb);
			mfbab     += omega * (-mfbab);
			mfbba     += omega * (-mfbba);
			//////////////////////////////////////////////////////////////////////////
			// linear combinations back
			mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-two*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - two* mxxMzz + mxxPyyPzz);

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
			//das ist der limiter
 			wadjust    = Oxyz+(one-Oxyz)*abs(mfbbb)/(abs(mfbbb)+qudricLimitD);
 			mfbbb     += wadjust * (-mfbbb);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimitP);
 			mxxyPyzz  += wadjust * (-mxxyPyzz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimitM);
 			mxxyMyzz  += wadjust * (-mxxyMyzz);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimitP);
 			mxxzPyyz  += wadjust * (-mxxzPyyz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimitM);
 			mxxzMyyz  += wadjust * (-mxxzMyyz);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimitP);
 			mxyyPxzz  += wadjust * (-mxyyPxzz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimitM);
 			mxyyMxzz  += wadjust * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
			//mfbbb     += OxyyMxzz * (-mfbbb);
			//mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
			//mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
			//mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
			//mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
			//mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
			//mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			//////////////////////////////////////////////////////////////////////////
			//mit limiter
 		//	wadjust    = O4+(one-O4)*abs(CUMacc)/(abs(CUMacc)+qudricLimit);
			//CUMacc    += wadjust * (-CUMacc);
 		//	wadjust    = O4+(one-O4)*abs(CUMcac)/(abs(CUMcac)+qudricLimit);
			//CUMcac    += wadjust * (-CUMcac); 
 		//	wadjust    = O4+(one-O4)*abs(CUMcca)/(abs(CUMcca)+qudricLimit);
			//CUMcca    += wadjust * (-CUMcca); 

 		//	wadjust    = O4+(one-O4)*abs(CUMbbc)/(abs(CUMbbc)+qudricLimit);
			//CUMbbc    += wadjust * (-CUMbbc); 
 		//	wadjust    = O4+(one-O4)*abs(CUMbcb)/(abs(CUMbcb)+qudricLimit);
			//CUMbcb    += wadjust * (-CUMbcb); 
 		//	wadjust    = O4+(one-O4)*abs(CUMcbb)/(abs(CUMcbb)+qudricLimit);
			//CUMcbb    += wadjust * (-CUMcbb); 
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
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
			mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho;
			mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho;
			mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho;
						   
			mfcca = CUMcca + (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			mfcac = CUMcac + (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			mfacc = CUMacc + (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			mfccc = CUMccc - ((-four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));

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
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + one* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - two* mfaab *  vvz         +  mfaaa                * (one- vz2)              - one* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + one* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - two* mfabb *  vvz         + mfaba * (one- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - two* mfacb *  vvz         +  mfaca                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - two* mfbab *  vvz         + mfbaa * (one- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - two* mfbbb *  vvz         + mfbba * (one- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - two* mfbcb *  vvz         + mfbca * (one- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - two* mfcab *  vvz         +  mfcaa                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - two* mfcbb *  vvz         + mfcba * (one- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - two* mfccb *  vvz         +  mfcca                  * (one- vz2)              - c1o9 * oMdrho * vz2; 
			m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcca = m0;
			mfccb = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfaca        - two* mfaba *  vvy         +  mfaaa                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - two* mfabb *  vvy         +  mfaab                  * (one- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - two* mfabc *  vvy         +  mfaac                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - two* mfbba *  vvy         + mfbaa * (one- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - two* mfbbb *  vvy         + mfbab * (one- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - two* mfbbc *  vvy         + mfbac * (one- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - two* mfcba *  vvy         +  mfcaa                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - two* mfcbb *  vvy         +  mfcab                  * (one- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - two* mfcbc *  vvy         +  mfcac                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcac = m0;
			mfcbc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcaa        - two* mfbaa *  vvx         +  mfaaa                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - two* mfbba *  vvx         +  mfaba                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - two* mfbca *  vvx         +  mfaca                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - two* mfbab *  vvx         +  mfaab                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - two* mfbbb *  vvx         +  mfabb                  * (one- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - two* mfbcb *  vvx         +  mfacb                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - two* mfbac *  vvx         +  mfaac                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - two* mfbbc *  vvx         +  mfabc                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - two* mfbcc *  vvx         +  mfacc                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dirE   ])[k   ] = mfabb;//(D.f[ dirE   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dirE   ])[k   ]                                                                     
			(D.f[ dirW   ])[kw  ] = mfcbb;//(D.f[ dirW   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dirW   ])[kw  ]                                                                   
			(D.f[ dirN   ])[k   ] = mfbab;//(D.f[ dirN   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ dirN   ])[k   ]
			(D.f[ dirS   ])[ks  ] = mfbcb;//(D.f[ dirS   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ dirS   ])[ks  ]
			(D.f[ dirT   ])[k   ] = mfbba;//(D.f[ dirT   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ dirT   ])[k   ]
			(D.f[ dirB   ])[kb  ] = mfbbc;//(D.f[ dirB   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ dirB   ])[kb  ]
			(D.f[ dirNE  ])[k   ] = mfaab;//(D.f[ dirNE  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dirNE  ])[k   ]
			(D.f[ dirSW  ])[ksw ] = mfccb;//(D.f[ dirSW  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dirSW  ])[ksw ]
			(D.f[ dirSE  ])[ks  ] = mfacb;//(D.f[ dirSE  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dirSE  ])[ks  ]
			(D.f[ dirNW  ])[kw  ] = mfcab;//(D.f[ dirNW  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dirNW  ])[kw  ]
			(D.f[ dirTE  ])[k   ] = mfaba;//(D.f[ dirTE  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dirTE  ])[k   ]
			(D.f[ dirBW  ])[kbw ] = mfcbc;//(D.f[ dirBW  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dirBW  ])[kbw ]
			(D.f[ dirBE  ])[kb  ] = mfabc;//(D.f[ dirBE  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dirBE  ])[kb  ]
			(D.f[ dirTW  ])[kw  ] = mfcba;//(D.f[ dirTW  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dirTW  ])[kw  ]
			(D.f[ dirTN  ])[k   ] = mfbaa;//(D.f[ dirTN  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ dirTN  ])[k   ]
			(D.f[ dirBS  ])[kbs ] = mfbcc;//(D.f[ dirBS  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ dirBS  ])[kbs ]
			(D.f[ dirBN  ])[kb  ] = mfbac;//(D.f[ dirBN  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ dirBN  ])[kb  ]
			(D.f[ dirTS  ])[ks  ] = mfbca;//(D.f[ dirTS  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ dirTS  ])[ks  ]
			(D.f[ dirZERO])[k   ] = mfbbb;//(D.f[ dirZERO])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ dirZERO])[k   ]
			(D.f[ dirTNE ])[k   ] = mfaaa;//(D.f[ dirTNE ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dirTNE ])[k   ]
			(D.f[ dirTSE ])[ks  ] = mfaca;//(D.f[ dirTSE ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dirTSE ])[ks  ]
			(D.f[ dirBNE ])[kb  ] = mfaac;//(D.f[ dirBNE ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dirBNE ])[kb  ]
			(D.f[ dirBSE ])[kbs ] = mfacc;//(D.f[ dirBSE ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dirBSE ])[kbs ]
			(D.f[ dirTNW ])[kw  ] = mfcaa;//(D.f[ dirTNW ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dirTNW ])[kw  ]
			(D.f[ dirTSW ])[ksw ] = mfcca;//(D.f[ dirTSW ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dirTSW ])[ksw ]
			(D.f[ dirBNW ])[kbw ] = mfcac;//(D.f[ dirBNW ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dirBNW ])[kbw ]
			(D.f[ dirBSW ])[kbsw] = mfccc;//(D.f[ dirBSW ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dirBSW ])[kbsw]
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Kum_New_Comp_Sponge_SP_27(real omegaIn,
															   unsigned int* bcMatD,
														       unsigned int* neighborX,
														       unsigned int* neighborY,
														       unsigned int* neighborZ,
													           real* coordX,
													           real* coordY,
													           real* coordZ,
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

	if(k<size_Mat)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			//unsigned int kzero= k;
			//unsigned int ke   = k;
			unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			unsigned int kbs  = neighborZ[ks];
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
			real mfcbb = (D.f[dirE   ])[k  ];
			real mfabb = (D.f[dirW   ])[kw ];
			real mfbcb = (D.f[dirN   ])[k  ];
			real mfbab = (D.f[dirS   ])[ks ];
			real mfbbc = (D.f[dirT   ])[k  ];
			real mfbba = (D.f[dirB   ])[kb ];
			real mfccb = (D.f[dirNE  ])[k  ];
			real mfaab = (D.f[dirSW  ])[ksw];
			real mfcab = (D.f[dirSE  ])[ks ];
			real mfacb = (D.f[dirNW  ])[kw ];
			real mfcbc = (D.f[dirTE  ])[k  ];
			real mfaba = (D.f[dirBW  ])[kbw];
			real mfcba = (D.f[dirBE  ])[kb ];
			real mfabc = (D.f[dirTW  ])[kw ];
			real mfbcc = (D.f[dirTN  ])[k  ];
			real mfbaa = (D.f[dirBS  ])[kbs];
			real mfbca = (D.f[dirBN  ])[kb ];
			real mfbac = (D.f[dirTS  ])[ks ];
			real mfbbb = (D.f[dirZERO])[k  ];
			real mfccc = (D.f[dirTNE ])[k  ];
			real mfaac = (D.f[dirTSW ])[ksw];
			real mfcac = (D.f[dirTSE ])[ks ];
			real mfacc = (D.f[dirTNW ])[kw ];
			real mfcca = (D.f[dirBNE ])[kb ];
			real mfaaa = (D.f[dirBSW ])[kbsw];
			real mfcaa = (D.f[dirBSE ])[kbs];
			real mfaca = (D.f[dirBNW ])[kbw];
			////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
							(((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
							((mfabb+mfcbb) + (mfbab+mfbcb) + (mfbba+mfbbc))) + mfbbb;

			real rho = one+drho;
			////////////////////////////////////////////////////////////////////////////////////
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) / rho;
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) / rho;
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			//the force be with you
			real fx = zero; //0.000000005;//(two/1600000.0) / 120.0; //zero;
			real fy = zero;
			real fz = zero;
			vvx += fx;
			vvy += fy;
			vvz += fz;
			////////////////////////////////////////////////////////////////////////////////////
			real oMdrho = one; // comp special
			////////////////////////////////////////////////////////////////////////////////////
			real m0, m1, m2;	
			real vx2 = vvx * vvx;
			real vy2 = vvy * vvy;
			real vz2 = vvz * vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real wadjust;
			real qudricLimit = 0.01f;
			////////////////////////////////////////////////////////////////////////////////////
			//sponge layer inflow
			real startXsponge = 1507.0f;//120.995703125;
			real endXsponge   = 1537.0f;//120.987890625;
			real sizeSponge = endXsponge - startXsponge; 
			real factor= one;
			real omega = factor * omegaIn;
			if(coordX[k] > startXsponge){
				factor = (((endXsponge - coordX[k]) / sizeSponge) * c1o2) + c1o2; 
				omega = factor * omegaIn;
			}
			////////////////////////////////////////////////////////////////////////////////////
			//sponge layer outflow
			endXsponge   = 30.0f;
			if(coordX[k] < endXsponge){
				factor = (((coordX[k]) / endXsponge) * c1o2) + c1o2; 
				omega = factor * omegaIn;
			}
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += one* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = one;

			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			real OxyyPxzz  = one;//two-omega;//(1000.*(-2. + omega))/(-1000. + 439.*omega);//(eight * (omega - two)) / (omega - eight);//two-omega;//eight*(two-omega)/(eight -omega);//one;//omega;//two-omega;//
			real OxyyMxzz  = one;//two-omega;//(1000.*(-2. + omega))/(-1000. + 439.*omega);//(eight * (omega - two)) / (omega - eight);//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
			//////////////////////////////
			//real OxyyPxzz  = two-omega;//
			//real OxyyMxzz  = two-omega;//
			//////////////////////////////
			//real OxyyPxzz  = (eight * (omega - two)) / (omega - eight);//Ginzburg
			//real OxyyMxzz  = (eight * (omega - two)) / (omega - eight);//Ginzburg
			//////////////////////////////
			//real OxyyPxzz  = omega;//BGK
			//real OxyyMxzz  = omega;//BGK
			//////////////////////////////
			//real OxyyPxzz  = (one + omega) / two;//1P5
			//real OxyyMxzz  = (one + omega) / two;//1P5
			//////////////////////////////
			//real OxyyPxzz  = (three - omega) / two;//0P5
			//real OxyyMxzz  = (three - omega) / two;//0P5
			//////////////////////////////
			//real OxyyPxzz  = (one + (eight * (omega - two)) / (omega - eight)) / two;//one + Ginzburg / two ... Car
			//real OxyyMxzz  = (one + (eight * (omega - two)) / (omega - eight)) / two;//one + Ginzburg / two ... Car
			////////////////////////////////////////////////////////////
			//4.
			//////////////////////////////
			real O4        = one;
			//////////////////////////////
			//real O4        = omega;//TRT
			////////////////////////////////////////////////////////////
			//5.
			//////////////////////////////
			real O5        = one;
			////////////////////////////////////////////////////////////
			//6.
			//////////////////////////////
			real O6        = one;
			////////////////////////////////////////////////////////////


			//central moments to cumulants
			//4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho;
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho;
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho;
				  	 		
			real CUMcca = mfcca - (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			real CUMcac = mfcac - (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			real CUMacc = mfacc - (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			real CUMccc = mfccc + ((-four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));
							//+ c1o27*(one -three/rho +two/(rho*rho)));





			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;
			
 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			//incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
 			{
 				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
 				real dyuy = dxux + omega * c3o2 * mxxMyy;
 				real dzuz = dxux + omega * c3o2 * mxxMzz;
 
 				//relax
 				mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 				mxxMyy    += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
 				mxxMzz    += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
 
 			}
 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			////no correction
 			//mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);//-magicBulk*OxxPyyPzz;
 			//mxxMyy    += -(-omega) * (-mxxMyy);
 			//mxxMzz    += -(-omega) * (-mxxMzz);
 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb     += omega * (-mfabb);
			mfbab     += omega * (-mfbab);
			mfbba     += omega * (-mfbba);

			// linear combinations back
			mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-two*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - two* mxxMzz + mxxPyyPzz);

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
			//das ist der limiter
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
 			mfbbb     += wadjust * (-mfbbb);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
 			mxxyPyzz  += wadjust * (-mxxyPyzz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
 			mxxyMyzz  += wadjust * (-mxxyMyzz);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
 			mxxzPyyz  += wadjust * (-mxxzPyyz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
 			mxxzMyyz  += wadjust * (-mxxzMyyz);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
 			mxyyPxzz  += wadjust * (-mxyyPxzz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
 			mxyyMxzz  += wadjust * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////
			mfbbb     += OxyyMxzz * (-mfbbb);
			mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
			mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
			mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
			mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
			mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
			mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);

			//// linear combinations back

			mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			//////////////////////////////////////////////////////////////////////////
			//mit limiter
 			wadjust    = O4+(one-O4)*abs(CUMacc)/(abs(CUMacc)+qudricLimit);
			CUMacc    += wadjust * (-CUMacc);
 			wadjust    = O4+(one-O4)*abs(CUMcac)/(abs(CUMcac)+qudricLimit);
			CUMcac    += wadjust * (-CUMcac); 
 			wadjust    = O4+(one-O4)*abs(CUMcca)/(abs(CUMcca)+qudricLimit);
			CUMcca    += wadjust * (-CUMcca); 

 			wadjust    = O4+(one-O4)*abs(CUMbbc)/(abs(CUMbbc)+qudricLimit);
			CUMbbc    += wadjust * (-CUMbbc); 
 			wadjust    = O4+(one-O4)*abs(CUMbcb)/(abs(CUMbcb)+qudricLimit);
			CUMbcb    += wadjust * (-CUMbcb); 
 			wadjust    = O4+(one-O4)*abs(CUMcbb)/(abs(CUMcbb)+qudricLimit);
			CUMcbb    += wadjust * (-CUMcbb); 
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
			//CUMacc += O4 * (-CUMacc); 
			//CUMcac += O4 * (-CUMcac); 
			//CUMcca += O4 * (-CUMcca); 

			//CUMbbc += O4 * (-CUMbbc); 
			//CUMbcb += O4 * (-CUMbcb); 
			//CUMcbb += O4 * (-CUMcbb); 
			//////////////////////////////////////////////////////////////////////////
			
					
			//5.
			CUMbcc += O5 * (-CUMbcc);
			CUMcbc += O5 * (-CUMcbc);
			CUMccb += O5 * (-CUMccb);

			//6.
			CUMccc += O6 * (-CUMccc);
			


			//back cumulants to central moments
			//4.
			mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho; 
			mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho; 
			mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho; 
						   
			mfcca = CUMcca + (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			mfcac = CUMcac + (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			mfacc = CUMacc + (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			mfccc = CUMccc - ((-four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));
							//+ c1o27*(one -three/rho +two/(rho*rho)));


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
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + one* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - two* mfaab *  vvz         +  mfaaa                * (one- vz2)              - one* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + one* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - two* mfabb *  vvz         + mfaba * (one- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - two* mfacb *  vvz         +  mfaca                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - two* mfbab *  vvz         + mfbaa * (one- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - two* mfbbb *  vvz         + mfbba * (one- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - two* mfbcb *  vvz         + mfbca * (one- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - two* mfcab *  vvz         +  mfcaa                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - two* mfcbb *  vvz         + mfcba * (one- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - two* mfccb *  vvz         +  mfcca                  * (one- vz2)              - c1o9 * oMdrho * vz2; 
			m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcca = m0;
			mfccb = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfaca        - two* mfaba *  vvy         +  mfaaa                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - two* mfabb *  vvy         +  mfaab                  * (one- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - two* mfabc *  vvy         +  mfaac                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - two* mfbba *  vvy         + mfbaa * (one- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - two* mfbbb *  vvy         + mfbab * (one- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - two* mfbbc *  vvy         + mfbac * (one- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - two* mfcba *  vvy         +  mfcaa                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - two* mfcbb *  vvy         +  mfcab                  * (one- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - two* mfcbc *  vvy         +  mfcac                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcac = m0;
			mfcbc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcaa        - two* mfbaa *  vvx         +  mfaaa                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - two* mfbba *  vvx         +  mfaba                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - two* mfbca *  vvx         +  mfaca                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - two* mfbab *  vvx         +  mfaab                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - two* mfbbb *  vvx         +  mfabb                  * (one- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - two* mfbcb *  vvx         +  mfacb                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - two* mfbac *  vvx         +  mfaac                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - two* mfbbc *  vvx         +  mfabc                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - two* mfbcc *  vvx         +  mfacc                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dirE   ])[k   ] = mfabb;                                                                     
			(D.f[ dirW   ])[kw  ] = mfcbb;                                                                   
			(D.f[ dirN   ])[k   ] = mfbab;
			(D.f[ dirS   ])[ks  ] = mfbcb;
			(D.f[ dirT   ])[k   ] = mfbba;
			(D.f[ dirB   ])[kb  ] = mfbbc;
			(D.f[ dirNE  ])[k   ] = mfaab;
			(D.f[ dirSW  ])[ksw ] = mfccb;
			(D.f[ dirSE  ])[ks  ] = mfacb;
			(D.f[ dirNW  ])[kw  ] = mfcab;
			(D.f[ dirTE  ])[k   ] = mfaba;
			(D.f[ dirBW  ])[kbw ] = mfcbc;
			(D.f[ dirBE  ])[kb  ] = mfabc;
			(D.f[ dirTW  ])[kw  ] = mfcba;
			(D.f[ dirTN  ])[k   ] = mfbaa;
			(D.f[ dirBS  ])[kbs ] = mfbcc;
			(D.f[ dirBN  ])[kb  ] = mfbac;
			(D.f[ dirTS  ])[ks  ] = mfbca;
			(D.f[ dirZERO])[k   ] = mfbbb;
			(D.f[ dirTNE ])[k   ] = mfaaa;
			(D.f[ dirTSE ])[ks  ] = mfaca;
			(D.f[ dirBNE ])[kb  ] = mfaac;
			(D.f[ dirBSE ])[kbs ] = mfacc;
			(D.f[ dirTNW ])[kw  ] = mfcaa;
			(D.f[ dirTSW ])[ksw ] = mfcca;
			(D.f[ dirBNW ])[kbw ] = mfcac;
			(D.f[ dirBSW ])[kbsw] = mfccc;
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Kum_IsoTest_SP_27( real omega,
														unsigned int* bcMatD,
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														real* DDStart,
														real* dxxUx,
														real* dyyUy,
														real* dzzUz,
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

	if(k<size_Mat)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			//unsigned int kzero= k;
			//unsigned int ke   = k;
			unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			unsigned int kbs  = neighborZ[ks];
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
			real mfcbb = (D.f[dirE   ])[k  ];//[ke   ];// +  c2over27 ;(D.f[dirE   ])[k  ];//ke
			real mfabb = (D.f[dirW   ])[kw ];//[kw   ];// +  c2over27 ;(D.f[dirW   ])[kw ];
			real mfbcb = (D.f[dirN   ])[k  ];//[kn   ];// +  c2over27 ;(D.f[dirN   ])[k  ];//kn
			real mfbab = (D.f[dirS   ])[ks ];//[ks   ];// +  c2over27 ;(D.f[dirS   ])[ks ];
			real mfbbc = (D.f[dirT   ])[k  ];//[kt   ];// +  c2over27 ;(D.f[dirT   ])[k  ];//kt
			real mfbba = (D.f[dirB   ])[kb ];//[kb   ];// +  c2over27 ;(D.f[dirB   ])[kb ];
			real mfccb = (D.f[dirNE  ])[k  ];//[kne  ];// +  c1over54 ;(D.f[dirNE  ])[k  ];//kne
			real mfaab = (D.f[dirSW  ])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dirSW  ])[ksw];
			real mfcab = (D.f[dirSE  ])[ks ];//[kse  ];// +  c1over54 ;(D.f[dirSE  ])[ks ];//kse
			real mfacb = (D.f[dirNW  ])[kw ];//[knw  ];// +  c1over54 ;(D.f[dirNW  ])[kw ];//knw
			real mfcbc = (D.f[dirTE  ])[k  ];//[kte  ];// +  c1over54 ;(D.f[dirTE  ])[k  ];//kte
			real mfaba = (D.f[dirBW  ])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dirBW  ])[kbw];
			real mfcba = (D.f[dirBE  ])[kb ];//[kbe  ];// +  c1over54 ;(D.f[dirBE  ])[kb ];//kbe
			real mfabc = (D.f[dirTW  ])[kw ];//[ktw  ];// +  c1over54 ;(D.f[dirTW  ])[kw ];//ktw
			real mfbcc = (D.f[dirTN  ])[k  ];//[ktn  ];// +  c1over54 ;(D.f[dirTN  ])[k  ];//ktn
			real mfbaa = (D.f[dirBS  ])[kbs];//[kbs  ];// +  c1over54 ;(D.f[dirBS  ])[kbs];
			real mfbca = (D.f[dirBN  ])[kb ];//[kbn  ];// +  c1over54 ;(D.f[dirBN  ])[kb ];//kbn
			real mfbac = (D.f[dirTS  ])[ks ];//[kts  ];// +  c1over54 ;(D.f[dirTS  ])[ks ];//kts
			real mfbbb = (D.f[dirZERO])[k  ];//[kzero];// +  c8over27 ;(D.f[dirZERO])[k  ];//kzero
			real mfccc = (D.f[dirTNE ])[k  ];//[ktne ];// +  c1over216;(D.f[dirTNE ])[k  ];//ktne
			real mfaac = (D.f[dirTSW ])[ksw];//[ktsw ];// +  c1over216;(D.f[dirTSW ])[ksw];//ktsw
			real mfcac = (D.f[dirTSE ])[ks ];//[ktse ];// +  c1over216;(D.f[dirTSE ])[ks ];//ktse
			real mfacc = (D.f[dirTNW ])[kw ];//[ktnw ];// +  c1over216;(D.f[dirTNW ])[kw ];//ktnw
			real mfcca = (D.f[dirBNE ])[kb ];//[kbne ];// +  c1over216;(D.f[dirBNE ])[kb ];//kbne
			real mfaaa = (D.f[dirBSW ])[kbsw];//[kbsw ];// +  c1over216;(D.f[dirBSW ])[kbsw];
			real mfcaa = (D.f[dirBSE ])[kbs];//[kbse ];// +  c1over216;(D.f[dirBSE ])[kbs];//kbse
			real mfaca = (D.f[dirBNW ])[kbw];//[kbnw ];// +  c1over216;(D.f[dirBNW ])[kbw];//kbnw
			////////////////////////////////////////////////////////////////////////////////////
			//slow
			//real oMdrho = one - ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
			//					   (((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
			//						((mfabb+mfcbb) + (mfbab+mfbcb)  +  (mfbba+mfbbc)));//fehlt mfbbb
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb));
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab));
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba));
			////////////////////////////////////////////////////////////////////////////////////
			//fast
			real oMdrho = one - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
								   mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
								   mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb);//fehlt mfbbb nicht mehr
			//real vvx    =mfccc-mfaaa + mfcac-mfaca + mfcaa-mfacc + mfcca-mfaac + 
			//				mfcba-mfabc + mfcbc-mfaba + mfcab-mfacb + mfccb-mfaab +
			//				mfcbb-mfabb;
			//real vvy    =mfccc-mfaaa + mfaca-mfcac + mfacc-mfcaa + mfcca-mfaac + 
			//				mfbca-mfbac + mfbcc-mfbaa + mfacb-mfcab + mfccb-mfaab +
			//				mfbcb-mfbab;
			//real vvz    =mfccc-mfaaa + mfcac-mfaca + mfacc-mfcaa + mfaac-mfcca + 
			//				mfbac-mfbca + mfbcc-mfbaa + mfabc-mfcba + mfcbc-mfaba +
			//				mfbbc-mfbba;
			////////////////////////////////////////////////////////////////////////////////////
			// oMdrho assembler style -------> faaaaaastaaaa
			// or much sloooowaaaa ... it depändssssss on sadaku
			real m0, m1, m2;	
			//real oMdrho;
			//{
			//	oMdrho=mfccc+mfaaa;
			//	m0=mfaca+mfcac;
			//	m1=mfacc+mfcaa;
			//	m2=mfaac+mfcca;
			//	oMdrho+=m0;
			//	m1+=m2;
			//	oMdrho+=m1;
			//	m0=mfbac+mfbca;
			//	m1=mfbaa+mfbcc;
			//	m0+=m1;
			//	m1=mfabc+mfcba;
			//	m2=mfaba+mfcbc;
			//	m1+=m2;
			//	m0+=m1;
			//	m1=mfacb+mfcab;
			//	m2=mfaab+mfccb;
			//	m1+=m2;
			//	m0+=m1;
			//	oMdrho+=m0;
			//	m0=mfabb+mfcbb;
			//	m1=mfbab+mfbcb;
			//	m2=mfbba+mfbbc;
			//	m0+=m1+m2;
			//	m0+=mfbbb; //hat gefehlt
			//	oMdrho = one - (oMdrho + m0);
			//}
			//real vvx;
			real vx2;
			//{
			//	vvx = mfccc-mfaaa;
			//	m0  = mfcac-mfaca;
			//	m1  = mfcaa-mfacc;
			//	m2  = mfcca-mfaac;
			//	vvx+= m0;
			//	m1 += m2;
			//	vvx+= m1;
			//	vx2 = mfcba-mfabc;
			//	m0  = mfcbc-mfaba;
			//	m1  = mfcab-mfacb;
			//	m2  = mfccb-mfaab;
			//	vx2+= m0;
			//	m1 += m2;
			//	vx2+= m1;
			//	vvx+= vx2;
			//	vx2 = mfcbb-mfabb;
			//	vvx+= vx2;
			//}
			//real vvy;
			real vy2;
			//{
			//	vvy = mfccc-mfaaa;
			//	m0  = mfaca-mfcac;
			//	m1  = mfacc-mfcaa;
			//	m2  = mfcca-mfaac;
			//	vvy+= m0;
			//	m1 += m2;
			//	vvy+= m1;
			//	vy2 = mfbca-mfbac;
			//	m0  = mfbcc-mfbaa;
			//	m1  = mfacb-mfcab;
			//	m2  = mfccb-mfaab;
			//	vy2+= m0;
			//	m1 += m2;
			//	vy2+= m1;
			//	vvy+= vy2;
			//	vy2 = mfbcb-mfbab;
			//	vvy+= vy2;
			//}
			//real vvz;
			real vz2;
			//{
			//	vvz = mfccc-mfaaa;
			//	m0  = mfcac-mfaca;
			//	m1  = mfacc-mfcaa;
			//	m2  = mfaac-mfcca;
			//	vvz+= m0;
			//	m1 += m2;
			//	vvz+= m1;
			//	vz2 = mfbac-mfbca;
			//	m0  = mfbcc-mfbaa;
			//	m1  = mfabc-mfcba;
			//	m2  = mfcbc-mfaba;
			//	vz2+= m0;
			//	m1 += m2;
			//	vz2+= m1;
			//	vvz+= vz2;
			//	vz2 = mfbbc-mfbba;
			//	vvz+= vz2;
			//}
			vx2=vvx*vvx;
			vy2=vvy*vvy;
			vz2=vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real wadjust;
			real qudricLimit = 0.01f;
			//real s9 = minusomega;
			//test
			//s9 = 0.;
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += one* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// BGK
			////////////////////////////////////////////////////////////////////////////////////
			////2.
			//mfabb += omega * (-mfabb);
			//mfbab += omega * (-mfbab);
			//mfbba += omega * (-mfbba);
			//
			//mfcaa += omega * (c1o3 * mfaaa - mfcaa);
			//mfaca += omega * (c1o3 * mfaaa - mfaca);
			//mfaac += omega * (c1o3 * mfaaa - mfaac);
			//
			////3.
			//mfabc += omega * (-mfabc);
			//mfbac += omega * (-mfbac);
			//
			//mfacb += omega * (-mfacb);
			//mfbca += omega * (-mfbca);

			//mfcab += omega * (-mfcab);
			//mfcba += omega * (-mfcba);

			//mfbbb += omega * (-mfbbb);

			////4.
			//mfacc += omega * (c1o9 * mfaaa - mfacc);
			//mfcac += omega * (c1o9 * mfaaa - mfcac);
			//mfcca += omega * (c1o9 * mfaaa - mfcca);

			//mfbbc += omega * (-mfbbc);
			//mfbcb += omega * (-mfbcb);
			//mfcbb += omega * (-mfcbb);

			////5.
			//mfbcc += omega * (-mfbcc);
			//mfcbc += omega * (-mfcbc);
			//mfccb += omega * (-mfccb);

			////6.
			//mfccc += omega * (c1o27 * mfaaa - mfccc);
			//////////////////////////////////////////////////////////////////////////////////////



			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = one;
			real OxyyPxzz  = one;//two-omega;//eight*(two-omega)/(eight -omega);//one;//omega;//two-omega;//
			real OxyyMxzz  = one;//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
			real O4        = one;
			real O5        = one;
			real O6        = one;

			//Cum 4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab);
			real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb);
			real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb);

			real CUMcca = mfcca - (mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho;
			real CUMcac = mfcac - (mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-one)*oMdrho;
			real CUMacc = mfacc - (mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho;

			//Cum 5.
			real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
			real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
			real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

			//Cum 6.
			real CUMccc = mfccc  +((-four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb)
							-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
							-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(one-two* oMdrho)- c1o27* oMdrho * oMdrho*(-two* oMdrho)
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) +c1o27*oMdrho;





			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Iso Test Part 1
			real dxuydyux = -three * omega * mfbba;
			real dxuzdzux = -three * omega * mfbab;
			real dyuzdzuy = -three * omega * mfabb;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			{
				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
				real dyuy = dxux + omega * c3o2 * mxxMyy;
				real dzuz = dxux + omega * c3o2 * mxxMzz;

				//relax
				mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
				mxxMyy    += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
				mxxMzz    += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

			}

			//mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);
			//mxxMyy    += -(-omega) * (-mxxMyy);
			//mxxMzz    += -(-omega) * (-mxxMzz);
			mfabb     += omega * (-mfabb);
			mfbab     += omega * (-mfbab);
			mfbba     += omega * (-mfbba);

			// linear combinations back
			mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-two*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - two* mxxMzz + mxxPyyPzz);

			//3.
			// linear combinations

			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Iso Test Part 2
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//precollision terms 3. moments
			real premxxyPyzz = mxxyPyzz;
			real premxxzPyyz = mxxzPyyz;
			real premxyyPxzz = mxyyPxzz;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			//relax
			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
			mfbbb     += wadjust * (-mfbbb);
			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
			mxxyPyzz  += wadjust * (-mxxyPyzz);
			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
			mxxyMyzz  += wadjust * (-mxxyMyzz);
			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
			mxxzPyyz  += wadjust * (-mxxzPyyz);
			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
			mxxzMyyz  += wadjust * (-mxxzMyyz);
			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
			mxyyPxzz  += wadjust * (-mxyyPxzz);
			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
			mxyyMxzz  += wadjust * (-mxyyMxzz);

			//// linear combinations back
			mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			CUMacc += O4 * (-CUMacc); 
			CUMcac += O4 * (-CUMcac); 
			CUMcca += O4 * (-CUMcca); 
			
			CUMbbc += O4 * (-CUMbbc); 
			CUMbcb += O4 * (-CUMbcb); 
			CUMcbb += O4 * (-CUMcbb); 
					
			//5.
			//CUMbcc += O5 * (-CUMbcc);
			//CUMcbc += O5 * (-CUMcbc);
			//CUMccb += O5 * (-CUMccb);

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Iso Test Part 3
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//precollision terms 5. moments
			real preCUMbcc = CUMbcc;
			real preCUMcbc = CUMcbc;
			real preCUMccb = CUMccb;
			//new calculation of 5. moments
			CUMbcc = c2o3 * (one - c1o2 * O5) * (vvy * dxuydyux + vvz * dxuzdzux) * zero + (one - O5) * CUMbcc;
			CUMcbc = c2o3 * (one - c1o2 * O5) * (vvx * dxuydyux + vvz * dyuzdzuy) * zero + (one - O5) * CUMcbc;
			CUMccb = c2o3 * (one - c1o2 * O5) * (vvx * dxuzdzux + vvy * dyuzdzuy) * zero + (one - O5) * CUMccb;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			//6.
			//CUMccc += O6 * (-CUMccc);
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Iso Test Part 4
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//new calculation of 6. moment
			CUMccc = O6 * c1o3 * (vx2 + vy2 + vz2) * zero + (one - O6) * CUMccc;
			// second derivation of ux, uy, uz
			real dxxux = three * (c3o2 * (CUMbcc - preCUMbcc) - zero * (mxyyPxzz - premxyyPxzz));
			real dyyuy = three * (c3o2 * (CUMcbc - preCUMcbc) - zero * (mxxyPyzz - premxxyPyzz));
			real dzzuz = three * (c3o2 * (CUMccb - preCUMccb) - zero * (mxxzPyyz - premxxzPyyz));
			// copy local values to global arrays for paraview files 
			dxxUx[k] = dxxux;
			dyyUy[k] = dyyuy;
			dzzUz[k] = dzzuz;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			//back cumulants to central moments
			//4.
			mfcbb = CUMcbb + ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab);
			mfbcb = CUMbcb + ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb);
			mfbbc = CUMbbc + ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb); 
						   
			mfcca = CUMcca + (mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho;
			mfcac = CUMcac + (mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-one)*oMdrho;
			mfacc = CUMacc + (mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho;

			//5.
			mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) * oMdrho;
			mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) * oMdrho;
			mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) * oMdrho;
			
			//6.
			mfccc = CUMccc  -(( -four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb)
							-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
							-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(one-two* oMdrho)- c1o27* oMdrho * oMdrho*(-two* oMdrho)
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) -c1o27*oMdrho;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			//back
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + one* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - two* mfaab *  vvz         +  mfaaa                * (one- vz2)              - one* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + one* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - two* mfabb *  vvz         + mfaba * (one- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - two* mfacb *  vvz         +  mfaca                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - two* mfbab *  vvz         + mfbaa * (one- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - two* mfbbb *  vvz         + mfbba * (one- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - two* mfbcb *  vvz         + mfbca * (one- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - two* mfcab *  vvz         +  mfcaa                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - two* mfcbb *  vvz         + mfcba * (one- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - two* mfccb *  vvz         +  mfcca                  * (one- vz2)              - c1o9 * oMdrho * vz2; 
			m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcca = m0;
			mfccb = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfaca        - two* mfaba *  vvy         +  mfaaa                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - two* mfabb *  vvy         +  mfaab                  * (one- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - two* mfabc *  vvy         +  mfaac                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - two* mfbba *  vvy         + mfbaa * (one- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - two* mfbbb *  vvy         + mfbab * (one- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - two* mfbbc *  vvy         + mfbac * (one- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - two* mfcba *  vvy         +  mfcaa                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - two* mfcbb *  vvy         +  mfcab                  * (one- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - two* mfcbc *  vvy         +  mfcac                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcac = m0;
			mfcbc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcaa        - two* mfbaa *  vvx         +  mfaaa                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - two* mfbba *  vvx         +  mfaba                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - two* mfbca *  vvx         +  mfaca                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - two* mfbab *  vvx         +  mfaab                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - two* mfbbb *  vvx         +  mfabb                  * (one- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - two* mfbcb *  vvx         +  mfacb                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - two* mfbac *  vvx         +  mfaac                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - two* mfbbc *  vvx         +  mfabc                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - two* mfbcc *  vvx         +  mfacc                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dirE   ])[k   ] = mfabb;//(D.f[ dirE   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dirE   ])[k   ]                                                                     
			(D.f[ dirW   ])[kw  ] = mfcbb;//(D.f[ dirW   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dirW   ])[kw  ]                                                                   
			(D.f[ dirN   ])[k   ] = mfbab;//(D.f[ dirN   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ dirN   ])[k   ]
			(D.f[ dirS   ])[ks  ] = mfbcb;//(D.f[ dirS   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ dirS   ])[ks  ]
			(D.f[ dirT   ])[k   ] = mfbba;//(D.f[ dirT   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ dirT   ])[k   ]
			(D.f[ dirB   ])[kb  ] = mfbbc;//(D.f[ dirB   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ dirB   ])[kb  ]
			(D.f[ dirNE  ])[k   ] = mfaab;//(D.f[ dirNE  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dirNE  ])[k   ]
			(D.f[ dirSW  ])[ksw ] = mfccb;//(D.f[ dirSW  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dirSW  ])[ksw ]
			(D.f[ dirSE  ])[ks  ] = mfacb;//(D.f[ dirSE  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dirSE  ])[ks  ]
			(D.f[ dirNW  ])[kw  ] = mfcab;//(D.f[ dirNW  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dirNW  ])[kw  ]
			(D.f[ dirTE  ])[k   ] = mfaba;//(D.f[ dirTE  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dirTE  ])[k   ]
			(D.f[ dirBW  ])[kbw ] = mfcbc;//(D.f[ dirBW  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dirBW  ])[kbw ]
			(D.f[ dirBE  ])[kb  ] = mfabc;//(D.f[ dirBE  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dirBE  ])[kb  ]
			(D.f[ dirTW  ])[kw  ] = mfcba;//(D.f[ dirTW  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dirTW  ])[kw  ]
			(D.f[ dirTN  ])[k   ] = mfbaa;//(D.f[ dirTN  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ dirTN  ])[k   ]
			(D.f[ dirBS  ])[kbs ] = mfbcc;//(D.f[ dirBS  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ dirBS  ])[kbs ]
			(D.f[ dirBN  ])[kb  ] = mfbac;//(D.f[ dirBN  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ dirBN  ])[kb  ]
			(D.f[ dirTS  ])[ks  ] = mfbca;//(D.f[ dirTS  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ dirTS  ])[ks  ]
			(D.f[ dirZERO])[k   ] = mfbbb;//(D.f[ dirZERO])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ dirZERO])[k   ]
			(D.f[ dirTNE ])[k   ] = mfaaa;//(D.f[ dirTNE ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dirTNE ])[k   ]
			(D.f[ dirTSE ])[ks  ] = mfaca;//(D.f[ dirTSE ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dirTSE ])[ks  ]
			(D.f[ dirBNE ])[kb  ] = mfaac;//(D.f[ dirBNE ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dirBNE ])[kb  ]
			(D.f[ dirBSE ])[kbs ] = mfacc;//(D.f[ dirBSE ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dirBSE ])[kbs ]
			(D.f[ dirTNW ])[kw  ] = mfcaa;//(D.f[ dirTNW ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dirTNW ])[kw  ]
			(D.f[ dirTSW ])[ksw ] = mfcca;//(D.f[ dirTSW ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dirTSW ])[ksw ]
			(D.f[ dirBNW ])[kbw ] = mfcac;//(D.f[ dirBNW ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dirBNW ])[kbw ]
			(D.f[ dirBSW ])[kbsw] = mfccc;//(D.f[ dirBSW ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dirBSW ])[kbsw]
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Kum_1h_SP_27(  real omega,
													real deltaPhi,
													real angularVelocity,
													unsigned int* bcMatD,
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													real* coordX,
													real* coordY,
													real* coordZ,
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

	if(k<size_Mat)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			//unsigned int kzero= k;
			//unsigned int ke   = k;
			unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			unsigned int kbs  = neighborZ[ks];
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

			//unsigned int kzero= k;
			//unsigned int ke   = k;
			//unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			//unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			//unsigned int kb   = neighborZ[k];
			//unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			//unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			//unsigned int kbs  = neighborZ[ks];
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
			//unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[dirE   ])[k  ];//[ke   ];// +  c2over27 ;(D.f[dirE   ])[k  ];//ke
			real mfabb = (D.f[dirW   ])[kw ];//[kw   ];// +  c2over27 ;(D.f[dirW   ])[kw ];
			real mfbcb = (D.f[dirN   ])[k  ];//[kn   ];// +  c2over27 ;(D.f[dirN   ])[k  ];//kn
			real mfbab = (D.f[dirS   ])[ks ];//[ks   ];// +  c2over27 ;(D.f[dirS   ])[ks ];
			real mfbbc = (D.f[dirT   ])[k  ];//[kt   ];// +  c2over27 ;(D.f[dirT   ])[k  ];//kt
			real mfbba = (D.f[dirB   ])[kb ];//[kb   ];// +  c2over27 ;(D.f[dirB   ])[kb ];
			real mfccb = (D.f[dirNE  ])[k  ];//[kne  ];// +  c1over54 ;(D.f[dirNE  ])[k  ];//kne
			real mfaab = (D.f[dirSW  ])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dirSW  ])[ksw];
			real mfcab = (D.f[dirSE  ])[ks ];//[kse  ];// +  c1over54 ;(D.f[dirSE  ])[ks ];//kse
			real mfacb = (D.f[dirNW  ])[kw ];//[knw  ];// +  c1over54 ;(D.f[dirNW  ])[kw ];//knw
			real mfcbc = (D.f[dirTE  ])[k  ];//[kte  ];// +  c1over54 ;(D.f[dirTE  ])[k  ];//kte
			real mfaba = (D.f[dirBW  ])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dirBW  ])[kbw];
			real mfcba = (D.f[dirBE  ])[kb ];//[kbe  ];// +  c1over54 ;(D.f[dirBE  ])[kb ];//kbe
			real mfabc = (D.f[dirTW  ])[kw ];//[ktw  ];// +  c1over54 ;(D.f[dirTW  ])[kw ];//ktw
			real mfbcc = (D.f[dirTN  ])[k  ];//[ktn  ];// +  c1over54 ;(D.f[dirTN  ])[k  ];//ktn
			real mfbaa = (D.f[dirBS  ])[kbs];//[kbs  ];// +  c1over54 ;(D.f[dirBS  ])[kbs];
			real mfbca = (D.f[dirBN  ])[kb ];//[kbn  ];// +  c1over54 ;(D.f[dirBN  ])[kb ];//kbn
			real mfbac = (D.f[dirTS  ])[ks ];//[kts  ];// +  c1over54 ;(D.f[dirTS  ])[ks ];//kts
			real mfbbb = (D.f[dirZERO])[k  ];//[kzero];// +  c8over27 ;(D.f[dirZERO])[k  ];//kzero
			real mfccc = (D.f[dirTNE ])[k  ];//[ktne ];// +  c1over216;(D.f[dirTNE ])[k  ];//ktne
			real mfaac = (D.f[dirTSW ])[ksw];//[ktsw ];// +  c1over216;(D.f[dirTSW ])[ksw];//ktsw
			real mfcac = (D.f[dirTSE ])[ks ];//[ktse ];// +  c1over216;(D.f[dirTSE ])[ks ];//ktse
			real mfacc = (D.f[dirTNW ])[kw ];//[ktnw ];// +  c1over216;(D.f[dirTNW ])[kw ];//ktnw
			real mfcca = (D.f[dirBNE ])[kb ];//[kbne ];// +  c1over216;(D.f[dirBNE ])[kb ];//kbne
			real mfaaa = (D.f[dirBSW ])[kbsw];//[kbsw ];// +  c1over216;(D.f[dirBSW ])[kbsw];
			real mfcaa = (D.f[dirBSE ])[kbs];//[kbse ];// +  c1over216;(D.f[dirBSE ])[kbs];//kbse
			real mfaca = (D.f[dirBNW ])[kbw];//[kbnw ];// +  c1over216;(D.f[dirBNW ])[kbw];//kbnw
			////////////////////////////////////////////////////////////////////////////////////
			//Ship
			real coord0X = 281.125f;//7.5f;
			real coord0Y = 388.125f;//7.5f;
			real ux = - angularVelocity * (coordY[k] - coord0Y);
			real uy =   angularVelocity * (coordX[k] - coord0X);
			////////////////////////////////////////////////////////////////////////////////////
			//slow
			//real oMdrho = one - ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
			//					   (((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
			//						((mfabb+mfcbb) + (mfbab+mfbcb)  +  (mfbba+mfbbc)));//fehlt mfbbb
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
				(((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
				(mfcbb-mfabb));
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				(((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				(mfbcb-mfbab));
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				(((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				(mfbbc-mfbba));
			////////////////////////////////////////////////////////////////////////////////////
			real vxNeu = cosf(deltaPhi) * vvx - sinf(deltaPhi) * vvy; 
			real vyNeu = sinf(deltaPhi) * vvx + cosf(deltaPhi) * vvy; 




			////////////////////////////////////////////////////////////////////////////////////
			//fast
			real oMdrho = one - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
				mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
				mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb);//fehlt mfbbb nicht mehr
			//real vvx    =mfccc-mfaaa + mfcac-mfaca + mfcaa-mfacc + mfcca-mfaac + 
			//				mfcba-mfabc + mfcbc-mfaba + mfcab-mfacb + mfccb-mfaab +
			//				mfcbb-mfabb;
			//real vvy    =mfccc-mfaaa + mfaca-mfcac + mfacc-mfcaa + mfcca-mfaac + 
			//				mfbca-mfbac + mfbcc-mfbaa + mfacb-mfcab + mfccb-mfaab +
			//				mfbcb-mfbab;
			//real vvz    =mfccc-mfaaa + mfcac-mfaca + mfacc-mfcaa + mfaac-mfcca + 
			//				mfbac-mfbca + mfbcc-mfbaa + mfabc-mfcba + mfcbc-mfaba +
			//				mfbbc-mfbba;
			////////////////////////////////////////////////////////////////////////////////////
			// oMdrho assembler style -------> faaaaaastaaaa
			// or much sloooowaaaa ... it depändssssss on sadaku
			real m0, m1, m2;	
			//real oMdrho;
			//{
			//	oMdrho=mfccc+mfaaa;
			//	m0=mfaca+mfcac;
			//	m1=mfacc+mfcaa;
			//	m2=mfaac+mfcca;
			//	oMdrho+=m0;
			//	m1+=m2;
			//	oMdrho+=m1;
			//	m0=mfbac+mfbca;
			//	m1=mfbaa+mfbcc;
			//	m0+=m1;
			//	m1=mfabc+mfcba;
			//	m2=mfaba+mfcbc;
			//	m1+=m2;
			//	m0+=m1;
			//	m1=mfacb+mfcab;
			//	m2=mfaab+mfccb;
			//	m1+=m2;
			//	m0+=m1;
			//	oMdrho+=m0;
			//	m0=mfabb+mfcbb;
			//	m1=mfbab+mfbcb;
			//	m2=mfbba+mfbbc;
			//	m0+=m1+m2;
			//	m0+=mfbbb; //hat gefehlt
			//	oMdrho = one - (oMdrho + m0);
			//}
			//real vvx;
			real vx2;
			//{
			//	vvx = mfccc-mfaaa;
			//	m0  = mfcac-mfaca;
			//	m1  = mfcaa-mfacc;
			//	m2  = mfcca-mfaac;
			//	vvx+= m0;
			//	m1 += m2;
			//	vvx+= m1;
			//	vx2 = mfcba-mfabc;
			//	m0  = mfcbc-mfaba;
			//	m1  = mfcab-mfacb;
			//	m2  = mfccb-mfaab;
			//	vx2+= m0;
			//	m1 += m2;
			//	vx2+= m1;
			//	vvx+= vx2;
			//	vx2 = mfcbb-mfabb;
			//	vvx+= vx2;
			//}
			//real vvy;
			real vy2;
			//{
			//	vvy = mfccc-mfaaa;
			//	m0  = mfaca-mfcac;
			//	m1  = mfacc-mfcaa;
			//	m2  = mfcca-mfaac;
			//	vvy+= m0;
			//	m1 += m2;
			//	vvy+= m1;
			//	vy2 = mfbca-mfbac;
			//	m0  = mfbcc-mfbaa;
			//	m1  = mfacb-mfcab;
			//	m2  = mfccb-mfaab;
			//	vy2+= m0;
			//	m1 += m2;
			//	vy2+= m1;
			//	vvy+= vy2;
			//	vy2 = mfbcb-mfbab;
			//	vvy+= vy2;
			//}
			//real vvz;
			real vz2;
			//{
			//	vvz = mfccc-mfaaa;
			//	m0  = mfcac-mfaca;
			//	m1  = mfacc-mfcaa;
			//	m2  = mfaac-mfcca;
			//	vvz+= m0;
			//	m1 += m2;
			//	vvz+= m1;
			//	vz2 = mfbac-mfbca;
			//	m0  = mfbcc-mfbaa;
			//	m1  = mfabc-mfcba;
			//	m2  = mfcbc-mfaba;
			//	vz2+= m0;
			//	m1 += m2;
			//	vz2+= m1;
			//	vvz+= vz2;
			//	vz2 = mfbbc-mfbba;
			//	vvz+= vz2;
			//}
			vx2=vvx*vvx;
			vy2=vvy*vvy;
			vz2=vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real wadjust;
			real qudricLimit = 0.01f;
			//real s9 = minusomega;
			//test
			//s9 = 0.;
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += one* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// BGK
			////////////////////////////////////////////////////////////////////////////////////
			////2.
			//mfabb += omega * (-mfabb);
			//mfbab += omega * (-mfbab);
			//mfbba += omega * (-mfbba);
			//
			//mfcaa += omega * (c1o3 * mfaaa - mfcaa);
			//mfaca += omega * (c1o3 * mfaaa - mfaca);
			//mfaac += omega * (c1o3 * mfaaa - mfaac);
			//
			////3.
			//mfabc += omega * (-mfabc);
			//mfbac += omega * (-mfbac);
			//
			//mfacb += omega * (-mfacb);
			//mfbca += omega * (-mfbca);

			//mfcab += omega * (-mfcab);
			//mfcba += omega * (-mfcba);

			//mfbbb += omega * (-mfbbb);

			////4.
			//mfacc += omega * (c1o9 * mfaaa - mfacc);
			//mfcac += omega * (c1o9 * mfaaa - mfcac);
			//mfcca += omega * (c1o9 * mfaaa - mfcca);

			//mfbbc += omega * (-mfbbc);
			//mfbcb += omega * (-mfbcb);
			//mfcbb += omega * (-mfcbb);

			////5.
			//mfbcc += omega * (-mfbcc);
			//mfcbc += omega * (-mfcbc);
			//mfccb += omega * (-mfccb);

			////6.
			//mfccc += omega * (c1o27 * mfaaa - mfccc);
			////////////////////////////////////////////////////////////////////////////////////



			//////////////////////////////////////////////////////////////////////////////////////////
			//////// Cumulants
			//////////////////////////////////////////////////////////////////////////////////////////
			//////real OxxPyyPzz = one;
			//////real OxyyPxzz  = one;//two-omega;//
			//////real OxyyMxzz  = one;//two-omega;//
			//////real O4        = one;
			//////real O5        = one;
			//////real O6        = one;

			////////Cum 4.
			//////real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + two* mfbba * mfbab);
			//////real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + two* mfbba * mfabb);
			//////real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + two* mfbab * mfabb); 

			//////real CUMcca = mfcca - (mfcaa * mfaca + two* mfbba * mfbba)- c1o3 * (mfcaa + mfaca);
			//////real CUMcac = mfcac - (mfcaa * mfaac + two* mfbab * mfbab)- c1o3 * (mfcaa + mfaac);
			//////real CUMacc = mfacc - (mfaac * mfaca + two* mfabb * mfabb)- c1o3 * (mfaac + mfaca);

			////////Cum 5.
			//////real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + four* mfabb * mfbbb + two* (mfbab * mfacb + mfbba * mfabc)) //O(eps^5) 
			//////				- c1o3 * (mfbca + mfbac); //O(eps^3)
			//////real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + four* mfbab * mfbbb + two* (mfabb * mfcab + mfbba * mfbac)) //O(eps^5)
			//////				- c1o3 * (mfcba + mfabc); //O(eps^3)
			//////real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + four* mfbba * mfbbb + two* (mfbab * mfbca + mfabb * mfcba)) //O(eps^5)
			//////				- c1o3 * (mfacb + mfcab);//O(eps^3)

			////////Cum 6.
			//////real CUMccc = mfccc +(-four*  mfbbb * mfbbb  //O(eps^6)
			//////	                   -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca) // O(eps^4)
			//////					   -  four* (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc) // O(eps^6) 
			//////					   -  two* (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) // O(esp^6)
			//////					   +( four* (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac) //O(eps^6)
			//////					   +  two* (mfcaa * mfaca * mfaac) //O(eps^6)
			//////					   + sixteen*  mfbba * mfbab * mfabb) //O(eps^6)
			//////					   - c1o3* (mfacc + mfcac + mfcca) //O(eps^2)
			//////					   + c1o9* (mfcaa + mfaca + mfaac) //O(eps^2)
			//////					   +( two* (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)//O(eps^4)
			//////					   +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3;//O(eps^4)


			////////2.
			//////// linear combinations
			//////real mxxPyyPzz = mfcaa + mfaca + mfaac;
			//////real mxxMyy    = mfcaa - mfaca;
			//////real mxxMzz	   = mfcaa - mfaac;

			//////{
			//////	real dxux = c1o2 * ((-omega) * (mxxMyy + mxxMzz) + (mfaaa - mxxPyyPzz));
			//////	real dyuy = dxux + omega * c3o2 * mxxMyy;
			//////	real dzuz = dxux + omega * c3o2 * mxxMzz;

			//////	//relax
			//////	mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three* (one- c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
			//////	mxxMyy    += omega * (-mxxMyy) - three* (one+ c1o2 * (-omega)) * (vx2 * dxux + vy2 * dyuy);
			//////	mxxMzz    += omega * (-mxxMzz) - three* (one+ c1o2 * (-omega)) * (vx2 * dxux + vz2 * dzuz);
			//////}
			//////mfabb     += omega * (-mfabb);
			//////mfbab     += omega * (-mfbab);
			//////mfbba     += omega * (-mfbba);

			//////// linear combinations back
			//////mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			//////mfaca = c1o3 * (-two*  mxxMyy +      mxxMzz + mxxPyyPzz);
			//////mfaac = c1o3 * (       mxxMyy - two* mxxMzz + mxxPyyPzz);

			////////3.
			//////// linear combinations
			//////real mxxyPyzz = mfcba + mfabc;
			//////real mxxyMyzz = mfcba - mfabc;

			//////real mxxzPyyz = mfcab + mfacb;
			//////real mxxzMyyz = mfcab - mfacb;

			//////real mxyyPxzz = mfbca + mfbac;
			//////real mxyyMxzz = mfbca - mfbac;

			////////relax
			//////wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
			//////mfbbb     += wadjust * (-mfbbb);
			//////wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
			//////mxxyPyzz  += wadjust * (-mxxyPyzz);
			//////wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
			//////mxxyMyzz  += wadjust * (-mxxyMyzz);
			//////wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
			//////mxxzPyyz  += wadjust * (-mxxzPyyz);
			//////wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
			//////mxxzMyyz  += wadjust * (-mxxzMyyz);
			//////wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
			//////mxyyPxzz  += wadjust * (-mxyyPxzz);
			//////wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
			//////mxyyMxzz  += wadjust * (-mxyyMxzz);

			//////// linear combinations back
			//////mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
			//////mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			//////mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
			//////mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			//////mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
			//////mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			////////4.
			//////CUMacc += O4 * (-CUMacc); 
			//////CUMcac += O4 * (-CUMcac); 
			//////CUMcca += O4 * (-CUMcca); 
			//////
			//////CUMbbc += O4 * (-CUMbbc); 
			//////CUMbcb += O4 * (-CUMbcb); 
			//////CUMcbb += O4 * (-CUMcbb); 
			//////		
			////////5.
			//////CUMbcc += O5 * (-CUMbcc);
			//////CUMcbc += O5 * (-CUMcbc);
			//////CUMccb += O5 * (-CUMccb);

			////////6.
			//////CUMccc += O6 * (-CUMccc);
			//////
			////////back cumulants to central moments
			////////4.
			//////mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + two* mfbba * mfbab);
			//////mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + two* mfbba * mfabb);
			//////mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + two* mfbab * mfabb); 
			//////			   
			//////mfcca = CUMcca + (mfcaa * mfaca + two* mfbba * mfbba) + c1o3 * (mfcaa + mfaca);
			//////mfcac = CUMcac + (mfcaa * mfaac + two* mfbab * mfbab) + c1o3 * (mfcaa + mfaac);
			//////mfacc = CUMacc + (mfaac * mfaca + two* mfabb * mfabb) + c1o3 * (mfaac + mfaca);

			////////5.
			//////mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + four* mfabb * mfbbb + two* (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac);
			//////mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + four* mfbab * mfbbb + two* (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc);
			//////mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + four* mfbba * mfbbb + two* (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab);
			//////
			////////6.
			//////mfccc = CUMccc  -((-four*  mfbbb * mfbbb  
			//////				-       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
			//////				-  four* (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
			//////				-  two* (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
			//////				+( four* (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
			//////				+  two* (mfcaa * mfaca * mfaac)
			//////				+ sixteen*  mfbba * mfbab * mfabb)
			//////				- c1o3* (mfacc + mfcac + mfcca)
			//////				+ c1o9* (mfcaa + mfaca + mfaac)
			//////				+( two* (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
			//////				+       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3);
			//////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = one;
			real OxyyPxzz  = one;//omega;//two-omega;//
			real OxyyMxzz  = two-omega;//one;// 
			real O4        = one;
			real O5        = one;
			real O6        = one;

			////Cum 4.
			//real CUMcbb;	real CUMbcb;	real CUMbbc;
			//real CUMcca;	real CUMcac;	real CUMacc;
			////Cum 5.
			//real CUMbcc;	real CUMcbc;	real CUMccb;
			////Cum 6.
			//real CUMccc;

			//Cum four
			real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab);
			real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb);
			real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb);

			real CUMcca = mfcca - (mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho;
			real CUMcac = mfcac - (mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-one)*oMdrho;
			real CUMacc = mfacc - (mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho;

			//Cum 5.
			//real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) //O(eps^5) 
			//				- c1o3 * (mfbca + mfbac); //O(eps^3)
			//real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) //O(eps^5)
			//				- c1o3 * (mfcba + mfabc); //O(eps^3)
			//real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) //O(eps^5)
			//				- c1o3 * (mfacb + mfcab);//O(eps^3)

			real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
			real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
			real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

			//Cum 6.
			real CUMccc = mfccc  +((-four *  mfbbb * mfbbb  
				-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
				+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
				+     two * (mfcaa * mfaca * mfaac)
				+ sixteen *  mfbba * mfbab * mfabb)
				-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
				-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(one-two* oMdrho)- c1o27* oMdrho * oMdrho*(-two* oMdrho)
				+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
				+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) +c1o27*oMdrho;





			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;

			{
				real dxux = zero;//c1o2 * (-omega) *(mxxMyy + mxxMzz) +  OxxPyyPzz * (mfaaa - mxxPyyPzz);
				real dyuy = zero;//dxux + omega * c3o2 * mxxMyy;
				real dzuz = zero;//dxux + omega * c3o2 * mxxMzz;

				//relax
				mxxPyyPzz += (OxxPyyPzz)*(mfaaa-two*(ux*vvx+uy*vvy)-ux*ux-uy*uy  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
				mxxMyy    += omega * (-two*(ux*vvx-uy*vvy)-(ux*ux-uy*uy)-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vy2 * dyuy);
				mxxMzz    += omega * (-two*ux*vvx -ux*ux                -mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vz2 * dzuz);
			}
			mfabb     += omega * (-uy*vvz-mfabb);
			mfbab     += omega * (-ux*vvz-mfbab);
			mfbba     += omega * (-ux*vvy-uy*vvx-ux*uy-mfbba);

			// linear combinations back
			mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-two*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - two* mxxMzz + mxxPyyPzz);

			//3.
			// linear combinations

			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

			//relax
			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
			mfbbb     += wadjust * (-mfbbb);
			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
			mxxyPyzz  += wadjust * (-mxxyPyzz);
			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
			mxxyMyzz  += wadjust * (-mxxyMyzz);
			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
			mxxzPyyz  += wadjust * (-mxxzPyyz);
			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
			mxxzMyyz  += wadjust * (-mxxzMyyz);
			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
			mxyyPxzz  += wadjust * (-mxyyPxzz);
			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
			mxyyMxzz  += wadjust * (-mxyyMxzz);

			//// linear combinations back
			////generic
			//mfcba =  zero;
			//mfabc =  zero;
			//mfcab =  zero;
			//mfacb =  zero;
			//mfbca =  zero;
			//mfbac =  zero;

			mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			//CUMacc =  zero; 
			//CUMcac =  zero; 
			//CUMcca =  zero; 
			//	   
			//CUMbbc =  zero; 
			//CUMbcb =  zero; 
			//CUMcbb =  zero; 
			//	  
			////5.   
			//CUMbcc =  zero;
			//CUMcbc =  zero;
			//CUMccb =  zero;
			//	   
			////6.   
			//CUMccc =  zero;

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
			mfcbb = CUMcbb + ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab);
			mfbcb = CUMbcb + ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb);
			mfbbc = CUMbbc + ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb); 

			mfcca = CUMcca + (mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho;
			mfcac = CUMcac + (mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-one)*oMdrho;
			mfacc = CUMacc + (mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho;

			//5.
			mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) * oMdrho;
			mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) * oMdrho;
			mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) * oMdrho;

			//6.
			mfccc = CUMccc  -(( -four *  mfbbb * mfbbb  
				-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
				+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
				+     two * (mfcaa * mfaca * mfaac)
				+ sixteen *  mfbba * mfbab * mfabb)
				-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
				-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(one-two* oMdrho)- c1o27* oMdrho * oMdrho*(-two* oMdrho)
				+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
				+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) -c1o27*oMdrho;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			//back

			//turning ship:

			vvx=vxNeu;
			vvy=vyNeu;

			////////////////////////////////////////////////////////////////////////////////////
			//mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + one* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - two* mfaab *  vvz         +  mfaaa                * (one- vz2)              - one* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + one* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - two* mfabb *  vvz         + mfaba * (one- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - two* mfacb *  vvz         +  mfaca                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - two* mfbab *  vvz         + mfbaa * (one- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - two* mfbbb *  vvz         + mfbba * (one- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - two* mfbcb *  vvz         + mfbca * (one- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - two* mfcab *  vvz         +  mfcaa                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - two* mfcbb *  vvz         + mfcba * (one- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - two* mfccb *  vvz         +  mfcca                  * (one- vz2)              - c1o9 * oMdrho * vz2; 
			m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcca = m0;
			mfccb = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfaca        - two* mfaba *  vvy         +  mfaaa                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - two* mfabb *  vvy         +  mfaab                  * (one- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - two* mfabc *  vvy         +  mfaac                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - two* mfbba *  vvy         + mfbaa * (one- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - two* mfbbb *  vvy         + mfbab * (one- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - two* mfbbc *  vvy         + mfbac * (one- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - two* mfcba *  vvy         +  mfcaa                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - two* mfcbb *  vvy         +  mfcab                  * (one- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - two* mfcbc *  vvy         +  mfcac                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcac = m0;
			mfcbc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcaa        - two* mfbaa *  vvx         +  mfaaa                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - two* mfbba *  vvx         +  mfaba                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - two* mfbca *  vvx         +  mfaca                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - two* mfbab *  vvx         +  mfaab                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - two* mfbbb *  vvx         +  mfabb                  * (one- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - two* mfbcb *  vvx         +  mfacb                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - two* mfbac *  vvx         +  mfaac                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - two* mfbbc *  vvx         +  mfabc                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - two* mfbcc *  vvx         +  mfacc                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dirE   ])[k   ] = mfabb;//(D.f[ dirE   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dirE   ])[k   ]                                                                     
			(D.f[ dirW   ])[kw  ] = mfcbb;//(D.f[ dirW   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dirW   ])[kw  ]                                                                   
			(D.f[ dirN   ])[k   ] = mfbab;//(D.f[ dirN   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ dirN   ])[k   ]
			(D.f[ dirS   ])[ks  ] = mfbcb;//(D.f[ dirS   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ dirS   ])[ks  ]
			(D.f[ dirT   ])[k   ] = mfbba;//(D.f[ dirT   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ dirT   ])[k   ]
			(D.f[ dirB   ])[kb  ] = mfbbc;//(D.f[ dirB   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ dirB   ])[kb  ]
			(D.f[ dirNE  ])[k   ] = mfaab;//(D.f[ dirNE  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dirNE  ])[k   ]
			(D.f[ dirSW  ])[ksw ] = mfccb;//(D.f[ dirSW  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dirSW  ])[ksw ]
			(D.f[ dirSE  ])[ks  ] = mfacb;//(D.f[ dirSE  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dirSE  ])[ks  ]
			(D.f[ dirNW  ])[kw  ] = mfcab;//(D.f[ dirNW  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dirNW  ])[kw  ]
			(D.f[ dirTE  ])[k   ] = mfaba;//(D.f[ dirTE  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dirTE  ])[k   ]
			(D.f[ dirBW  ])[kbw ] = mfcbc;//(D.f[ dirBW  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dirBW  ])[kbw ]
			(D.f[ dirBE  ])[kb  ] = mfabc;//(D.f[ dirBE  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dirBE  ])[kb  ]
			(D.f[ dirTW  ])[kw  ] = mfcba;//(D.f[ dirTW  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dirTW  ])[kw  ]
			(D.f[ dirTN  ])[k   ] = mfbaa;//(D.f[ dirTN  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ dirTN  ])[k   ]
			(D.f[ dirBS  ])[kbs ] = mfbcc;//(D.f[ dirBS  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ dirBS  ])[kbs ]
			(D.f[ dirBN  ])[kb  ] = mfbac;//(D.f[ dirBN  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ dirBN  ])[kb  ]
			(D.f[ dirTS  ])[ks  ] = mfbca;//(D.f[ dirTS  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ dirTS  ])[ks  ]
			(D.f[ dirZERO])[k   ] = mfbbb;//(D.f[ dirZERO])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ dirZERO])[k   ]
			(D.f[ dirTNE ])[k   ] = mfaaa;//(D.f[ dirTNE ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dirTNE ])[k   ]
			(D.f[ dirTSE ])[ks  ] = mfaca;//(D.f[ dirTSE ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dirTSE ])[ks  ]
			(D.f[ dirBNE ])[kb  ] = mfaac;//(D.f[ dirBNE ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dirBNE ])[kb  ]
			(D.f[ dirBSE ])[kbs ] = mfacc;//(D.f[ dirBSE ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dirBSE ])[kbs ]
			(D.f[ dirTNW ])[kw  ] = mfcaa;//(D.f[ dirTNW ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dirTNW ])[kw  ]
			(D.f[ dirTSW ])[ksw ] = mfcca;//(D.f[ dirTSW ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dirTSW ])[ksw ]
			(D.f[ dirBNW ])[kbw ] = mfcac;//(D.f[ dirBNW ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dirBNW ])[kbw ]
			(D.f[ dirBSW ])[kbsw] = mfccc;//(D.f[ dirBSW ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dirBSW ])[kbsw]
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Kum_New_SP_27(     real omega,
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

	if(k<size_Mat)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			//unsigned int kzero= k;
			//unsigned int ke   = k;
			unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			unsigned int kbs  = neighborZ[ks];
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

			//unsigned int kzero= k;
			//unsigned int ke   = k;
			//unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			//unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			//unsigned int kb   = neighborZ[k];
			//unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			//unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			//unsigned int kbs  = neighborZ[ks];
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
			//unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[dirE   ])[k  ];//[ke   ];// +  c2over27 ;(D.f[dirE   ])[k  ];//ke
			real mfabb = (D.f[dirW   ])[kw ];//[kw   ];// +  c2over27 ;(D.f[dirW   ])[kw ];
			real mfbcb = (D.f[dirN   ])[k  ];//[kn   ];// +  c2over27 ;(D.f[dirN   ])[k  ];//kn
			real mfbab = (D.f[dirS   ])[ks ];//[ks   ];// +  c2over27 ;(D.f[dirS   ])[ks ];
			real mfbbc = (D.f[dirT   ])[k  ];//[kt   ];// +  c2over27 ;(D.f[dirT   ])[k  ];//kt
			real mfbba = (D.f[dirB   ])[kb ];//[kb   ];// +  c2over27 ;(D.f[dirB   ])[kb ];
			real mfccb = (D.f[dirNE  ])[k  ];//[kne  ];// +  c1over54 ;(D.f[dirNE  ])[k  ];//kne
			real mfaab = (D.f[dirSW  ])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dirSW  ])[ksw];
			real mfcab = (D.f[dirSE  ])[ks ];//[kse  ];// +  c1over54 ;(D.f[dirSE  ])[ks ];//kse
			real mfacb = (D.f[dirNW  ])[kw ];//[knw  ];// +  c1over54 ;(D.f[dirNW  ])[kw ];//knw
			real mfcbc = (D.f[dirTE  ])[k  ];//[kte  ];// +  c1over54 ;(D.f[dirTE  ])[k  ];//kte
			real mfaba = (D.f[dirBW  ])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dirBW  ])[kbw];
			real mfcba = (D.f[dirBE  ])[kb ];//[kbe  ];// +  c1over54 ;(D.f[dirBE  ])[kb ];//kbe
			real mfabc = (D.f[dirTW  ])[kw ];//[ktw  ];// +  c1over54 ;(D.f[dirTW  ])[kw ];//ktw
			real mfbcc = (D.f[dirTN  ])[k  ];//[ktn  ];// +  c1over54 ;(D.f[dirTN  ])[k  ];//ktn
			real mfbaa = (D.f[dirBS  ])[kbs];//[kbs  ];// +  c1over54 ;(D.f[dirBS  ])[kbs];
			real mfbca = (D.f[dirBN  ])[kb ];//[kbn  ];// +  c1over54 ;(D.f[dirBN  ])[kb ];//kbn
			real mfbac = (D.f[dirTS  ])[ks ];//[kts  ];// +  c1over54 ;(D.f[dirTS  ])[ks ];//kts
			real mfbbb = (D.f[dirZERO])[k  ];//[kzero];// +  c8over27 ;(D.f[dirZERO])[k  ];//kzero
			real mfccc = (D.f[dirTNE ])[k  ];//[ktne ];// +  c1over216;(D.f[dirTNE ])[k  ];//ktne
			real mfaac = (D.f[dirTSW ])[ksw];//[ktsw ];// +  c1over216;(D.f[dirTSW ])[ksw];//ktsw
			real mfcac = (D.f[dirTSE ])[ks ];//[ktse ];// +  c1over216;(D.f[dirTSE ])[ks ];//ktse
			real mfacc = (D.f[dirTNW ])[kw ];//[ktnw ];// +  c1over216;(D.f[dirTNW ])[kw ];//ktnw
			real mfcca = (D.f[dirBNE ])[kb ];//[kbne ];// +  c1over216;(D.f[dirBNE ])[kb ];//kbne
			real mfaaa = (D.f[dirBSW ])[kbsw];//[kbsw ];// +  c1over216;(D.f[dirBSW ])[kbsw];
			real mfcaa = (D.f[dirBSE ])[kbs];//[kbse ];// +  c1over216;(D.f[dirBSE ])[kbs];//kbse
			real mfaca = (D.f[dirBNW ])[kbw];//[kbnw ];// +  c1over216;(D.f[dirBNW ])[kbw];//kbnw
			////////////////////////////////////////////////////////////////////////////////////
			//slow
			//real oMdrho = one - ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
			//					   (((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
			//						((mfabb+mfcbb) + (mfbab+mfbcb)  +  (mfbba+mfbbc)));//fehlt mfbbb
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb));
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab));
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba));
			////////////////////////////////////////////////////////////////////////////////////
			//fast
			real oMdrho = one - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
								   mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
								   mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb);//fehlt mfbbb nicht mehr
			//real vvx    =mfccc-mfaaa + mfcac-mfaca + mfcaa-mfacc + mfcca-mfaac + 
			//				mfcba-mfabc + mfcbc-mfaba + mfcab-mfacb + mfccb-mfaab +
			//				mfcbb-mfabb;
			//real vvy    =mfccc-mfaaa + mfaca-mfcac + mfacc-mfcaa + mfcca-mfaac + 
			//				mfbca-mfbac + mfbcc-mfbaa + mfacb-mfcab + mfccb-mfaab +
			//				mfbcb-mfbab;
			//real vvz    =mfccc-mfaaa + mfcac-mfaca + mfacc-mfcaa + mfaac-mfcca + 
			//				mfbac-mfbca + mfbcc-mfbaa + mfabc-mfcba + mfcbc-mfaba +
			//				mfbbc-mfbba;
			////////////////////////////////////////////////////////////////////////////////////
			// oMdrho assembler style -------> faaaaaastaaaa
			// or much sloooowaaaa ... it depändssssss on sadaku
			real m0, m1, m2;	
			//real oMdrho;
			//{
			//	oMdrho=mfccc+mfaaa;
			//	m0=mfaca+mfcac;
			//	m1=mfacc+mfcaa;
			//	m2=mfaac+mfcca;
			//	oMdrho+=m0;
			//	m1+=m2;
			//	oMdrho+=m1;
			//	m0=mfbac+mfbca;
			//	m1=mfbaa+mfbcc;
			//	m0+=m1;
			//	m1=mfabc+mfcba;
			//	m2=mfaba+mfcbc;
			//	m1+=m2;
			//	m0+=m1;
			//	m1=mfacb+mfcab;
			//	m2=mfaab+mfccb;
			//	m1+=m2;
			//	m0+=m1;
			//	oMdrho+=m0;
			//	m0=mfabb+mfcbb;
			//	m1=mfbab+mfbcb;
			//	m2=mfbba+mfbbc;
			//	m0+=m1+m2;
			//	m0+=mfbbb; //hat gefehlt
			//	oMdrho = one - (oMdrho + m0);
			//}
			//real vvx;
			real vx2;
			//{
			//	vvx = mfccc-mfaaa;
			//	m0  = mfcac-mfaca;
			//	m1  = mfcaa-mfacc;
			//	m2  = mfcca-mfaac;
			//	vvx+= m0;
			//	m1 += m2;
			//	vvx+= m1;
			//	vx2 = mfcba-mfabc;
			//	m0  = mfcbc-mfaba;
			//	m1  = mfcab-mfacb;
			//	m2  = mfccb-mfaab;
			//	vx2+= m0;
			//	m1 += m2;
			//	vx2+= m1;
			//	vvx+= vx2;
			//	vx2 = mfcbb-mfabb;
			//	vvx+= vx2;
			//}
			//real vvy;
			real vy2;
			//{
			//	vvy = mfccc-mfaaa;
			//	m0  = mfaca-mfcac;
			//	m1  = mfacc-mfcaa;
			//	m2  = mfcca-mfaac;
			//	vvy+= m0;
			//	m1 += m2;
			//	vvy+= m1;
			//	vy2 = mfbca-mfbac;
			//	m0  = mfbcc-mfbaa;
			//	m1  = mfacb-mfcab;
			//	m2  = mfccb-mfaab;
			//	vy2+= m0;
			//	m1 += m2;
			//	vy2+= m1;
			//	vvy+= vy2;
			//	vy2 = mfbcb-mfbab;
			//	vvy+= vy2;
			//}
			//real vvz;
			real vz2;
			//{
			//	vvz = mfccc-mfaaa;
			//	m0  = mfcac-mfaca;
			//	m1  = mfacc-mfcaa;
			//	m2  = mfaac-mfcca;
			//	vvz+= m0;
			//	m1 += m2;
			//	vvz+= m1;
			//	vz2 = mfbac-mfbca;
			//	m0  = mfbcc-mfbaa;
			//	m1  = mfabc-mfcba;
			//	m2  = mfcbc-mfaba;
			//	vz2+= m0;
			//	m1 += m2;
			//	vz2+= m1;
			//	vvz+= vz2;
			//	vz2 = mfbbc-mfbba;
			//	vvz+= vz2;
			//}
			vx2=vvx*vvx;
			vy2=vvy*vvy;
			vz2=vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real wadjust;
			real qudricLimit = 0.01f;
			//real s9 = minusomega;
			//test
			//s9 = 0.;
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += one* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// BGK
			////////////////////////////////////////////////////////////////////////////////////
			////2.
			//mfabb += omega * (-mfabb);
			//mfbab += omega * (-mfbab);
			//mfbba += omega * (-mfbba);
			//
			//mfcaa += omega * (c1o3 * mfaaa - mfcaa);
			//mfaca += omega * (c1o3 * mfaaa - mfaca);
			//mfaac += omega * (c1o3 * mfaaa - mfaac);
			//
			////3.
			//mfabc += omega * (-mfabc);
			//mfbac += omega * (-mfbac);
			//
			//mfacb += omega * (-mfacb);
			//mfbca += omega * (-mfbca);

			//mfcab += omega * (-mfcab);
			//mfcba += omega * (-mfcba);

			//mfbbb += omega * (-mfbbb);

			////4.
			//mfacc += omega * (c1o9 * mfaaa - mfacc);
			//mfcac += omega * (c1o9 * mfaaa - mfcac);
			//mfcca += omega * (c1o9 * mfaaa - mfcca);

			//mfbbc += omega * (-mfbbc);
			//mfbcb += omega * (-mfbcb);
			//mfcbb += omega * (-mfcbb);

			////5.
			//mfbcc += omega * (-mfbcc);
			//mfcbc += omega * (-mfcbc);
			//mfccb += omega * (-mfccb);

			////6.
			//mfccc += omega * (c1o27 * mfaaa - mfccc);
			//////////////////////////////////////////////////////////////////////////////////////



			//////////////////////////////////////////////////////////////////////////////////////////
			//////// Cumulants
			//////////////////////////////////////////////////////////////////////////////////////////
			//////real OxxPyyPzz = one;
			//////real OxyyPxzz  = one;//two-omega;//
			//////real OxyyMxzz  = one;//two-omega;//
			//////real O4        = one;
			//////real O5        = one;
			//////real O6        = one;

			////////Cum 4.
			//////real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + two* mfbba * mfbab);
			//////real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + two* mfbba * mfabb);
			//////real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + two* mfbab * mfabb); 

			//////real CUMcca = mfcca - (mfcaa * mfaca + two* mfbba * mfbba)- c1o3 * (mfcaa + mfaca);
			//////real CUMcac = mfcac - (mfcaa * mfaac + two* mfbab * mfbab)- c1o3 * (mfcaa + mfaac);
			//////real CUMacc = mfacc - (mfaac * mfaca + two* mfabb * mfabb)- c1o3 * (mfaac + mfaca);

			////////Cum 5.
			//////real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + four* mfabb * mfbbb + two* (mfbab * mfacb + mfbba * mfabc)) //O(eps^5) 
			//////				- c1o3 * (mfbca + mfbac); //O(eps^3)
			//////real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + four* mfbab * mfbbb + two* (mfabb * mfcab + mfbba * mfbac)) //O(eps^5)
			//////				- c1o3 * (mfcba + mfabc); //O(eps^3)
			//////real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + four* mfbba * mfbbb + two* (mfbab * mfbca + mfabb * mfcba)) //O(eps^5)
			//////				- c1o3 * (mfacb + mfcab);//O(eps^3)

			////////Cum 6.
			//////real CUMccc = mfccc +(-four*  mfbbb * mfbbb  //O(eps^6)
			//////	                   -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca) // O(eps^4)
			//////					   -  four* (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc) // O(eps^6) 
			//////					   -  two* (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) // O(esp^6)
			//////					   +( four* (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac) //O(eps^6)
			//////					   +  two* (mfcaa * mfaca * mfaac) //O(eps^6)
			//////					   + sixteen*  mfbba * mfbab * mfabb) //O(eps^6)
			//////					   - c1o3* (mfacc + mfcac + mfcca) //O(eps^2)
			//////					   + c1o9* (mfcaa + mfaca + mfaac) //O(eps^2)
			//////					   +( two* (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)//O(eps^4)
			//////					   +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3;//O(eps^4)


			////////2.
			//////// linear combinations
			//////real mxxPyyPzz = mfcaa + mfaca + mfaac;
			//////real mxxMyy    = mfcaa - mfaca;
			//////real mxxMzz	   = mfcaa - mfaac;

			//////{
			//////	real dxux = c1o2 * ((-omega) * (mxxMyy + mxxMzz) + (mfaaa - mxxPyyPzz));
			//////	real dyuy = dxux + omega * c3o2 * mxxMyy;
			//////	real dzuz = dxux + omega * c3o2 * mxxMzz;

			//////	//relax
			//////	mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three* (one- c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
			//////	mxxMyy    += omega * (-mxxMyy) - three* (one+ c1o2 * (-omega)) * (vx2 * dxux + vy2 * dyuy);
			//////	mxxMzz    += omega * (-mxxMzz) - three* (one+ c1o2 * (-omega)) * (vx2 * dxux + vz2 * dzuz);
			//////}
			//////mfabb     += omega * (-mfabb);
			//////mfbab     += omega * (-mfbab);
			//////mfbba     += omega * (-mfbba);

			//////// linear combinations back
			//////mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			//////mfaca = c1o3 * (-two*  mxxMyy +      mxxMzz + mxxPyyPzz);
			//////mfaac = c1o3 * (       mxxMyy - two* mxxMzz + mxxPyyPzz);

			////////3.
			//////// linear combinations
			//////real mxxyPyzz = mfcba + mfabc;
			//////real mxxyMyzz = mfcba - mfabc;

			//////real mxxzPyyz = mfcab + mfacb;
			//////real mxxzMyyz = mfcab - mfacb;

			//////real mxyyPxzz = mfbca + mfbac;
			//////real mxyyMxzz = mfbca - mfbac;

			////////relax
			//////wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
			//////mfbbb     += wadjust * (-mfbbb);
			//////wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
			//////mxxyPyzz  += wadjust * (-mxxyPyzz);
			//////wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
			//////mxxyMyzz  += wadjust * (-mxxyMyzz);
			//////wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
			//////mxxzPyyz  += wadjust * (-mxxzPyyz);
			//////wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
			//////mxxzMyyz  += wadjust * (-mxxzMyyz);
			//////wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
			//////mxyyPxzz  += wadjust * (-mxyyPxzz);
			//////wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
			//////mxyyMxzz  += wadjust * (-mxyyMxzz);

			//////// linear combinations back
			//////mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
			//////mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			//////mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
			//////mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			//////mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
			//////mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			////////4.
			//////CUMacc += O4 * (-CUMacc); 
			//////CUMcac += O4 * (-CUMcac); 
			//////CUMcca += O4 * (-CUMcca); 
			//////
			//////CUMbbc += O4 * (-CUMbbc); 
			//////CUMbcb += O4 * (-CUMbcb); 
			//////CUMcbb += O4 * (-CUMcbb); 
			//////		
			////////5.
			//////CUMbcc += O5 * (-CUMbcc);
			//////CUMcbc += O5 * (-CUMcbc);
			//////CUMccb += O5 * (-CUMccb);

			////////6.
			//////CUMccc += O6 * (-CUMccc);
			//////
			////////back cumulants to central moments
			////////4.
			//////mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + two* mfbba * mfbab);
			//////mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + two* mfbba * mfabb);
			//////mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + two* mfbab * mfabb); 
			//////			   
			//////mfcca = CUMcca + (mfcaa * mfaca + two* mfbba * mfbba) + c1o3 * (mfcaa + mfaca);
			//////mfcac = CUMcac + (mfcaa * mfaac + two* mfbab * mfbab) + c1o3 * (mfcaa + mfaac);
			//////mfacc = CUMacc + (mfaac * mfaca + two* mfabb * mfabb) + c1o3 * (mfaac + mfaca);

			////////5.
			//////mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + four* mfabb * mfbbb + two* (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac);
			//////mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + four* mfbab * mfbbb + two* (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc);
			//////mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + four* mfbba * mfbbb + two* (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab);
			//////
			////////6.
			//////mfccc = CUMccc  -((-four*  mfbbb * mfbbb  
			//////				-       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
			//////				-  four* (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
			//////				-  two* (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
			//////				+( four* (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
			//////				+  two* (mfcaa * mfaca * mfaac)
			//////				+ sixteen*  mfbba * mfbab * mfabb)
			//////				- c1o3* (mfacc + mfcac + mfcca)
			//////				+ c1o9* (mfcaa + mfaca + mfaac)
			//////				+( two* (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
			//////				+       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3);
			//////////////////////////////////////////////////////////////////////////////////////////



			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = one;
			real OxyyPxzz  = eight*(two-omega)/(eight -omega);//one;//two-omega;//eight*(two-omega)/(eight -omega);//one;//omega;//two-omega;//
			real OxyyMxzz  = eight*(two-omega)/(eight -omega);//one;//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
			real O4        = one;
			real O5        = one;
			real O6        = one;

			////Cum 4.
			//real CUMcbb;	real CUMbcb;	real CUMbbc;
			//real CUMcca;	real CUMcac;	real CUMacc;
			////Cum 5.
			//real CUMbcc;	real CUMcbc;	real CUMccb;
			////Cum 6.
			//real CUMccc;

			//Cum 4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab); // /rho
			real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb); // /rho
			real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb); // /rho

			real CUMcca = mfcca - ((mfcaa * mfaca + two * mfbba * mfbba) /* /rho*/ + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho);
			real CUMcac = mfcac - ((mfcaa * mfaac + two * mfbab * mfbab) /* /rho*/ + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-one)*oMdrho);
			real CUMacc = mfacc - ((mfaac * mfaca + two * mfabb * mfabb) /* /rho*/ + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho);

			//Cum 5.
			//real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) //O(eps^5) 
			//				- c1o3 * (mfbca + mfbac); //O(eps^3)
			//real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) //O(eps^5)
			//				- c1o3 * (mfcba + mfabc); //O(eps^3)
			//real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) //O(eps^5)
			//				- c1o3 * (mfacb + mfcab);//O(eps^3)

			real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
			real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
			real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

			//Cum 6.
			real CUMccc = mfccc  +((-four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb)
							-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
							-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(one-two* oMdrho)- c1o27* oMdrho * oMdrho*(-two* oMdrho)
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) +c1o27*oMdrho;





			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;
			
			//////////////////////////////////////////////////////////////////////////
			real magicBulk=(CUMacc+CUMcac+CUMcca)*(one/OxxPyyPzz-c1o2)*c3o2*8.;

			//////////////////////////////////////////////////////////////////////////
			//limiter-Scheise Teil 1
			//real oxxyy,oxxzz,oxy,oxz,oyz;
			//real smag=0.001;
			//oxxyy    = omega+(one-omega)*abs(mxxMyy)/(abs(mxxMyy)+smag);
			//oxxzz    = omega+(one-omega)*abs(mxxMzz)/(abs(mxxMzz)+smag);
			//oxy      = omega+(one-omega)*abs(mfbba)/(abs(mfbba)+smag);
			//oxz      = omega+(one-omega)*abs(mfbab)/(abs(mfbab)+smag);
			//oyz      = omega+(one-omega)*abs(mfabb)/(abs(mfabb)+smag);

			////////////////////////////////////////////////////////////////////////////
			////Teil 1b
			//real constante = 1000.0;
			//real nuEddi = constante * abs(mxxPyyPzz);
			//real omegaLimit = one / (one / omega + three * nuEddi);

			//{
			//	real dxux = c1o2 * (-omegaLimit) *(mxxMyy + mxxMzz) +  OxxPyyPzz * (mfaaa - mxxPyyPzz);
			//	real dyuy = dxux + omegaLimit * c3o2 * mxxMyy;
			//	real dzuz = dxux + omegaLimit * c3o2 * mxxMzz;

				////relax
				//mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
				//mxxMyy    += omegaLimit * (-mxxMyy) - three * (one + c1o2 * (-omegaLimit)) * (vx2 * dxux + vy2 * dyuy);
				//mxxMzz    += omegaLimit * (-mxxMzz) - three * (one + c1o2 * (-omegaLimit)) * (vx2 * dxux + vz2 * dzuz);

			//}
			//mfabb     += omegaLimit * (-mfabb);
			//mfbab     += omegaLimit * (-mfbab);
			//mfbba     += omegaLimit * (-mfbba);
			////////////////////////////////////////////////////////////////////////////

 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			////incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
 			//{
 			//	real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
 			//	real dyuy = dxux + omega * c3o2 * mxxMyy;
 			//	real dzuz = dxux + omega * c3o2 * mxxMzz;
 
 			//	//relax
 			//	mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 			//	mxxMyy    += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
 			//	mxxMzz    += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
 
 			//	//////////////////////////////////////////////////////////////////////////
 			//	//limiter-Scheise Teil 2
 			//	//mxxMyy    += oxxyy * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vy2 * dyuy);
 			//	//mxxMzz    += oxxzz * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vz2 * dzuz);
 			//	//////////////////////////////////////////////////////////////////////////
 
 			//}
 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//no correction
			mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);//-magicBulk*OxxPyyPzz;
			mxxMyy    += -(-omega) * (-mxxMyy);
			mxxMzz    += -(-omega) * (-mxxMzz);
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb     += omega * (-mfabb);
			mfbab     += omega * (-mfbab);
			mfbba     += omega * (-mfbba);

			//////////////////////////////////////////////////////////////////////////
			//limiter-Scheise Teil 3
			//mfabb     += oyz * (-mfabb);
			//mfbab     += oxz * (-mfbab);
			//mfbba     += oxy * (-mfbba);
			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-two*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - two* mxxMzz + mxxPyyPzz);

			//3.
			// linear combinations

			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

			//relax
// 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
// 			mfbbb     += wadjust * (-mfbbb);
// 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
// 			mxxyPyzz  += wadjust * (-mxxyPyzz);
// 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
// 			mxxyMyzz  += wadjust * (-mxxyMyzz);
// 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
// 			mxxzPyyz  += wadjust * (-mxxzPyyz);
// 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
// 			mxxzMyyz  += wadjust * (-mxxzMyyz);
// 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
// 			mxyyPxzz  += wadjust * (-mxyyPxzz);
// 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
// 			mxyyMxzz  += wadjust * (-mxyyMxzz);
			mfbbb     += OxyyMxzz * (-mfbbb);
			mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
			mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
			mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
			mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
			mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
			mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);

			//// linear combinations back
			////generic
			//mfcba =  zero;
			//mfabc =  zero;
			//mfcab =  zero;
			//mfacb =  zero;
			//mfbca =  zero;
			//mfbac =  zero;

			mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			//CUMacc =  zero; 
			//CUMcac =  zero; 
			//CUMcca =  zero; 
			//	   
			//CUMbbc =  zero; 
			//CUMbcb =  zero; 
			//CUMcbb =  zero; 
			//	  
			////5.   
			//CUMbcc =  zero;
			//CUMcbc =  zero;
			//CUMccb =  zero;
			//	   
			////6.   
			//CUMccc =  zero;

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
			mfcbb = CUMcbb + ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab);
			mfbcb = CUMbcb + ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb);
			mfbbc = CUMbbc + ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb); 
						   
			mfcca = CUMcca + (mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho;
			mfcac = CUMcac + (mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-one)*oMdrho;
			mfacc = CUMacc + (mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho;

			//5.
			mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) * oMdrho;
			mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) * oMdrho;
			mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) * oMdrho;
			
			//6.
			mfccc = CUMccc  -(( -four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb)
							-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
							-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(one-two* oMdrho)- c1o27* oMdrho * oMdrho*(-two* oMdrho)
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) -c1o27*oMdrho;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			//back
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + one* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - two* mfaab *  vvz         +  mfaaa                * (one- vz2)              - one* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + one* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - two* mfabb *  vvz         + mfaba * (one- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - two* mfacb *  vvz         +  mfaca                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - two* mfbab *  vvz         + mfbaa * (one- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - two* mfbbb *  vvz         + mfbba * (one- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - two* mfbcb *  vvz         + mfbca * (one- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - two* mfcab *  vvz         +  mfcaa                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - two* mfcbb *  vvz         + mfcba * (one- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - two* mfccb *  vvz         +  mfcca                  * (one- vz2)              - c1o9 * oMdrho * vz2; 
			m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcca = m0;
			mfccb = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfaca        - two* mfaba *  vvy         +  mfaaa                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - two* mfabb *  vvy         +  mfaab                  * (one- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - two* mfabc *  vvy         +  mfaac                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - two* mfbba *  vvy         + mfbaa * (one- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - two* mfbbb *  vvy         + mfbab * (one- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - two* mfbbc *  vvy         + mfbac * (one- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - two* mfcba *  vvy         +  mfcaa                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - two* mfcbb *  vvy         +  mfcab                  * (one- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - two* mfcbc *  vvy         +  mfcac                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcac = m0;
			mfcbc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcaa        - two* mfbaa *  vvx         +  mfaaa                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - two* mfbba *  vvx         +  mfaba                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - two* mfbca *  vvx         +  mfaca                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - two* mfbab *  vvx         +  mfaab                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - two* mfbbb *  vvx         +  mfabb                  * (one- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - two* mfbcb *  vvx         +  mfacb                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - two* mfbac *  vvx         +  mfaac                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - two* mfbbc *  vvx         +  mfabc                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - two* mfbcc *  vvx         +  mfacc                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dirE   ])[k   ] = mfabb;//(D.f[ dirE   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dirE   ])[k   ]                                                                     
			(D.f[ dirW   ])[kw  ] = mfcbb;//(D.f[ dirW   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dirW   ])[kw  ]                                                                   
			(D.f[ dirN   ])[k   ] = mfbab;//(D.f[ dirN   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ dirN   ])[k   ]
			(D.f[ dirS   ])[ks  ] = mfbcb;//(D.f[ dirS   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ dirS   ])[ks  ]
			(D.f[ dirT   ])[k   ] = mfbba;//(D.f[ dirT   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ dirT   ])[k   ]
			(D.f[ dirB   ])[kb  ] = mfbbc;//(D.f[ dirB   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ dirB   ])[kb  ]
			(D.f[ dirNE  ])[k   ] = mfaab;//(D.f[ dirNE  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dirNE  ])[k   ]
			(D.f[ dirSW  ])[ksw ] = mfccb;//(D.f[ dirSW  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dirSW  ])[ksw ]
			(D.f[ dirSE  ])[ks  ] = mfacb;//(D.f[ dirSE  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dirSE  ])[ks  ]
			(D.f[ dirNW  ])[kw  ] = mfcab;//(D.f[ dirNW  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dirNW  ])[kw  ]
			(D.f[ dirTE  ])[k   ] = mfaba;//(D.f[ dirTE  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dirTE  ])[k   ]
			(D.f[ dirBW  ])[kbw ] = mfcbc;//(D.f[ dirBW  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dirBW  ])[kbw ]
			(D.f[ dirBE  ])[kb  ] = mfabc;//(D.f[ dirBE  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dirBE  ])[kb  ]
			(D.f[ dirTW  ])[kw  ] = mfcba;//(D.f[ dirTW  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dirTW  ])[kw  ]
			(D.f[ dirTN  ])[k   ] = mfbaa;//(D.f[ dirTN  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ dirTN  ])[k   ]
			(D.f[ dirBS  ])[kbs ] = mfbcc;//(D.f[ dirBS  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ dirBS  ])[kbs ]
			(D.f[ dirBN  ])[kb  ] = mfbac;//(D.f[ dirBN  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ dirBN  ])[kb  ]
			(D.f[ dirTS  ])[ks  ] = mfbca;//(D.f[ dirTS  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ dirTS  ])[ks  ]
			(D.f[ dirZERO])[k   ] = mfbbb;//(D.f[ dirZERO])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ dirZERO])[k   ]
			(D.f[ dirTNE ])[k   ] = mfaaa;//(D.f[ dirTNE ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dirTNE ])[k   ]
			(D.f[ dirTSE ])[ks  ] = mfaca;//(D.f[ dirTSE ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dirTSE ])[ks  ]
			(D.f[ dirBNE ])[kb  ] = mfaac;//(D.f[ dirBNE ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dirBNE ])[kb  ]
			(D.f[ dirBSE ])[kbs ] = mfacc;//(D.f[ dirBSE ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dirBSE ])[kbs ]
			(D.f[ dirTNW ])[kw  ] = mfcaa;//(D.f[ dirTNW ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dirTNW ])[kw  ]
			(D.f[ dirTSW ])[ksw ] = mfcca;//(D.f[ dirTSW ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dirTSW ])[ksw ]
			(D.f[ dirBNW ])[kbw ] = mfcac;//(D.f[ dirBNW ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dirBNW ])[kbw ]
			(D.f[ dirBSW ])[kbsw] = mfccc;//(D.f[ dirBSW ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dirBSW ])[kbsw]
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Kum_New_Comp_SP_27(real omega,
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

	if(k<size_Mat)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( BC >= GEO_FLUID/*(BC != GEO_SOLID) && (BC != GEO_VOID)*/ )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			//unsigned int kzero= k;
			//unsigned int ke   = k;
			unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			unsigned int kbs  = neighborZ[ks];
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

			//unsigned int kzero= k;
			//unsigned int ke   = k;
			//unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			//unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			//unsigned int kb   = neighborZ[k];
			//unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			//unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			//unsigned int kbs  = neighborZ[ks];
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
			//unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[dirE   ])[k  ];//[ke   ];// +  c2over27 ;(D.f[dirE   ])[k  ];//ke
			real mfabb = (D.f[dirW   ])[kw ];//[kw   ];// +  c2over27 ;(D.f[dirW   ])[kw ];
			real mfbcb = (D.f[dirN   ])[k  ];//[kn   ];// +  c2over27 ;(D.f[dirN   ])[k  ];//kn
			real mfbab = (D.f[dirS   ])[ks ];//[ks   ];// +  c2over27 ;(D.f[dirS   ])[ks ];
			real mfbbc = (D.f[dirT   ])[k  ];//[kt   ];// +  c2over27 ;(D.f[dirT   ])[k  ];//kt
			real mfbba = (D.f[dirB   ])[kb ];//[kb   ];// +  c2over27 ;(D.f[dirB   ])[kb ];
			real mfccb = (D.f[dirNE  ])[k  ];//[kne  ];// +  c1over54 ;(D.f[dirNE  ])[k  ];//kne
			real mfaab = (D.f[dirSW  ])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dirSW  ])[ksw];
			real mfcab = (D.f[dirSE  ])[ks ];//[kse  ];// +  c1over54 ;(D.f[dirSE  ])[ks ];//kse
			real mfacb = (D.f[dirNW  ])[kw ];//[knw  ];// +  c1over54 ;(D.f[dirNW  ])[kw ];//knw
			real mfcbc = (D.f[dirTE  ])[k  ];//[kte  ];// +  c1over54 ;(D.f[dirTE  ])[k  ];//kte
			real mfaba = (D.f[dirBW  ])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dirBW  ])[kbw];
			real mfcba = (D.f[dirBE  ])[kb ];//[kbe  ];// +  c1over54 ;(D.f[dirBE  ])[kb ];//kbe
			real mfabc = (D.f[dirTW  ])[kw ];//[ktw  ];// +  c1over54 ;(D.f[dirTW  ])[kw ];//ktw
			real mfbcc = (D.f[dirTN  ])[k  ];//[ktn  ];// +  c1over54 ;(D.f[dirTN  ])[k  ];//ktn
			real mfbaa = (D.f[dirBS  ])[kbs];//[kbs  ];// +  c1over54 ;(D.f[dirBS  ])[kbs];
			real mfbca = (D.f[dirBN  ])[kb ];//[kbn  ];// +  c1over54 ;(D.f[dirBN  ])[kb ];//kbn
			real mfbac = (D.f[dirTS  ])[ks ];//[kts  ];// +  c1over54 ;(D.f[dirTS  ])[ks ];//kts
			real mfbbb = (D.f[dirZERO])[k  ];//[kzero];// +  c8over27 ;(D.f[dirZERO])[k  ];//kzero
			real mfccc = (D.f[dirTNE ])[k  ];//[ktne ];// +  c1over216;(D.f[dirTNE ])[k  ];//ktne
			real mfaac = (D.f[dirTSW ])[ksw];//[ktsw ];// +  c1over216;(D.f[dirTSW ])[ksw];//ktsw
			real mfcac = (D.f[dirTSE ])[ks ];//[ktse ];// +  c1over216;(D.f[dirTSE ])[ks ];//ktse
			real mfacc = (D.f[dirTNW ])[kw ];//[ktnw ];// +  c1over216;(D.f[dirTNW ])[kw ];//ktnw
			real mfcca = (D.f[dirBNE ])[kb ];//[kbne ];// +  c1over216;(D.f[dirBNE ])[kb ];//kbne
			real mfaaa = (D.f[dirBSW ])[kbsw];//[kbsw ];// +  c1over216;(D.f[dirBSW ])[kbsw];
			real mfcaa = (D.f[dirBSE ])[kbs];//[kbse ];// +  c1over216;(D.f[dirBSE ])[kbs];//kbse
			real mfaca = (D.f[dirBNW ])[kbw];//[kbnw ];// +  c1over216;(D.f[dirBNW ])[kbw];//kbnw
			////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
							(((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
							((mfabb+mfcbb) + (mfbab+mfbcb) + (mfbba+mfbbc))) + mfbbb;

			real rho = one+drho;
			////////////////////////////////////////////////////////////////////////////////////
			//slow
			//real oMdrho = one - ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
			//					   (((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
			//						((mfabb+mfcbb) + (mfbab+mfbcb)  +  (mfbba+mfbbc)));//fehlt mfbbb
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) / rho;
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) / rho;
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			//the force be with you
			real fx = forces[0]/(pow(two,level)); //zero;//0.0032653/(pow(two,level)); //0.000000005;//(two/1600000.0) / 120.0; //
			real fy = forces[1]/(pow(two,level)); //zero;
			real fz = forces[2]/(pow(two,level)); //zero;
			vvx += fx*c1o2;
			vvy += fy*c1o2;
			vvz += fz*c1o2;
			////////////////////////////////////////////////////////////////////////////////////
			//real omega = omega_in;
			////////////////////////////////////////////////////////////////////////////////////
			//fast
			real oMdrho = one; // comp special
			//real oMdrho = one - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
			//					   mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
			//					   mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb);//fehlt mfbbb nicht mehr
			//real vvx    =mfccc-mfaaa + mfcac-mfaca + mfcaa-mfacc + mfcca-mfaac + 
			//				mfcba-mfabc + mfcbc-mfaba + mfcab-mfacb + mfccb-mfaab +
			//				mfcbb-mfabb;
			//real vvy    =mfccc-mfaaa + mfaca-mfcac + mfacc-mfcaa + mfcca-mfaac + 
			//				mfbca-mfbac + mfbcc-mfbaa + mfacb-mfcab + mfccb-mfaab +
			//				mfbcb-mfbab;
			//real vvz    =mfccc-mfaaa + mfcac-mfaca + mfacc-mfcaa + mfaac-mfcca + 
			//				mfbac-mfbca + mfbcc-mfbaa + mfabc-mfcba + mfcbc-mfaba +
			//				mfbbc-mfbba;
			////////////////////////////////////////////////////////////////////////////////////
			// oMdrho assembler style -------> faaaaaastaaaa
			// or much sloooowaaaa ... it depändssssss on sadaku
			real m0, m1, m2;	
			//real oMdrho;
			//{
			//	oMdrho=mfccc+mfaaa;
			//	m0=mfaca+mfcac;
			//	m1=mfacc+mfcaa;
			//	m2=mfaac+mfcca;
			//	oMdrho+=m0;
			//	m1+=m2;
			//	oMdrho+=m1;
			//	m0=mfbac+mfbca;
			//	m1=mfbaa+mfbcc;
			//	m0+=m1;
			//	m1=mfabc+mfcba;
			//	m2=mfaba+mfcbc;
			//	m1+=m2;
			//	m0+=m1;
			//	m1=mfacb+mfcab;
			//	m2=mfaab+mfccb;
			//	m1+=m2;
			//	m0+=m1;
			//	oMdrho+=m0;
			//	m0=mfabb+mfcbb;
			//	m1=mfbab+mfbcb;
			//	m2=mfbba+mfbbc;
			//	m0+=m1+m2;
			//	m0+=mfbbb; //hat gefehlt
			//	oMdrho = one - (oMdrho + m0);
			//}
			//real vvx;
			real vx2;
			//{
			//	vvx = mfccc-mfaaa;
			//	m0  = mfcac-mfaca;
			//	m1  = mfcaa-mfacc;
			//	m2  = mfcca-mfaac;
			//	vvx+= m0;
			//	m1 += m2;
			//	vvx+= m1;
			//	vx2 = mfcba-mfabc;
			//	m0  = mfcbc-mfaba;
			//	m1  = mfcab-mfacb;
			//	m2  = mfccb-mfaab;
			//	vx2+= m0;
			//	m1 += m2;
			//	vx2+= m1;
			//	vvx+= vx2;
			//	vx2 = mfcbb-mfabb;
			//	vvx+= vx2;
			//}
			//real vvy;
			real vy2;
			//{
			//	vvy = mfccc-mfaaa;
			//	m0  = mfaca-mfcac;
			//	m1  = mfacc-mfcaa;
			//	m2  = mfcca-mfaac;
			//	vvy+= m0;
			//	m1 += m2;
			//	vvy+= m1;
			//	vy2 = mfbca-mfbac;
			//	m0  = mfbcc-mfbaa;
			//	m1  = mfacb-mfcab;
			//	m2  = mfccb-mfaab;
			//	vy2+= m0;
			//	m1 += m2;
			//	vy2+= m1;
			//	vvy+= vy2;
			//	vy2 = mfbcb-mfbab;
			//	vvy+= vy2;
			//}
			//real vvz;
			real vz2;
			//{
			//	vvz = mfccc-mfaaa;
			//	m0  = mfcac-mfaca;
			//	m1  = mfacc-mfcaa;
			//	m2  = mfaac-mfcca;
			//	vvz+= m0;
			//	m1 += m2;
			//	vvz+= m1;
			//	vz2 = mfbac-mfbca;
			//	m0  = mfbcc-mfbaa;
			//	m1  = mfabc-mfcba;
			//	m2  = mfcbc-mfaba;
			//	vz2+= m0;
			//	m1 += m2;
			//	vz2+= m1;
			//	vvz+= vz2;
			//	vz2 = mfbbc-mfbba;
			//	vvz+= vz2;
			//}
			vx2=vvx*vvx;
			vy2=vvy*vvy;
			vz2=vvz*vvz;
			//////////////////////////////////////////////////////////////////////////////////////
			//// test rundungsfehler....
			//mfcbb = (mfcbb - c2over27* (rho-one))/(rho);
			//mfabb = (mfabb - c2over27* (rho-one))/(rho);
			//mfbcb = (mfbcb - c2over27* (rho-one))/(rho);
			//mfbab = (mfbab - c2over27* (rho-one))/(rho);
			//mfbbc = (mfbbc - c2over27* (rho-one))/(rho);
			//mfbba = (mfbba - c2over27* (rho-one))/(rho);
			//mfccb = (mfccb - c1over54* (rho-one))/(rho);
			//mfaab = (mfaab - c1over54* (rho-one))/(rho);
			//mfcab = (mfcab - c1over54* (rho-one))/(rho);
			//mfacb = (mfacb - c1over54* (rho-one))/(rho);
			//mfcbc = (mfcbc - c1over54* (rho-one))/(rho);
			//mfaba = (mfaba - c1over54* (rho-one))/(rho);
			//mfcba = (mfcba - c1over54* (rho-one))/(rho);
			//mfabc = (mfabc - c1over54* (rho-one))/(rho);
			//mfbcc = (mfbcc - c1over54* (rho-one))/(rho);
			//mfbaa = (mfbaa - c1over54* (rho-one))/(rho);
			//mfbca = (mfbca - c1over54* (rho-one))/(rho);
			//mfbac = (mfbac - c1over54* (rho-one))/(rho);
			//mfbbb = (mfbbb - c8over27* (rho-one))/(rho);
			//mfccc = (mfccc - c1over216*(rho-one))/(rho);
			//mfaac = (mfaac - c1over216*(rho-one))/(rho);
			//mfcac = (mfcac - c1over216*(rho-one))/(rho);
			//mfacc = (mfacc - c1over216*(rho-one))/(rho);
			//mfcca = (mfcca - c1over216*(rho-one))/(rho);
			//mfaaa = (mfaaa - c1over216*(rho-one))/(rho);
			//mfcaa = (mfcaa - c1over216*(rho-one))/(rho);
			//mfaca = (mfaca - c1over216*(rho-one))/(rho);
			////////////////////////////////////////////////////////////////////////////////////
			real wadjust;
			real qudricLimitP = c1o100;// * 0.0001f;
			real qudricLimitM = c1o100;// * 0.0001f;
			real qudricLimitD = c1o100;// * 0.001f;
			//real s9 = minusomega;
			//test
			//s9 = 0.;
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += one* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = one; //omega; // one;	//set the bulk viscosity one is high / two is very low and zero is (too) high

			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			real OxyyPxzz  = one;//three  * (two - omega) / (three  - omega);//one;//(1000.*(-2. + omega))/(-1000. + 439.*omega);//(eight * (omega - two)) / (omega - eight);//two-omega;//eight*(two-omega)/(eight -omega);//one;//omega;//two-omega;//
			real OxyyMxzz  = one;//six    * (two - omega) / (six    - omega);//one;//two-omega;//(1000.*(-2. + omega))/(-1000. + 439.*omega);//(eight * (omega - two)) / (omega - eight);//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
			real Oxyz      = one;//twelve * (two - omega) / (twelve + omega);//one;//two-omega;//(1000.*(-2. + omega))/(-1000. + 439.*omega);//(eight * (omega - two)) / (omega - eight);//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
			//////////////////////////////
			//real OxyyPxzz  = two-omega;//
			//real OxyyMxzz  = two-omega;//
			//////////////////////////////
			//real OxyyPxzz  = (eight * (omega - two)) / (omega - eight);//Ginzburg
			//real OxyyMxzz  = (eight * (omega - two)) / (omega - eight);//Ginzburg
			//////////////////////////////
			//real OxyyPxzz  = omega;//BGK
			//real OxyyMxzz  = omega;//BGK
			//////////////////////////////
			//real OxyyPxzz  = (one + omega) / two;//1P5
			//real OxyyMxzz  = (one + omega) / two;//1P5
			//////////////////////////////
			//real OxyyPxzz  = (three - omega) / two;//0P5
			//real OxyyMxzz  = (three - omega) / two;//0P5
			//////////////////////////////
			//real OxyyPxzz  = (one + (eight * (omega - two)) / (omega - eight)) / two;//one + Ginzburg / two ... Car
			//real OxyyMxzz  = (one + (eight * (omega - two)) / (omega - eight)) / two;//one + Ginzburg / two ... Car
			////////////////////////////////////////////////////////////
			//4.
			//////////////////////////////
			real O4        = one;
			//////////////////////////////
			//real O4        = omega;//TRT
			////////////////////////////////////////////////////////////
			//5.
			//////////////////////////////
			real O5        = one;
			////////////////////////////////////////////////////////////
			//6.
			//////////////////////////////
			real O6        = one;
			////////////////////////////////////////////////////////////


			//central moments to cumulants
			//4.
			//real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab) / rho;  //bis 15.05.2015 verwendet
			//real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb) / rho;  //bis 15.05.2015 verwendet
			//real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb) / rho;  //bis 15.05.2015 verwendet
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho;	//ab 15.05.2015 verwendet
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho; //ab 15.05.2015 verwendet
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho; //ab 15.05.2015 verwendet
				  	 		
			real CUMcca = mfcca - (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			real CUMcac = mfcac - (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			real CUMacc = mfacc - (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			real CUMccc = mfccc + ((-four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));
							//+ c1o27*(one -three/rho +two/(rho*rho)));




			////Cum 4.
			//real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab)/rho; // 
			//real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb)/rho; // 
			//real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb)/rho; // 

			//real CUMcca = mfcca - ((mfcaa * mfaca + two * mfbba * mfbba) /rho + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho);
			//real CUMcac = mfcac - ((mfcaa * mfaac + two * mfbab * mfbab) /rho + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-one)*oMdrho);
			//real CUMacc = mfacc - ((mfaac * mfaca + two * mfabb * mfabb) /rho + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho);

			////Cum 5.
			//real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) / rho - c1o3 * (mfbca + mfbac) * oMdrho;
			//real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) / rho - c1o3 * (mfcba + mfabc) * oMdrho;
			//real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) / rho - c1o3 * (mfacb + mfcab) * oMdrho;

			////Cum 6.
			//real CUMccc = mfccc  +((-four *  mfbbb * mfbbb  
			//				-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
			//				-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
			//				-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
			//				+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
			//				+     two * (mfcaa * mfaca * mfaac)
			//				+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
			//				-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
			//				-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(one-two* oMdrho)- c1o27* oMdrho * oMdrho*(-two* oMdrho)
			//				+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
			//				+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) / rho * c2o3*oMdrho) +c1o27*oMdrho;





			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;
			
			//////////////////////////////////////////////////////////////////////////
// 			real magicBulk=(CUMacc+CUMcac+CUMcca)*(one/OxxPyyPzz-c1o2)*c3o2*8.;

			//////////////////////////////////////////////////////////////////////////
			//limiter-Scheise Teil 1
			//real oxxyy,oxxzz,oxy,oxz,oyz;
			//real smag=0.001;
			//oxxyy    = omega+(one-omega)*abs(mxxMyy)/(abs(mxxMyy)+smag);
			//oxxzz    = omega+(one-omega)*abs(mxxMzz)/(abs(mxxMzz)+smag);
			//oxy      = omega+(one-omega)*abs(mfbba)/(abs(mfbba)+smag);
			//oxz      = omega+(one-omega)*abs(mfbab)/(abs(mfbab)+smag);
			//oyz      = omega+(one-omega)*abs(mfabb)/(abs(mfabb)+smag);

			////////////////////////////////////////////////////////////////////////////
			////Teil 1b
			//real constante = 1000.0;
			//real nuEddi = constante * abs(mxxPyyPzz);
			//real omegaLimit = one / (one / omega + three * nuEddi);

			//{
			//	real dxux = c1o2 * (-omegaLimit) *(mxxMyy + mxxMzz) +  OxxPyyPzz * (mfaaa - mxxPyyPzz);
			//	real dyuy = dxux + omegaLimit * c3o2 * mxxMyy;
			//	real dzuz = dxux + omegaLimit * c3o2 * mxxMzz;

				////relax
				//mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
				//mxxMyy    += omegaLimit * (-mxxMyy) - three * (one + c1o2 * (-omegaLimit)) * (vx2 * dxux + vy2 * dyuy);
				//mxxMzz    += omegaLimit * (-mxxMzz) - three * (one + c1o2 * (-omegaLimit)) * (vx2 * dxux + vz2 * dzuz);

			//}
			//mfabb     += omegaLimit * (-mfabb);
			//mfbab     += omegaLimit * (-mfbab);
			//mfbba     += omegaLimit * (-mfbba);
			////////////////////////////////////////////////////////////////////////////

 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			//incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
 			{
 				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
 				real dyuy = dxux + omega * c3o2 * mxxMyy;
 				real dzuz = dxux + omega * c3o2 * mxxMzz;
 
 				//relax
 				mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 				mxxMyy    += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
 				mxxMzz    += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
 
 				//////////////////////////////////////////////////////////////////////////
 				//limiter-Scheise Teil 2
 				//mxxMyy    += oxxyy * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vy2 * dyuy);
 				//mxxMzz    += oxxzz * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vz2 * dzuz);
 				//////////////////////////////////////////////////////////////////////////
 
 			}
 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			////no correction
 			//mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);//-magicBulk*OxxPyyPzz;
 			//mxxMyy    += -(-omega) * (-mxxMyy);
 			//mxxMzz    += -(-omega) * (-mxxMzz);
 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb     += omega * (-mfabb);
			mfbab     += omega * (-mfbab);
			mfbba     += omega * (-mfbba);

			//////////////////////////////////////////////////////////////////////////
			//limiter-Scheise Teil 3
			//mfabb     += oyz * (-mfabb);
			//mfbab     += oxz * (-mfbab);
			//mfbba     += oxy * (-mfbba);
			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-two*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - two* mxxMzz + mxxPyyPzz);

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
			//das ist der limiter
 			wadjust    = Oxyz+(one-Oxyz)*abs(mfbbb)/(abs(mfbbb)+qudricLimitD);
 			mfbbb     += wadjust * (-mfbbb);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimitP);
 			mxxyPyzz  += wadjust * (-mxxyPyzz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimitM);
 			mxxyMyzz  += wadjust * (-mxxyMyzz);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimitP);
 			mxxzPyyz  += wadjust * (-mxxzPyyz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimitM);
 			mxxzMyyz  += wadjust * (-mxxzMyyz);
 			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimitP);
 			mxyyPxzz  += wadjust * (-mxyyPxzz);
 			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimitM);
 			mxyyMxzz  += wadjust * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
			//mfbbb     += OxyyMxzz * (-mfbbb);
			//mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
			//mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
			//mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
			//mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
			//mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
			//mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////

			//// linear combinations back
			////generic
			//mfcba =  zero;
			//mfabc =  zero;
			//mfcab =  zero;
			//mfacb =  zero;
			//mfbca =  zero;
			//mfbac =  zero;

			mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			//CUMacc =  zero; 
			//CUMcac =  zero; 
			//CUMcca =  zero; 
			//	   
			//CUMbbc =  zero; 
			//CUMbcb =  zero; 
			//CUMcbb =  zero; 
			//	  
			////5.   
			//CUMbcc =  zero;
			//CUMcbc =  zero;
			//CUMccb =  zero;
			//	   
			////6.   
			//CUMccc =  zero;

			//4.
			//////////////////////////////////////////////////////////////////////////
			//mit limiter
 		//	wadjust    = O4+(one-O4)*abs(CUMacc)/(abs(CUMacc)+qudricLimit);
			//CUMacc    += wadjust * (-CUMacc);
 		//	wadjust    = O4+(one-O4)*abs(CUMcac)/(abs(CUMcac)+qudricLimit);
			//CUMcac    += wadjust * (-CUMcac); 
 		//	wadjust    = O4+(one-O4)*abs(CUMcca)/(abs(CUMcca)+qudricLimit);
			//CUMcca    += wadjust * (-CUMcca); 

 		//	wadjust    = O4+(one-O4)*abs(CUMbbc)/(abs(CUMbbc)+qudricLimit);
			//CUMbbc    += wadjust * (-CUMbbc); 
 		//	wadjust    = O4+(one-O4)*abs(CUMbcb)/(abs(CUMbcb)+qudricLimit);
			//CUMbcb    += wadjust * (-CUMbcb); 
 		//	wadjust    = O4+(one-O4)*abs(CUMcbb)/(abs(CUMcbb)+qudricLimit);
			//CUMcbb    += wadjust * (-CUMcbb); 
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
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
			//mfcbb = CUMcbb + ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab) / rho; //bis 15.05.2015 verwendet
			//mfbcb = CUMbcb + ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb) / rho; //bis 15.05.2015 verwendet
			//mfbbc = CUMbbc + ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb) / rho; //bis 15.05.2015 verwendet
			mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho; //ab 15.05.2015 verwendet
			mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho; //ab 15.05.2015 verwendet
			mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho; //ab 15.05.2015 verwendet
						   
			mfcca = CUMcca + (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));//(one/rho-one));
			mfcac = CUMcac + (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));//(one/rho-one));
			mfacc = CUMacc + (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));//(one/rho-one));

			//5.
			mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			mfccc = CUMccc - ((-four *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two * (mfcaa * mfaca * mfaac)
							+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));
							//+ c1o27*(one -three/rho +two/(rho*rho)));


			//mfccc = CUMccc - ((-four *  mfbbb * mfbbb  
			//				-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
			//				-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
			//				-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
			//				+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
			//				+     two * (mfcaa * mfaca * mfaac)
			//				+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
			//				-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
			//				-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(one-two* oMdrho)- c1o27* oMdrho * oMdrho*(-two* oMdrho)
			//				+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
			//				+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) / rho * c2o3*oMdrho) -c1o27*oMdrho;





			//mfccc = CUMccc  -(( -four *  mfbbb * mfbbb  
			//				-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
			//				-    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
			//				-     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
			//				+(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
			//				+     two * (mfcaa * mfaca * mfaac)
			//				+ sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
			//				+    (c1o3 * (mfacc + mfcac + mfcca) )/rho
			//				+	 (c1o9 * (mfcaa + mfaca + mfaac) )/rho
			//				- (four*c1o3*(mfabb * mfabb + mfbab * mfbab + mfbba * mfbba) + two * ( c1o3 *(mfcaa * mfaca + mfcaa * mfaac + mfaca * mfaac) + c1o9* (mfcaa + mfaca + mfaac )))/(rho*rho) 
			//				+ c1o27*(three/rho - two/(rho * rho) - one);
			//				//-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(one-two* oMdrho)- c1o27* oMdrho * oMdrho*(-two* oMdrho) //????
			//				//+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
			//				//+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) / rho * c2o3*oMdrho) -c1o27*oMdrho;

			//Test Cascade
			//4.
			//mfacc = mfaaa * c1o9 * O4 + (one - O4) * mfacc; 
			//mfcac = mfaaa * c1o9 * O4 + (one - O4) * mfcac; 
			//mfcca = mfaaa * c1o9 * O4 + (one - O4) * mfcca; 
			//
			//mfbbc += O4 * (-mfbbc); 
			//mfbcb += O4 * (-mfbcb); 
			//mfcbb += O4 * (-mfcbb); 
					
			////5.
			//mfbcc += O5 * (-mfbcc);
			//mfcbc += O5 * (-mfcbc);
			//mfccb += O5 * (-mfccb);

			//6.
			//mfccc = mfaaa * c1o27 * O6 + (one - O6) * mfccc;

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
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + one* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - two* mfaab *  vvz         +  mfaaa                * (one- vz2)              - one* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + one* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - two* mfabb *  vvz         + mfaba * (one- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - two* mfacb *  vvz         +  mfaca                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - two* mfbab *  vvz         + mfbaa * (one- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - two* mfbbb *  vvz         + mfbba * (one- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - two* mfbcb *  vvz         + mfbca * (one- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - two* mfcab *  vvz         +  mfcaa                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - two* mfcbb *  vvz         + mfcba * (one- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - two* mfccb *  vvz         +  mfcca                  * (one- vz2)              - c1o9 * oMdrho * vz2; 
			m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcca = m0;
			mfccb = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfaca        - two* mfaba *  vvy         +  mfaaa                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - two* mfabb *  vvy         +  mfaab                  * (one- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - two* mfabc *  vvy         +  mfaac                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - two* mfbba *  vvy         + mfbaa * (one- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - two* mfbbb *  vvy         + mfbab * (one- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - two* mfbbc *  vvy         + mfbac * (one- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - two* mfcba *  vvy         +  mfcaa                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - two* mfcbb *  vvy         +  mfcab                  * (one- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - two* mfcbc *  vvy         +  mfcac                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcac = m0;
			mfcbc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcaa        - two* mfbaa *  vvx         +  mfaaa                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - two* mfbba *  vvx         +  mfaba                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - two* mfbca *  vvx         +  mfaca                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - two* mfbab *  vvx         +  mfaab                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - two* mfbbb *  vvx         +  mfabb                  * (one- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - two* mfbcb *  vvx         +  mfacb                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - two* mfbac *  vvx         +  mfaac                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - two* mfbbc *  vvx         +  mfabc                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - two* mfbcc *  vvx         +  mfacc                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////

			//////////////////////////////////////////////////////////////////////////////////////
			//// test rundungsfehler....
			//mfcbb = (mfcbb * rho + c2over27* (rho-one));
			//mfabb = (mfabb * rho + c2over27* (rho-one));
			//mfbcb = (mfbcb * rho + c2over27* (rho-one));
			//mfbab = (mfbab * rho + c2over27* (rho-one));
			//mfbbc = (mfbbc * rho + c2over27* (rho-one));
			//mfbba = (mfbba * rho + c2over27* (rho-one));
			//mfccb = (mfccb * rho + c1over54* (rho-one));
			//mfaab = (mfaab * rho + c1over54* (rho-one));
			//mfcab = (mfcab * rho + c1over54* (rho-one));
			//mfacb = (mfacb * rho + c1over54* (rho-one));
			//mfcbc = (mfcbc * rho + c1over54* (rho-one));
			//mfaba = (mfaba * rho + c1over54* (rho-one));
			//mfcba = (mfcba * rho + c1over54* (rho-one));
			//mfabc = (mfabc * rho + c1over54* (rho-one));
			//mfbcc = (mfbcc * rho + c1over54* (rho-one));
			//mfbaa = (mfbaa * rho + c1over54* (rho-one));
			//mfbca = (mfbca * rho + c1over54* (rho-one));
			//mfbac = (mfbac * rho + c1over54* (rho-one));
			//mfbbb = (mfbbb * rho + c8over27* (rho-one));
			//mfccc = (mfccc * rho + c1over216*(rho-one));
			//mfaac = (mfaac * rho + c1over216*(rho-one));
			//mfcac = (mfcac * rho + c1over216*(rho-one));
			//mfacc = (mfacc * rho + c1over216*(rho-one));
			//mfcca = (mfcca * rho + c1over216*(rho-one));
			//mfaaa = (mfaaa * rho + c1over216*(rho-one));
			//mfcaa = (mfcaa * rho + c1over216*(rho-one));
			//mfaca = (mfaca * rho + c1over216*(rho-one));
			real drhoPost = 
				((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
				((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
			mfbbb += drho - drhoPost;
			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dirE   ])[k   ] = mfabb;//(D.f[ dirE   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dirE   ])[k   ]                                                                     
			(D.f[ dirW   ])[kw  ] = mfcbb;//(D.f[ dirW   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dirW   ])[kw  ]                                                                   
			(D.f[ dirN   ])[k   ] = mfbab;//(D.f[ dirN   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ dirN   ])[k   ]
			(D.f[ dirS   ])[ks  ] = mfbcb;//(D.f[ dirS   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ dirS   ])[ks  ]
			(D.f[ dirT   ])[k   ] = mfbba;//(D.f[ dirT   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ dirT   ])[k   ]
			(D.f[ dirB   ])[kb  ] = mfbbc;//(D.f[ dirB   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ dirB   ])[kb  ]
			(D.f[ dirNE  ])[k   ] = mfaab;//(D.f[ dirNE  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dirNE  ])[k   ]
			(D.f[ dirSW  ])[ksw ] = mfccb;//(D.f[ dirSW  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dirSW  ])[ksw ]
			(D.f[ dirSE  ])[ks  ] = mfacb;//(D.f[ dirSE  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dirSE  ])[ks  ]
			(D.f[ dirNW  ])[kw  ] = mfcab;//(D.f[ dirNW  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dirNW  ])[kw  ]
			(D.f[ dirTE  ])[k   ] = mfaba;//(D.f[ dirTE  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dirTE  ])[k   ]
			(D.f[ dirBW  ])[kbw ] = mfcbc;//(D.f[ dirBW  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dirBW  ])[kbw ]
			(D.f[ dirBE  ])[kb  ] = mfabc;//(D.f[ dirBE  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dirBE  ])[kb  ]
			(D.f[ dirTW  ])[kw  ] = mfcba;//(D.f[ dirTW  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dirTW  ])[kw  ]
			(D.f[ dirTN  ])[k   ] = mfbaa;//(D.f[ dirTN  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ dirTN  ])[k   ]
			(D.f[ dirBS  ])[kbs ] = mfbcc;//(D.f[ dirBS  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ dirBS  ])[kbs ]
			(D.f[ dirBN  ])[kb  ] = mfbac;//(D.f[ dirBN  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ dirBN  ])[kb  ]
			(D.f[ dirTS  ])[ks  ] = mfbca;//(D.f[ dirTS  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ dirTS  ])[ks  ]
			(D.f[ dirZERO])[k   ] = mfbbb;//(D.f[ dirZERO])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ dirZERO])[k   ]
			(D.f[ dirTNE ])[k   ] = mfaaa;//(D.f[ dirTNE ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dirTNE ])[k   ]
			(D.f[ dirTSE ])[ks  ] = mfaca;//(D.f[ dirTSE ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dirTSE ])[ks  ]
			(D.f[ dirBNE ])[kb  ] = mfaac;//(D.f[ dirBNE ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dirBNE ])[kb  ]
			(D.f[ dirBSE ])[kbs ] = mfacc;//(D.f[ dirBSE ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dirBSE ])[kbs ]
			(D.f[ dirTNW ])[kw  ] = mfcaa;//(D.f[ dirTNW ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dirTNW ])[kw  ]
			(D.f[ dirTSW ])[ksw ] = mfcca;//(D.f[ dirTSW ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dirTSW ])[ksw ]
			(D.f[ dirBNW ])[kbw ] = mfcac;//(D.f[ dirBNW ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dirBNW ])[kbw ]
			(D.f[ dirBSW ])[kbsw] = mfccc;//(D.f[ dirBSW ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dirBSW ])[kbsw]
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Kum_Comp_SP_27(    real omega,
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

	if(k<size_Mat)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			unsigned int kzero= k;
			unsigned int ke   = k;
			unsigned int kw   = neighborX[k];
			unsigned int kn   = k;
			unsigned int ks   = neighborY[k];
			unsigned int kt   = k;
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			unsigned int kne  = k;
			unsigned int kse  = ks;
			unsigned int knw  = kw;
			unsigned int kbw  = neighborZ[kw];
			unsigned int kte  = k;
			unsigned int kbe  = kb;
			unsigned int ktw  = kw;
			unsigned int kbs  = neighborZ[ks];
			unsigned int ktn  = k;
			unsigned int kbn  = kb;
			unsigned int kts  = ks;
			unsigned int ktse = ks;
			unsigned int kbnw = kbw;
			unsigned int ktnw = kw;
			unsigned int kbse = kbs;
			unsigned int ktsw = ksw;
			unsigned int kbne = kb;
			unsigned int ktne = k;
			unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real E     = (D.f[dirE   ])[ke   ];// +  c2over27 ;
			real W     = (D.f[dirW   ])[kw   ];// +  c2over27 ;
			real N     = (D.f[dirN   ])[kn   ];// +  c2over27 ;
			real S     = (D.f[dirS   ])[ks   ];// +  c2over27 ;
			real F     = (D.f[dirT   ])[kt   ];// +  c2over27 ;
			real B     = (D.f[dirB   ])[kb   ];// +  c2over27 ;
			real Ne    = (D.f[dirNE  ])[kne  ];// +  c1over54 ;
			real Sw    = (D.f[dirSW  ])[ksw  ];// +  c1over54 ;
			real Se    = (D.f[dirSE  ])[kse  ];// +  c1over54 ;
			real Nw    = (D.f[dirNW  ])[knw  ];// +  c1over54 ;
			real Ef    = (D.f[dirTE  ])[kte  ];// +  c1over54 ;
			real Wb    = (D.f[dirBW  ])[kbw  ];// +  c1over54 ;
			real Eb    = (D.f[dirBE  ])[kbe  ];// +  c1over54 ;
			real Wf    = (D.f[dirTW  ])[ktw  ];// +  c1over54 ;
			real Nf    = (D.f[dirTN  ])[ktn  ];// +  c1over54 ;
			real Sb    = (D.f[dirBS  ])[kbs  ];// +  c1over54 ;
			real Nb    = (D.f[dirBN  ])[kbn  ];// +  c1over54 ;
			real Sf    = (D.f[dirTS  ])[kts  ];// +  c1over54 ;
			real R     = (D.f[dirZERO])[kzero];// +  c8over27 ;
			real Nef   = (D.f[dirTNE ])[ktne ];// +  c1over216;
			real Swf   = (D.f[dirTSW ])[ktsw ];// +  c1over216;
			real Sef   = (D.f[dirTSE ])[ktse ];// +  c1over216;
			real Nwf   = (D.f[dirTNW ])[ktnw ];// +  c1over216;
			real Neb   = (D.f[dirBNE ])[kbne ];// +  c1over216;
			real Swb   = (D.f[dirBSW ])[kbsw ];// +  c1over216;
			real Seb   = (D.f[dirBSE ])[kbse ];// +  c1over216;
			real Nwb   = (D.f[dirBNW ])[kbnw ];// +  c1over216;
			////////////////////////////////////////////////////////////////////////////////////
			real fx = zero;
			real fy = zero;
			real fz = zero;
			////////////////////////////////////////////////////////////////////////////////////
			real rho=Nw+W+Sw+S+Se+E+Ne+N+R+Nf+Nb+Sf+Sb+Ef+Eb+Wf+Wb+Nwf+Nwb+Nef+Neb+Swf+Swb+Sef+Seb+F+B+one;// ACHTUNG ne EINS !!!!!!!!
			real pix=(Ne+E+Se+Ef+Eb-Nw-W-Sw-Wf-Wb+Nef+Neb+Sef+Seb-Nwf-Nwb-Swf-Swb);
			real piy=(Ne+N+Nw+Nf+Nb-Se-S-Sw-Sf-Sb+Nef+Neb+Nwf+Nwb-Sef-Seb-Swf-Swb);
			real piz=(Nf+Sf+Wf+Ef+F-Nb-Sb-Wb-Eb-B+Nef+Nwf+Sef+Swf-Neb-Nwb-Seb-Swb);
			real vvx=pix/rho + fx;
			real vvy=piy/rho + fy;
			real vvz=piz/rho + fz;
			real vx2=vvx*vvx;
			real vy2=vvy*vvy;
			real vz2=vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real mfaaa = Swb;
			real mfaab = Sw;
			real mfaac = Swf;
			real mfaba = Wb;
			real mfabb = W;
			real mfabc = Wf;
			real mfbaa = Sb;
			real mfbab = S;
			real mfbac = Sf;
			real mfbba = B;
			real mfbbb = R;
			real mfbbc = F;
			real mfaca = Nwb;
			real mfacb = Nw;
			real mfacc = Nwf;
			real mfcaa = Seb;
			real mfcab = Se;
			real mfcac = Sef;
			real mfcca = Neb;
			real mfccb = Ne;
			real mfccc = Nef;
			real mfbca = Nb;
			real mfbcb = N;
			real mfbcc = Nf;
			real mfcba = Eb;
			real mfcbb = E;
			real mfcbc = Ef;
			real m0, m1, m2;
			real wadjust;
			real qudricLimit = c1o100;
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - two*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - two*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += one;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - two*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - two*	m1 * vvx + vx2 * m0;
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
			real OxxPyyPzz = one;
			real OxyyPxzz  = one;//two+(-omega);//one;
			real OxyyMxzz  = one;//two+(-omega);//one;
			real O4        = one;
			real O5        = one;
			real O6        = one;

			//Cum 4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3 * rho) * mfabb + two* mfbba * mfbab) / rho;
			real CUMbcb = mfbcb - ((mfaca + c1o3 * rho) * mfbab + two* mfbba * mfabb) / rho;
			real CUMbbc = mfbbc - ((mfaac + c1o3 * rho) * mfbba + two* mfbab * mfabb) / rho; 

			real CUMcca = mfcca - (mfcaa * mfaca + two* mfbba * mfbba) / rho - c1o3 * (mfcaa + mfaca);
			real CUMcac = mfcac - (mfcaa * mfaac + two* mfbab * mfbab) / rho - c1o3 * (mfcaa + mfaac);
			real CUMacc = mfacc - (mfaac * mfaca + two* mfabb * mfabb) / rho - c1o3 * (mfaac + mfaca);

			//Cum 5.
			real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + four* mfabb * mfbbb + two* (mfbab * mfacb + mfbba * mfabc)) / rho - c1o3 * (mfbca + mfbac);
			real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + four* mfbab * mfbbb + two* (mfabb * mfcab + mfbba * mfbac)) / rho - c1o3 * (mfcba + mfabc);
			real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + four* mfbba * mfbbb + two* (mfbab * mfbca + mfabb * mfcba)) / rho - c1o3 * (mfacb + mfcab);

			//Cum 6.
			real CUMccc = mfccc  +(  -four*  mfbbb * mfbbb  
									-          (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
									-    four* (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
									-     two* (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
									+(   four* (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
									+     two* (mfcaa * mfaca * mfaac)
									+ sixteen*  mfbba * mfbab * mfabb) / (rho * rho)
									-    c1o3* (mfacc + mfcac + mfcca)
									+    c1o9* (mfcaa + mfaca + mfaac)
									+(    two* (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
									+          (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3 / rho;


			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;

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
 				mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 				mxxMyy    += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
 				mxxMzz    += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
 
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
			mfabb     += -(-omega) * (-mfabb);
			mfbab     += -(-omega) * (-mfbab);
			mfbba     += -(-omega) * (-mfbba);

			// linear combinations back
			mfcaa = c1o3 * (       mxxMyy +       mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-two * mxxMyy +       mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - two * mxxMzz + mxxPyyPzz);

			//3.
			// linear combinations
			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

			//relax
			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
			mfbbb     += wadjust * (-mfbbb);
			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
			mxxyPyzz  += wadjust * (-mxxyPyzz);
			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
			mxxyMyzz  += wadjust * (-mxxyMyzz);
			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
			mxxzPyyz  += wadjust * (-mxxzPyyz);
			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
			mxxzMyyz  += wadjust * (-mxxzMyyz);
			wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
			mxyyPxzz  += wadjust * (-mxyyPxzz);
			wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
			mxyyMxzz  += wadjust * (-mxyyMxzz);

			// linear combinations back
			mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
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
			mfcbb = CUMcbb + ((mfcaa + c1o3 * rho) * mfabb + two* mfbba * mfbab) / rho;
			mfbcb = CUMbcb + ((mfaca + c1o3 * rho) * mfbab + two* mfbba * mfabb) / rho;
			mfbbc = CUMbbc + ((mfaac + c1o3 * rho) * mfbba + two* mfbab * mfabb) / rho; 

			mfcca = CUMcca + (mfcaa * mfaca + two* mfbba * mfbba) / rho + c1o3 * (mfcaa + mfaca);
			mfcac = CUMcac + (mfcaa * mfaac + two* mfbab * mfbab) / rho + c1o3 * (mfcaa + mfaac);
			mfacc = CUMacc + (mfaac * mfaca + two* mfabb * mfabb) / rho + c1o3 * (mfaac + mfaca);

			//5.
			mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + four* mfabb * mfbbb + two* (mfbab * mfacb + mfbba * mfabc)) / rho + c1o3 * (mfbca + mfbac);
			mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + four* mfbab * mfbbb + two* (mfabb * mfcab + mfbba * mfbac)) / rho + c1o3 * (mfcba + mfabc);
			mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + four* mfbba * mfbbb + two* (mfbab * mfbca + mfabb * mfcba)) / rho + c1o3 * (mfacb + mfcab);

			//6.
			mfccc = CUMccc  -(( -four*  mfbbb * mfbbb  
							-          (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    four* (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     two* (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   four* (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     two* (mfcaa * mfaca * mfaac)
							+ sixteen*  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3* (mfacc + mfcac + mfcca)
							+    c1o9* (mfcaa + mfaca + mfaac)
							+(    two* (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
							+          (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3 / rho);
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
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + 1.) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - two* mfaab *  vvz         +  mfaaa       * (one- vz2)              - one* vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + 1.) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - two* mfabb *  vvz         + mfaba * (one- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - two* mfacb *  vvz         +  mfaca         * (one- vz2)              - c1o3 * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - two* mfbab *  vvz         + mfbaa * (one- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - two* mfbbb *  vvz         + mfbba * (one- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - two* mfbcb *  vvz         + mfbca * (one- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - two* mfcab *  vvz         +  mfcaa         * (one- vz2)              - c1o3 * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - two* mfcbb *  vvz         + mfcba * (one- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - two* mfccb *  vvz         +  mfcca         * (one- vz2)              - c1o9 * vz2; 
			m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9) * (     vz2 + vvz) * c1o2;
			mfcca = m0;
			mfccb = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6) * (     vy2 - vvy) * c1o2; 
			m1 = -mfaca        - two* mfaba *  vvy         +  mfaaa         * (one- vy2)              - c1o6 * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - two* mfabb *  vvy         +  mfaab         * (one- vy2)              - c2o3 * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - two* mfabc *  vvy         +  mfaac         * (one- vy2)              - c1o6 * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - two* mfbba *  vvy         + mfbaa * (one- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - two* mfbbb *  vvy         + mfbab * (one- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - two* mfbbc *  vvy         + mfbac * (one- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - two* mfcba *  vvy         +  mfcaa          * (one- vy2)              - c1o18 * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - two* mfcbb *  vvy         +  mfcab         * (one- vy2)              - c2o9 * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - two* mfcbc *  vvy         +  mfcac          * (one- vy2)              - c1o18 * vy2; 
			m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18) * (     vy2 + vvy) * c1o2;
			mfcac = m0;
			mfcbc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcaa        - two* mfbaa *  vvx         +  mfaaa          * (one- vx2)              - c1o36 * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - two* mfbba *  vvx         +  mfaba         * (one- vx2)              - c1o9 * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - two* mfbca *  vvx         +  mfaca          * (one- vx2)              - c1o36 * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - two* mfbab *  vvx         +  mfaab         * (one- vx2)              - c1o9 * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - two* mfbbb *  vvx         +  mfabb         * (one- vx2)              - c4o9 * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - two* mfbcb *  vvx         +  mfacb         * (one- vx2)              - c1o9 * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - two* mfbac *  vvx         +  mfaac          * (one- vx2)              - c1o36 * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - two* mfbbc *  vvx         +  mfabc         * (one- vx2)              - c1o9 * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - two* mfbcc *  vvx         +  mfacc          * (one- vx2)              - c1o36 * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dirE   ])[ke   ] = mfabb;// -  c2over27 ;//                                                                     
			(D.f[ dirW   ])[kw   ] = mfcbb;// -  c2over27 ;                                                                     
			(D.f[ dirN   ])[kn   ] = mfbab;// -  c2over27 ;
			(D.f[ dirS   ])[ks   ] = mfbcb;// -  c2over27 ;
			(D.f[ dirT   ])[kt   ] = mfbba;// -  c2over27 ;
			(D.f[ dirB   ])[kb   ] = mfbbc;// -  c2over27 ;
			(D.f[ dirNE  ])[kne  ] = mfaab;// -  c1over54 ;
			(D.f[ dirSW  ])[ksw  ] = mfccb;// -  c1over54 ;
			(D.f[ dirSE  ])[kse  ] = mfacb;// -  c1over54 ;
			(D.f[ dirNW  ])[knw  ] = mfcab;// -  c1over54 ;
			(D.f[ dirTE  ])[kte  ] = mfaba;// -  c1over54 ;
			(D.f[ dirBW  ])[kbw  ] = mfcbc;// -  c1over54 ;
			(D.f[ dirBE  ])[kbe  ] = mfabc;// -  c1over54 ;
			(D.f[ dirTW  ])[ktw  ] = mfcba;// -  c1over54 ;
			(D.f[ dirTN  ])[ktn  ] = mfbaa;// -  c1over54 ;
			(D.f[ dirBS  ])[kbs  ] = mfbcc;// -  c1over54 ;
			(D.f[ dirBN  ])[kbn  ] = mfbac;// -  c1over54 ;
			(D.f[ dirTS  ])[kts  ] = mfbca;// -  c1over54 ;
			(D.f[ dirZERO])[kzero] = mfbbb;// -  c8over27 ;
			(D.f[ dirTNE ])[ktne ] = mfaaa;// -  c1over216;
			(D.f[ dirTSE ])[ktse ] = mfaca;// -  c1over216;
			(D.f[ dirBNE ])[kbne ] = mfaac;// -  c1over216;
			(D.f[ dirBSE ])[kbse ] = mfacc;// -  c1over216;
			(D.f[ dirTNW ])[ktnw ] = mfcaa;// -  c1over216;
			(D.f[ dirTSW ])[ktsw ] = mfcca;// -  c1over216;
			(D.f[ dirBNW ])[kbnw ] = mfcac;// -  c1over216;
			(D.f[ dirBSW ])[kbsw ] = mfccc;// -  c1over216;
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Casc_Kum_SP_27(   real omega,
                                                       unsigned int* bcMatD,
                                                       unsigned int* neighborX,
                                                       unsigned int* neighborY,
                                                       unsigned int* neighborZ,
                                                       real* DDStart,
                                                       int size_Mat,
                                                       bool EvenOrOdd)
{
  // ////////////////////////////////////////////////////////////////////////////////
  // const unsigned  x = threadIdx.x;  // Globaler x-Index 
  // const unsigned  y = blockIdx.x;   // Globaler y-Index 
  // const unsigned  z = blockIdx.y;   // Globaler z-Index 

  // const unsigned nx = blockDim.x;
  // const unsigned ny = gridDim.x;

  // const unsigned k = nx*(ny*z + y) + x;
  // //////////////////////////////////////////////////////////////////////////

  // if(k<size_Mat)
  // {
  //    ////////////////////////////////////////////////////////////////////////////////
  //    unsigned int BC;
  //    BC        =   bcMatD[k];

  //    if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
  //    {
  //       Distributions27 D;
  //       if (EvenOrOdd==true)
  //       {
  //          D.f[dirE   ] = &DDStart[dirE   *size_Mat];
  //          D.f[dirW   ] = &DDStart[dirW   *size_Mat];
  //          D.f[dirN   ] = &DDStart[dirN   *size_Mat];
  //          D.f[dirS   ] = &DDStart[dirS   *size_Mat];
  //          D.f[dirT   ] = &DDStart[dirT   *size_Mat];
  //          D.f[dirB   ] = &DDStart[dirB   *size_Mat];
  //          D.f[dirNE  ] = &DDStart[dirNE  *size_Mat];
  //          D.f[dirSW  ] = &DDStart[dirSW  *size_Mat];
  //          D.f[dirSE  ] = &DDStart[dirSE  *size_Mat];
  //          D.f[dirNW  ] = &DDStart[dirNW  *size_Mat];
  //          D.f[dirTE  ] = &DDStart[dirTE  *size_Mat];
  //          D.f[dirBW  ] = &DDStart[dirBW  *size_Mat];
  //          D.f[dirBE  ] = &DDStart[dirBE  *size_Mat];
  //          D.f[dirTW  ] = &DDStart[dirTW  *size_Mat];
  //          D.f[dirTN  ] = &DDStart[dirTN  *size_Mat];
  //          D.f[dirBS  ] = &DDStart[dirBS  *size_Mat];
  //          D.f[dirBN  ] = &DDStart[dirBN  *size_Mat];
  //          D.f[dirTS  ] = &DDStart[dirTS  *size_Mat];
  //          D.f[dirZERO] = &DDStart[dirZERO*size_Mat];
  //          D.f[dirTNE ] = &DDStart[dirTNE *size_Mat];
  //          D.f[dirTSW ] = &DDStart[dirTSW *size_Mat];
  //          D.f[dirTSE ] = &DDStart[dirTSE *size_Mat];
  //          D.f[dirTNW ] = &DDStart[dirTNW *size_Mat];
  //          D.f[dirBNE ] = &DDStart[dirBNE *size_Mat];
  //          D.f[dirBSW ] = &DDStart[dirBSW *size_Mat];
  //          D.f[dirBSE ] = &DDStart[dirBSE *size_Mat];
  //          D.f[dirBNW ] = &DDStart[dirBNW *size_Mat];
  //       }
  //       else
  //       {
  //          D.f[dirW   ] = &DDStart[dirE   *size_Mat];
  //          D.f[dirE   ] = &DDStart[dirW   *size_Mat];
  //          D.f[dirS   ] = &DDStart[dirN   *size_Mat];
  //          D.f[dirN   ] = &DDStart[dirS   *size_Mat];
  //          D.f[dirB   ] = &DDStart[dirT   *size_Mat];
  //          D.f[dirT   ] = &DDStart[dirB   *size_Mat];
  //          D.f[dirSW  ] = &DDStart[dirNE  *size_Mat];
  //          D.f[dirNE  ] = &DDStart[dirSW  *size_Mat];
  //          D.f[dirNW  ] = &DDStart[dirSE  *size_Mat];
  //          D.f[dirSE  ] = &DDStart[dirNW  *size_Mat];
  //          D.f[dirBW  ] = &DDStart[dirTE  *size_Mat];
  //          D.f[dirTE  ] = &DDStart[dirBW  *size_Mat];
  //          D.f[dirTW  ] = &DDStart[dirBE  *size_Mat];
  //          D.f[dirBE  ] = &DDStart[dirTW  *size_Mat];
  //          D.f[dirBS  ] = &DDStart[dirTN  *size_Mat];
  //          D.f[dirTN  ] = &DDStart[dirBS  *size_Mat];
  //          D.f[dirTS  ] = &DDStart[dirBN  *size_Mat];
  //          D.f[dirBN  ] = &DDStart[dirTS  *size_Mat];
  //          D.f[dirZERO] = &DDStart[dirZERO*size_Mat];
  //          D.f[dirBSW ] = &DDStart[dirTNE *size_Mat];
  //          D.f[dirBNE ] = &DDStart[dirTSW *size_Mat];
  //          D.f[dirBNW ] = &DDStart[dirTSE *size_Mat];
  //          D.f[dirBSE ] = &DDStart[dirTNW *size_Mat];
  //          D.f[dirTSW ] = &DDStart[dirBNE *size_Mat];
  //          D.f[dirTNE ] = &DDStart[dirBSW *size_Mat];
  //          D.f[dirTNW ] = &DDStart[dirBSE *size_Mat];
  //          D.f[dirTSE ] = &DDStart[dirBNW *size_Mat];
  //       }

  //       ////////////////////////////////////////////////////////////////////////////////
  //       //index
  //       unsigned int kzero= k;
  //       unsigned int ke   = k;
  //       unsigned int kw   = neighborX[k];
  //       unsigned int kn   = k;
  //       unsigned int ks   = neighborY[k];
  //       unsigned int kt   = k;
  //       unsigned int kb   = neighborZ[k];
  //       unsigned int ksw  = neighborY[kw];
  //       unsigned int kne  = k;
  //       unsigned int kse  = ks;
  //       unsigned int knw  = kw;
  //       unsigned int kbw  = neighborZ[kw];
  //       unsigned int kte  = k;
  //       unsigned int kbe  = kb;
  //       unsigned int ktw  = kw;
  //       unsigned int kbs  = neighborZ[ks];
  //       unsigned int ktn  = k;
  //       unsigned int kbn  = kb;
  //       unsigned int kts  = ks;
  //       unsigned int ktse = ks;
  //       unsigned int kbnw = kbw;
  //       unsigned int ktnw = kw;
  //       unsigned int kbse = kbs;
  //       unsigned int ktsw = ksw;
  //       unsigned int kbne = kb;
  //       unsigned int ktne = k;
  //       unsigned int kbsw = neighborZ[ksw];
  //       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //       real E     = (D.f[dirE   ])[ke   ] +  c2over27 ;
  //       real W     = (D.f[dirW   ])[kw   ] +  c2over27 ;
  //       real N     = (D.f[dirN   ])[kn   ] +  c2over27 ;
  //       real S     = (D.f[dirS   ])[ks   ] +  c2over27 ;
  //       real F     = (D.f[dirT   ])[kt   ] +  c2over27 ;
  //       real B     = (D.f[dirB   ])[kb   ] +  c2over27 ;
  //       real Ne    = (D.f[dirNE  ])[kne  ] +  c1over54 ;
  //       real Sw    = (D.f[dirSW  ])[ksw  ] +  c1over54 ;
  //       real Se    = (D.f[dirSE  ])[kse  ] +  c1over54 ;
  //       real Nw    = (D.f[dirNW  ])[knw  ] +  c1over54 ;
  //       real Ef    = (D.f[dirTE  ])[kte  ] +  c1over54 ;
  //       real Wb    = (D.f[dirBW  ])[kbw  ] +  c1over54 ;
  //       real Eb    = (D.f[dirBE  ])[kbe  ] +  c1over54 ;
  //       real Wf    = (D.f[dirTW  ])[ktw  ] +  c1over54 ;
  //       real Nf    = (D.f[dirTN  ])[ktn  ] +  c1over54 ;
  //       real Sb    = (D.f[dirBS  ])[kbs  ] +  c1over54 ;
  //       real Nb    = (D.f[dirBN  ])[kbn  ] +  c1over54 ;
  //       real Sf    = (D.f[dirTS  ])[kts  ] +  c1over54 ;
  //       real R     = (D.f[dirZERO])[kzero] +  c8over27 ;
  //       real Nef   = (D.f[dirTNE ])[ktne ] +  c1over216;
  //       real Swf   = (D.f[dirTSW ])[ktsw ] +  c1over216;
  //       real Sef   = (D.f[dirTSE ])[ktse ] +  c1over216;
  //       real Nwf   = (D.f[dirTNW ])[ktnw ] +  c1over216;
  //       real Neb   = (D.f[dirBNE ])[kbne ] +  c1over216;
  //       real Swb   = (D.f[dirBSW ])[kbsw ] +  c1over216;
  //       real Seb   = (D.f[dirBSE ])[kbse ] +  c1over216;
  //       real Nwb   = (D.f[dirBNW ])[kbnw ] +  c1over216;
  //       ////////////////////////////////////////////////////////////////////////////////////
		// real o1   = omega;
		// real om3  = two - omega;
		// real om4  = two - omega;
		// real om5  = one;
		// real om6  = one;
		// real om7  = one;
		// real om8  = two - omega;
		// real om9  = one;
		// real om10 = one;
		// real om11 = one;
		// real QuadricLimiter = c1o100;
		// real vvx, vvy, vvz;
		// real pix, piy, piz;
		// real vx2,vy2,vz2;
  //       real rho=Nw+W+Sw+S+Se+E+Ne+N+R+Nf+Nb+Sf+Sb+Ef+Eb+Wf+Wb+Nwf+Nwb+Nef+Neb+Swf+Swb+Sef+Seb+F+B;
		// pix=(Ne+E+Se+Ef+Eb-Nw-W-Sw-Wf-Wb+Nef+Neb+Sef+Seb-Nwf-Nwb-Swf-Swb);
		// piy=(Ne+N+Nw+Nf+Nb-Se-S-Sw-Sf-Sb+Nef+Neb+Nwf+Nwb-Sef-Seb-Swf-Swb);
		// piz=(Nf+Sf+Wf+Ef+F-Nb-Sb-Wb-Eb-B+Nef+Nwf+Sef+Swf-Neb-Nwb-Seb-Swb);
  //       vvx=pix/rho;
  //       vvy=piy/rho;
  //       vvz=piz/rho;
  //       
  //       vx2=vvx*vvx;
  //       vy2=vvy*vvy;
  //       vz2=vvz*vvz;

  //       real x,y,z,xxx,yyy,zzz,uxy,uxz,uyz,pe,vxy,p,vxz,a,c,xxyyzz,xyyzz,xxyzz,xxyyz,uxxyz,uxyyz,uxyzz,xyz;
  //       real UXY,UXZ,UYZ,UXYZZ,UXYYZ,UXXYZ,XYZ,EPa,EPb,EPc,A,C,CA,XpXXX,XnXXX,YpYYY,YnYYY,ZpZZZ,ZnZZZ,X_,Y_,Z_;
  //
		// uxy=(o1*(Ne+Neb+Nef-Nw-Nwb-Nwf-Se-Seb-Sef+Sw+Swb+Swf+vvx*vvy*rho-pix*vvy-piy*vvx)*c1o4);
		// uxz=(o1*(-Eb+Ef-Neb+Nef+Nwb-Nwf-Seb+Sef+Swb-Swf+Wb-Wf+vvx*vvz*rho-pix*vvz-piz*vvx)*c1o4);
		// uyz=(o1*(-Nb-Neb+Nef+Nf-Nwb+Nwf+Sb+Seb-Sef-Sf+Swb-Swf+vvz*vvy*rho-piz*vvy-piy*vvz)*c1o4);

		// vxy=(o1*(c1o6)*(-B-E-F+Nb+Ne+Nf+Nw+Sb+Se+Sf+Sw-W+two*(-Eb-Ef+N+S-Wb-Wf)+two*(-two*piy*vvy+pix*vvx+piz*vvz)-rho*(vz2-two*vy2+vx2)));
		// vxz=(o1*(c1o6)*(-E+Eb+Ef-N+Nb+Nf-S+Sb+Sf-W+Wb+Wf+two*(B+F-Ne-Nw-Se-Sw)+two*(-two*piz*vvz+pix*vvx+piy*vvy)-rho*(vx2+vy2-two*vz2)));
		// pe=(om5*(c1over126)*(-B-E-F-N-S-W+two*(-Eb-Ef-Nb-Ne-Nf-Nw-Sb-Se-Sf-Sw-Wb-Wf)+three*(-Neb-Nef-Nwb-Nwf-Seb-Sef-Swb-Swf)+(rho+vvy*(two*piy-vvy*rho)+vvz*(two*piz-vvz*rho)+vvx*(two*pix-vvx*rho))));

		// UXY=(-Ne-Neb-Nef+Nw+Nwb+Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*uxy);
		// UXZ=(Eb-Ef+Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-Wb+Wf+four*uxz);
		// UYZ=(Nb+Neb-Nef-Nf+Nwb-Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*uyz);

		// ///dirty hack2
		// x=abs(((c1o8*(Nw+Sw+Wf+Wb-Eb-Ef-Ne-Se+vvx*(B+Eb+Ef+F+N+Ne+Nw+S+Se+Sw+Wb+Wf+two*(Nb+Neb+Nef+Nf+Nwb+Nwf+Sb+Seb+Sef+Sf+Swb+Swf)))
		//	+c1o4*(-Neb-Nef+Nwb+Nwf-Seb-Sef+Swb+Swf-vvy*(UXY-four*uxy)-vvz*(UXZ-four*uxz))
  //          -c1o4*vvx*(vvy*piy+vvz*piz)
  //          +c1o8*(vy2+vz2)*(vvx*rho-pix))));//+c1o8*vvx*(eightyfour*pe+two*(-vxy-vxz))-vvy*uxy-vvz*uxz;
		// y=abs(((c1o8*(Sb+Sf+Se+Sw-Nw-Ne-Nf-Nb+vvy*(B+E+F+Nb+Ne+Nf+Nw+Sb+Se+Sf+Sw+W+two*(Eb+Ef+Neb+Nef+Nwb+Nwf+Seb+Sef+Swb+Swf+Wb+Wf )))
		//	+c1o4*(-Neb-Nef-Nwb-Nwf+Seb+Sef+Swb+Swf-vvx*(UXY-four*uxy)-vvz*(UYZ-four*uyz))
  //          -c1o4*vvy*(vvz*piz+vvx*pix)
  //          +c1o8*(vx2+vz2)*(vvy*rho-piy))));//+c1o8*vvy*(eightyfour*pe+two*vxy)-vvx*uxy-vvz*uyz;
		// z=abs(((c1o8*(Eb+Nb+Sb+Wb-Ef-Nf-Sf-Wf+vvz*(E+Eb+Ef+N+Nb+Nf+S+Sb+Sf+W+Wb+Wf+two*(Ne+Neb+Nef+Nw+Nwb+Nwf+Se+Seb+Sef+Sw+Swb+Swf )))
		//	+c1o4*(-Nef+Neb+Nwb-Nwf+Seb-Sef+Swb-Swf-vvx*(UXZ-four*uxz)-vvy*(UYZ-four*uyz))
  //          -c1o4*vvz*(vvx*pix+vvy*piy)
  //          +c1o8*(vx2+vy2)*(vvz*rho-piz))));//+c1o8*vvz*(eightyfour*pe+two*vxz)-vvx*uxz-vvy*uyz;

		//xxx=abs(((c1o8*(Eb+Ef-Ne+Nw-Se+Sw-Wb-Wf+vvx*(-B-Eb-Ef-F+N+Ne+Nw+S+Se+Sw-Wb-Wf))+c1o4*(vvz*(UXZ-four*uxz)-vvy*(UXY-four*uxy)+vvx*(vvz*piz-vvy*piy))+c1o8*((pix-vvx*rho)*(vz2-vy2)))));//+c1o4*(vxz-vxy)*vvx +uxz*vvz-uxy*vvy);
		//yyy=abs(((c1o8*(Ne-Nb-Nf+Nw+Sb-Se+Sf-Sw+vvy*(-E+B+F+Nb-Ne+Nf-Nw+Sb-Se+Sf-Sw-W))+c1o4*(vvx*(UXY-four*uxy)-vvz*(UYZ-four*uyz)+vvy*(vvx*pix-vvz*piz))+c1o8*((piy-vvy*rho)*(vx2-vz2)))));//+c1o4*(-vxy-two*vxz)*vvy+uxy*vvx-uyz*vvz);
  //      zzz=abs(((c1o8*(Eb-Ef-Nb+Nf-Sb+Sf+Wb-Wf+vvz*( E+Eb+Ef-N-Nb-Nf-S-Sb-Sf+W+Wb+Wf))+c1o4*(vvy*(UYZ-four*uyz)-vvx*(UXZ-four*uxz)+vvz*(vvy*piy-vvx*pix))+c1o8*((piz-vvz*rho)*(vy2-vx2)))));//+c1o4*(vxz+two*vxy)*vvz +uyz*vvy-uxz*vvx);

		//xyz=abs(c1o8*(Neb-Nef-Nwb+Nwf-Seb+Sef+Swb-Swf-vvx*(UYZ-four*uyz)-vvy*(UXZ-four*uxz)-vvz*(UXY-four*uxy)+rho*vvx*vvy*vvz-(vvx*vvy*piz+vvx*piy*vvz+pix*vvy*vvz)));//-c1o2*(vvx*uyz+vvy*uxz+vvz*uxy);

  //      real oo;
		//
		//oo=om3+(one-om3)*(x/(x+QuadricLimiter/rho));
		//x=(oo*(c1o8*(Nw+Sw+Wf+Wb-Eb-Ef-Ne-Se+vvx*(B+Eb+Ef+F+N+Ne+Nw+S+Se+Sw+Wb+Wf+two*(Nb+Neb+Nef+Nf+Nwb+Nwf+Sb+Seb+Sef+Sf+Swb+Swf)))
		//	+c1o4*(-Neb-Nef+Nwb+Nwf-Seb-Sef+Swb+Swf-vvy*(UXY-four*uxy)-vvz*(UXZ-four*uxz))
  //          -c1o4*vvx*(vvy*piy+vvz*piz)
  //          +c1o8*(vy2+vz2)*(vvx*rho-pix)))+c1o8*vvx*(eightyfour*pe+two*(-vxy-vxz))-vvy*uxy-vvz*uxz;
		//
		//oo=om3+(one-om3)*(y/(y+QuadricLimiter/rho));
		//y=(oo*(c1o8*(Sb+Sf+Se+Sw-Nw-Ne-Nf-Nb+vvy*(B+E+F+Nb+Ne+Nf+Nw+Sb+Se+Sf+Sw+W+two*(Eb+Ef+Neb+Nef+Nwb+Nwf+Seb+Sef+Swb+Swf+Wb+Wf )))
		//	+c1o4*(-Neb-Nef-Nwb-Nwf+Seb+Sef+Swb+Swf-vvx*(UXY-four*uxy)-vvz*(UYZ-four*uyz))
  //          -c1o4*vvy*(vvz*piz+vvx*pix)
  //          +c1o8*(vx2+vz2)*(vvy*rho-piy)))+c1o8*vvy*(eightyfour*pe+two*vxy)-vvx*uxy-vvz*uyz;
  //      
		//oo=om3+(one-om3)*(z/(z+QuadricLimiter/rho));
		//z=(oo*(c1o8*(Eb+Nb+Sb+Wb-Ef-Nf-Sf-Wf+vvz*(E+Eb+Ef+N+Nb+Nf+S+Sb+Sf+W+Wb+Wf+two*(Ne+Neb+Nef+Nw+Nwb+Nwf+Se+Seb+Sef+Sw+Swb+Swf )))
		//	+c1o4*(-Nef+Neb+Nwb-Nwf+Seb-Sef+Swb-Swf-vvx*(UXZ-four*uxz)-vvy*(UYZ-four*uyz))
  //          -c1o4*vvz*(vvx*pix+vvy*piy)
  //          +c1o8*(vx2+vy2)*(vvz*rho-piz)))+c1o8*vvz*(eightyfour*pe+two*vxz)-vvx*uxz-vvy*uyz;

  //      oo=om4+(one-om4)*(xxx/(xxx+QuadricLimiter/rho));
		//xxx=(oo*(c1o8*(Eb+Ef-Ne+Nw-Se+Sw-Wb-Wf+vvx*(-B-Eb-Ef-F+N+Ne+Nw+S+Se+Sw-Wb-Wf))+c1o4*(vvz*(UXZ-four*uxz)-vvy*(UXY-four*uxy)+vvx*(vvz*piz-vvy*piy))+c1o8*((pix-vvx*rho)*(vz2-vy2))))+c1o4*(vxz-vxy)*vvx +uxz*vvz-uxy*vvy;
  //      
		//oo=om4+(one-om4)*(yyy/(yyy+QuadricLimiter/rho));
		//yyy=(oo*(c1o8*(Ne-Nb-Nf+Nw+Sb-Se+Sf-Sw+vvy*(-E+B+F+Nb-Ne+Nf-Nw+Sb-Se+Sf-Sw-W))+c1o4*(vvx*(UXY-four*uxy)-vvz*(UYZ-four*uyz)+vvy*(vvx*pix-vvz*piz))+c1o8*((piy-vvy*rho)*(vx2-vz2))))+c1o4*(-vxy-two*vxz)*vvy+uxy*vvx-uyz*vvz;
  //      
		//oo=om4+(one-om4)*(zzz/(zzz+QuadricLimiter/rho));
  //      zzz=(oo*(c1o8*(Eb-Ef-Nb+Nf-Sb+Sf+Wb-Wf+vvz*( E+Eb+Ef-N-Nb-Nf-S-Sb-Sf+W+Wb+Wf))+c1o4*(vvy*(UYZ-four*uyz)-vvx*(UXZ-four*uxz)+vvz*(vvy*piy-vvx*pix))+c1o8*((piz-vvz*rho)*(vy2-vx2))))+c1o4*(vxz+two*vxy)*vvz +uyz*vvy-uxz*vvx;

  //      oo=om8+(one-om8)*(xyz/(xyz+QuadricLimiter/rho));
		//xyz=(oo*c1o8*(Neb-Nef-Nwb+Nwf-Seb+Sef+Swb-Swf-vvx*(UYZ-four*uyz)-vvy*(UXZ-four*uxz)-vvz*(UXY-four*uxy)+rho*vvx*vvy*vvz-(vvx*vvy*piz+vvx*piy*vvz+pix*vvy*vvz)))-c1o2*(vvx*uyz+vvy*uxz+vvz*uxy);

  //              ///!dirty hack2

		//X_=(Eb+Ef+Ne-Nw+Se-Sw-Wb-Wf+two*(Neb+Nef-Nwb-Nwf+Seb+Sef-Swb-Swf)+eight*x);
		//Y_=(Nb+Ne+Nf+Nw-Sb-Se-Sf-Sw+two*(Neb+Nef+Nwb+Nwf-Seb-Sef-Swb-Swf)+eight*y);
		//Z_=(Ef-Eb-Nb+Nf-Sb+Sf-Wb+Wf+two*(Nef-Neb-Nwb+Nwf-Seb+Sef-Swb+Swf)+eight*z);

		//XnXXX=(Eb+Ef+Neb+Nef-Nwb-Nwf+Seb+Sef-Swb-Swf-Wb-Wf+four*(x-xxx));
		//XpXXX=(Ne+Neb+Nef-Nw-Nwb-Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*(x+xxx));
		//YnYYY=(Ne+Neb+Nef+Nw+Nwb+Nwf-Se-Seb-Sef-Sw-Swb-Swf+four*(y-yyy));
		//YpYYY=(Nb+Neb+Nef+Nf+Nwb+Nwf-Sb-Seb-Sef-Sf-Swb-Swf+four*(y+yyy));
		//ZnZZZ=(Nef-Nb-Neb+Nf-Nwb+Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*(z-zzz));
		//ZpZZZ=(Nef-Eb+Ef-Neb-Nwb+Nwf-Seb+Sef-Swb+Swf-Wb+Wf+four*(z+zzz));

		//XYZ=(-Neb+Nef+Nwb-Nwf+Seb-Sef-Swb+Swf+eight*xyz);
		//EPa=(-fortytwo*pe-E-Eb-Ef-Ne-Neb-Nef-Nw-Nwb-Nwf-Se-Seb-Sef-Sw-Swb-Swf-W-Wb-Wf-two*(vxy+vxz));
		//EPb=(-fortytwo*pe-B-Eb-Ef-F-Nb-Neb-Nef-Nf-Nwb-Nwf-Sb-Seb-Sef-Sf-Swb-Swf-Wb-Wf+two*vxz);
		//EPc=(-fortytwo*pe-N-Nb-Ne-Neb-Nef-Nf-Nw-Nwb-Nwf-S-Sb-Se-Seb-Sef-Sf-Sw-Swb-Swf+two*vxy);


  //      real kxx=(-two*vvx*pix+vx2*rho+(E+Ne+Se+Ef+Eb+Nef+Neb+Sef+Seb+W+Nw+Sw+Wf+Wb+Nwf+Nwb+Swf+Swb)+(fortytwo*pe+two*vxy+two*vxz)/om7)/rho;
  //      real kyy=(-two*vvy*piy+vy2*rho+(N+Ne+Nw+Nwf+Nwb+Nef+Neb+Nf+Nb+S+Se+Sw+Swf+Swb+Sef+Seb+Sf+Sb)+(fortytwo*pe-two*vxy)/om7)/rho;
  //      real kzz=(-two*vvz*piz+vz2*rho+(F+Nf+Wf+Ef+Sf+Nwf+Nef+Sef+Swf+B+Nb+Eb+Sb+Wb+Nwb+Neb+Seb+Swb)+(fortytwo*pe-two*vxz)/om7)/rho;

  //      real kxy=(-vvy*pix-vvx*piy+vvx*vvy*rho+(Ne+Neb+Nef+Sw+Swb+Swf-Nw-Nwb-Nwf-Se-Seb-Sef)-four*uxy/om7)/rho;
  //      real kxz=(-vvz*pix-vvx*piz+vvx*vvz*rho+(Ef+Nef+Nwb+Swb+Wb+Sef-Eb-Neb-Nwf-Seb-Swf-Wf)-four*uxz/om7)/rho;
  //      real kyz=(-vvy*piz-vvz*piy+vvy*vvz*rho+(Nef+Nf+Nwf+Sb+Seb+Swb-Nb-Neb-Nwb-Sef-Sf-Swf)-four*uyz/om7)/rho;

		//p=(om7*(c1o12)*((kxx*(kyy+kzz)+kyy*kzz+two*(kxy*kxy+kxz*kxz+kyz*kyz))*rho-ninetysix*pe/om7-Eb-Ef-Nb-Ne-Nf-Nw-Sb-Se-Sf-Sw-Wb-Wf-three*(Nwf+Nwb+Nef+Neb+Swf+Swb+Sef+Seb)
		//	+two*(+vvx*(Eb+Ef+Ne-Nw+Se-Sw-Wb-Wf+two*(Neb+Nef-Nwb-Nwf+Seb+Sef-Swb-Swf)+eight*x/om7)
		//	+vvy*(Nb+Ne+Nf+Nw-Sb-Se-Sf-Sw+two*(Neb+Nef+Nwb+Nwf-Seb-Sef-Swb-Swf)+eight*y/om7)+vvz*(Ef-Eb-Nb+Nf-Sb+Sf-Wb+Wf+two*(Nef-Neb-Nwb+Nwf-Seb+Sef-Swb+Swf)+eight*z/om7))
		//	+vx2*(-eightyfour*pe/om7-B-Eb-Ef-F-N-Ne-Nw-S-Se-Sw-Wb-Wf+two*((vxy+vxz)/om7-Nb-Neb-Nef-Nf-Nwb-Nwf-Sb-Seb-Sef-Sf-Swb-Swf))
		//	+vy2*(-eightyfour*pe/om7-B-E-F-Nb-Ne-Nf-Nw-Sb-Se-Sf-Sw-W+two*(-vxy/om7 -Eb-Ef-Neb-Nef-Nwb-Nwf-Seb-Sef-Swb-Swf-Wb-Wf))
		//	+vz2*(-eightyfour*pe/om7-E-Eb-Ef-N-Nb-Nf-S-Sb-Sf-W-Wb-Wf+two*(-vxz/om7 -Ne-Neb-Nef-Nw-Nwb-Nwf-Se-Seb-Sef-Sw-Swb-Swf))
		//	+four*(vvx*vvy*(-Ne-Neb-Nef+Nw+Nwb+Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*uxy/om7)
		//	+vvx*vvz*(Eb-Ef+Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-Wb+Wf+four*uxz/om7)+vvy*vvz*(Nb+Neb-Nef-Nf+Nwb-Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*uyz/om7))
  //          -rho*(vx2*vy2+vx2*vz2+vy2*vz2)+two*(piy*vvy*(vx2+vz2)+pix*vvx*(vy2+vz2)+piz*vvz*(vx2+vy2))));

		////brauchen wir nicht solange om6 == om7
	 // //kxx=(-two*vvx*pix+vx2*rho+(E+Ne+Se+Ef+Eb+Nef+Neb+Sef+Seb+W+Nw+Sw+Wf+Wb+Nwf+Nwb+Swf+Swb)+(fortytwo.*pe+two*vxy+two*vxz)/om6)/rho;
	 // //kyy=(-two*vvy*piy+vy2*rho+(N+Ne+Nw+Nwf+Nwb+Nef+Neb+Nf+Nb+S+Se+Sw+Swf+Swb+Sef+Seb+Sf+Sb)+(fortytwo.*pe-two*vxy)/om6)/rho;
	 // //kzz=(-two*vvz*piz+vz2*rho+(F+Nf+Wf+Ef+Sf+Nwf+Nef+Sef+Swf+B+Nb+Eb+Sb+Wb+Nwb+Neb+Seb+Swb)+(fortytwo.*pe-two*vxz)/om6)/rho;

	 // //kxy=(-vvy*pix-vvx*piy+vvx*vvy*rho+(Ne+Neb+Nef+Sw+Swb+Swf-Nw-Nwb-Nwf-Se-Seb-Sef)-four*uxy/om6)/rho;
	 // //kxz=(-vvz*pix-vvx*piz+vvx*vvz*rho+(Ef+Nef+Nwb+Swb+Wb+Sef-Eb-Neb-Nwf-Seb-Swf-Wf)-four*uxz/om6)/rho;
	 // //kyz=(-vvy*piz-vvz*piy+vvy*vvz*rho+(Nef+Nf+Nwf+Sb+Seb+Swb-Nb-Neb-Nwb-Sef-Sf-Swf)-four*uyz/om6)/rho;

	 // a=(om6*(c1o12)*((kxx*(kyy+kzz)-two*kyy*kzz+two*(kxy*kxy+kxz*kxz-two*kyz*kyz))*rho-Eb-Ef-Ne-Nw-Se-Sw-Wb-Wf+two*(Nb+Nf+Sb+Sf
		//+vvx*(Eb+Ef+Ne-Nw+Se-Sw-Wb-Wf+two*(Neb+Nef-Nwb-Nwf+Seb+Sef-Swb-Swf)+eight*x/om6)
		//+vvy*(Ne-Neb-Nef+Nw-Nwb-Nwf-Se+Seb+Sef-Sw+Swb+Swf+two*(Sf+Sb-Nb-Nf-two*y/om6-six*yyy/om6))
		//+vvz*(Ef-Eb+Neb-Nef+Nwb-Nwf+Seb-Sef+Swb-Swf-Wb+Wf+two*(Nb-Nf+Sb-Sf-two*z/om6+six*zzz/om6)))
		//+vx2*(-B-Eb-Ef-F-N-Ne-Nw-S-Se-Sw-Wb-Wf-eightyfour*pe/om6-two*(Nb+Neb+Nef+Nf+Nwb+Nwf+Sb+Seb+Sef+Sf+Swb+Swf-vxy/om6-vxz/om6))
		//+vy2*(Eb-E+Ef-Ne+Neb+Nef-Nw+Nwb+Nwf-Se+Seb+Sef-Sw+Swb+Swf-W+Wb+Wf+fortytwo*pe/om6+two*(B+F+Nb+Nf+Sb+Sf-vxy/om6-three*vxz/om6))
		//+vz2*(Ne-E-Eb-Ef+Neb+Nef+Nw+Nwb+Nwf+Se+Seb+Sef+Sw+Swb+Swf-W-Wb-Wf+fortytwo*pe/om6+two*(N+Nb+Nf+S+Sb+Sf-vxz/om6-three*vxy/om6))
		//+four*(vvx*(vvy*(-Ne-Neb-Nef+Nw+Nwb+Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*uxy/om6)+vvz*(Eb-Ef+Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-Wb+Wf+four*uxz/om6)))
		//-eight*vvy*vvz*(Nb+Neb-Nef-Nf+Nwb-Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*uyz/om6)+(two*vy2*vz2-vx2*(vy2+vz2))*rho
		//+two*(piy*vvy*(vx2-two*vz2)+piz*vvz*(vx2-two*vy2)+pix*vvx*(vy2+vz2))));
  //
	 // c=(om6*(c1o12)*((kyy*(kxx+kzz)-two*kxx*kzz+two*(kxy*kxy+kyz*kyz-two*kxz*kxz))*rho-Nb-Ne-Nf-Nw-Sb-Se-Sf-Sw+two*(Eb+Ef+Wb+Wf
		//  +vvx*(Ne-Neb-Nef-Nw+Nwb+Nwf+Se-Seb-Sef-Sw+Swb+Swf+two*(Wb+Wf-Eb-Ef-two*x/om6+six*xxx/om6))
		//  +vvy*(Nb+Ne+Nf+Nw-Sb-Se-Sf-Sw+two*(Neb+Nef+Nwb+Nwf-Seb-Sef-Swb-Swf)+eight*y/om6)
		//  +vvz*(Nf-Nb+Neb-Nef+Nwb-Nwf-Sb+Seb-Sef+Sf+Swb-Swf+two*(Eb-Ef+Wb-Wf-two*z/om6-six*zzz/om6)))
		//  +vx2*(Nb-N-Ne+Neb+Nef+Nf-Nw+Nwb+Nwf-S+Sb-Se+Seb+Sef+Sf-Sw+Swb+Swf+fortytwo*pe/om6+two*(B+Eb+Ef+F+Wb+Wf+vxy/om6-two*vxz/om6))
		//  -vy2*(B+E+F+Nb+Ne+Nf+Nw+Sb+Se+Sf+Sw+W+eightyfour*pe/om6+two*(Eb+Ef+Neb+Nef+Nwb+Nwf+Seb+Sef+Swb+Swf+Wb+Wf+vxy/om6))
		//  +vz2*(Ne-N-Nb+Neb+Nef-Nf+Nw+Nwb+Nwf-S-Sb+Se+Seb+Sef-Sf+Sw+Swb+Swf+fortytwo*pe/om6+two*(E+Eb+Ef+W+Wb+Wf+three*vxy/om6+two*vxz/om6))
		//  +four*vvy*(vvx*(-Ne-Neb-Nef+Nw+Nwb+Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*uxy/om6)+vvz*(Nb+Neb-Nef-Nf+Nwb-Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*uyz/om6))
		//  -eight*vvx*vvz*(Eb-Ef+Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-Wb+Wf+four*uxz/om6)+(two*vx2*vz2-(vx2+vz2)*vy2)*rho
		//  +two*(pix*vvx*(vy2-two*vz2)+piz*vvz*(vy2-two*vx2)+piy*vvy*(vx2+vz2))));

		////brauchen wir nicht solange om9 == om6
		////kxx=(-two*vvx*pix+vx2*rho+(E+Ne+Se+Ef+Eb+Nef+Neb+Sef+Seb+W+Nw+Sw+Wf+Wb+Nwf+Nwb+Swf+Swb)+(fortytwo.*pe+two*vxy+two*vxz)/om9)/rho;
		////kyy=(-two*vvy*piy+vy2*rho+(N+Ne+Nw+Nwf+Nwb+Nef+Neb+Nf+Nb+S+Se+Sw+Swf+Swb+Sef+Seb+Sf+Sb)+(fortytwo.*pe-two*vxy)/om9)/rho;
		////kzz=(-two*vvz*piz+vz2*rho+(F+Nf+Wf+Ef+Sf+Nwf+Nef+Sef+Swf+B+Nb+Eb+Sb+Wb+Nwb+Neb+Seb+Swb)+(fortytwo.*pe-two*vxz)/om9)/rho;

		////kxy=(-vvy*pix-vvx*piy+vvx*vvy*rho+(Ne+Neb+Nef+Sw+Swb+Swf-Nw-Nwb-Nwf-Se-Seb-Sef)-4.*uxy/om9)/rho;
		////kxz=(-vvz*pix-vvx*piz+vvx*vvz*rho+(Ef+Nef+Nwb+Swb+Wb+Sef-Eb-Neb-Nwf-Seb-Swf-Wf)-4.*uxz/om9)/rho;
		////kyz=(-vvy*piz-vvz*piy+vvy*vvz*rho+(Nef+Nf+Nwf+Sb+Seb+Swb-Nb-Neb-Nwb-Sef-Sf-Swf)-4.*uyz/om9)/rho;

		//uxxyz=(om9*(c1o8)*((kxx*kyz+two*kxy*kxz)*rho+Neb-Nef+Nwb-Nwf-Seb+Sef-Swb+Swf
		//	+vvy*(Nef-Eb+Ef-Neb-Nwb+Nwf-Seb+Sef-Swb+Swf-Wb+Wf+four*(z+zzz)/om9)+vvz*(Ne+Neb+Nef+Nw+Nwb+Nwf-Se-Seb-Sef-Sw-Swb-Swf+four*(y-yyy)/om9)
		//	+two*vvx*((-Neb+Nef+Nwb-Nwf+Seb-Sef-Swb+Swf+eight*xyz/om9)+vvy*(Eb-Ef+Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-Wb+Wf+four*uxz/om9)+vvz*(-Ne-Neb-Nef+Nw+Nwb+Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*uxy/om9))+vx2*(Nb+Neb-Nef-Nf+Nwb-Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*uyz/om9)
		//	+vvy*vvz*(-fortytwo*pe/om9-E-Eb-Ef-Ne-Neb-Nef-Nw-Nwb-Nwf-Se-Seb-Sef-Sw-Swb-Swf-W-Wb-Wf-two*(vxy+vxz)/om9)
  //          +(two*pix-rho*vvx)*vvx*vvy*vvz+vx2*(piy*vvz+piz*vvy)));
  //
		//uxyyz=(om9*(c1o8)*((kxz*kyy+two*kxy*kyz)*rho+Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf
		//	  +vvx*(Nef-Nb-Neb+Nf-Nwb+Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*(z-zzz)/om9)/*YnYYY*/+vvz*(Ne+Neb+Nef-Nw-Nwb-Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*(x+xxx)/om9)
		//	  +two*vvy*((-Neb+Nef+Nwb-Nwf+Seb-Sef-Swb+Swf+eight*xyz/om9)+vvx*(Nb+Neb-Nef-Nf+Nwb-Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*uyz/om9)+vvz*(-Ne-Neb-Nef+Nw+Nwb+Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*uxy/om9))+vy2*(Eb-Ef+Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-Wb+Wf+four*uxz/om9)
		//	  +vvx*vvz*(-fortytwo*pe/om9-N-Nb-Ne-Neb-Nef-Nf-Nw-Nwb-Nwf-S-Sb-Se-Seb-Sef-Sf-Sw-Swb-Swf+two*vxy/om9)
  //            //+rho*3*vy2*vvx*vvz
  //            +(two*piy-vvy*rho)*vvy*vvx*vvz
  //            +vy2*(pix*vvz+piz*vvx)));
		//uxyzz=(om9*(c1o8)*((kxy*kzz+two*kxz*kyz)*rho+Nwb-Neb-Nef+Nwf+Seb+Sef-Swb-Swf
		//	+vvx*(Nb+Neb+Nef+Nf+Nwb+Nwf-Sb-Seb-Sef-Sf-Swb-Swf+four*(y+yyy)/om9)+vvy*(Eb+Ef+Neb+Nef-Nwb-Nwf+Seb+Sef-Swb-Swf-Wb-Wf+four*(x-xxx)/om9)
		//	+two*vvz*((-Neb+Nef+Nwb-Nwf+Seb-Sef-Swb+Swf+eight*xyz/om9)+vvy*(Eb-Ef+Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-Wb+Wf+four*uxz/om9)+vvx*(Nb+Neb-Nef-Nf+Nwb-Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*uyz/om9))+vz2*(-Ne-Neb-Nef+Nw+Nwb+Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*uxy/om9)
		//	+vvx*vvy*(-fortytwo*pe/om9-B-Eb-Ef-F-Nb-Neb-Nef-Nf-Nwb-Nwf-Sb-Seb-Sef-Sf-Swb-Swf-Wb-Wf+two*vxz/om9)
  //          //+rho*3*vz2*vvx*vvy
  //          +(two*piz-vvz*rho)*vvx*vvy*vvz
  //          +vz2*(piy*vvx+pix*vvy)));

		////brauchen wir nicht solange om10 == om9
		////kxx=(-two*vvx*pix+vx2*rho+(E+Ne+Se+Ef+Eb+Nef+Neb+Sef+Seb+W+Nw+Sw+Wf+Wb+Nwf+Nwb+Swf+Swb)+(fortytwo.*pe+two*vxy+two*vxz)/om10)/rho;
		////kyy=(-two*vvy*piy+vy2*rho+(N+Ne+Nw+Nwf+Nwb+Nef+Neb+Nf+Nb+S+Se+Sw+Swf+Swb+Sef+Seb+Sf+Sb)+(fortytwo.*pe-two*vxy)/om10)/rho;
		////kzz=(-two*vvz*piz+vz2*rho+(F+Nf+Wf+Ef+Sf+Nwf+Nef+Sef+Swf+B+Nb+Eb+Sb+Wb+Nwb+Neb+Seb+Swb)+(fortytwo.*pe-two*vxz)/om10)/rho;

		////kxy=(-vvy*pix-vvx*piy+vvx*vvy*rho+(Ne+Neb+Nef+Sw+Swb+Swf-Nw-Nwb-Nwf-Se-Seb-Sef)-four*uxy/om10)/rho;
		////kxz=(-vvz*pix-vvx*piz+vvx*vvz*rho+(Ef+Nef+Nwb+Swb+Wb+Sef-Eb-Neb-Nwf-Seb-Swf-Wf)-four*uxz/om10)/rho;
		////kyz=(-vvy*piz-vvz*piy+vvy*vvz*rho+(Nef+Nf+Nwf+Sb+Seb+Swb-Nb-Neb-Nwb-Sef-Sf-Swf)-four*uyz/om10)/rho;

  //      real kxyy=(Ne + Neb + Nef - Nw - Nwb - Nwf + Se + Seb + Sef - Sw -Swb - Swf + (-two*Ne - two*Neb - two*Nef + two*Nw + two*Nwb + 
		//				two*Nwf + two*Se + two*Seb + two*Sef - two*Sw - two*Swb - two*Swf + eight*uxy/om10)*vvy + (two*N + two*Nb + two*Ne + two*Neb + two*Nef + two*Nf + 
		//				two*Nw + two*Nwb + two*Nwf - two*S - two*Sb - two*Se - two*Seb - two*Sef - two*Sf - two*Sw - two*Swb - two*Swf)*vvx*vvy + vvx*(-N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb 
		//				- Nwf - fortytwo*pe/om10 - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf + two*vxy/om10 -rho*vvy*vvy) + vvy*vvy*(E + Eb + Ef + Ne + 
		//				Neb + Nef - Nw - Nwb - Nwf + Se + Seb + Sef - Sw - Swb - Swf - W - Wb - Wf)+ four*x/om10 + four*xxx/om10)/rho;

  //      real kxxy=(Ne + Neb + Nef + Nw + Nwb + Nwf - Se - Seb - Sef - Sw - Swb - Swf + (-two*Ne - two*Neb - two*Nef + two*Nw + two*Nwb 
		//				+ two*Nwf + two*Se + two*Seb + two*Sef - two*Sw - two*Swb - two*Swf + eight*uxy/om10)*vvx + vvx*vvx*(N + Nb + Ne + Neb + Nef + Nf + Nw 
		//				+ Nwb + Nwf - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf -rho*vvy) + vvx*vvy*(two*E + two*Eb + two*Ef + two*Ne + two*Neb + 
		//				two*Nef - two*Nw - two*Nwb - two*Nwf + two*Se + two*Seb + two*Sef - two*Sw - two*Swb - two*Swf - two*W - two*Wb - two*Wf)+ vvy*(-E - Eb - Ef - Ne - Neb - Nef - Nw - 
		//				Nwb - Nwf - fortytwo*pe/om10 - Se - Seb - Sef - Sw - Swb - Swf - two*vxy/om10 - two*vxz/om10 - W - Wb - Wf) + four*y/om10 - four*yyy/om10)/rho;

  //     real kxxz=(-Eb + Ef - Neb + Nef - Nwb + Nwf - Seb + Sef - Swb + Swf - Wb + vvx*vvz*(two*E + two*Eb + two*Ef + two*Ne + two*Neb 
		//			+ two*Nef - two*Nw - two*Nwb - two*Nwf + two*Se + two*Seb + two*Sef - two*Sw - two*Swb - two*Swf - two*W - two*Wb - two*Wf)
  //                  + vvz*(-E - Eb - Ef - Ne - Neb - Nef - Nw - Nwb - Nwf - fortytwo*pe/om10 - Se - Seb - Sef - Sw - Swb - Swf - two*vxy/om10 - two*vxz/om10 - W - Wb - Wf)
  //                  + Wf + vvx*vvx*(-B - Eb + Ef + F - Nb - Neb + Nef + Nf - Nwb + Nwf - Sb - Seb + Sef + Sf - Swb + Swf - Wb
  //                  + vvz*(-B - E - Eb - Ef - F - N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - R - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf - W - Wb - Wf) + Wf)
  //                  + vvx*(two*Eb - two*Ef + two*Neb - two*Nef - two*Nwb + two*Nwf + two*Seb - two*Sef - two*Swb + two*Swf + eight*uxz/om10 - two*Wb + two*Wf) + four*z/om10 + four*zzz/om10)/rho;

  //     real kxzz=(Eb + Ef + Neb + Nef - Nwb - Nwf + Seb + Sef - Swb - Swf - Wb + vvz*vvz*(E + Eb + Ef + Ne + Neb + Nef - Nw - Nwb - Nwf + Se + Seb + Sef - Sw - Swb - Swf - W - Wb - Wf)
  //                  + vvx*(-B - Eb - Ef - F - Nb - Neb - Nef - Nf - Nwb - Nwf - fortytwo*pe/om10 - Sb - Seb - Sef - Sf - Swb - Swf + two*vxz/om10 - Wb
  //                  + vvz*vvz*(-B - E - Eb - Ef - F - N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - R - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf - W - Wb - Wf) - Wf) - Wf
  //                  + vvx*vvz*(-two*B - two*Eb + two*Ef + two*F - two*Nb - two*Neb + two*Nef + two*Nf - two*Nwb + two*Nwf - two*Sb - two*Seb + two*Sef + two*Sf - two*Swb + two*Swf - two*Wb + two*Wf)
  //                  + vvz*(two*Eb - two*Ef + two*Neb - two*Nef - two*Nwb + two*Nwf + two*Seb - two*Sef - two*Swb + two*Swf + eight*uxz/om10 - two*Wb + two*Wf) + four*x/om10 - four*xxx/om10)/rho;

  //     real kyzz=(Nb + Neb + Nef + Nf + Nwb + Nwf - Sb - Seb - Sef - Sf - Swb - Swf + (two*Nb + two*Neb - two*Nef - two*Nf + two*Nwb - two*Nwf - two*Sb - two*Seb + two*Sef + two*Sf - two*Swb + two*Swf + eight*uyz/om10)*vvz
  //                  + (N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf)*vvz*vvz + vvy*(-B - Eb - Ef - F - Nb - Neb - Nef - Nf - Nwb - Nwf - fortytwo*pe/om10 - Sb - Seb - Sef - Sf - Swb - Swf + two*vxz/om10 - Wb
  //                  + vvz*vvz*(-B - E - Eb - Ef - F - N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - R - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf - W - Wb - Wf) - Wf)
  //                  + vvy*vvz*(-two*B - two*Eb + two*Ef + two*F - two*Nb - two*Neb + two*Nef + two*Nf - two*Nwb + two*Nwf - two*Sb - two*Seb + two*Sef + two*Sf - two*Swb + two*Swf - two*Wb + two*Wf) + four*y/om10 + four*yyy/om10)/rho;

  //     real kyyz=(-Nb - Neb + Nef + Nf - Nwb + Nwf - Sb - Seb + Sef + Sf - Swb + Swf + (two*Nb + two*Neb - two*Nef - two*Nf + two*Nwb - two*Nwf - two*Sb - two*Seb + two*Sef + two*Sf - two*Swb + two*Swf + eight*uyz/om10)*vvy
  //                   + (-N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - fortytwo*pe/om10 - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf + two*vxy/om10)*vvz
  //                   + (two*N + two*Nb + two*Ne + two*Neb + two*Nef + two*Nf + two*Nw + two*Nwb + two*Nwf - two*S - two*Sb - two*Se - two*Seb - two*Sef - two*Sf - two*Sw - two*Swb - two*Swf)*vvy*vvz
  //                   + vvy*vvy*(-B - Eb + Ef + F - Nb - Neb + Nef + Nf - Nwb + Nwf - Sb - Seb + Sef + Sf - Swb + Swf - Wb
		//			 + vvz*(-B - E - Eb - Ef - F - N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - R - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf - W - Wb - Wf) + Wf) + four*z/om10 - four*zzz/om10)/rho;

  //     real kxyz=(-Neb + Nef + Nwb - Nwf + Seb - Sef - Swb + Swf + (Nb + Neb - Nef - Nf + Nwb - Nwf - Sb - Seb + Sef + Sf - Swb + Swf + four*uyz/om10)*vvx
  //                  + (N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf)*vvx*vvz + vvz*(-Ne - Neb - Nef + Nw + Nwb + Nwf + Se + Seb + Sef - Sw - Swb - Swf + four*uxy/om10
  //                  + vvx*vvy*(-B - E - Eb - Ef - F - N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - R - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf - W - Wb - Wf)) + vvy*vvz*(E + Eb + Ef + Ne + Neb + Nef - Nw - Nwb - Nwf + Se + Seb + Sef - Sw - Swb - Swf - W - Wb - Wf)
  //                  + vvx*vvy*(-B - Eb + Ef + F - Nb - Neb + Nef + Nf - Nwb + Nwf - Sb - Seb + Sef + Sf - Swb + Swf - Wb + Wf) + vvy*(Eb - Ef + Neb - Nef - Nwb + Nwf + Seb - Sef - Swb + Swf + four*uxz/om10 - Wb + Wf) + eight*xyz/om10)/rho;

		//xyyzz=(om10*(c1o8)*((kzz*kxyy+kyy*kxzz+two*(kxz*kyyz+kxy*kyzz)+four*kyz*kxyz)*rho+Nwb-Neb-Nef+Nwf-Seb-Sef+Swb+Swf
		//	  -vvx*(four*a/om10-thirtytwo*pe/om10-Nb-Neb-Nef-Nf-Nwb-Nwf-four*p/om10-Sb-Seb-Sef-Sf-Swb-Swf)-two*(vvy*(Nwb-Neb-Nef+Nwf+Seb+Sef-Swb-Swf-eight*uxyzz/om10)+vvz*(Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-eight*uxyyz/om10))
		//	  -vy2*(Eb+Ef+Neb+Nef-Nwb-Nwf+Seb+Sef-Swb-Swf-Wb-Wf+four*(x-xxx)/om10)-vz2*(Ne+Neb+Nef-Nw-Nwb-Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*(x+xxx)/om10)+two*vvx*(-vvy*(Nb+Neb+Nef+Nf+Nwb+Nwf-Sb-Seb-Sef-Sf-Swb-Swf+four*(y+yyy)/om10)-vvz*(Nef-Nb-Neb+Nf-Nwb+Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*(z-zzz)/om10))-four*vvy*vvz*((-Neb+Nef+Nwb-Nwf+Seb-Sef-Swb+Swf+eight*xyz/om10)+vvx*(Nb+Neb-Nef-Nf+Nwb-Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*uyz/om10))
		//	  -vvx*(vz2*(-fortytwo*pe/om10-N-Nb-Ne-Neb-Nef-Nf-Nw-Nwb-Nwf-S-Sb-Se-Seb-Sef-Sf-Sw-Swb-Swf+two*vxy/om10)+vy2*(-fortytwo*pe/om10-B-Eb-Ef-F-Nb-Neb-Nef-Nf-Nwb-Nwf-Sb-Seb-Sef-Sf-Swb-Swf-Wb-Wf+two*vxz/om10))
		//	  -two*(vy2*vvz*(Eb-Ef+Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-Wb+Wf+four*uxz/om10)+vvy*vz2*(-Ne-Neb-Nef+Nw+Nwb+Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*uxy/om10))
  //                       //-four*vvx*vy2*vz2*rho
		//	  -two*vvx*vvy*vvz*(piy*vvz+piz*vvy)+vy2*vz2*(vvx*rho-pix)));

		//xxyzz=(om10*(c1o8)*((kxx*kyzz+kzz*kxxy+two*(kxy*kxzz+kyz*kxxz)+four*kxz*kxyz)*rho+Seb-Neb-Nef-Nwb-Nwf+Sef+Swb+Swf
		//	  -vvy*(four*(c-p)/om10-Eb-Ef-thirtytwo*pe/om10-Neb-Nef-Nwb-Nwf-Seb-Sef-Swb-Swf-Wb-Wf)-two*(vvx*(Nwb-Neb-Nef+Nwf+Seb+Sef-Swb-Swf-eight*uxyzz/om10)+vvz*(Neb-Nef+Nwb-Nwf-Seb+Sef-Swb+Swf-eight*uxxyz/om10))
		//	  -vx2*(Nb+Neb+Nef+Nf+Nwb+Nwf-Sb-Seb-Sef-Sf-Swb-Swf+four*(y+yyy)/om10)-vz2*(Ne+Neb+Nef+Nw+Nwb+Nwf-Se-Seb-Sef-Sw-Swb-Swf+four*(y-yyy)/om10)+two*vvy*(-vvx*(Eb+Ef+Neb+Nef-Nwb-Nwf+Seb+Sef-Swb-Swf-Wb-Wf+four*(x-xxx)/om10)-vvz*(Nef-Eb+Ef-Neb-Nwb+Nwf-Seb+Sef-Swb+Swf-Wb+Wf+four*(z+zzz)/om10))-four*vvx*vvz*((-Neb+Nef+Nwb-Nwf+Seb-Sef-Swb+Swf+eight*xyz/om10)+vvy*(Eb-Ef+Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-Wb+Wf+four*uxz/om10))
		//	  -vvy*(vx2*(-fortytwo*pe/om10-B-Eb-Ef-F-Nb-Neb-Nef-Nf-Nwb-Nwf-Sb-Seb-Sef-Sf-Swb-Swf-Wb-Wf+two*vxz/om10)+vz2*(-fortytwo*pe/om10-E-Eb-Ef-Ne-Neb-Nef-Nw-Nwb-Nwf-Se-Seb-Sef-Sw-Swb-Swf-W-Wb-Wf-two*(vxy+vxz)/om10))
		//	  -two*(vx2*vvz*(Nb+Neb-Nef-Nf+Nwb-Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*uyz/om10)+vvx*vz2*(-Ne-Neb-Nef+Nw+Nwb+Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*uxy/om10))
		//							 //-four*vvy*vx2*vz2*rho
		//	  -two*vvx*vvy*vvz*(pix*vvz+piz*vvx)+vx2*vz2*(vvy*rho-piy)));
		//xxyyz=(om10*(c1o8)*((kxx*kyyz+kyy*kxxz+two*(kxz*kxyy+kyz*kxxy)+four*kxy*kxyz)*rho+Neb-Nef+Nwb-Nwf+Seb-Sef+Swb-Swf
		//	  -vvz*(-four*(a+c+p)/om10-thirtytwo*pe/om10-Ne-Neb-Nef-Nw-Nwb-Nwf-Se-Seb-Sef-Sw-Swb-Swf)-two*(vvx*(Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-eight*uxyyz/om10)+vvy*(Neb-Nef+Nwb-Nwf-Seb+Sef-Swb+Swf-eight*uxxyz/om10))
		//	  -vx2*(Nef-Nb-Neb+Nf-Nwb+Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*(z-zzz)/om10)-vy2*(Nef-Eb+Ef-Neb-Nwb+Nwf-Seb+Sef-Swb+Swf-Wb+Wf+four*(z+zzz)/om10)+two*vvz*(-vvx*(Ne+Neb+Nef-Nw-Nwb-Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*(x+xxx)/om10)-vvy*(Ne+Neb+Nef+Nw+Nwb+Nwf-Se-Seb-Sef-Sw-Swb-Swf+four*(y-yyy)/om10))-four*vvx*vvy*((-Neb+Nef+Nwb-Nwf+Seb-Sef-Swb+Swf+eight*xyz/om10)+vvz*(-Ne-Neb-Nef+Nw+Nwb+Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*uxy/om10))
		//	  -vvz*(vx2*(-fortytwo*pe/om10-N-Nb-Ne-Neb-Nef-Nf-Nw-Nwb-Nwf-S-Sb-Se-Seb-Sef-Sf-Sw-Swb-Swf+two*vxy/om10)+vy2*(-fortytwo*pe/om10-E-Eb-Ef-Ne-Neb-Nef-Nw-Nwb-Nwf-Se-Seb-Sef-Sw-Swb-Swf-W-Wb-Wf-two*(vxy+vxz)/om10))
		//	  -two*(vx2*vvy*(Nb+Neb-Nef-Nf+Nwb-Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*uyz/om10)+vvx*vy2*(Eb-Ef+Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-Wb+Wf+four*uxz/om10))
		//							 //-four*vx2*vy2*vvz*rho
		//	  -two*vvx*vvy*vvz*(pix*vvy+piy*vvx)+vx2*vy2*(vvz*rho-piz)));


		//	//brauchen wir nicht solange om11 == om10
		//  //kxx=(-two*vvx*pix+vx2*rho+(E+Ne+Se+Ef+Eb+Nef+Neb+Sef+Seb+W+Nw+Sw+Wf+Wb+Nwf+Nwb+Swf+Swb)+(fortytwo.*pe+two*vxy+two*vxz)/om11)/rho;
		//  //kyy=(-two*vvy*piy+vy2*rho+(N+Ne+Nw+Nwf+Nwb+Nef+Neb+Nf+Nb+S+Se+Sw+Swf+Swb+Sef+Seb+Sf+Sb)+(fortytwo.*pe-two*vxy)/om11)/rho;
		//  //kzz=(-two*vvz*piz+vz2*rho+(F+Nf+Wf+Ef+Sf+Nwf+Nef+Sef+Swf+B+Nb+Eb+Sb+Wb+Nwb+Neb+Seb+Swb)+(fortytwo.*pe-two*vxz)/om11)/rho;

		//  //kxy=(-vvy*pix-vvx*piy+vvx*vvy*rho+(Ne+Neb+Nef+Sw+Swb+Swf-Nw-Nwb-Nwf-Se-Seb-Sef)-four*uxy/om11)/rho;
		//  //kxz=(-vvz*pix-vvx*piz+vvx*vvz*rho+(Ef+Nef+Nwb+Swb+Wb+Sef-Eb-Neb-Nwf-Seb-Swf-Wf)-four*uxz/om11)/rho;
		//  //kyz=(-vvy*piz-vvz*piy+vvy*vvz*rho+(Nef+Nf+Nwf+Sb+Seb+Swb-Nb-Neb-Nwb-Sef-Sf-Swf)-four*uyz/om11)/rho;

  //  //      kxyy=(Ne + Neb + Nef - Nw - Nwb - Nwf + Se + Seb + Sef - Sw -Swb - Swf + (-two*Ne - two*Neb - two*Nef + two*Nw + two*Nwb + two*Nwf + two*Se + two*Seb + two*Sef - two*Sw - two*Swb - two*Swf + eight*uxy/om11)*vvy + (two*N + two*Nb + two*Ne + two*Neb + two*Nef + two*Nf + 
		//		//two*Nw + two*Nwb + two*Nwf - two*S - two*Sb - two*Se - two*Seb - two*Sef - two*Sf - two*Sw - two*Swb - two*Swf)*vvx*vvy + vvx*(-N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - fortytwo*pe/om11 - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf + two*vxy/om11
  //  //            -rho*vvy*vvy) + vvy*vvy*(E + Eb + Ef + Ne + Neb + Nef - Nw - Nwb - Nwf + Se + Seb + Sef - Sw - Swb - Swf - W - Wb - Wf)+ four*x/om11 + four*xxx/om11)/rho;

  //  //      kxxy=(Ne + Neb + Nef + Nw + Nwb + Nwf - Se - Seb - Sef - Sw - Swb - Swf + (-two*Ne - two*Neb - two*Nef + two*Nw + two*Nwb + two*Nwf + two*Se + two*Seb + two*Sef - two*Sw - two*Swb - two*Swf + eight*uxy/om11)*vvx
  //  //           + vvx*vvx*(N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf -rho*vvy) + vvx*vvy*(two*E + two*Eb + two*Ef + two*Ne + two*Neb + 
		//		//two*Nef - two*Nw - two*Nwb - two*Nwf + two*Se + two*Seb + two*Sef - two*Sw - two*Swb - two*Swf - two*W - two*Wb - two*Wf)
  //  //           + vvy*(-E - Eb - Ef - Ne - Neb - Nef - Nw - Nwb - Nwf - fortytwo*pe/om11 - Se - Seb - Sef - Sw - Swb - Swf - two*vxy/om11 - two*vxz/om11 - W - Wb - Wf) + four*y/om11 - four*yyy/om11)/rho;

  //  //      kxxz=(-Eb + Ef - Neb + Nef - Nwb + Nwf - Seb + Sef - Swb + Swf - Wb + vvx*vvz*(two*E + two*Eb + two*Ef + two*Ne + two*Neb + two*Nef - two*Nw - two*Nwb - two*Nwf + two*Se + two*Seb + two*Sef - two*Sw - two*Swb - two*Swf - two*W - two*Wb - two*Wf)
  //  //           + vvz*(-E - Eb - Ef - Ne - Neb - Nef - Nw - Nwb - Nwf - fortytwo*pe/om11 - Se - Seb - Sef - Sw - Swb - Swf - two*vxy/om11 - two*vxz/om11 - W - Wb - Wf) + Wf + vvx*vvx*(-B - Eb + Ef + F - Nb - Neb + Nef + Nf - Nwb + Nwf - Sb - Seb + Sef + Sf - Swb + Swf - Wb
  //  //           + vvz*(-B - E - Eb - Ef - F - N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - R - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf - W - Wb - Wf) + Wf)
  //  //           + vvx*(two*Eb - two*Ef + two*Neb - two*Nef - two*Nwb + two*Nwf + two*Seb - two*Sef - two*Swb + two*Swf + eight*uxz/om11 - two*Wb + two*Wf) + four*z/om11 + four*zzz/om11)/rho;

  //  //      kxzz=(Eb + Ef + Neb + Nef - Nwb - Nwf + Seb + Sef - Swb - Swf - Wb + vvz*vvz*(E + Eb + Ef + Ne + Neb + Nef - Nw - Nwb - Nwf + Se + Seb + Sef - Sw - Swb - Swf - W - Wb - Wf)
  //  //           + vvx*(-B - Eb - Ef - F - Nb - Neb - Nef - Nf - Nwb - Nwf - fortytwo*pe/om11 - Sb - Seb - Sef - Sf - Swb - Swf + two*vxz/om11 - Wb
  //  //           + vvz*vvz*(-B - E - Eb - Ef - F - N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - R - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf - W - Wb - Wf) - Wf) - Wf
  //  //           + vvx*vvz*(-two*B - two*Eb + two*Ef + two*F - two*Nb - two*Neb + two*Nef + two*Nf - two*Nwb + two*Nwf - two*Sb - two*Seb + two*Sef + two*Sf - two*Swb + two*Swf - two*Wb + two*Wf)
  //  //           + vvz*(two*Eb - two*Ef + two*Neb - two*Nef - two*Nwb + two*Nwf + two*Seb - two*Sef - two*Swb + two*Swf + eight*uxz/om11 - two*Wb + two*Wf) + four*x/om11 - four*xxx/om11)/rho;

  //  //      kyzz=(Nb + Neb + Nef + Nf + Nwb + Nwf - Sb - Seb - Sef - Sf - Swb - Swf + (two*Nb + two*Neb - two*Nef - two*Nf + two*Nwb - two*Nwf - two*Sb - two*Seb + two*Sef + two*Sf - two*Swb + two*Swf + eight*uyz/om11)*vvz
  //  //           + (N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf)*vvz*vvz
  //  //           + vvy*(-B - Eb - Ef - F - Nb - Neb - Nef - Nf - Nwb - Nwf - fortytwo*pe/om11 - Sb - Seb - Sef - Sf - Swb - Swf + two*vxz/om11 - Wb
  //  //           + vvz*vvz*(-B - E - Eb - Ef - F - N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - R - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf - W - Wb - Wf) - Wf)
  //  //           + vvy*vvz*(-two*B - two*Eb + two*Ef + two*F - two*Nb - two*Neb + two*Nef + two*Nf - two*Nwb + two*Nwf - two*Sb - two*Seb + two*Sef + two*Sf - two*Swb + two*Swf - two*Wb + two*Wf) + four*y/om11 + four*yyy/om11)/rho;

  //  //      kyyz=(-Nb - Neb + Nef + Nf - Nwb + Nwf - Sb - Seb + Sef + Sf - Swb + Swf + (two*Nb + two*Neb - two*Nef - two*Nf + two*Nwb - two*Nwf - two*Sb - two*Seb + two*Sef + two*Sf - two*Swb + two*Swf + eight*uyz/om11)*vvy
  //  //           + (-N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - fortytwo*pe/om11 - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf + two*vxy/om11)*vvz
  //  //           + (two*N + two*Nb + two*Ne + two*Neb + two*Nef + two*Nf + two*Nw + two*Nwb + two*Nwf - two*S - two*Sb - two*Se - two*Seb - two*Sef - two*Sf - two*Sw - two*Swb - two*Swf)*vvy*vvz
  //  //           + vvy*vvy*(-B - Eb + Ef + F - Nb - Neb + Nef + Nf - Nwb + Nwf - Sb - Seb + Sef + Sf - Swb + Swf - Wb                         
		//	 //  + vvz*(-B - E - Eb - Ef - F - N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - R - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf - W - Wb - Wf) + Wf) + four*z/om11 - four*zzz/om11)/rho;

  //  //      kxyz=(-Neb + Nef + Nwb - Nwf + Seb - Sef - Swb + Swf + (Nb + Neb - Nef - Nf + Nwb - Nwf - Sb - Seb + Sef + Sf - Swb + Swf + four*uyz/om11)*vvx
  //  //           + (N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf)*vvx*vvz + vvz*(-Ne - Neb - Nef + Nw + Nwb + Nwf + Se + Seb + Sef - Sw - Swb - Swf + four*uxy/om11
  //  //           + vvx*vvy*(-B - E - Eb - Ef - F - N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - R - S - Sb - Se - Seb - Sef - Sf - Sw - Swb - Swf - W - Wb - Wf)) + vvy*vvz*(E + Eb + Ef + Ne + Neb + Nef - Nw - Nwb - Nwf + Se + Seb + Sef - Sw - Swb - Swf - W - Wb - Wf)
  //  //           + vvx*vvy*(-B - Eb + Ef + F - Nb - Neb + Nef + Nf - Nwb + Nwf - Sb - Seb + Sef + Sf - Swb + Swf - Wb + Wf) + vvy*(Eb - Ef + Neb - Nef - Nwb + Nwf + Seb - Sef - Swb + Swf + four*uxz/om11 - Wb + Wf) + eight*xyz/om11)/rho;


  //        real kxxyy=(four*a/om11 + four*c/om11 + Ne + Neb + Nef + Nw + Nwb + Nwf + four*p/om11 + thirtytwo*pe/om11 + Se + Seb + Sef + Sw + Swb + Swf + (four*Ne + four*Neb + four*Nef - four*Nw - four*Nwb - four*Nwf - four*Se - four*Seb - four*Sef + four*Sw + four*Swb + four*Swf - sixteen*uxy/om11)*vvx*vvy
  //                       + vvy*vvy*(E + Eb + Ef + Ne + Neb + Nef + Nw + Nwb + Nwf + fortytwo*pe/om11 + Se + Seb + Sef + Sw + Swb + Swf + two*vxy/om11 + two*vxz/om11 + W + Wb + Wf)
  //                       + vvx*vvx*(N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf + fortytwo*pe/om11 + S + Sb + Se + Seb + Sef + Sf + Sw + Swb + Swf - two*vxy/om11
  //                       + (-two*N - two*Nb - two*Ne - two*Neb - two*Nef - two*Nf - two*Nw - two*Nwb - two*Nwf + two*S + two*Sb + two*Se + two*Seb + two*Sef + two*Sf + two*Sw + two*Swb + two*Swf)*vvy
  //                       + vvy*vvy*(B + E + Eb + Ef + F + N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf + R + S + Sb + Se + Seb + Sef + Sf + Sw + Swb + Swf + W + Wb + Wf))
  //                       + vvx*(-two*Ne - two*Neb - two*Nef + two*Nw + two*Nwb + two*Nwf - two*Se - two*Seb - two*Sef + two*Sw + two*Swb + two*Swf
  //                       + vvy*vvy*(-two*E - two*Eb - two*Ef - two*Ne - two*Neb - two*Nef + two*Nw + two*Nwb + two*Nwf - two*Se - two*Seb - two*Sef + two*Sw + two*Swb + two*Swf + two*W + two*Wb + two*Wf) - eight*x/om11 - eight*xxx/om11)
  //                       + vvy*(-two*Ne - two*Neb - two*Nef - two*Nw - two*Nwb - two*Nwf + two*Se + two*Seb + two*Sef + two*Sw + two*Swb + two*Swf - eight*y/om11 + eight*yyy/om11))/rho;

  //        real kxxzz=(-four*c/om11 + Eb + Ef + Neb + Nef + Nwb + Nwf + four*p/om11 + thirtytwo*pe/om11 + Seb + Sef + Swb + Swf + Wb + vvx*vvz*(-four*Eb + four*Ef - four*Neb + four*Nef + four*Nwb - four*Nwf - four*Seb + four*Sef + four*Swb - four*Swf - sixteen*uxz/om11 + four*Wb - four*Wf) + Wf
  //                       + vvz*vvz*(E + Eb + Ef + Ne + Neb + Nef + Nw + Nwb + Nwf + fortytwo*pe/om11 + Se + Seb + Sef + Sw + Swb + Swf + two*vxy/om11 + two*vxz/om11 + W + Wb + Wf)
  //                       + vvx*vvx*(B + Eb + Ef + F + Nb + Neb + Nef + Nf + Nwb + Nwf + fortytwo*pe/om11 + Sb + Seb + Sef + Sf + Swb + Swf - two*vxz/om11 + Wb
  //                       + vvz*(two*B + two*Eb - two*Ef - two*F + two*Nb + two*Neb - two*Nef - two*Nf + two*Nwb - two*Nwf + two*Sb + two*Seb - two*Sef - two*Sf + two*Swb - two*Swf + two*Wb - two*Wf) + Wf
  //                       + vvz*vvz*(B + E + Eb + Ef + F + N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf + R + S + Sb + Se + Seb + Sef + Sf + Sw + Swb + Swf + W + Wb + Wf))
  //                       + vvx*(-two*Eb - two*Ef - two*Neb - two*Nef + two*Nwb + two*Nwf - two*Seb - two*Sef + two*Swb + two*Swf + two*Wb + two*Wf  
		//				 + vvz*vvz*(-two*E - two*Eb - two*Ef - two*Ne - two*Neb - two*Nef + two*Nw + two*Nwb + two*Nwf - two*Se - two*Seb - two*Sef + two*Sw + two*Swb + two*Swf + two*W + two*Wb + two*Wf) - eight*x/om11 + eight*xxx/om11)
  //                       + vvz*(two*Eb - two*Ef + two*Neb - two*Nef + two*Nwb - two*Nwf + two*Seb - two*Sef + two*Swb - two*Swf + two*Wb - two*Wf - eight*z/om11 - eight*zzz/om11))/rho;

  //        real kyyzz=(-four*a/om11 + Nb + Neb + Nef + Nf + Nwb + Nwf + four*p/om11 + thirtytwo*pe/om11 + Sb + Seb + Sef + Sf + Swb + Swf + (-four*Nb - four*Neb + four*Nef + four*Nf - four*Nwb + four*Nwf + four*Sb + four*Seb - four*Sef - four*Sf + four*Swb - four*Swf - sixteen*uyz/om11)*vvy*vvz
  //                       + (N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf + fortytwo*pe/om11 + S + Sb + Se + Seb + Sef + Sf + Sw + Swb + Swf - two*vxy/om11)*vvz*vvz
  //                       + vvy*vvy*(B + Eb + Ef + F + Nb + Neb + Nef + Nf + Nwb + Nwf + fortytwo*pe/om11 + Sb + Seb + Sef + Sf + Swb + Swf - two*vxz/om11 + Wb
  //                       + vvz*(two*B + two*Eb - two*Ef - two*F + two*Nb + two*Neb - two*Nef - two*Nf + two*Nwb - two*Nwf + two*Sb + two*Seb - two*Sef - two*Sf + two*Swb - two*Swf + two*Wb - two*Wf) + Wf
  //                       + vvz*vvz*(B + E + Eb + Ef + F + N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf + R + S + Sb + Se + Seb + Sef + Sf + Sw + Swb + Swf + W + Wb + Wf))
  //                       + vvy*(-two*Nb - two*Neb - two*Nef - two*Nf - two*Nwb - two*Nwf + two*Sb + two*Seb + two*Sef + two*Sf + two*Swb + two*Swf
  //                       + (-two*N - two*Nb - two*Ne - two*Neb - two*Nef - two*Nf - two*Nw - two*Nwb - two*Nwf + two*S + two*Sb + two*Se + two*Seb + two*Sef + two*Sf + two*Sw + two*Swb + two*Swf)*vvz*vvz - eight*y/om11 - eight*yyy/om11)
  //                       + vvz*(two*Nb + two*Neb - two*Nef - two*Nf + two*Nwb - two*Nwf + two*Sb + two*Seb - two*Sef - two*Sf + two*Swb - two*Swf - eight*z/om11 + eight*zzz/om11))/rho;

  //        real kxxyz=(-Neb + Nef - Nwb + Nwf + Seb - Sef + Swb - Swf + eight*uxxyz/om11 + (two*Ne + two*Neb + two*Nef - two*Nw - two*Nwb - two*Nwf - two*Se - two*Seb - two*Sef + two*Sw + two*Swb + two*Swf - eight*uxy/om11)*vvx*vvz
  //                       + vvx*vvy*(-two*Eb + two*Ef - two*Neb + two*Nef + two*Nwb - two*Nwf - two*Seb + two*Sef + two*Swb - two*Swf - eight*uxz/om11 + two*Wb - two*Wf)
  //                       + vvy*vvz*(E + Eb + Ef + Ne + Neb + Nef + Nw + Nwb + Nwf + fortytwo*pe/om11 + Se + Seb + Sef + Sw + Swb + Swf + two*vxy/om11 + two*vxz/om11 + W + Wb + Wf)
  //                       + vvx*vvx*(-Nb - Neb + Nef + Nf - Nwb + Nwf + Sb + Seb - Sef - Sf + Swb - Swf - four*uyz/om11 + (-N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf + S + Sb + Se + Seb + Sef + Sf + Sw + Swb + Swf)*vvz
  //                       + vvy*(B + Eb - Ef - F + Nb + Neb - Nef - Nf + Nwb - Nwf + Sb + Seb - Sef - Sf + Swb - Swf + Wb - Wf)
  //                       + vvy*vvz*(B + E + Eb + Ef + F + N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf + R + S + Sb + Se + Seb + Sef + Sf + Sw + Swb + Swf + W + Wb + Wf))
  //                       + vvx*(two*Neb - two*Nef - two*Nwb + two*Nwf - two*Seb + two*Sef + two*Swb - two*Swf - sixteen*xyz/om11)
  //                       + vvz*(-Ne - Neb - Nef - Nw - Nwb - Nwf + Se + Seb + Sef + Sw + Swb + Swf
  //                       + vvx*vvy*(-two*E - two*Eb - two*Ef - two*Ne - two*Neb - two*Nef + two*Nw + two*Nwb + two*Nwf - two*Se - two*Seb - two*Sef + two*Sw + two*Swb + two*Swf + two*W + two*Wb + two*Wf) - four*y/om11 + four*yyy/om11)
  //                       + vvy*(Eb - Ef + Neb - Nef + Nwb - Nwf + Seb - Sef + Swb - Swf + Wb - Wf - four*z/om11 - four*zzz/om11))/rho;

  //        real kxyyz=(-Neb + Nef + Nwb - Nwf - Seb + Sef + Swb - Swf + eight*uxyyz/om11 + (-two*Nb - two*Neb + two*Nef + two*Nf - two*Nwb + two*Nwf + two*Sb + two*Seb - two*Sef - two*Sf + two*Swb - two*Swf - eight*uyz/om11)*vvx*vvy
  //                       + vvx*(N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf + fortytwo*pe/om11 + S + Sb + Se + Seb + Sef + Sf + Sw + Swb + Swf - two*vxy/om11)*vvz
  //                       + (two*Ne + two*Neb + two*Nef - two*Nw - two*Nwb - two*Nwf - two*Se - two*Seb - two*Sef + two*Sw + two*Swb + two*Swf - eight*uxy/om11)*vvy*vvz
  //                       + vvy*vvy*(-Eb + Ef - Neb + Nef + Nwb - Nwf - Seb + Sef + Swb - Swf - four*uxz/om11 + Wb - Wf
  //                       + vvz*(-E - Eb - Ef - Ne - Neb - Nef + Nw + Nwb + Nwf - Se - Seb - Sef + Sw + Swb + Swf + W + Wb + Wf)
  //                       + vvx*vvz*(B + E + Eb + Ef + F + N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf + R + S + Sb + Se + Seb + Sef + Sf + Sw + Swb + Swf + W + Wb + Wf))
  //                       + vvz*(-Ne - Neb - Nef + Nw + Nwb + Nwf - Se - Seb - Sef + Sw + Swb + Swf
  //                       + (-two*N - two*Nb - two*Ne - two*Neb - two*Nef - two*Nf - two*Nw - two*Nwb - two*Nwf + two*S + two*Sb + two*Se + two*Seb + two*Sef + two*Sf + two*Sw + two*Swb + two*Swf)*vvx*vvy - four*x/om11 - four*xxx/om11)
  //                       + vvy*(two*Neb - two*Nef - two*Nwb + two*Nwf - two*Seb + two*Sef + two*Swb - two*Swf - sixteen*xyz/om11) + vvx*(Nb + Neb - Nef - Nf + Nwb - Nwf + Sb + Seb - Sef - Sf + Swb - Swf
  //                       + vvy*vvy*(B + Eb - Ef - F + Nb + Neb - Nef - Nf + Nwb - Nwf + Sb + Seb - Sef - Sf + Swb - Swf + Wb - Wf) - four*z/om11 + four*zzz/om11))/rho;

  //        real kxyzz=(Neb + Nef - Nwb - Nwf - Seb - Sef + Swb + Swf + eight*uxyzz/om11 + (-two*Nb - two*Neb + two*Nef + two*Nf - two*Nwb + two*Nwf + two*Sb + two*Seb - two*Sef - two*Sf + two*Swb - two*Swf - eight*uyz/om11)*vvx*vvz
  //                       + vvy*vvz*(-two*Eb + two*Ef - two*Neb + two*Nef + two*Nwb - two*Nwf - two*Seb + two*Sef + two*Swb - two*Swf - eight*uxz/om11 + two*Wb - two*Wf)
  //                       + vvx*vvy*(B + Eb + Ef + F + Nb + Neb + Nef + Nf + Nwb + Nwf + fortytwo*pe/om11 + Sb + Seb + Sef + Sf + Swb + Swf - two*vxz/om11 + Wb + Wf)
  //                       + vvz*vvz*(Ne + Neb + Nef - Nw - Nwb - Nwf - Se - Seb - Sef + Sw + Swb + Swf - four*uxy/om11
  //                       + vvx*vvy*(B + E + Eb + Ef + F + N + Nb + Ne + Neb + Nef + Nf + Nw + Nwb + Nwf + R + S + Sb + Se + Seb + Sef + Sf + Sw + Swb + Swf + W + Wb + Wf))
  //                       + vvy*(-Eb - Ef - Neb - Nef + Nwb + Nwf - Seb - Sef + Swb + Swf + Wb + Wf
  //                       + vvz*vvz*(-E - Eb - Ef - Ne - Neb - Nef + Nw + Nwb + Nwf - Se - Seb - Sef + Sw + Swb + Swf + W + Wb + Wf) - four*x/om11 + four*xxx/om11)
  //                       + vvz*(two*Neb - two*Nef - two*Nwb + two*Nwf - two*Seb + two*Sef + two*Swb - two*Swf + vvx*vvy*(two*B + two*Eb - two*Ef - two*F + two*Nb + two*Neb - two*Nef - two*Nf + two*Nwb - two*Nwf + two*Sb + two*Seb - two*Sef - two*Sf + two*Swb - two*Swf + two*Wb - two*Wf) - sixteen*xyz/om11)
  //                       + vvx*(-Nb - Neb - Nef - Nf - Nwb - Nwf + Sb + Seb + Sef + Sf + Swb + Swf + (-N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf + S + Sb + Se + Seb + Sef + Sf + Sw + Swb + Swf)*vvz*vvz - four*y/om11 - four*yyy/om11))/rho;


		// xxyyzz=(om11*(c1o8)*((kxx*kyyzz+kyy*kxxzz+kzz*kxxyy+four*(kyz*kxxyz+kxz*kxyyz+kxy*kxyzz)+four*kxyz*kxyz+two*(kxyy*kxzz+kxxy*kyzz+kxxz*kyyz)-two*(kxx*kyy*kzz+two*(kxz*kxz*kyy+kyz*kyz*kxx+kxy*kxy*kzz)+eight*kxy*kxz*kyz))*rho-Neb-Nef-Nwb-Nwf-Seb-Sef-Swb-Swf // mit oder ohne rho???
		//			+two*(
		//		  +vvx*(Neb+Nef-Nwb-Nwf+Seb+Sef-Swb-Swf+eight*xyyzz/om11)
		//		  +vvy*(Neb+Nef+Nwb+Nwf-Seb-Sef-Swb-Swf+eight*xxyzz/om11)
		//		  +vvz*(-Neb+Nef-Nwb+Nwf-Seb+Sef-Swb+Swf+eight*xxyyz/om11)
		//								 )
		//		  +vx2*(four*a/om11-thirtytwo*pe/om11-Nb-Neb-Nef-Nf-Nwb-Nwf-four*p/om11-Sb-Seb-Sef-Sf-Swb-Swf)+vy2*(four*(c-p)/om11-Eb-Ef-thirtytwo*pe/om11-Neb-Nef-Nwb-Nwf-Seb-Sef-Swb-Swf-Wb-Wf)+vz2*(-four*(a+c+p)/om11-thirtytwo*pe/om11-Ne-Neb-Nef-Nw-Nwb-Nwf-Se-Seb-Sef-Sw-Swb-Swf)
		//		  +four*(vvx*vvy*(Nwb-Neb-Nef+Nwf+Seb+Sef-Swb-Swf-eight*uxyzz/om11)+vvx*vvz*(Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-eight*uxyyz/om11)+vvy*vvz*(Neb-Nef+Nwb-Nwf-Seb+Sef-Swb+Swf-eight*uxxyz/om11))
		//		  +two*(vx2*vvy*(Nb+Neb+Nef+Nf+Nwb+Nwf-Sb-Seb-Sef-Sf-Swb-Swf+four*(y+yyy)/om11)+vx2*vvz*(Nef-Nb-Neb+Nf-Nwb+Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*(z-zzz)/om11)+vy2*vvz*(Nef-Eb+Ef-Neb-Nwb+Nwf-Seb+Sef-Swb+Swf-Wb+Wf+four*(z+zzz)/om11)+vvx*vy2*(Eb+Ef+Neb+Nef-Nwb-Nwf+Seb+Sef-Swb-Swf-Wb-Wf+four*(x-xxx)/om11)+vvx*vz2*(Ne+Neb+Nef-Nw-Nwb-Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*(x+xxx)/om11)+vz2*vvy*(Ne+Neb+Nef+Nw+Nwb+Nwf-Se-Seb-Sef-Sw-Swb-Swf+four*(y-yyy)/om11))
		//		  +eight*vvx*vvy*vvz*(-Neb+Nef+Nwb-Nwf+Seb-Sef-Swb+Swf+eight*xyz/om11)
		//		  +vx2*vz2*(-fortytwo*pe/om11-N-Nb-Ne-Neb-Nef-Nf-Nw-Nwb-Nwf-S-Sb-Se-Seb-Sef-Sf-Sw-Swb-Swf+two*vxy/om11)+vx2*vy2*(-fortytwo*pe/om11-B-Eb-Ef-F-Nb-Neb-Nef-Nf-Nwb-Nwf-Sb-Seb-Sef-Sf-Swb-Swf-Wb-Wf+two*vxz/om11)+vy2*vz2*(-fortytwo*pe/om11-E-Eb-Ef-Ne-Neb-Nef-Nw-Nwb-Nwf-Se-Seb-Sef-Sw-Swb-Swf-W-Wb-Wf-two*(vxy+vxz)/om11)
		//		  +four*vvx*vvy*vvz*(vvz*(-Ne-Neb-Nef+Nw+Nwb+Nwf+Se+Seb+Sef-Sw-Swb-Swf+four*uxy/om11)+vvy*(Eb-Ef+Neb-Nef-Nwb+Nwf+Seb-Sef-Swb+Swf-Wb+Wf+four*uxz/om11)+vvx*(Nb+Neb-Nef-Nf+Nwb-Nwf-Sb-Seb+Sef+Sf-Swb+Swf+four*uyz/om11))
		//		  +two*(vx2*vy2*vvz*piz+vx2*vvy*piy*vz2+vvx*pix*vy2*vz2)-rho*vx2*vy2*vz2));

  //       (D.f[ dirE   ])[ke   ] = W-two*a+four*x-eleven*pe+vxy+vxz-four*p+four*(-xyyzz+xxyyzz)                -  c2over27  ;                                                                     
  //       (D.f[ dirW   ])[kw   ] = E-two*a-four*x-eleven*pe+vxy+vxz-four*p+four*(xyyzz+xxyyzz)                 -  c2over27  ;                                                                     
  //       (D.f[ dirN   ])[kn   ] = S-two*c+four*y-eleven*pe-vxy-four*p+four*(-xxyzz+xxyyzz)                    -  c2over27  ;
  //       (D.f[ dirS   ])[ks   ] = N-two*c-four*y-eleven*pe-vxy-four*p+four*(xxyzz+xxyyzz)                     -  c2over27  ;
  //       (D.f[ dirT   ])[kt   ] = B+two*a+two*c+four*z-eleven*pe-vxz-four*p+four*(-xxyyz+xxyyzz)                -  c2over27  ;
  //       (D.f[ dirB   ])[kb   ] = F+two*a+two*c-four*z-eleven*pe-vxz-four*p+four*(xxyyzz+xxyyz)                 -  c2over27  ;
  //       (D.f[ dirNE  ])[kne  ] = Sw+a+c-xxx+yyy-x-y-uxy+p+eight*pe+two*(-xxyyzz+xxyzz+xyyzz-uxyzz)  -  c1over54  ;
  //       (D.f[ dirSW  ])[ksw  ] = Ne+c+a+x+y+xxx-yyy-uxy+p+eight*pe+two*(-xxyyzz-xxyzz-xyyzz-uxyzz)  -  c1over54  ;
  //       (D.f[ dirSE  ])[kse  ] = Nw+c+a-xxx-yyy-x+y+uxy+p+eight*pe+two*(-xxyyzz-xxyzz+xyyzz+uxyzz)  -  c1over54  ;
  //       (D.f[ dirNW  ])[knw  ] = Se+a+c+xxx+yyy+x-y+uxy+p+eight*pe+two*(-xxyyzz+xxyzz-xyyzz+uxyzz)  -  c1over54  ;
  //       (D.f[ dirTE  ])[kte  ] = Wb-c+xxx-zzz-x-z-uxz+p+eight*pe+two*(-xxyyzz+xyyzz+xxyyz-uxyyz)    -  c1over54  ;
  //       (D.f[ dirBW  ])[kbw  ] = Ef-c-xxx+zzz+x+z-uxz+p+eight*pe+two*(-xxyyzz-xyyzz-xxyyz-uxyyz)    -  c1over54  ;
  //       (D.f[ dirBE  ])[kbe  ] = Wf-c+xxx+zzz-x+z+uxz+p+eight*pe+two*(-xxyyzz+xyyzz-xxyyz+uxyyz)    -  c1over54  ;
  //       (D.f[ dirTW  ])[ktw  ] = Eb-c-xxx-zzz+x-z+uxz+p+eight*pe+two*(-xxyyzz-xyyzz+xxyyz+uxyyz)    -  c1over54  ;
  //       (D.f[ dirTN  ])[ktn  ] = Sb-a-yyy+zzz-y-z-uyz+p+eight*pe+two*(-xxyyzz+xxyzz+xxyyz-uxxyz)    -  c1over54  ;
  //       (D.f[ dirBS  ])[kbs  ] = Nf-a+yyy-zzz+y+z-uyz+p+eight*pe+two*(-xxyyzz-xxyzz-xxyyz-uxxyz)    -  c1over54  ;
  //       (D.f[ dirBN  ])[kbn  ] = Sf-a-yyy-zzz-y+z+uyz+p+eight*pe+two*(-xxyyzz+xxyzz-xxyyz+uxxyz)    -  c1over54  ;
  //       (D.f[ dirTS  ])[kts  ] = Nb-a+yyy+zzz+y-z+uyz+p+eight*pe+two*(-xxyyzz-xxyzz+xxyyz+uxxyz)    -  c1over54  ;
  //       (D.f[ dirZERO])[kzero] = R+(twelve*p-thirty*pe-eight*xxyyzz)                                      -  c8over27  ;
  //       (D.f[ dirTNE ])[ktne ] = Swb-xyz+uxxyz+uxyyz+uxyzz-(xxyyz+xyyzz+xxyzz)+xxyyzz         -  c1over216 ;
  //       (D.f[ dirTSE ])[ktse ] = Nwb+xyz-uxxyz+uxyyz-uxyzz-(xxyyz+xyyzz-xxyzz)+xxyyzz         -  c1over216 ;
  //       (D.f[ dirBNE ])[kbne ] = Swf+xyz-uxxyz-uxyyz+uxyzz-(-xxyyz+xyyzz+xxyzz)+xxyyzz        -  c1over216 ;
  //       (D.f[ dirBSE ])[kbse ] = Nwf-xyz+uxxyz-uxyyz-uxyzz-(-xxyyz+xyyzz-xxyzz)+xxyyzz        -  c1over216 ;
  //       (D.f[ dirTNW ])[ktnw ] = Seb+xyz+uxxyz-uxyyz-uxyzz-(xxyyz-xyyzz+xxyzz)+xxyyzz         -  c1over216 ;
  //       (D.f[ dirTSW ])[ktsw ] = Neb-xyz-uxxyz-uxyyz+uxyzz-(xxyyz-xyyzz-xxyzz)+xxyyzz         -  c1over216 ;
  //       (D.f[ dirBNW ])[kbnw ] = Sef-xyz-uxxyz+uxyyz-uxyzz-(-xxyyz-xyyzz+xxyzz)+xxyyzz        -  c1over216 ;
  //       (D.f[ dirBSW ])[kbsw ] = Nef+xyz+uxxyz+uxyyz+uxyzz-(-xxyyz-xyyzz-xxyzz)+xxyyzz        -  c1over216 ;




		// 
		// //////////////////////////////////////////////////////////////////////////////////
  //       ////index
  //       ////unsigned int kzero= k;
  //       ////unsigned int ke   = k;
  //       //unsigned int kw   = neighborX[k];
  //       ////unsigned int kn   = k;
  //       //unsigned int ks   = neighborY[k];
  //       ////unsigned int kt   = k;
  //       //unsigned int kb   = neighborZ[k];
  //       //unsigned int ksw  = neighborY[kw];
  //       ////unsigned int kne  = k;
  //       ////unsigned int kse  = ks;
  //       ////unsigned int knw  = kw;
  //       //unsigned int kbw  = neighborZ[kw];
  //       ////unsigned int kte  = k;
  //       ////unsigned int kbe  = kb;
  //       ////unsigned int ktw  = kw;
  //       //unsigned int kbs  = neighborZ[ks];
  //       ////unsigned int ktn  = k;
  //       ////unsigned int kbn  = kb;
  //       ////unsigned int kts  = ks;
  //       ////unsigned int ktse = ks;
  //       ////unsigned int kbnw = kbw;
  //       ////unsigned int ktnw = kw;
  //       ////unsigned int kbse = kbs;
  //       ////unsigned int ktsw = ksw;
  //       ////unsigned int kbne = kb;
  //       ////unsigned int ktne = k;
  //       //unsigned int kbsw = neighborZ[ksw];
  //       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //       //real fE    =  (D.f[dirE   ])[k  ];//ke
  //       //real fW    =  (D.f[dirW   ])[kw ];
  //       //real fN    =  (D.f[dirN   ])[k  ];//kn
  //       //real fS    =  (D.f[dirS   ])[ks ];
  //       //real fT    =  (D.f[dirT   ])[k  ];//kt
  //       //real fB    =  (D.f[dirB   ])[kb ];
  //       //real fNE   =  (D.f[dirNE  ])[k  ];//kne
  //       //real fSW   =  (D.f[dirSW  ])[ksw];
  //       //real fSE   =  (D.f[dirSE  ])[ks ];//kse
  //       //real fNW   =  (D.f[dirNW  ])[kw ];//knw
  //       //real fTE   =  (D.f[dirTE  ])[k  ];//kte
  //       //real fBW   =  (D.f[dirBW  ])[kbw];
  //       //real fBE   =  (D.f[dirBE  ])[kb ];//kbe
  //       //real fTW   =  (D.f[dirTW  ])[kw ];//ktw
  //       //real fTN   =  (D.f[dirTN  ])[k  ];//ktn
  //       //real fBS   =  (D.f[dirBS  ])[kbs];
  //       //real fBN   =  (D.f[dirBN  ])[kb ];//kbn
  //       //real fTS   =  (D.f[dirTS  ])[ks ];//kts
  //       //real fZERO =  (D.f[dirZERO])[k  ];//kzero
  //       //real fTNE   = (D.f[dirTNE ])[k  ];//ktne
  //       //real fTSW   = (D.f[dirTSW ])[ksw];//ktsw
  //       //real fTSE   = (D.f[dirTSE ])[ks ];//ktse
  //       //real fTNW   = (D.f[dirTNW ])[kw ];//ktnw
  //       //real fBNE   = (D.f[dirBNE ])[kb ];//kbne
  //       //real fBSW   = (D.f[dirBSW ])[kbsw];
  //       //real fBSE   = (D.f[dirBSE ])[kbs];//kbse
  //       //real fBNW   = (D.f[dirBNW ])[kbw];//kbnw
  //       ////////////////////////////////////////////////////////////////////////////////////
  //       //real rho0   =  (fTNE+fBSW)+(fTSW+fBNE)+(fTSE+fBNW)+(fTNW+fBSE)+(fNE+fSW)+(fNW+fSE)+(fTE+fBW)+(fBE+fTW)+(fTN+fBS)+(fBN+fTS)+(fE+fW)+(fN+fS)+(fT+fB)+fZERO;
  //       //real rho    =  rho0 + one;
  //       //real OORho  =  one/rho;
  //       //real vx     =  OORho*((fTNE-fBSW)+(fBNE-fTSW)+(fTSE-fBNW)+(fBSE-fTNW) +(fNE-fSW)+(fSE-fNW)+(fTE-fBW)+(fBE-fTW)+(fE-fW));
  //       //real vy     =  OORho*((fTNE-fBSW)+(fBNE-fTSW)+(fBNW-fTSE)+(fTNW-fBSE) +(fNE-fSW)+(fNW-fSE)+(fTN-fBS)+(fBN-fTS)+(fN-fS));
  //       //real vz     =  OORho*((fTNE-fBSW)+(fTSW-fBNE)+(fTSE-fBNW)+(fTNW-fBSE) +(fTE-fBW)+(fTW-fBE)+(fTN-fBS)+(fTS-fBN)+(fT-fB));
  //       ////////////////////////////////////////////////////////////////////////////////////



  //       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //       //BGK
  //       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //       //real drho    =  fZERO+fE+fW+fN+fS+fT+fB+fNE+fSW+fSE+fNW+fTE+fBW+fBE+fTW+fTN+fBS+fBN+fTS+fTNE+fTSW+fTSE+fTNW+fBNE+fBSW+fBSE+fBNW;
  //       //real vx1     =  (fE -fW +fNE-fSW+fSE-fNW+fTE-fBW+fBE-fTW+ fTNE-fTSW+fTSE-fTNW+ fBNE-fBSW+fBSE-fBNW);
  //       //real vx2     =  (fN -fS +fNE-fSW-fSE+fNW+fTN-fBS+fBN-fTS+ fTNE-fTSW-fTSE+fTNW+ fBNE-fBSW-fBSE+fBNW);
  //       //real vx3     =  (fT -fB +fTE-fBW-fBE+fTW+fTN-fBS-fBN+fTS+ fTNE+fTSW+fTSE+fTNW- fBNE-fBSW-fBSE-fBNW);
  //       //real cusq    =  c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

  //       //fZERO = fZERO *(one+(-omega))-(-omega)*   c8over27*  (drho-cusq);
  //       //fE    = fE    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq);
  //       //fW    = fW    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq);
  //       //fN    = fN    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq);
  //       //fS    = fS    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq);
  //       //fT    = fT    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq);
  //       //fB    = fB    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq);
  //       //fNE   = fNE   *(one+(-omega))-(-omega)*   c1over54*  (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq);
  //       //fSW   = fSW   *(one+(-omega))-(-omega)*   c1over54*  (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq);
  //       //fSE   = fSE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq);
  //       //fNW   = fNW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq);
  //       //fTE   = fTE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq);
  //       //fBW   = fBW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq);
  //       //fBE   = fBE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq);
  //       //fTW   = fTW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq);
  //       //fTN   = fTN   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq);
  //       //fBS   = fBS   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq);
  //       //fBN   = fBN   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq);
  //       //fTS   = fTS   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq);
  //       //fTNE  = fTNE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq);
  //       //fBSW  = fBSW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq);
  //       //fBNE  = fBNE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq);
  //       //fTSW  = fTSW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq);
  //       //fTSE  = fTSE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq);
  //       //fBNW  = fBNW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq);
  //       //fBSE  = fBSE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq);
  //       //fTNW  = fTNW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq);
  //       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  // //      //////////////////////////////////////////////////////////////////////////////////
  // //      real vx2    = vx*vx;
  // //      real vy2    = vy*vy;
  // //      real vz2    = vz*vz;
  // //      real vxy    = vx*vy;
  // //      real vxz    = vx*vz;
  // //      real vyz    = vy*vz;
  // //      real vx2y   = vx*vx*vy;
  // //      real vx2z   = vx*vx*vz;
  // //      real vy2z   = vy*vy*vz;
  // //      real vxy2   = vx*vy*vy;
  // //      real vyz2   = vy*vz*vz;
  // //      real vxz2   = vx*vz*vz;
  // //      real vxyz   = vx*vy*vz;
		// //////nur für 6.
  // ////      real vx2y2  = vx*vx*vy*vy;
  // ////      real vx2z2  = vx*vx*vz*vz;
  // ////      real vy2z2  = vy*vy*vz*vz;
  // ////      real vx2yz  = vx*vx*vy*vz;
  // ////      real vxyz2  = vx*vy*vz*vz;
  // ////      real vxy2z  = vx*vy*vy*vz;
  // ////      real vx2y2z = vx*vx*vy*vy*vz;
  // ////      real vx2yz2 = vx*vx*vy*vz*vz;
  // ////      real vxy2z2 = vx*vy*vy*vz*vz;
  // ////      real vx2y2z2= vx*vx*vy*vy*vz*vz;
  // //      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // //      real mu200   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE))) + ((fNE+fSW) + (fNW+fSE)) + ((fTE+fBW) + (fTW+fBE))                           + (fE+fW)                    );
  // //      real mu020   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE))) + ((fNE+fSW) + (fNW+fSE))                           + ((fTN+fBS) + (fTS+fBN))           + (fN+fS)          );
  // //      real mu002   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE)))                           + ((fTE+fBW) + (fTW+fBE)) + ((fTN+fBS) + (fTS+fBN))                     + (fT+fB));
  // //      real mu110   = OORho * ((((fTNE+fBSW) - (fTNW+fBSE)) - ((fTSE+fBNW) - (fTSW+fBNE))) + ((fNE+fSW) - (fNW+fSE))                                                                                  );
  // //      real mu101   = OORho * ((((fTNE+fBSW) - (fTNW+fBSE)) + ((fTSE+fBNW) - (fTSW+fBNE)))                           + ((fTE+fBW) - (fTW+fBE))                                                        );
  // //      real mu011   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) - ((fTSE+fBNW) + (fTSW+fBNE)))                                                     + ((fTN+fBS) - (fTS+fBN))                              );
  // //      real mu210   = OORho * ((((fTNE-fBSW) + (fTNW-fBSE)) - ((fTSE-fBNW) + (fTSW-fBNE))) + ((fNE-fSW) + (fNW-fSE))                                                                                  );
  // //      real mu120   = OORho * ((((fTNE-fBSW) - (fTNW-fBSE)) + ((fTSE-fBNW) - (fTSW-fBNE))) + ((fNE-fSW) - (fNW-fSE))                                                                                  );
  // //      real mu102   = OORho * ((((fTNE-fBSW) - (fTNW-fBSE)) + ((fTSE-fBNW) - (fTSW-fBNE)))                           + ((fTE-fBW) - (fTW-fBE))                                                        );
  // //      real mu111   = OORho * ((((fTNE-fBSW) - (fTNW-fBSE)) - ((fTSE-fBNW) - (fTSW-fBNE)))                                                                                                            );
  // //      real mu201   = OORho * ((((fTNE-fBSW) + (fTNW-fBSE)) + ((fTSE-fBNW) + (fTSW-fBNE)))                           + ((fTE-fBW) + (fTW-fBE))                                                        );
  // //      real mu021   = OORho * ((((fTNE-fBSW) + (fTNW-fBSE)) + ((fTSE-fBNW) + (fTSW-fBNE)))                                                     + ((fTN-fBS) + (fTS-fBN))                              );
  // //      real mu012   = OORho * ((((fTNE-fBSW) + (fTNW-fBSE)) - ((fTSE-fBNW) + (fTSW-fBNE)))                                                     + ((fTN-fBS) - (fTS-fBN))                              );
  // //      real mu220   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE))) + ((fNE+fSW) + (fNW+fSE))                                                                                  );
  // //      real mu121   = OORho * ((((fTNE+fBSW) - (fTNW+fBSE)) + ((fTSE+fBNW) - (fTSW+fBNE)))                                                                                                            );
  // //      real mu202   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE)))                           + ((fTE+fBW) + (fTW+fBE))                                                        );
  // //      real mu211   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) - ((fTSE+fBNW) + (fTSW+fBNE)))                                                                                                            );
  // //      real mu112   = OORho * ((((fTNE+fBSW) - (fTNW+fBSE)) - ((fTSE+fBNW) - (fTSW+fBNE)))                                                                                                            );
  // //      real mu022   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE)))                                                     + ((fTN+fBS) + (fTS+fBN))                              );
  // //      real mu221   = OORho * ((((fTNE-fBSW) + (fTNW-fBSE)) + ((fTSE-fBNW) + (fTSW-fBNE)))                                                                                                            );
  // //      real mu122   = OORho * ((((fTNE-fBSW) - (fTNW-fBSE)) + ((fTSE-fBNW) - (fTSW-fBNE)))                                                                                                            );
  // //      real mu212   = OORho * ((((fTNE-fBSW) + (fTNW-fBSE)) - ((fTSE-fBNW) + (fTSW-fBNE)))                                                                                                            );
  // //      real mu222   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE)))                                                                                                            );
  // //      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // //      real MzXX,MzYY,MzZZ,MzXY,MzXZ,MzYZ,MzXXY,MzXYY,MzXXZ,MzXZZ,MzYYZ,MzYZZ,MzXYZ,MzXXYY,MzXXZZ,MzYYZZ,MzXXYZ,MzXYYZ,MzXYZZ,MzXXYYZ,MzXXYZZ,MzXYYZZ,MzXXYYZZ;
  // //      {
  // //         //2.
  // //         MzXX      =   mu200-vx2;
  // //         MzXY      =   mu110-vxy;
  // //         MzXZ      =   mu101-vxz;
  // //         MzYY      =   mu020-vy2;
  // //         MzYZ      =   mu011-vyz;
  // //         MzZZ      =   mu002-vz2;

  // //         real pimpmu200 = mu200 + c1o3 * OORho;
  // //         real pimpmu020 = mu020 + c1o3 * OORho;
  // //         real pimpmu002 = mu002 + c1o3 * OORho;

  // //         //3.
  // //         MzXXY     =    two*vx2y - two*vx*mu110 - vy*pimpmu200 + mu210; 
  // //         MzXXZ     =    two*vx2z - two*vx*mu101 - vz*pimpmu200 + mu201; 
  // //         MzXYY     =    two*vxy2 - pimpmu020*vx - two*vy*mu110 + mu120; 
  // //         MzXYZ     =    two*vxyz - mu011*vx-vy*mu101 -vz*mu110 + mu111;
  // //         MzXZZ     =    two*vxz2 - pimpmu002*vx - two*vz*mu101 + mu102; 
  // //         MzYYZ     =    two*vy2z - two*vy*mu011 - vz*pimpmu020 + mu021; 
  // //         MzYZZ     =    two*vyz2 - pimpmu002*vy - two*vz*mu011 + mu012; 

  // //         //4.
  // //         //MzXXYY    =   -three*vx2y2+pimpmu020*vx2+four*vxy*mu110-two*vx*mu120+vy2*pimpmu200-two*vy*mu210+mu220; 
  // //         //MzXXYZ    =   -three*vx2yz+mu011*vx2+two*vxy*mu101+two*vxz*mu110-two*vx*mu111+vyz*pimpmu200-vy*mu201-vz*mu210+mu211; 
  // //         //MzXXZZ    =   -three*vx2z2+pimpmu002*vx2+four*vxz*mu101-two*vx*mu102+vz2*pimpmu200-two*vz*mu201+mu202; 
  // //         //MzXYYZ    =   -three*vxy2z+two*vxy*mu011+vxz*pimpmu020-mu021*vx+vy2*mu101+two*vyz*mu110-two*vy*mu111-vz*mu120+mu121; 
  // //         //MzXYZZ    =   -three*vxyz2+pimpmu002*vxy+two*vxz*mu011-mu012*vx+two*vyz*mu101-vy*mu102+vz2*mu110-two*vz*mu111+mu112; 
  // //         //MzYYZZ    =   -three*vy2z2+pimpmu002*vy2+four*vyz*mu011-two*vy*mu012+vz2*pimpmu020-two*vz*mu021+mu022; 
		//	////schnell
		//	//MzXXYY    =   /*-three*vx2y2+*/pimpmu020*vx2/*+four*vxy*mu110*/-two*vx*mu120+vy2*pimpmu200-two*vy*mu210+mu220; 
  // //         MzXXYZ    =   /*-three*vx2yz+mu011*vx2+two*vxy*mu101+two*vxz*mu110-two*vx*mu111*/+vyz*pimpmu200-vy*mu201-vz*mu210+mu211; 
  // //         MzXXZZ    =   /*-three*vx2z2+*/pimpmu002*vx2/*+four*vxz*mu101*/-two*vx*mu102+vz2*pimpmu200-two*vz*mu201+mu202; 
  // //         MzXYYZ    =   /*-three*vxy2z+two*vxy*mu011*/+vxz*pimpmu020-mu021*vx/*+vy2*mu101+two*vyz*mu110-two*vy*mu111*/-vz*mu120+mu121; 
  // //         MzXYZZ    =   /*-three*vxyz2+*/pimpmu002*vxy/*+two*vxz*mu011*/-mu012*vx/*+two*vyz*mu101*/-vy*mu102/*+vz2*mu110-two*vz*mu111*/+mu112; 
  // //         MzYYZZ    =   /*-three*vy2z2+*/pimpmu002*vy2/*+four*vyz*mu011*/-two*vy*mu012+vz2*pimpmu020-two*vz*mu021+mu022; 

  // //         real pimpmu220 = mu220 + c1o9 * OORho;
  // //         real pimpmu202 = mu202 + c1o9 * OORho;
  // //         real pimpmu022 = mu022 + c1o9 * OORho;

  // //         //5.
  // //         //MzXXYYZ   =    four*(vx2y2z-vxyz*mu110+vxy*mu111)+ 
  // //         //               two*(vxz*mu120-vxy2*mu101-vx2y*mu011+vyz*mu210-vy*mu211-vx*mu121)+
  // //         //               vy2*mu201-vx2z*pimpmu020+mu021*vx2-vz*pimpmu220-vy2z*pimpmu200+mu221; 
  // //         //MzXXYZZ   =    four*(vx2yz2-vxyz*mu101+vxz*mu111)+
  // //         //               two*(vxy*mu102-vxz2*mu110-vx2z*mu011-vx*mu112-vz*mu211+vyz*mu201)+
  // //         //               vz2*mu210-vy*pimpmu202-pimpmu002*vx2y+mu012*vx2-vyz2*pimpmu200+mu212;
  // //         //MzXYYZZ   =    four*(vxy2z2-vxyz*mu011+vyz*mu111)+
  // //         //               two*(vxy*mu012-vyz2*mu110-vy*mu112+vxz*mu021-vz*mu121-vy2z*mu101)+
  // //         //               vy2*mu102+vz2*mu120-vxz2*pimpmu020-pimpmu022*vx-pimpmu002*vxy2+mu122;    
		//	////schnell
		//	//MzXXYYZ   =    /*four*(vx2y2z-vxyz*mu110+vxy*mu111)+*/ 
  // //                        two*(vxz*mu120/*-vxy2*mu101-vx2y*mu011*/+vyz*mu210-vy*mu211-vx*mu121)+
  // //                        vy2*mu201-vx2z*pimpmu020+mu021*vx2-vz*pimpmu220-vy2z*pimpmu200+mu221; 
  // //         MzXXYZZ   =    /*four*(vx2yz2-vxyz*mu101+vxz*mu111)+*/
  // //                        two*(vxy*mu102/*-vxz2*mu110-vx2z*mu011*/-vx*mu112-vz*mu211+vyz*mu201)+
  // //                        vz2*mu210-vy*pimpmu202-pimpmu002*vx2y+mu012*vx2-vyz2*pimpmu200+mu212;
  // //         MzXYYZZ   =    /*four*(vxy2z2-vxyz*mu011+vyz*mu111)+*/
  // //                        two*(vxy*mu012/*-vyz2*mu110*/-vy*mu112+vxz*mu021-vz*mu121/*-vy2z*mu101*/)+
  // //                        vy2*mu102+vz2*mu120-vxz2*pimpmu020-pimpmu022*vx-pimpmu002*vxy2+mu122;    

  // //         //6.
  // //         //MzXXYYZZ  =   mu222 + pimpmu022*vx2 - two*mu122*vx + pimpmu202*vy2 + pimpmu002*vx2*vy2 - two*mu102*vx*vy2 - two*mu212*vy - two*mu012*vx2*vy + four*mu112*vx*vy + pimpmu220*vz2 + pimpmu020*vx2*vz2 - two*mu120*vx*vz2 + 
  // //         //              pimpmu200*vy2*vz2 - five*vx2*vy2*vz2 - two*mu210*vy*vz2 + four*mu110*vx*vy*vz2 - two*mu221*vz - two*mu021*vx2*vz + four*mu121*vx*vz - two*mu201*vy2*vz + four*mu101*vx*vy2*vz + 
  // //         //              four*mu211*vy*vz + four*mu011*vx2*vy*vz - eight*mu111*vx*vy*vz;
  // //         //schnell
		//	//MzXXYYZZ  =   mu222 + pimpmu022*vx2 - two*mu122*vx + pimpmu202*vy2 /*+ pimpmu002*vx2*vy2 - two*mu102*vx*vy2*/ - two*mu212*vy /*- two*mu012*vx2*vy + four*mu112*vx*vy*/ + pimpmu220*vz2 /*+ pimpmu020*vx2*vz2 - two*mu120*vx*vz2 + 
  // //                       pimpmu200*vy2*vz2 - five*vx2*vy2*vz2 - two*mu210*vy*vz2 + four*mu110*vx*vy*vz2*/ - two*mu221*vz /*- two*mu021*vx2*vz + four*mu121*vx*vz - two*mu201*vy2*vz + four*mu101*vx*vy2*vz + 
  // //                       four*mu211*vy*vz + four*mu011*vx2*vy*vz - eight*mu111*vx*vy*vz*/;

  // //         //MzXXYYZZ  =   /*-five*vx2y2z2-eight*vxyz*mu111+*/vy2*pimpmu202/*+pimpmu002*vx2y2+vy2z2*pimpmu200+vx2z2*pimpmu020*/+vz2*pimpmu220+mu222+pimpmu022*vx2//+
  // //         //               /*four*(vx2yz*mu011+vxy2z*mu101+vxyz2*mu110+vxy*mu112+vxz*mu121+vyz*mu211)*/- 
  // //         //               two*(vy*mu212/*-vy2z*mu201-vxz2*mu120-vyz2*mu210-vxy2*mu102-vx2z*mu021-vx2y*mu012*/-vx*mu122-vz*mu221); 
  // //      }
  // //      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // //      real MXXpMYYpMZZ,MXXmMYY, MXXmMZZ, MXXYpMYZZ,MXXZpMYYZ,MXYYpMXZZ, MXXYmMYZZ,MXXZmMYYZ,MXYYmMXZZ;
  // //      {
  // //         //coll faktoren:
  // //         real w1 =  (-omega);
  // //         real w2 = -one;//(-omega);
  // //         real w3 = -(two+(-omega));//-one;
  // //         real w4 = -(two+(-omega));//-one;
  // //         real w5 = -(two+(-omega));//-one;
  // //         real w6 = -one;
  // //         real w7 = -one;
  // //         real w8 = -one;
  // //         real w9 = -one;
  // //         real w10= -one;

		//	//real wadjust;
		//	//real qudricLimit = c1o100;

  // //         ////////lin kombi bilden:
  // //         MXXpMYYpMZZ =  MzXX + MzYY + MzZZ;
  // //         MXXmMYY     =  MzXX - MzYY;
  // //         MXXmMZZ     =  MzXX - MzZZ;

  // //         MXXYpMYZZ   =  MzXXY+MzYZZ;
  // //         MXXYmMYZZ   =  MzXXY-MzYZZ;
  // //         MXXZpMYYZ   =  MzXXZ+MzYYZ;
  // //         MXXZmMYYZ   =  MzXXZ-MzYYZ;
  // //         MXYYpMXZZ   =  MzXYY+MzXZZ;
  // //         MXYYmMXZZ   =  MzXYY-MzXZZ;

  // //         real MXXYYppp    = MzXXYY + MzXXZZ + MzYYZZ;
  // //         real MXXYYpm2p   = MzXXYY - two*MzXXZZ + MzYYZZ;
  // //         real MXXYYppm2   = MzXXYY + MzXXZZ - two*MzYYZZ;

  // //         //relaxation:
  // //         MXXpMYYpMZZ -= w2*(one-OORho-MXXpMYYpMZZ);

  // //         MzXZ       *=  one+w1;
  // //         MzYZ       *=  one+w1;
  // //         MzXY       *=  one+w1;

  // //         MXXmMYY    *=  onef+w1;
  // //         MXXmMZZ    *=  onef+w1;

		//	//wadjust     =  w5-(one+w5)*abs(MzXYZ)/(abs(MzXYZ)+qudricLimit);
		//	//MzXYZ      *=  one+wadjust;
		//	////MzXYZ      *=  one;

		//	//wadjust     =  w3-(one+w3)*abs(MXXYpMYZZ)/(abs(MXXYpMYZZ)+qudricLimit);
  // //         MXXYpMYZZ  *=  one+wadjust;
  // //         //MXXYpMYZZ  *=  one;
		//	//wadjust     =  w3-(one+w3)*abs(MXXZpMYYZ)/(abs(MXXZpMYYZ)+qudricLimit);
  // //         MXXZpMYYZ  *=  one+wadjust;
  // //         //MXXZpMYYZ  *=  one;
		//	//wadjust     =  w3-(one+w3)*abs(MXYYpMXZZ)/(abs(MXYYpMXZZ)+qudricLimit);
  // //         MXYYpMXZZ  *=  one+wadjust;
  // //         //MXYYpMXZZ  *=  one;
		//	//wadjust     =  w4-(one+w4)*abs(MXXYmMYZZ)/(abs(MXXYmMYZZ)+qudricLimit);
  // //         MXXYmMYZZ  *=  one+wadjust;
  // //         //MXXYmMYZZ  *=  one;
		//	//wadjust     =  w4-(one+w4)*abs(MXXZmMYYZ)/(abs(MXXZmMYYZ)+qudricLimit);
  // //         MXXZmMYYZ  *=  one+wadjust;
  // //         //MXXZmMYYZ  *=  one;
		//	//wadjust     =  w4-(one+w4)*abs(MXYYmMXZZ)/(abs(MXYYmMXZZ)+qudricLimit);
  // //         MXYYmMXZZ  *=  one+wadjust;
  // //         //MXYYmMXZZ  *=  one;

  // //         ////////von Lin Kombis zurueck:
  // //         MzXX  =  c1o3 * (      MXXmMYY +      MXXmMZZ + MXXpMYYpMZZ + OORho);//verdachtsmomente
  // //         MzYY  =  c1o3 * (-two * MXXmMYY +      MXXmMZZ + MXXpMYYpMZZ + OORho);
  // //         MzZZ  =  c1o3 * (      MXXmMYY - two * MXXmMZZ + MXXpMYYpMZZ + OORho);

  // //         MzXXY = (MXXYmMYZZ + MXXYpMYZZ)*c1o2;
  // //         MzYZZ = (MXXYpMYZZ - MXXYmMYZZ)*c1o2;
  // //         MzXYY = (MXYYmMXZZ + MXYYpMXZZ)*c1o2;
  // //         MzXZZ = (MXYYpMXZZ - MXYYmMXZZ)*c1o2;
  // //         MzXXZ = (MXXZmMYYZ + MXXZpMYYZ)*c1o2;
  // //         MzYYZ = (MXXZpMYYZ - MXXZmMYYZ)*c1o2;

  // //         //faktorisierte atraktoren:
  // //         MXXYYppp  -=  w7*MzXX*MzYY + w7*MzXX*MzZZ     + w7*    MzZZ*MzYY - w7*MXXYYppp - w7*c1o3*OORho;//verdachtsmoment nicht erhaertet
  // //         MXXYYpm2p -=  w6*MzXX*MzYY - w6*two*MzXX*MzZZ + w6*    MzZZ*MzYY - w6*MXXYYpm2p;
  // //         MXXYYppm2 -=  w6*MzXX*MzYY + w6*    MzXX*MzZZ - w6*two*MzZZ*MzYY - w6*MXXYYppm2;
  // //         MzXXYYZZ  -= w10*MzXX*MzYY*MzZZ - w10*c1o27*OORho - w10*MzXXYYZZ; // verdachtsmoment ---10.12.12
  // //         MzXYYZ    -=  w8*MzYY*MzXZ - w8*MzXYYZ;
  // //         MzXYZZ    -=  w8*MzZZ*MzXY - w8*MzXYZZ;
  // //         MzXXYZ    -=  w8*MzXX*MzYZ - w8*MzXXYZ;

  // //         MzXXYYZ *= one+w9;
  // //         MzXXYZZ *= one+w9;
  // //         MzXYYZZ *= one+w9;

  // //         MzXXYY =  c1o3 * (MXXYYpm2p + MXXYYppm2 + MXXYYppp);
  // //         MzXXZZ =  c1o3 * (MXXYYppp - MXXYYpm2p);
  // //         MzYYZZ =  c1o3 * (MXXYYppp-MXXYYppm2);
  // //      }

  // //      //2.
  // //      mu200 = vx2 + c1o3 * (MXXmMYY + MXXmMZZ + MXXpMYYpMZZ);
  // //      mu020 = vy2 + c1o3 * (MXXmMZZ +  MXXpMYYpMZZ - 2. * MXXmMYY);
  // //      mu002 = vz2 + c1o3 * (MXXmMYY - 2. * MXXmMZZ + MXXpMYYpMZZ);
  // //      mu110 = vxy + MzXY;
  // //      mu101 = vxz + MzXZ;
  // //      mu011 = vyz + MzYZ;

  // //      //3.
  // //      mu111 = vxyz +     vx*MzYZ +     vy*MzXZ + vz*MzXY + MzXYZ;
  // //      mu210 = vx2y + two*vx*MzXY +     vy*MzXX + MzXXY;
  // //      mu120 = vxy2 +     vx*MzYY + two*vy*MzXY + MzXYY;
  // //      mu102 = vxz2 +     vx*MzZZ + two*vz*MzXZ + MzXZZ;
  // //      mu201 = vx2z + two*vx*MzXZ +     vz*MzXX + MzXXZ;
  // //      mu021 = vy2z + two*vy*MzYZ +     vz*MzYY + MzYYZ;
  // //      mu012 = vyz2 +     vy*MzZZ + two*vz*MzYZ + MzYZZ;

  // //      //4.
  // //      //mu211 = vx2yz +     vx2*MzYZ + two*vxy*MzXZ + two*vxz*MzXY + two*vx*MzXYZ +     vyz*MzXX +     vy*MzXXZ +     vz*MzXXY + MzXXYZ;
  // //      //mu121 = vxy2z + two*vxy*MzYZ +     vxz*MzYY +     vx*MzYYZ +     vy2*MzXZ + two*vyz*MzXY + two*vy*MzXYZ +     vz*MzXYY + MzXYYZ;
  // //      //mu112 = vxyz2 +     vxy*MzZZ + two*vxz*MzYZ +     vx*MzYZZ + two*vyz*MzXZ +     vy*MzXZZ +     vz2*MzXY + two*vz*MzXYZ + MzXYZZ;
  // //      //mu220 = vx2y2 +     vx2*MzYY + four*vxy*MzXY + two*vx*MzXYY +     vy2*MzXX + two*vy*MzXXY + MzXXYY ;
  // //      //mu202 = vx2z2 +     vx2*MzZZ + four*vxz*MzXZ + two*vx*MzXZZ +     vz2*MzXX + two*vz*MzXXZ + MzXXZZ;
  // //      //mu022 = vy2z2 +     vy2*MzZZ + four*vyz*MzYZ + two*vy*MzYZZ +     vz2*MzYY + two*vz*MzYYZ + MzYYZZ;
		// ////schnell
		// //mu211 = /*vx2yz +     vx2*MzYZ + two*vxy*MzXZ + two*vxz*MzXY + two*vx*MzXYZ +*/     vyz*MzXX +     vy*MzXXZ +     vz*MzXXY + MzXXYZ;
  // //      mu121 = /*vxy2z + two*vxy*MzYZ +*/     vxz*MzYY +     vx*MzYYZ /*+     vy2*MzXZ + two*vyz*MzXY + two*vy*MzXYZ*/ +     vz*MzXYY + MzXYYZ;
  // //      mu112 = /*vxyz2 +*/     vxy*MzZZ /*+ two*vxz*MzYZ*/ +     vx*MzYZZ /*+ two*vyz*MzXZ*/ +     vy*MzXZZ /*+     vz2*MzXY + two*vz*MzXYZ*/ + MzXYZZ;
  // //      mu220 = /*vx2y2 +*/     vx2*MzYY /*+ four*vxy*MzXY*/ + two*vx*MzXYY +     vy2*MzXX + two*vy*MzXXY + MzXXYY ;
  // //      mu202 = /*vx2z2 +*/     vx2*MzZZ /*+ four*vxz*MzXZ*/ + two*vx*MzXZZ +     vz2*MzXX + two*vz*MzXXZ + MzXXZZ;
  // //      mu022 = /*vy2z2 +*/     vy2*MzZZ /*+ four*vyz*MzYZ*/ + two*vy*MzYZZ +     vz2*MzYY + two*vz*MzYYZ + MzYYZZ;

  // //      MzXXYY += c1o9 * OORho;//verdachtsmomente
  // //      MzXXZZ += c1o9 * OORho;
  // //      MzYYZZ += c1o9 * OORho;

  // //      //5.
  // //      //mu221 = vx2y2z+vx2z*MzYY+vx2*MzYYZ+vy2z*MzXX+vy2*MzXXZ+vz*MzXXYY+MzXXYYZ+
  // //      //         four*(vxyz*MzXY+vxy*MzXYZ)+
  // //      //         two*(vxz*MzXYY+vx*MzXYYZ+vx2y*MzYZ+vyz*MzXXY+vy*MzXXYZ+vxy2*MzXZ);
  // //      //mu212 = vx2yz2 + vx2y*MzZZ + vx2*MzYZZ + vy*MzXXZZ + vz2*MzXXY + vyz2*MzXX+MzXXYZZ+
  // //      //         four*(vxyz*MzXZ + vxz*MzXYZ)+
  // //      //         two*(vxy*MzXZZ + vxz2*MzXY + vx*MzXYZZ + vyz*MzXXZ + vx2z*MzYZ + vz*MzXXYZ);
  // //      //mu122 = vxy2z2 + vxy2*MzZZ + MzXYYZZ + vz2*MzXYY + vy2*MzXZZ + vx*MzYYZZ + vxz2*MzYY+
  // //      //         four*(vxyz*MzYZ + vyz*MzXYZ)+
  // //      //         two*(vxy*MzYZZ + vxz*MzYYZ + vy2z*MzXZ + vyz2*MzXY + vy*MzXYZZ + vz*MzXYYZ);
		// ////schnell
		// //mu221 = /*vx2y2z+*/vx2z*MzYY+vx2*MzYYZ+vy2z*MzXX+vy2*MzXXZ+vz*MzXXYY+MzXXYYZ+
  // //               /*four*(vxyz*MzXY+vxy*MzXYZ)+*/
  // //               two*(vxz*MzXYY+vx*MzXYYZ/*+vx2y*MzYZ*/+vyz*MzXXY+vy*MzXXYZ/*+vxy2*MzXZ*/);
  // //      mu212 = /*vx2yz2 +*/ vx2y*MzZZ + vx2*MzYZZ + vy*MzXXZZ + vz2*MzXXY + vyz2*MzXX+MzXXYZZ+
  // //               /*four*(vxyz*MzXZ + vxz*MzXYZ)+*/
  // //               two*(vxy*MzXZZ /*+ vxz2*MzXY*/ + vx*MzXYZZ + vyz*MzXXZ /*+ vx2z*MzYZ*/ + vz*MzXXYZ);
  // //      mu122 = /*vxy2z2*/ + vxy2*MzZZ + MzXYYZZ + vz2*MzXYY + vy2*MzXZZ + vx*MzYYZZ + vxz2*MzYY+
  // //               /*four*(vxyz*MzYZ + vyz*MzXYZ)+*/
  // //               two*(vxy*MzYZZ + vxz*MzYYZ/* + vy2z*MzXZ + vyz2*MzXY*/ + vy*MzXYZZ + vz*MzXYYZ);

  // //      //6.
  // //      //mu222 = vx2y2z2 + vx2y2*MzZZ + vx2z2*MzYY + vx2*MzYYZZ + vy2z2*MzXX + vy2*MzXXZZ + vz2*MzXXYY + MzXXYYZZ + eight*vxyz*MzXYZ +
  // //      //   two*(vx2y*MzYZZ +  vx*MzXYYZZ + vz*MzXXYYZ + vyz2*MzXXY + vx2z*MzYYZ + vxy2*MzXZZ + vxz2*MzXYY + vy*MzXXYZZ + vy2z*MzXXZ)+ 
  // //      //   four*(vxy2z*MzXZ + vx2yz*MzYZ + vxyz2*MzXY + vxy*MzXYZZ + vxz*MzXYYZ + vyz*MzXXYZ);
  // //      //schnell
		// //mu222 = /*vx2y2z2 + vx2y2*MzZZ + vx2z2*MzYY */+ vx2*MzYYZZ /*+ vy2z2*MzXX*/ + vy2*MzXXZZ + vz2*MzXXYY + MzXXYYZZ /*+ eight*vxyz*MzXYZ*/ +
  // //         two*(/*vx2y*MzYZZ +*/  vx*MzXYYZZ + vz*MzXXYYZ /*+ vyz2*MzXXY + vx2z*MzYYZ + vxy2*MzXZZ + vxz2*MzXYY*/ + vy*MzXXYZZ/* + vy2z*MzXXZ*/);//+ 
  // //         //four*(/*vxy2z*MzXZ + vx2yz*MzYZ + vxyz2*MzXY +*/ vxy*MzXYZZ + vxz*MzXYYZ + vyz*MzXXYZ);

  // //      //(D.f[ dirE   ])[k   ] =   c1o2*rho*( mu200  - mu220 + mu222 - mu202 +  mu120 - mu122 + mu102 - vx   );   //ke
  // //      //(D.f[ dirW   ])[kw  ] =   c1o2*rho*( mu200  - mu220 + mu222 - mu202 -  mu120 + mu122 - mu102 + vx   );   
  // //      //(D.f[ dirN   ])[k   ] =   c1o2*rho*( mu210  - mu220 + mu222 - mu212 +  mu020 - mu022 + mu012 - vy   );   //kn
  // //      //(D.f[ dirS   ])[ks  ] =   c1o2*rho*(-mu210  - mu220 + mu222 + mu212 +  mu020 - mu022 - mu012 + vy   );   
  // //      //(D.f[ dirT   ])[k   ] =   c1o2*rho*(-mu221  + mu222 + mu201 - mu202 +  mu021 - mu022 + mu002 - vz   );   //kt
  // //      //(D.f[ dirB   ])[kb  ] =   c1o2*rho*( mu221  + mu222 - mu201 - mu202 -  mu021 - mu022 + mu002 + vz   );   
  // //      //(D.f[ dirNE  ])[k   ] =  c1o4*rho*(-mu210  + mu220 - mu222 + mu212 +  mu110 - mu120 + mu122 - mu112);   //kne
  // //      //(D.f[ dirSW  ])[ksw ] =  c1o4*rho*( mu210  + mu220 - mu222 - mu212 +  mu110 + mu120 - mu122 - mu112);   
  // //      //(D.f[ dirSE  ])[ks  ] =  c1o4*rho*( mu210  + mu220 - mu222 - mu212 -  mu110 - mu120 + mu122 + mu112);   //kse
  // //      //(D.f[ dirNW  ])[kw  ] =  c1o4*rho*(-mu210  + mu220 - mu222 + mu212 -  mu110 + mu120 - mu122 + mu112);   //knw
  // //      //(D.f[ dirTE  ])[k   ] =  c1o4*rho*( mu221  - mu222 - mu201 + mu202 -  mu121 + mu122 + mu101 - mu102);   //kte
  // //      //(D.f[ dirBW  ])[kbw ] =  c1o4*rho*(-mu221  - mu222 + mu201 + mu202 -  mu121 - mu122 + mu101 + mu102);   
  // //      //(D.f[ dirBE  ])[kb  ] =  c1o4*rho*(-mu221  - mu222 + mu201 + mu202 +  mu121 + mu122 - mu101 - mu102);   //kbe
  // //      //(D.f[ dirTW  ])[kw  ] =  c1o4*rho*( mu221  - mu222 - mu201 + mu202 +  mu121 - mu122 - mu101 + mu102);   //ktw
  // //      //(D.f[ dirTN  ])[k   ] =  c1o4*rho*( mu221  - mu222 - mu211 + mu212 -  mu021 + mu022 + mu011 - mu012);   //ktn
  // //      //(D.f[ dirBS  ])[kbs ] =  c1o4*rho*(-mu221  - mu222 - mu211 - mu212 +  mu021 + mu022 + mu011 + mu012);   
  // //      //(D.f[ dirBN  ])[kb  ] =  c1o4*rho*(-mu221  - mu222 + mu211 + mu212 +  mu021 + mu022 - mu011 - mu012);   //kbn
  // //      //(D.f[ dirTS  ])[ks  ] =  c1o4*rho*( mu221  - mu222 + mu211 - mu212 -  mu021 + mu022 - mu011 + mu012);   //kts
  // //      //(D.f[ dirZERO])[k   ] =       rho*(-mu200  + mu220 - mu222 + mu202 -  mu020 + mu022 - mu002        )+rho0;   //kzero
  // //      //(D.f[ dirTNE ])[k   ] = c1o8*rho*(-mu221  + mu222 + mu211 - mu212 +  mu121 - mu122 - mu111 + mu112);   //ktne
  // //      //(D.f[ dirTSE ])[ks  ] = c1o8*rho*(-mu221  + mu222 - mu211 + mu212 +  mu121 - mu122 + mu111 - mu112);   //ktse
  // //      //(D.f[ dirBNE ])[kb  ] = c1o8*rho*( mu221  + mu222 - mu211 - mu212 -  mu121 - mu122 + mu111 + mu112);   //kbne
  // //      //(D.f[ dirBSE ])[kbs ] = c1o8*rho*( mu221  + mu222 + mu211 + mu212 -  mu121 - mu122 - mu111 - mu112);   //kbse
  // //      //(D.f[ dirTNW ])[kw  ] = c1o8*rho*(-mu221  + mu222 + mu211 - mu212 -  mu121 + mu122 + mu111 - mu112);   //ktnw
  // //      //(D.f[ dirTSW ])[ksw ] = c1o8*rho*(-mu221  + mu222 - mu211 + mu212 -  mu121 + mu122 - mu111 + mu112);   //ktsw
  // //      //(D.f[ dirBNW ])[kbw ] = c1o8*rho*( mu221  + mu222 - mu211 - mu212 +  mu121 + mu122 - mu111 - mu112);   //kbnw
  // //      //(D.f[ dirBSW ])[kbsw] = c1o8*rho*( mu221  + mu222 + mu211 + mu212 +  mu121 + mu122 + mu111 + mu112);   
  // //      (D.f[ dirE   ])[k   ] =   c1o2*rho*(+ mu222 + (                          - mu220 - mu202        ) + (                         + mu200                ) + (                 - mu122) + (                                          + mu102 + mu120) + ( - vx            ) );   //ke
  // //      (D.f[ dirW   ])[kw  ] =   c1o2*rho*(+ mu222 + (                          - mu220 - mu202        ) + (                         + mu200                ) + (                 + mu122) + (                                          - mu102 - mu120) + ( + vx            ) );   
  // //      (D.f[ dirN   ])[k   ] =   c1o2*rho*(+ mu222 + (                          - mu220         - mu022) + (                                 + mu020        ) + (         - mu212        ) + (                          + mu012 + mu210                ) + (      - vy       ) );   //kn
  // //      (D.f[ dirS   ])[ks  ] =   c1o2*rho*(+ mu222 + (                          - mu220         - mu022) + (                                 + mu020        ) + (         + mu212        ) + (                          - mu012 - mu210                ) + (      + vy       ) );   
  // //      (D.f[ dirT   ])[k   ] =   c1o2*rho*(+ mu222 + (                                  - mu202 - mu022) + (                                         + mu002) + ( - mu221                ) + (         + mu201 +  mu021                                ) + (           - vz  ) );   //kt
  // //      (D.f[ dirB   ])[kb  ] =   c1o2*rho*(+ mu222 + (                                  - mu202 - mu022) + (                                         + mu002) + ( + mu221                ) + (         - mu201 -  mu021                                ) + (           + vz  ) );   
  // //      (D.f[ dirNE  ])[k   ] =  c1o4*rho*(- mu222 + (                  - mu112 + mu220                ) + (+  mu110                                        ) + (         + mu212 + mu122) + (                                  - mu210         - mu120)                       );   //kne
  // //      (D.f[ dirSW  ])[ksw ] =  c1o4*rho*(- mu222 + (                  - mu112 + mu220                ) + (+  mu110                                        ) + (         - mu212 - mu122) + (                                  + mu210         + mu120)                       );   
  // //      (D.f[ dirSE  ])[ks  ] =  c1o4*rho*(- mu222 + (                  + mu112 + mu220                ) + (-  mu110                                        ) + (         - mu212 + mu122) + (                                  + mu210         - mu120)                       );   //kse
  // //      (D.f[ dirNW  ])[kw  ] =  c1o4*rho*(- mu222 + (                  + mu112 + mu220                ) + (-  mu110                                        ) + (         + mu212 - mu122) + (                                  - mu210         + mu120)                       );   //knw
  // //      (D.f[ dirTE  ])[k   ] =  c1o4*rho*(- mu222 + (        -  mu121                  + mu202        ) + (         + mu101                                ) + ( + mu221         + mu122) + (         - mu201                          - mu102        )                       );   //kte
  // //      (D.f[ dirBW  ])[kbw ] =  c1o4*rho*(- mu222 + (        -  mu121                  + mu202        ) + (         + mu101                                ) + ( - mu221         - mu122) + (         + mu201                          + mu102        )                       );   
  // //      (D.f[ dirBE  ])[kb  ] =  c1o4*rho*(- mu222 + (        +  mu121                  + mu202        ) + (         - mu101                                ) + ( - mu221         + mu122) + (         + mu201                          - mu102        )                       );   //kbe
  // //      (D.f[ dirTW  ])[kw  ] =  c1o4*rho*(- mu222 + (        +  mu121                  + mu202        ) + (         - mu101                                ) + ( + mu221         - mu122) + (         - mu201                          + mu102        )                       );   //ktw
  // //      (D.f[ dirTN  ])[k   ] =  c1o4*rho*(- mu222 + (- mu211                                   + mu022) + (                 + mu011                        ) + ( + mu221 + mu212        ) + (                 -  mu021 - mu012                        )                       );   //ktn
  // //      (D.f[ dirBS  ])[kbs ] =  c1o4*rho*(- mu222 + (- mu211                                   + mu022) + (                 + mu011                        ) + ( - mu221 - mu212        ) + (                 +  mu021 + mu012                        )                       );   
  // //      (D.f[ dirBN  ])[kb  ] =  c1o4*rho*(- mu222 + (+ mu211                                   + mu022) + (                 - mu011                        ) + ( - mu221 + mu212        ) + (                 +  mu021 - mu012                        )                       );   //kbn
  // //      (D.f[ dirTS  ])[ks  ] =  c1o4*rho*(- mu222 + (+ mu211                                   + mu022) + (                 - mu011                        ) + ( + mu221 - mu212        ) + (                 -  mu021 + mu012                        )                       );   //kts
  // //      (D.f[ dirZERO])[k   ] =       rho*(- mu222 + (                          + mu220 + mu202 + mu022) + (                         - mu200 - mu020 - mu002)                                                                                                                  )+rho0;   //kzero
  // //      (D.f[ dirTNE ])[k   ] = c1o8*rho*(+ mu222 + (+ mu211 +  mu121 + mu112                         )                                                      + ( - mu221 - mu212 - mu122) + ( - mu111                                                 )                       );   //ktne
  // //      (D.f[ dirTSE ])[ks  ] = c1o8*rho*(+ mu222 + (- mu211 +  mu121 - mu112                         )                                                      + ( - mu221 + mu212 - mu122) + ( + mu111                                                 )                       );   //ktse
  // //      (D.f[ dirBNE ])[kb  ] = c1o8*rho*(+ mu222 + (- mu211 -  mu121 + mu112                         )                                                      + ( + mu221 - mu212 - mu122) + ( + mu111                                                 )                       );   //kbne
  // //      (D.f[ dirBSE ])[kbs ] = c1o8*rho*(+ mu222 + (+ mu211 -  mu121 - mu112                         )                                                      + ( + mu221 + mu212 - mu122) + ( - mu111                                                 )                       );   //kbse
  // //      (D.f[ dirTNW ])[kw  ] = c1o8*rho*(+ mu222 + (+ mu211 -  mu121 - mu112                         )                                                      + ( - mu221 - mu212 + mu122) + ( + mu111                                                 )                       );   //ktnw
  // //      (D.f[ dirTSW ])[ksw ] = c1o8*rho*(+ mu222 + (- mu211 -  mu121 + mu112                         )                                                      + ( - mu221 + mu212 + mu122) + ( - mu111                                                 )                       );   //ktsw
  // //      (D.f[ dirBNW ])[kbw ] = c1o8*rho*(+ mu222 + (- mu211 +  mu121 - mu112                         )                                                      + ( + mu221 - mu212 + mu122) + ( - mu111                                                 )                       );   //kbnw
  // //      (D.f[ dirBSW ])[kbsw] = c1o8*rho*(+ mu222 + (+ mu211 +  mu121 + mu112                         )                                                      + ( + mu221 + mu212 + mu122) + ( + mu111                                                 )                       );   
  // //                                                                                                                                                                                                                                                             
  //                                                                                                                                                                                                                                                              
  //       //////////////////////////////////////////////////////////////////////////                                                                                                                                                                             
  //       //BGK                                                                                                 
  //       //////////////////////////////////////////////////////////////////////////                            
  //       //(D.f[ dirE   ])[k   ] = fW    ;                                                                     
  //       //(D.f[ dirW   ])[kw  ] = fE    ;                                                                     
  //       //(D.f[ dirN   ])[k   ] = fS    ;
  //       //(D.f[ dirS   ])[ks  ] = fN    ;
  //       //(D.f[ dirT   ])[k   ] = fB    ;
  //       //(D.f[ dirB   ])[kb  ] = fT    ;
  //       //(D.f[ dirNE  ])[k   ] = fSW   ;
  //       //(D.f[ dirSW  ])[ksw ] = fNE   ;
  //       //(D.f[ dirSE  ])[ks  ] = fNW   ;
  //       //(D.f[ dirNW  ])[kw  ] = fSE   ;
  //       //(D.f[ dirTE  ])[k   ] = fBW   ;
  //       //(D.f[ dirBW  ])[kbw ] = fTE   ;
  //       //(D.f[ dirBE  ])[kb  ] = fTW   ;
  //       //(D.f[ dirTW  ])[kw  ] = fBE   ;
  //       //(D.f[ dirTN  ])[k   ] = fBS   ;
  //       //(D.f[ dirBS  ])[kbs ] = fTN   ;
  //       //(D.f[ dirBN  ])[kb  ] = fTS   ;
  //       //(D.f[ dirTS  ])[ks  ] = fBN   ;
  //       //(D.f[ dirZERO])[k   ] = fZERO ;
  //       //(D.f[ dirTNE ])[k   ] = fBSW  ;
  //       //(D.f[ dirTSE ])[ks  ] = fBNW  ;
  //       //(D.f[ dirBNE ])[kb  ] = fTSW  ;
  //       //(D.f[ dirBSE ])[kbs ] = fTNW  ;
  //       //(D.f[ dirTNW ])[kw  ] = fBSE  ;
  //       //(D.f[ dirTSW ])[ksw ] = fBNE  ;
  //       //(D.f[ dirBNW ])[kbw ] = fTSE  ;
  //       //(D.f[ dirBSW ])[kbsw] = fTNE  ;
  //    }                                                                                                                    
  // }
}
////////////////////////////////////////////////////////////////////////////////

