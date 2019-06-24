/* Device code */
#include "LBM/D3Q27.h"
#include "math.h"
#include "GPU/constant.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Cascade_SP_27(     real omega,
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
			real qudricLimit = 0.01f;
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
			// Cascade
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = omega;
			real OxyyPxzz  = one;//two-omega;//eight*(two-omega)/(eight -omega);//one;//omega;//two-omega;//
			real OxyyMxzz  = one;//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
			real O4        = one;
			real O5        = one;
			real O6        = one;

			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//incl. correction
			{
				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
				real dyuy = dxux + omega * c3o2 * mxxMyy;
				real dzuz = dxux + omega * c3o2 * mxxMzz;

				//relax
				mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
				mxxMyy    += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
				mxxMzz    += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
// 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 			//no correction
// 			mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);
// 			mxxMyy    += -(-omega) * (-mxxMyy);
// 			mxxMzz    += -(-omega) * (-mxxMzz);
// 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
			mfacc = mfaaa * c1o9 * O4 + (one - O4) * mfacc; 
			mfcac = mfaaa * c1o9 * O4 + (one - O4) * mfcac; 
			mfcca = mfaaa * c1o9 * O4 + (one - O4) * mfcca; 
			
			mfbbc += O4 * (-mfbbc); 
			mfbcb += O4 * (-mfbcb); 
			mfcbb += O4 * (-mfcbb); 
					
			//5.
			mfbcc += O5 * (-mfbcc);
			mfcbc += O5 * (-mfcbc);
			mfccb += O5 * (-mfccb);

			//6.
			mfccc = mfaaa * c1o27 * O6 + (one - O6) * mfccc;
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
extern "C" __global__ void LB_Kernel_Casc_Comp_SP_27(      real omega,
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
		 real rho=Nw+W+Sw+S+Se+E+Ne+N+R+Nf+Nb+Sf+Sb+Ef+Eb+Wf+Wb+Nwf+Nwb+Nef+Neb+Swf+Swb+Sef+Seb+F+B+one;// ACHTUNG ne EINS !!!!!!!!
		 real pix=(Ne+E+Se+Ef+Eb-Nw-W-Sw-Wf-Wb+Nef+Neb+Sef+Seb-Nwf-Nwb-Swf-Swb);
		 real piy=(Ne+N+Nw+Nf+Nb-Se-S-Sw-Sf-Sb+Nef+Neb+Nwf+Nwb-Sef-Seb-Swf-Swb);
		 real piz=(Nf+Sf+Wf+Ef+F-Nb-Sb-Wb-Eb-B+Nef+Nwf+Sef+Swf-Neb-Nwb-Seb-Swb);
		 real vvx=pix/rho;
		 real vvy=piy/rho;
		 real vvz=piz/rho;
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
		 // Cascaded simple
		 ////////////////////////////////////////////////////////////////////////////////////
		 real OxxPyyPzz = one;
		 real OxyyPxzz  = one;
		 real OxyyMxzz  = one;
		 real O4        = one;
		 real O5        = one;
		 real O6        = one;


		 //2.
		 // linear combinations
		 real mxxPyyPzz = mfcaa + mfaca + mfaac;
		 real mxxMyy    = mfcaa - mfaca;
		 real mxxMzz	   = mfcaa - mfaac;

		 //relax
		 mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);
		 mxxMyy    += -(-omega) * (-mxxMyy);
		 mxxMzz    += -(-omega) * (-mxxMzz);
		 mfabb     += -(-omega) * (-mfabb);
		 mfbab     += -(-omega) * (-mfbab);
		 mfbba     += -(-omega) * (-mfbba);

		 // linear combinations back
		 mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
		 mfaca = c1o3 * (-two*  mxxMyy +      mxxMzz + mxxPyyPzz);
		 mfaac = c1o3 * (       mxxMyy - two* mxxMzz + mxxPyyPzz);

		 ////Test
		 //mfabb += -(-omega) * (-mfabb);
		 //mfbab += -(-omega) * (-mfbab);
		 //mfbba += -(-omega) * (-mfbba);

		 //mfcaa += -(-omega) * (c1o3 * mfaaa - mfcaa);
		 //mfaca += -(-omega) * (c1o3 * mfaaa - mfaca);
		 //mfaac += -(-omega) * (c1o3 * mfaaa - mfaac);



		 //3.
		 // linear combinations
		 real mxxyPyzz = mfcba + mfabc;
		 real mxxyMyzz = mfcba - mfabc;
		 
		 real mxxzPyyz = mfcab + mfacb;
		 real mxxzMyyz = mfcab - mfacb;

		 real mxyyPxzz = mfbca + mfbac;
		 real mxyyMxzz = mfbca - mfbac;

		 //relax
		 mfbbb     += OxyyMxzz * (-mfbbb);
		 mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
		 mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
		 mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
		 mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
		 mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
		 mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);

		 // linear combinations back
		 mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
		 mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
		 mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
		 mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
		 mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
		 mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

		 ////Test
		 //mfabc += -(-omega) * (-mfabc);
		 //mfbac += -(-omega) * (-mfbac);
		 //
		 //mfacb += -(-omega) * (-mfacb);
		 //mfbca += -(-omega) * (-mfbca);

		 //mfcab += -(-omega) * (-mfcab);
		 //mfcba += -(-omega) * (-mfcba);

		 //mfbbb += -(-omega) * (-mfbbb);


		 //4.
		 mfacc += O4 * (c1o9 * mfaaa - mfacc);
		 mfcac += O4 * (c1o9 * mfaaa - mfcac);
		 mfcca += O4 * (c1o9 * mfaaa - mfcca);

		 mfbbc += O4 * (-mfbbc);
		 mfbcb += O4 * (-mfbcb);
		 mfcbb += O4 * (-mfcbb);

		 //5.
		 mfbcc += O5 * (-mfbcc);
		 mfcbc += O5 * (-mfcbc);
		 mfccb += O5 * (-mfccb);

		 //6.
		 mfccc += O6 * (c1o27 * mfaaa - mfccc);
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
extern "C" __global__ void LB_Kernel_Casc_SP_MS_OHM_27(  real omega,
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
      BC        =   bcMatD[k];

      if( (BC != GEO_SOLID) && (BC != GEO_VOID))
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
         real fE    =  (D.f[dirE   ])[k  ];//ke
         real fW    =  (D.f[dirW   ])[kw ];
         real fN    =  (D.f[dirN   ])[k  ];//kn
         real fS    =  (D.f[dirS   ])[ks ];
         real fT    =  (D.f[dirT   ])[k  ];//kt
         real fB    =  (D.f[dirB   ])[kb ];
         real fNE   =  (D.f[dirNE  ])[k  ];//kne
         real fSW   =  (D.f[dirSW  ])[ksw];
         real fSE   =  (D.f[dirSE  ])[ks ];//kse
         real fNW   =  (D.f[dirNW  ])[kw ];//knw
         real fTE   =  (D.f[dirTE  ])[k  ];//kte
         real fBW   =  (D.f[dirBW  ])[kbw];
         real fBE   =  (D.f[dirBE  ])[kb ];//kbe
         real fTW   =  (D.f[dirTW  ])[kw ];//ktw
         real fTN   =  (D.f[dirTN  ])[k  ];//ktn
         real fBS   =  (D.f[dirBS  ])[kbs];
         real fBN   =  (D.f[dirBN  ])[kb ];//kbn
         real fTS   =  (D.f[dirTS  ])[ks ];//kts
         real fZERO =  (D.f[dirZERO])[k  ];//kzero
         real fTNE   = (D.f[dirTNE ])[k  ];//ktne
         real fTSW   = (D.f[dirTSW ])[ksw];//ktsw
         real fTSE   = (D.f[dirTSE ])[ks ];//ktse
         real fTNW   = (D.f[dirTNW ])[kw ];//ktnw
         real fBNE   = (D.f[dirBNE ])[kb ];//kbne
         real fBSW   = (D.f[dirBSW ])[kbsw];
         real fBSE   = (D.f[dirBSE ])[kbs];//kbse
         real fBNW   = (D.f[dirBNW ])[kbw];//kbnw
         ////////////////////////////////////////////////////////////////////////////////
         real rho0   =  (fTNE+fBSW)+(fTSW+fBNE)+(fTSE+fBNW)+(fTNW+fBSE)+(fNE+fSW)+(fNW+fSE)+(fTE+fBW)+(fBE+fTW)+(fTN+fBS)+(fBN+fTS)+(fE+fW)+(fN+fS)+(fT+fB)+fZERO;
         real rho    =  rho0 + one;
         real OORho  =  one/rho;
         real vx     =  OORho*((fTNE-fBSW)+(fBNE-fTSW)+(fTSE-fBNW)+(fBSE-fTNW) +(fNE-fSW)+(fSE-fNW)+(fTE-fBW)+(fBE-fTW)+(fE-fW));
         real vy     =  OORho*((fTNE-fBSW)+(fBNE-fTSW)+(fBNW-fTSE)+(fTNW-fBSE) +(fNE-fSW)+(fNW-fSE)+(fTN-fBS)+(fBN-fTS)+(fN-fS));
         real vz     =  OORho*((fTNE-fBSW)+(fTSW-fBNE)+(fTSE-fBNW)+(fTNW-fBSE) +(fTE-fBW)+(fTW-fBE)+(fTN-fBS)+(fTS-fBN)+(fT-fB));
         ////////////////////////////////////////////////////////////////////////////////



         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         ////BGK incomp
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //real drho    =  fZERO+fE+fW+fN+fS+fT+fB+fNE+fSW+fSE+fNW+fTE+fBW+fBE+fTW+fTN+fBS+fBN+fTS+fTNE+fTSW+fTSE+fTNW+fBNE+fBSW+fBSE+fBNW;
         //real vx1     =  (fE -fW +fNE-fSW+fSE-fNW+fTE-fBW+fBE-fTW+ fTNE-fTSW+fTSE-fTNW+ fBNE-fBSW+fBSE-fBNW);
         //real vx2     =  (fN -fS +fNE-fSW-fSE+fNW+fTN-fBS+fBN-fTS+ fTNE-fTSW-fTSE+fTNW+ fBNE-fBSW-fBSE+fBNW);
         //real vx3     =  (fT -fB +fTE-fBW-fBE+fTW+fTN-fBS-fBN+fTS+ fTNE+fTSW+fTSE+fTNW- fBNE-fBSW-fBSE-fBNW);
         //real cusq    =  c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         //fZERO = fZERO *(one+(-omega))-(-omega)*   c8over27*  (drho-cusq);
         //fE    = fE    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq);
         //fW    = fW    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq);
         //fN    = fN    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq);
         //fS    = fS    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq);
         //fT    = fT    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq);
         //fB    = fB    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq);
         //fNE   = fNE   *(one+(-omega))-(-omega)*   c1over54*  (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq);
         //fSW   = fSW   *(one+(-omega))-(-omega)*   c1over54*  (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq);
         //fSE   = fSE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq);
         //fNW   = fNW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq);
         //fTE   = fTE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq);
         //fBW   = fBW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq);
         //fBE   = fBE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq);
         //fTW   = fTW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq);
         //fTN   = fTN   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq);
         //fBS   = fBS   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq);
         //fBN   = fBN   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq);
         //fTS   = fTS   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq);
         //fTNE  = fTNE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq);
         //fBSW  = fBSW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq);
         //fBNE  = fBNE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq);
         //fTSW  = fTSW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq);
         //fTSE  = fTSE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq);
         //fBNW  = fBNW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq);
         //fBSE  = fBSE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq);
         //fTNW  = fTNW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq);
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //      //BGK comp
   //      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //      real drho    =  fZERO+fE+fW+fN+fS+fT+fB+fNE+fSW+fSE+fNW+fTE+fBW+fBE+fTW+fTN+fBS+fBN+fTS+fTNE+fTSW+fTSE+fTNW+fBNE+fBSW+fBSE+fBNW;
   //      real rho     =  one + drho;
   //      real OoRho   =  one / rho;
   //      real vx1     =  OoRho * (fE -fW +fNE-fSW+fSE-fNW+fTE-fBW+fBE-fTW+ fTNE-fTSW+fTSE-fTNW+ fBNE-fBSW+fBSE-fBNW);
   //      real vx2     =  OoRho * (fN -fS +fNE-fSW-fSE+fNW+fTN-fBS+fBN-fTS+ fTNE-fTSW-fTSE+fTNW+ fBNE-fBSW-fBSE+fBNW);
   //      real vx3     =  OoRho * (fT -fB +fTE-fBW-fBE+fTW+fTN-fBS-fBN+fTS+ fTNE+fTSW+fTSE+fTNW- fBNE-fBSW-fBSE-fBNW);
   //      real cusq    =  c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

   //      fZERO = fZERO *(one+(-omega))-(-omega)*   c8over27*  (drho-rho * cusq);
   //      fE    = fE    *(one+(-omega))-(-omega)*   c2over27*  (drho+rho * (three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq));
   //      fW    = fW    *(one+(-omega))-(-omega)*   c2over27*  (drho+rho * (three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq));
   //      fN    = fN    *(one+(-omega))-(-omega)*   c2over27*  (drho+rho * (three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq));
   //      fS    = fS    *(one+(-omega))-(-omega)*   c2over27*  (drho+rho * (three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq));
   //      fT    = fT    *(one+(-omega))-(-omega)*   c2over27*  (drho+rho * (three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq));
   //      fB    = fB    *(one+(-omega))-(-omega)*   c2over27*  (drho+rho * (three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq));
   //      fNE   = fNE   *(one+(-omega))-(-omega)*   c1over54*  (drho+rho * (three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq));
   //      fSW   = fSW   *(one+(-omega))-(-omega)*   c1over54*  (drho+rho * (three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
   //      fSE   = fSE   *(one+(-omega))-(-omega)*    c1over54* (drho+rho * (three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq));
   //      fNW   = fNW   *(one+(-omega))-(-omega)*    c1over54* (drho+rho * (three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
   //      fTE   = fTE   *(one+(-omega))-(-omega)*    c1over54* (drho+rho * (three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq));
   //      fBW   = fBW   *(one+(-omega))-(-omega)*    c1over54* (drho+rho * (three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
   //      fBE   = fBE   *(one+(-omega))-(-omega)*    c1over54* (drho+rho * (three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq));
   //      fTW   = fTW   *(one+(-omega))-(-omega)*    c1over54* (drho+rho * (three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
   //      fTN   = fTN   *(one+(-omega))-(-omega)*    c1over54* (drho+rho * (three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq));
   //      fBS   = fBS   *(one+(-omega))-(-omega)*    c1over54* (drho+rho * (three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
   //      fBN   = fBN   *(one+(-omega))-(-omega)*    c1over54* (drho+rho * (three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq));
   //      fTS   = fTS   *(one+(-omega))-(-omega)*    c1over54* (drho+rho * (three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
   //      fTNE  = fTNE  *(one+(-omega))-(-omega)*    c1over216*(drho+rho * (three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
   //      fBSW  = fBSW  *(one+(-omega))-(-omega)*    c1over216*(drho+rho * (three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
   //      fBNE  = fBNE  *(one+(-omega))-(-omega)*    c1over216*(drho+rho * (three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
   //      fTSW  = fTSW  *(one+(-omega))-(-omega)*    c1over216*(drho+rho * (three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
   //      fTSE  = fTSE  *(one+(-omega))-(-omega)*    c1over216*(drho+rho * (three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
   //      fBNW  = fBNW  *(one+(-omega))-(-omega)*    c1over216*(drho+rho * (three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
   //      fBSE  = fBSE  *(one+(-omega))-(-omega)*    c1over216*(drho+rho * (three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
   //      fTNW  = fTNW  *(one+(-omega))-(-omega)*    c1over216*(drho+rho * (three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));
   //      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



         ////////////////////////////////////////////////////////////////////////////////
         real vx2    = vx*vx;
         real vy2    = vy*vy;
         real vz2    = vz*vz;
         real vxy    = vx*vy;
         real vxz    = vx*vz;
         real vyz    = vy*vz;
         real vx2y   = vx*vx*vy;
         real vx2z   = vx*vx*vz;
         real vy2z   = vy*vy*vz;
         real vxy2   = vx*vy*vy;
         real vyz2   = vy*vz*vz;
         real vxz2   = vx*vz*vz;
         real vxyz   = vx*vy*vz;
         //real vx2y2  = vx*vx*vy*vy;
         //real vx2z2  = vx*vx*vz*vz;
         //real vy2z2  = vy*vy*vz*vz;
         //real vx2yz  = vx*vx*vy*vz;
         //real vxyz2  = vx*vy*vz*vz;
         //real vxy2z  = vx*vy*vy*vz;
         //real vx2y2z = vx*vx*vy*vy*vz;
         //real vx2yz2 = vx*vx*vy*vz*vz;
         //real vxy2z2 = vx*vy*vy*vz*vz;
         //real vx2y2z2= vx*vx*vy*vy*vz*vz;
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //real mu200   =(fE+fW+fNE+fSW+fSE+fNW+fTE+fBW+fBE+fTW+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         //real mu020   =(fN+fS+fNE+fSW+fSE+fNW+fTN+fBS+fBN+fTS+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         //real mu002   =(fT+fB+fTE+fBW+fBE+fTW+fTN+fBS+fBN+fTS+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         //real mu110   =(fNE-fSE+fSW-fNW+fTNE-fTSE+fBNE-fBSE+fTSW-fTNW+fBSW-fBNW) * OORho;
         //real mu101   =(fTE+fBW-fBE-fTW+fTNE-fBNE+fTSE-fBSE-fTNW+fBNW-fTSW+fBSW) * OORho;
         //real mu011   =(fTN+fBS-fBN-fTS+fTNE-fBNE-fTSE+fBSE+fTNW-fBNW-fTSW+fBSW) * OORho;
         //real mu210   =(fNE-fSW-fSE+fNW+fTNE+fBNE-fTSE-fBSE+fTNW+fBNW-fTSW-fBSW) * OORho;
         //real mu120   =(fNE-fSW+fSE-fNW+fTNE-fTNW+fBNE-fBNW+fTSE-fTSW+fBSE-fBSW) * OORho;
         //real mu102   =(fTE-fBW+fBE-fTW+fTNE-fTNW+fBNE-fBNW+fTSE-fTSW+fBSE-fBSW) * OORho;
         //real mu111   =(fTNE-fBNE-fTSE+fBSE-fTNW+fBNW+fTSW-fBSW) * OORho;
         //real mu201   =(fTE-fBW-fBE+fTW+fTNE-fBNE+fTSE-fBSE+fTNW-fBNW+fTSW-fBSW) * OORho;
         //real mu021   =(fTN-fBS-fBN+fTS+fTNE-fBNE+fTSE-fBSE+fTNW-fBNW+fTSW-fBSW) * OORho;
         //real mu012   =(fTN-fBS+fBN-fTS+fTNE+fBNE-fTSE-fBSE+fTNW+fBNW-fTSW-fBSW) * OORho;
         //real mu220   =(fNE+fSW+fSE+fNW+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         //real mu121   =(fTNE-fBNE+fTSE-fBSE-fTNW+fBNW-fTSW+fBSW) * OORho;
         //real mu202   =(fTE+fBW+fBE+fTW+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         //real mu211   =(fTNE-fBNE-fTSE+fBSE+fTNW-fBNW-fTSW+fBSW) * OORho;
         //real mu112   =(fTNE+fBNE-fTSE-fBSE-fTNW-fBNW+fTSW+fBSW) * OORho;
         //real mu022   =(fTN+fBS+fBN+fTS+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         //real mu221   =(fTNE-fBNE+fTSE-fBSE+fTNW-fBNW+fTSW-fBSW) * OORho;
         //real mu122   =(fTNE-fTNW+fBNE-fBNW+fTSE-fTSW+fBSE-fBSW) * OORho;
         //real mu212   =(fTNE+fBNE-fTSE-fBSE+fTNW+fBNW-fTSW-fBSW) * OORho;
         //real mu222   =(fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         real mu200   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE))) + ((fNE+fSW) + (fNW+fSE)) + ((fTE+fBW) + (fTW+fBE))                           + (fE+fW)                    );
         real mu020   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE))) + ((fNE+fSW) + (fNW+fSE))                           + ((fTN+fBS) + (fTS+fBN))           + (fN+fS)          );
         real mu002   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE)))                           + ((fTE+fBW) + (fTW+fBE)) + ((fTN+fBS) + (fTS+fBN))                     + (fT+fB));
         real mu110   = OORho * ((((fTNE+fBSW) - (fTNW+fBSE)) - ((fTSE+fBNW) - (fTSW+fBNE))) + ((fNE+fSW) - (fNW+fSE))                                                                                  );
         real mu101   = OORho * ((((fTNE+fBSW) - (fTNW+fBSE)) + ((fTSE+fBNW) - (fTSW+fBNE)))                           + ((fTE+fBW) - (fTW+fBE))                                                        );
         real mu011   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) - ((fTSE+fBNW) + (fTSW+fBNE)))                                                     + ((fTN+fBS) - (fTS+fBN))                              );
         real mu210   = OORho * ((((fTNE-fBSW) + (fTNW-fBSE)) - ((fTSE-fBNW) + (fTSW-fBNE))) + ((fNE-fSW) + (fNW-fSE))                                                                                  );
         real mu120   = OORho * ((((fTNE-fBSW) - (fTNW-fBSE)) + ((fTSE-fBNW) - (fTSW-fBNE))) + ((fNE-fSW) + (fNW-fSE))                                                                                  );
         real mu102   = OORho * ((((fTNE-fBSW) - (fTNW-fBSE)) + ((fTSE-fBNW) - (fTSW-fBNE)))                           + ((fTE-fBW) - (fTW-fBE))                                                        );
         real mu111   = OORho * ((((fTNE-fBSW) - (fTNW-fBSE)) - ((fTSE-fBNW) - (fTSW-fBNE)))                                                                                                            );
         real mu201   = OORho * ((((fTNE-fBSW) + (fTNW-fBSE)) + ((fTSE-fBNW) + (fTSW-fBNE)))                           + ((fTE-fBW) + (fTW-fBE))                                                        );
         real mu021   = OORho * ((((fTNE-fBSW) + (fTNW-fBSE)) + ((fTSE-fBNW) + (fTSW-fBNE)))                                                     + ((fTN-fBS) + (fTS-fBN))                              );
         real mu012   = OORho * ((((fTNE-fBSW) + (fTNW-fBSE)) - ((fTSE-fBNW) + (fTSW-fBNE)))                                                     + ((fTN-fBS) - (fTS-fBN))                              );
         real mu220   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE))) + ((fNE+fSW) + (fNW+fSE))                                                                                  );
         real mu121   = OORho * ((((fTNE+fBSW) - (fTNW+fBSE)) + ((fTSE+fBNW) - (fTSW+fBNE)))                                                                                                            );
         real mu202   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE)))                           + ((fTE+fBW) + (fTW+fBE))                                                        );
         real mu211   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) - ((fTSE+fBNW) + (fTSW+fBNE)))                                                                                                            );
         real mu112   = OORho * ((((fTNE+fBSW) - (fTNW+fBSE)) - ((fTSE+fBNW) - (fTSW+fBNE)))                                                                                                            );
         real mu022   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE)))                                                     + ((fTN+fBS) + (fTS+fBN))                              );
         real mu221   = OORho * ((((fTNE-fBSW) + (fTNW-fBSE)) + ((fTSE-fBNW) + (fTSW-fBNE)))                                                                                                            );
         real mu122   = OORho * ((((fTNE-fBSW) - (fTNW-fBSE)) + ((fTSE-fBNW) - (fTSW-fBNE)))                                                                                                            );
         real mu212   = OORho * ((((fTNE-fBSW) + (fTNW-fBSE)) - ((fTSE-fBNW) + (fTSW-fBNE)))                                                                                                            );
         real mu222   = OORho * ((((fTNE+fBSW) + (fTNW+fBSE)) + ((fTSE+fBNW) + (fTSW+fBNE)))                                                                                                            );
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         real MzXX,MzYY,MzZZ,MzXY,MzXZ,MzYZ,MzXXY,MzXYY,MzXXZ,MzXZZ,MzYYZ,MzYZZ,MzXYZ,MzXXYY,MzXXZZ,MzYYZZ,MzXXYZ,MzXYYZ,MzXYZZ,MzXXYYZ,MzXXYZZ,MzXYYZZ,MzXXYYZZ;
         {
            //2.
            MzXX      =   mu200-vx2;
            MzXY      =   mu110-vxy;
            MzXZ      =   mu101-vxz;
            MzYY      =   mu020-vy2;
            MzYZ      =   mu011-vyz;
            MzZZ      =   mu002-vz2;

            real pimpmu200 = mu200 + c1o3 * OORho;
            real pimpmu020 = mu020 + c1o3 * OORho;
            real pimpmu002 = mu002 + c1o3 * OORho;

            //3.
            MzXXY     =    two*vx2y - two*vx*mu110 - vy*pimpmu200 + mu210; 
            MzXXZ     =    two*vx2z - two*vx*mu101 - vz*pimpmu200 + mu201; 
            MzXYY     =    two*vxy2 - pimpmu020*vx - two*vy*mu110 + mu120; 
            MzXYZ     =    two*vxyz - mu011*vx-vy*mu101 -vz*mu110 + mu111;
            MzXZZ     =    two*vxz2 - pimpmu002*vx - two*vz*mu101 + mu102; 
            MzYYZ     =    two*vy2z - two*vy*mu011 - vz*pimpmu020 + mu021; 
            MzYZZ     =    two*vyz2 - pimpmu002*vy - two*vz*mu011 + mu012; 

            //4.
            MzXXYY    =   /*-three*vx2y2+*/pimpmu020*vx2/*+four*vxy*mu110*/-two*vx*mu120+vy2*pimpmu200-two*vy*mu210+mu220; 
            MzXXYZ    =   /*-three*vx2yz+mu011*vx2+two*vxy*mu101+two*vxz*mu110-two*vx*mu111*/+vyz*pimpmu200-vy*mu201-vz*mu210+mu211; 
            MzXXZZ    =   /*-three*vx2z2+*/pimpmu002*vx2/*+four*vxz*mu101*/-two*vx*mu102+vz2*pimpmu200-two*vz*mu201+mu202; 
            MzXYYZ    =   /*-three*vxy2z+two*vxy*mu011*/+vxz*pimpmu020-mu021*vx/*+vy2*mu101+two*vyz*mu110-two*vy*mu111*/-vz*mu120+mu121; 
            MzXYZZ    =   /*-three*vxyz2+*/pimpmu002*vxy/*+two*vxz*mu011*/-mu012*vx/*+two*vyz*mu101*/-vy*mu102/*+vz2*mu110-two*vz*mu111*/+mu112; 
            MzYYZZ    =   /*-three*vy2z2+*/pimpmu002*vy2/*+four*vyz*mu011*/-two*vy*mu012+vz2*pimpmu020-two*vz*mu021+mu022; 

            real pimpmu220 = mu220 + c1o9 * OORho;
            real pimpmu202 = mu202 + c1o9 * OORho;
            real pimpmu022 = mu022 + c1o9 * OORho;

            //5.
            MzXXYYZ   =    /*four*(vx2y2z-vxyz*mu110+vxy*mu111)+*/ 
                           two*(vxz*mu120/*-vxy2*mu101-vx2y*mu011*/+vyz*mu210-vy*mu211-vx*mu121)+
                           vy2*mu201-vx2z*pimpmu020+mu021*vx2-vz*pimpmu220-vy2z*pimpmu200+mu221; 
            MzXXYZZ   =    /*four*(vx2yz2-vxyz*mu101+vxz*mu111)+*/
                           two*(vxy*mu102/*-vxz2*mu110-vx2z*mu011*/-vx*mu112-vz*mu211+vyz*mu201)+
                           vz2*mu210-vy*pimpmu202-pimpmu002*vx2y+mu012*vx2-vyz2*pimpmu200+mu212;
            MzXYYZZ   =    /*four*(vxy2z2-vxyz*mu011+vyz*mu111)+*/
                           two*(vxy*mu012/*-vyz2*mu110*/-vy*mu112+vxz*mu021-vz*mu121/*-vy2z*mu101*/)+
                           vy2*mu102+vz2*mu120-vxz2*pimpmu020-pimpmu022*vx-pimpmu002*vxy2+mu122;    

            //6.
            MzXXYYZZ  =   mu222 + pimpmu022*vx2 - two*mu122*vx + pimpmu202*vy2 /*+ pimpmu002*vx2*vy2 - two*mu102*vx*vy2*/ - two*mu212*vy /*- two*mu012*vx2*vy + four*mu112*vx*vy*/ + pimpmu220*vz2 /*+ pimpmu020*vx2*vz2 - two*mu120*vx*vz2 + 
                          pimpmu200*vy2*vz2 - five*vx2*vy2*vz2 - two*mu210*vy*vz2 + 4.*mu110*vx*vy*vz2*/ - two*mu221*vz /*- two*mu021*vx2*vz + four*mu121*vx*vz - two*mu201*vy2*vz + four*mu101*vx*vy2*vz + 
                          four*mu211*vy*vz + four*mu011*vx2*vy*vz - eight*mu111*vx*vy*vz*/;

            //MzXXYYZZ  =   /*-five*vx2y2z2-eight*vxyz*mu111+*/vy2*pimpmu202/*+pimpmu002*vx2y2+vy2z2*pimpmu200+vx2z2*pimpmu020*/+vz2*pimpmu220+mu222+pimpmu022*vx2//+
            //               /*four*(vx2yz*mu011+vxy2z*mu101+vxyz2*mu110+vxy*mu112+vxz*mu121+vyz*mu211)*/- 
            //               two*(vy*mu212/*-vy2z*mu201-vxz2*mu120-vyz2*mu210-vxy2*mu102-vx2z*mu021-vx2y*mu012*/-vx*mu122-vz*mu221); 
         }
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         real MXXpMYYpMZZ,MXXmMYY, MXXmMZZ, MXXYpMYZZ,MXXZpMYYZ,MXYYpMXZZ, MXXYmMYZZ,MXXZmMYYZ,MXYYmMXZZ;
         {
            //coll faktoren:
            real w1 =   (-omega);
            real w2 = -one;//(-omega);
            real w3 = -(two+(-omega));//-one;
            real w4 = -(two+(-omega));//-one;
            real w5 = -(two+(-omega));//-one;
            real w6 = -one;
            real w7 = -one;
            real w8 = -one;
            real w9 = -one;
            real w10= -one;

			real wadjust;
			real qudricLimit = c1o100;

            ////////lin kombi bilden:
            MXXpMYYpMZZ =  MzXX + MzYY + MzZZ;
            MXXmMYY     =  MzXX - MzYY;
            MXXmMZZ     =  MzXX - MzZZ;

            MXXYpMYZZ   =  MzXXY+MzYZZ;
            MXXYmMYZZ   =  MzXXY-MzYZZ;
            MXXZpMYYZ   =  MzXXZ+MzYYZ;
            MXXZmMYYZ   =  MzXXZ-MzYYZ;
            MXYYpMXZZ   =  MzXYY+MzXZZ;
            MXYYmMXZZ   =  MzXYY-MzXZZ;

            real MXXYYppp    = MzXXYY + MzXXZZ + MzYYZZ;
            real MXXYYpm2p   = MzXXYY - two*MzXXZZ + MzYYZZ;
            real MXXYYppm2   = MzXXYY + MzXXZZ - two*MzYYZZ;

            //relaxation:
            MXXpMYYpMZZ -= w2*(one-OORho-MXXpMYYpMZZ);

            MzXZ       *=  one+w1;
            MzYZ       *=  one+w1;
            MzXY       *=  one+w1;

			
            MXXmMYY    *=  one+w1;
            MXXmMZZ    *=  one+w1;

			wadjust     =  w5-(one+w5)*abs(MzXYZ)/(abs(MzXYZ)+qudricLimit);
			MzXYZ      *=  one+wadjust;

			wadjust     =  w3-(one+w3)*abs(MXXYpMYZZ)/(abs(MXXYpMYZZ)+qudricLimit);
            MXXYpMYZZ  *=  one+wadjust;
			wadjust     =  w3-(one+w3)*abs(MXXZpMYYZ)/(abs(MXXZpMYYZ)+qudricLimit);
            MXXZpMYYZ  *=  one+wadjust;
			wadjust     =  w3-(one+w3)*abs(MXYYpMXZZ)/(abs(MXYYpMXZZ)+qudricLimit);
            MXYYpMXZZ  *=  one+wadjust;
			wadjust     =  w4-(one+w4)*abs(MXXYmMYZZ)/(abs(MXXYmMYZZ)+qudricLimit);
            MXXYmMYZZ  *=  one+wadjust;
			wadjust     =  w4-(one+w4)*abs(MXXZmMYYZ)/(abs(MXXZmMYYZ)+qudricLimit);
            MXXZmMYYZ  *=  one+wadjust;
			wadjust     =  w4-(one+w4)*abs(MXYYmMXZZ)/(abs(MXYYmMXZZ)+qudricLimit);
            MXYYmMXZZ  *=  one+wadjust;

            ////////von Lin Kombis zurueck:
            MzXX  =  c1o3 * (       MXXmMYY +       MXXmMZZ + MXXpMYYpMZZ + OORho);
            MzYY  =  c1o3 * (-two * MXXmMYY +       MXXmMZZ + MXXpMYYpMZZ + OORho);
            MzZZ  =  c1o3 * (       MXXmMYY - two * MXXmMZZ + MXXpMYYpMZZ + OORho);

            MzXXY = (MXXYmMYZZ + MXXYpMYZZ)*c1o2;
            MzYZZ = (MXXYpMYZZ - MXXYmMYZZ)*c1o2;
            MzXYY = (MXYYmMXZZ + MXYYpMXZZ)*c1o2;
            MzXZZ = (MXYYpMXZZ - MXYYmMXZZ)*c1o2;
            MzXXZ = (MXXZmMYYZ + MXXZpMYYZ)*c1o2;
            MzYYZ = (MXXZpMYYZ - MXXZmMYYZ)*c1o2;

            //faktorisierte atraktoren:
            MXXYYppp  -=  w7*MzXX*MzYY + w7*MzXX*MzZZ     + w7*    MzZZ*MzYY - w7*MXXYYppp - w7*c1o3*OORho;
            MXXYYpm2p -=  w6*MzXX*MzYY - w6*two*MzXX*MzZZ + w6*    MzZZ*MzYY - w6*MXXYYpm2p;
            MXXYYppm2 -=  w6*MzXX*MzYY + w6*    MzXX*MzZZ - w6*two*MzZZ*MzYY - w6*MXXYYppm2;
            MzXXYYZZ  -= w10*MzXX*MzYY*MzZZ - w10*c1o27 - w10*MzXXYYZZ;
            MzXYYZ    -=  w8*MzYY*MzXZ - w8*MzXYYZ;
            MzXYZZ    -=  w8*MzZZ*MzXY - w8*MzXYZZ;
            MzXXYZ    -=  w8*MzXX*MzYZ - w8*MzXXYZ;

            MzXXYYZ *= one+w9;
            MzXXYZZ *= one+w9;
            MzXYYZZ *= one+w9;

            MzXXYY =  c1o3 * (MXXYYpm2p + MXXYYppm2 + MXXYYppp);
            MzXXZZ =  c1o3 * (MXXYYppp - MXXYYpm2p);
            MzYYZZ =  c1o3 * (MXXYYppp-MXXYYppm2);
         }

         //2.
         mu200 = vx2 + c1o3 * (MXXmMYY + MXXmMZZ + MXXpMYYpMZZ);
         mu020 = vy2 + c1o3 * (MXXmMZZ +  MXXpMYYpMZZ - two * MXXmMYY);
         mu002 = vz2 + c1o3 * (MXXmMYY - two * MXXmMZZ + MXXpMYYpMZZ);
         mu110 = vxy + MzXY;
         mu101 = vxz + MzXZ;
         mu011 = vyz + MzYZ;

         //3.
         mu111 = vxyz +     vx*MzYZ +     vy*MzXZ + vz*MzXY + MzXYZ;
         mu210 = vx2y + two*vx*MzXY +     vy*MzXX + MzXXY;
         mu120 = vxy2 +     vx*MzYY + two*vy*MzXY + MzXYY;
         mu102 = vxz2 +     vx*MzZZ + two*vz*MzXZ + MzXZZ;
         mu201 = vx2z + two*vx*MzXZ +     vz*MzXX + MzXXZ;
         mu021 = vy2z + two*vy*MzYZ +     vz*MzYY + MzYYZ;
         mu012 = vyz2 +     vy*MzZZ + two*vz*MzYZ + MzYZZ;

         //4.
         mu211 = /*vx2yz +     vx2*MzYZ + two*vxy*MzXZ + two*vxz*MzXY + two*vx*MzXYZ +*/     vyz*MzXX +     vy*MzXXZ +     vz*MzXXY + MzXXYZ;
         mu121 = /*vxy2z + two*vxy*MzYZ +*/     vxz*MzYY +     vx*MzYYZ /*+     vy2*MzXZ + two*vyz*MzXY + two*vy*MzXYZ*/ +     vz*MzXYY + MzXYYZ;
         mu112 = /*vxyz2 +*/     vxy*MzZZ /*+ two*vxz*MzYZ*/ +     vx*MzYZZ /*+ two*vyz*MzXZ*/ +     vy*MzXZZ /*+     vz2*MzXY + two*vz*MzXYZ*/ + MzXYZZ;
         mu220 = /*vx2y2 +*/     vx2*MzYY /*+ 4.f*vxy*MzXY*/ + two*vx*MzXYY +     vy2*MzXX + two*vy*MzXXY + MzXXYY ;
         mu202 = /*vx2z2 +*/     vx2*MzZZ /*+ 4.f*vxz*MzXZ*/ + two*vx*MzXZZ +     vz2*MzXX + two*vz*MzXXZ + MzXXZZ;
         mu022 = /*vy2z2 +*/     vy2*MzZZ /*+ 4.f*vyz*MzYZ*/ + two*vy*MzYZZ +     vz2*MzYY + two*vz*MzYYZ + MzYYZZ;

         MzXXYY += c1o9 * OORho;
         MzXXZZ += c1o9 * OORho;
         MzYYZZ += c1o9 * OORho;

         //5.
         mu221 = /*vx2y2z+*/vx2z*MzYY+vx2*MzYYZ+vy2z*MzXX+vy2*MzXXZ+vz*MzXXYY+MzXXYYZ+
                  /*four*(vxyz*MzXY+vxy*MzXYZ)+*/
                  two*(vxz*MzXYY+vx*MzXYYZ/*+vx2y*MzYZ*/+vyz*MzXXY+vy*MzXXYZ/*+vxy2*MzXZ*/);
         mu212 = /*vx2yz2 +*/ vx2y*MzZZ + vx2*MzYZZ + vy*MzXXZZ + vz2*MzXXY + vyz2*MzXX+MzXXYZZ+
                  /*four*(vxyz*MzXZ + vxz*MzXYZ)+*/
                  two*(vxy*MzXZZ /*+ vxz2*MzXY*/ + vx*MzXYZZ + vyz*MzXXZ /*+ vx2z*MzYZ*/ + vz*MzXXYZ);
         mu122 = /*vxy2z2*/ + vxy2*MzZZ + MzXYYZZ + vz2*MzXYY + vy2*MzXZZ + vx*MzYYZZ + vxz2*MzYY+
                  /*four*(vxyz*MzYZ + vyz*MzXYZ)+*/
                  two*(vxy*MzYZZ + vxz*MzYYZ/* + vy2z*MzXZ + vyz2*MzXY*/ + vy*MzXYZZ + vz*MzXYYZ);

         //6.
         mu222 = /*vx2y2z2 + vx2y2*MzZZ + vx2z2*MzYY */+ vx2*MzYYZZ /*+ vy2z2*MzXX*/ + vy2*MzXXZZ + vz2*MzXXYY + MzXXYYZZ /*+ eight*vxyz*MzXYZ*/ +
            two*(/*vx2y*MzYZZ +*/  vx*MzXYYZZ + vz*MzXXYYZ /*+ vyz2*MzXXY + vx2z*MzYYZ + vxy2*MzXZZ + vxz2*MzXYY*/ + vy*MzXXYZZ/* + vy2z*MzXXZ*/);//+ 
            //four*(/*vxy2z*MzXZ + vx2yz*MzYZ + vxyz2*MzXY +*/ vxy*MzXYZZ + vxz*MzXYYZ + vyz*MzXXYZ);

         //(D.f[ dirE   ])[k   ] =   c1o2*rho*( mu200  - mu220 + mu222 - mu202 +  mu120 - mu122 + mu102 - vx   );   //ke
         //(D.f[ dirW   ])[kw  ] =   c1o2*rho*( mu200  - mu220 + mu222 - mu202 -  mu120 + mu122 - mu102 + vx   );   
         //(D.f[ dirN   ])[k   ] =   c1o2*rho*( mu210  - mu220 + mu222 - mu212 +  mu020 - mu022 + mu012 - vy   );   //kn
         //(D.f[ dirS   ])[ks  ] =   c1o2*rho*(-mu210  - mu220 + mu222 + mu212 +  mu020 - mu022 - mu012 + vy   );   
         //(D.f[ dirT   ])[k   ] =   c1o2*rho*(-mu221  + mu222 + mu201 - mu202 +  mu021 - mu022 + mu002 - vz   );   //kt
         //(D.f[ dirB   ])[kb  ] =   c1o2*rho*( mu221  + mu222 - mu201 - mu202 -  mu021 - mu022 + mu002 + vz   );   
         //(D.f[ dirNE  ])[k   ] =  c1o4*rho*(-mu210  + mu220 - mu222 + mu212 +  mu110 - mu120 + mu122 - mu112);   //kne
         //(D.f[ dirSW  ])[ksw ] =  c1o4*rho*( mu210  + mu220 - mu222 - mu212 +  mu110 + mu120 - mu122 - mu112);   
         //(D.f[ dirSE  ])[ks  ] =  c1o4*rho*( mu210  + mu220 - mu222 - mu212 -  mu110 - mu120 + mu122 + mu112);   //kse
         //(D.f[ dirNW  ])[kw  ] =  c1o4*rho*(-mu210  + mu220 - mu222 + mu212 -  mu110 + mu120 - mu122 + mu112);   //knw
         //(D.f[ dirTE  ])[k   ] =  c1o4*rho*( mu221  - mu222 - mu201 + mu202 -  mu121 + mu122 + mu101 - mu102);   //kte
         //(D.f[ dirBW  ])[kbw ] =  c1o4*rho*(-mu221  - mu222 + mu201 + mu202 -  mu121 - mu122 + mu101 + mu102);   
         //(D.f[ dirBE  ])[kb  ] =  c1o4*rho*(-mu221  - mu222 + mu201 + mu202 +  mu121 + mu122 - mu101 - mu102);   //kbe
         //(D.f[ dirTW  ])[kw  ] =  c1o4*rho*( mu221  - mu222 - mu201 + mu202 +  mu121 - mu122 - mu101 + mu102);   //ktw
         //(D.f[ dirTN  ])[k   ] =  c1o4*rho*( mu221  - mu222 - mu211 + mu212 -  mu021 + mu022 + mu011 - mu012);   //ktn
         //(D.f[ dirBS  ])[kbs ] =  c1o4*rho*(-mu221  - mu222 - mu211 - mu212 +  mu021 + mu022 + mu011 + mu012);   
         //(D.f[ dirBN  ])[kb  ] =  c1o4*rho*(-mu221  - mu222 + mu211 + mu212 +  mu021 + mu022 - mu011 - mu012);   //kbn
         //(D.f[ dirTS  ])[ks  ] =  c1o4*rho*( mu221  - mu222 + mu211 - mu212 -  mu021 + mu022 - mu011 + mu012);   //kts
         //(D.f[ dirZERO])[k   ] =       rho*(-mu200  + mu220 - mu222 + mu202 -  mu020 + mu022 - mu002        )+rho0;   //kzero
         //(D.f[ dirTNE ])[k   ] = c1o8*rho*(-mu221  + mu222 + mu211 - mu212 +  mu121 - mu122 - mu111 + mu112);   //ktne
         //(D.f[ dirTSE ])[ks  ] = c1o8*rho*(-mu221  + mu222 - mu211 + mu212 +  mu121 - mu122 + mu111 - mu112);   //ktse
         //(D.f[ dirBNE ])[kb  ] = c1o8*rho*( mu221  + mu222 - mu211 - mu212 -  mu121 - mu122 + mu111 + mu112);   //kbne
         //(D.f[ dirBSE ])[kbs ] = c1o8*rho*( mu221  + mu222 + mu211 + mu212 -  mu121 - mu122 - mu111 - mu112);   //kbse
         //(D.f[ dirTNW ])[kw  ] = c1o8*rho*(-mu221  + mu222 + mu211 - mu212 -  mu121 + mu122 + mu111 - mu112);   //ktnw
         //(D.f[ dirTSW ])[ksw ] = c1o8*rho*(-mu221  + mu222 - mu211 + mu212 -  mu121 + mu122 - mu111 + mu112);   //ktsw
         //(D.f[ dirBNW ])[kbw ] = c1o8*rho*( mu221  + mu222 - mu211 - mu212 +  mu121 + mu122 - mu111 - mu112);   //kbnw
         //(D.f[ dirBSW ])[kbsw] = c1o8*rho*( mu221  + mu222 + mu211 + mu212 +  mu121 + mu122 + mu111 + mu112);   
         (D.f[ dirE   ])[k   ] =   c1o2*rho*(+ mu222 + (                          - mu220 - mu202        ) + (                         + mu200                ) + (                 - mu122) + (                                          + mu102 + mu120) + ( - vx            ) );   //ke
         (D.f[ dirW   ])[kw  ] =   c1o2*rho*(+ mu222 + (                          - mu220 - mu202        ) + (                         + mu200                ) + (                 + mu122) + (                                          - mu102 - mu120) + ( + vx            ) );   
         (D.f[ dirN   ])[k   ] =   c1o2*rho*(+ mu222 + (                          - mu220         - mu022) + (                                 + mu020        ) + (         - mu212        ) + (                          + mu012 + mu210                ) + (      - vy       ) );   //kn
         (D.f[ dirS   ])[ks  ] =   c1o2*rho*(+ mu222 + (                          - mu220         - mu022) + (                                 + mu020        ) + (         + mu212        ) + (                          - mu012 - mu210                ) + (      + vy       ) );   
         (D.f[ dirT   ])[k   ] =   c1o2*rho*(+ mu222 + (                                  - mu202 - mu022) + (                                         + mu002) + ( - mu221                ) + (         + mu201 +  mu021                                ) + (           - vz  ) );   //kt
         (D.f[ dirB   ])[kb  ] =   c1o2*rho*(+ mu222 + (                                  - mu202 - mu022) + (                                         + mu002) + ( + mu221                ) + (         - mu201 -  mu021                                ) + (           + vz  ) );   
         (D.f[ dirNE  ])[k   ] =  c1o4*rho*(- mu222 + (                  - mu112 + mu220                ) + (+  mu110                                        ) + (         + mu212 + mu122) + (                                  - mu210         - mu120)                       );   //kne
         (D.f[ dirSW  ])[ksw ] =  c1o4*rho*(- mu222 + (                  - mu112 + mu220                ) + (+  mu110                                        ) + (         - mu212 - mu122) + (                                  + mu210         + mu120)                       );   
         (D.f[ dirSE  ])[ks  ] =  c1o4*rho*(- mu222 + (                  + mu112 + mu220                ) + (-  mu110                                        ) + (         - mu212 + mu122) + (                                  + mu210         - mu120)                       );   //kse
         (D.f[ dirNW  ])[kw  ] =  c1o4*rho*(- mu222 + (                  + mu112 + mu220                ) + (-  mu110                                        ) + (         + mu212 - mu122) + (                                  - mu210         + mu120)                       );   //knw
         (D.f[ dirTE  ])[k   ] =  c1o4*rho*(- mu222 + (        -  mu121                  + mu202        ) + (         + mu101                                ) + ( + mu221         + mu122) + (         - mu201                          - mu102        )                       );   //kte
         (D.f[ dirBW  ])[kbw ] =  c1o4*rho*(- mu222 + (        -  mu121                  + mu202        ) + (         + mu101                                ) + ( - mu221         - mu122) + (         + mu201                          + mu102        )                       );   
         (D.f[ dirBE  ])[kb  ] =  c1o4*rho*(- mu222 + (        +  mu121                  + mu202        ) + (         - mu101                                ) + ( - mu221         + mu122) + (         + mu201                          - mu102        )                       );   //kbe
         (D.f[ dirTW  ])[kw  ] =  c1o4*rho*(- mu222 + (        +  mu121                  + mu202        ) + (         - mu101                                ) + ( + mu221         - mu122) + (         - mu201                          + mu102        )                       );   //ktw
         (D.f[ dirTN  ])[k   ] =  c1o4*rho*(- mu222 + (- mu211                                   + mu022) + (                 + mu011                        ) + ( + mu221 + mu212        ) + (                 -  mu021 - mu012                        )                       );   //ktn
         (D.f[ dirBS  ])[kbs ] =  c1o4*rho*(- mu222 + (- mu211                                   + mu022) + (                 + mu011                        ) + ( - mu221 - mu212        ) + (                 +  mu021 + mu012                        )                       );   
         (D.f[ dirBN  ])[kb  ] =  c1o4*rho*(- mu222 + (+ mu211                                   + mu022) + (                 - mu011                        ) + ( - mu221 + mu212        ) + (                 +  mu021 - mu012                        )                       );   //kbn
         (D.f[ dirTS  ])[ks  ] =  c1o4*rho*(- mu222 + (+ mu211                                   + mu022) + (                 - mu011                        ) + ( + mu221 - mu212        ) + (                 -  mu021 + mu012                        )                       );   //kts
         (D.f[ dirZERO])[k   ] =       rho*(- mu222 + (                          + mu220 + mu202 + mu022) + (                         - mu200 - mu020 - mu002)                                                                                                                  )+rho0;   //kzero
         (D.f[ dirTNE ])[k   ] = c1o8*rho*(+ mu222 + (+ mu211 +  mu121 + mu112                         )                                                      + ( - mu221 - mu212 - mu122) + ( - mu111                                                 )                       );   //ktne
         (D.f[ dirTSE ])[ks  ] = c1o8*rho*(+ mu222 + (- mu211 +  mu121 - mu112                         )                                                      + ( - mu221 + mu212 - mu122) + ( + mu111                                                 )                       );   //ktse
         (D.f[ dirBNE ])[kb  ] = c1o8*rho*(+ mu222 + (- mu211 -  mu121 + mu112                         )                                                      + ( + mu221 - mu212 - mu122) + ( + mu111                                                 )                       );   //kbne
         (D.f[ dirBSE ])[kbs ] = c1o8*rho*(+ mu222 + (+ mu211 -  mu121 - mu112                         )                                                      + ( + mu221 + mu212 - mu122) + ( - mu111                                                 )                       );   //kbse
         (D.f[ dirTNW ])[kw  ] = c1o8*rho*(+ mu222 + (+ mu211 -  mu121 - mu112                         )                                                      + ( - mu221 - mu212 + mu122) + ( + mu111                                                 )                       );   //ktnw
         (D.f[ dirTSW ])[ksw ] = c1o8*rho*(+ mu222 + (- mu211 -  mu121 + mu112                         )                                                      + ( - mu221 + mu212 + mu122) + ( - mu111                                                 )                       );   //ktsw
         (D.f[ dirBNW ])[kbw ] = c1o8*rho*(+ mu222 + (- mu211 +  mu121 - mu112                         )                                                      + ( + mu221 - mu212 + mu122) + ( - mu111                                                 )                       );   //kbnw
         (D.f[ dirBSW ])[kbsw] = c1o8*rho*(+ mu222 + (+ mu211 +  mu121 + mu112                         )                                                      + ( + mu221 + mu212 + mu122) + ( + mu111                                                 )                       );   
                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                
         //////////////////////////////////////////////////////////////////////////                                                                                                                                                                             
         //BGK                                                                                                 
         //////////////////////////////////////////////////////////////////////////                            
         //(D.f[ dirE   ])[k   ] = fW    ;                                                                     
         //(D.f[ dirW   ])[kw  ] = fE    ;                                                                     
         //(D.f[ dirN   ])[k   ] = fS    ;
         //(D.f[ dirS   ])[ks  ] = fN    ;
         //(D.f[ dirT   ])[k   ] = fB    ;
         //(D.f[ dirB   ])[kb  ] = fT    ;
         //(D.f[ dirNE  ])[k   ] = fSW   ;
         //(D.f[ dirSW  ])[ksw ] = fNE   ;
         //(D.f[ dirSE  ])[ks  ] = fNW   ;
         //(D.f[ dirNW  ])[kw  ] = fSE   ;
         //(D.f[ dirTE  ])[k   ] = fBW   ;
         //(D.f[ dirBW  ])[kbw ] = fTE   ;
         //(D.f[ dirBE  ])[kb  ] = fTW   ;
         //(D.f[ dirTW  ])[kw  ] = fBE   ;
         //(D.f[ dirTN  ])[k   ] = fBS   ;
         //(D.f[ dirBS  ])[kbs ] = fTN   ;
         //(D.f[ dirBN  ])[kb  ] = fTS   ;
         //(D.f[ dirTS  ])[ks  ] = fBN   ;
         //(D.f[ dirZERO])[k   ] = fZERO ;
         //(D.f[ dirTNE ])[k   ] = fBSW  ;
         //(D.f[ dirTSE ])[ks  ] = fBNW  ;
         //(D.f[ dirBNE ])[kb  ] = fTSW  ;
         //(D.f[ dirBSE ])[kbs ] = fTNW  ;
         //(D.f[ dirTNW ])[kw  ] = fBSE  ;
         //(D.f[ dirTSW ])[ksw ] = fBNE  ;
         //(D.f[ dirBNW ])[kbw ] = fTSE  ;
         //(D.f[ dirBSW ])[kbsw] = fTNE  ;
      }                                                                                                                    
   }
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Casc_SP_MS_27(   real omega,
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
      BC        =   bcMatD[k];

      if( (BC != GEO_SOLID) && (BC != GEO_VOID))
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
         real fE    =  (D.f[dirE   ])[k  ];//ke
         real fW    =  (D.f[dirW   ])[kw ];
         real fN    =  (D.f[dirN   ])[k  ];//kn
         real fS    =  (D.f[dirS   ])[ks ];
         real fT    =  (D.f[dirT   ])[k  ];//kt
         real fB    =  (D.f[dirB   ])[kb ];
         real fNE   =  (D.f[dirNE  ])[k  ];//kne
         real fSW   =  (D.f[dirSW  ])[ksw];
         real fSE   =  (D.f[dirSE  ])[ks ];//kse
         real fNW   =  (D.f[dirNW  ])[kw ];//knw
         real fTE   =  (D.f[dirTE  ])[k  ];//kte
         real fBW   =  (D.f[dirBW  ])[kbw];
         real fBE   =  (D.f[dirBE  ])[kb ];//kbe
         real fTW   =  (D.f[dirTW  ])[kw ];//ktw
         real fTN   =  (D.f[dirTN  ])[k  ];//ktn
         real fBS   =  (D.f[dirBS  ])[kbs];
         real fBN   =  (D.f[dirBN  ])[kb ];//kbn
         real fTS   =  (D.f[dirTS  ])[ks ];//kts
         real fZERO =  (D.f[dirZERO])[k  ];//kzero
         real fTNE   = (D.f[dirTNE ])[k  ];//ktne
         real fTSW   = (D.f[dirTSW ])[ksw];//ktsw
         real fTSE   = (D.f[dirTSE ])[ks ];//ktse
         real fTNW   = (D.f[dirTNW ])[kw ];//ktnw
         real fBNE   = (D.f[dirBNE ])[kb ];//kbne
         real fBSW   = (D.f[dirBSW ])[kbsw];
         real fBSE   = (D.f[dirBSE ])[kbs];//kbse
         real fBNW   = (D.f[dirBNW ])[kbw];//kbnw
         ////////////////////////////////////////////////////////////////////////////////
         real rho0   =  fZERO+fE+fW+fN+fS+fT+fB+fNE+fSW+fSE+fNW+fTE+fBW+fBE+fTW+fTN+fBS+fBN+fTS+fTNE+fTSW+fTSE+fTNW+fBNE+fBSW+fBSE+fBNW;
         real rho    =  rho0 + one;
         real OORho  =  one/rho;
         real vx     =  (fE -fW +fNE-fSW+fSE-fNW+fTE-fBW+fBE-fTW+ fTNE-fTSW+fTSE-fTNW+ fBNE-fBSW+fBSE-fBNW) * OORho;
         real vy     =  (fN -fS +fNE-fSW-fSE+fNW+fTN-fBS+fBN-fTS+ fTNE-fTSW-fTSE+fTNW+ fBNE-fBSW-fBSE+fBNW) * OORho;
         real vz     =  (fT -fB +fTE-fBW-fBE+fTW+fTN-fBS-fBN+fTS+ fTNE+fTSW+fTSE+fTNW- fBNE-fBSW-fBSE-fBNW) * OORho;

         //////////////////////////////////////////////////////////////////////////
         //BGK
         //////////////////////////////////////////////////////////////////////////
         //real drho   =  fZERO+fE+fW+fN+fS+fT+fB+fNE+fSW+fSE+fNW+fTE+fBW+fBE+fTW+fTN+fBS+fBN+fTS+fTNE+fTSW+fTSE+fTNW+fBNE+fBSW+fBSE+fBNW;
         //real vx1     =  (fE -fW +fNE-fSW+fSE-fNW+fTE-fBW+fBE-fTW+ fTNE-fTSW+fTSE-fTNW+ fBNE-fBSW+fBSE-fBNW);
         //real vx2     =  (fN -fS +fNE-fSW-fSE+fNW+fTN-fBS+fBN-fTS+ fTNE-fTSW-fTSE+fTNW+ fBNE-fBSW-fBSE+fBNW);
         //real vx3     =  (fT -fB +fTE-fBW-fBE+fTW+fTN-fBS-fBN+fTS+ fTNE+fTSW+fTSE+fTNW- fBNE-fBSW-fBSE-fBNW);


         //real cusq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         //fZERO = fZERO *(one+(-omega))-(-omega)*   c8over27*  (drho-cusq);
         //fE    = fE    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq);
         //fW    = fW    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq);
         //fN    = fN    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq);
         //fS    = fS    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq);
         //fT    = fT    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq);
         //fB    = fB    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq);
         //fNE   = fNE   *(one+(-omega))-(-omega)*   c1over54*  (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq);
         //fSW   = fSW   *(one+(-omega))-(-omega)*   c1over54*  (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq);
         //fSE   = fSE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq);
         //fNW   = fNW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq);
         //fTE   = fTE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq);
         //fBW   = fBW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq);
         //fBE   = fBE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq);
         //fTW   = fTW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq);
         //fTN   = fTN   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq);
         //fBS   = fBS   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq);
         //fBN   = fBN   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq);
         //fTS   = fTS   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq);
         //fTNE  = fTNE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq);
         //fBSW  = fBSW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq);
         //fBNE  = fBNE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq);
         //fTSW  = fTSW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq);
         //fTSE  = fTSE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq);
         //fBNW  = fBNW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq);
         //fBSE  = fBSE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq);
         //fTNW  = fTNW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq);


         real vx2    = vx*vx;
         real vy2    = vy*vy;
         real vz2    = vz*vz;
         real vxy    = vx*vy;
         real vxz    = vx*vz;
         real vyz    = vy*vz;
         real vx2y   = vx*vx*vy;
         real vx2z   = vx*vx*vz;
         real vy2z   = vy*vy*vz;
         real vxy2   = vx*vy*vy;
         real vyz2   = vy*vz*vz;
         real vxz2   = vx*vz*vz;
         real vxyz   = vx*vy*vz;
         real vx2y2  = vx*vx*vy*vy;
         real vx2z2  = vx*vx*vz*vz;
         real vy2z2  = vy*vy*vz*vz;
         real vx2yz  = vx*vx*vy*vz;
         real vxyz2  = vx*vy*vz*vz;
         real vxy2z  = vx*vy*vy*vz;
         real vx2y2z = vx*vx*vy*vy*vz;
         real vx2yz2 = vx*vx*vy*vz*vz;
         real vxy2z2 = vx*vy*vy*vz*vz;
         real vx2y2z2= vx*vx*vy*vy*vz*vz;
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         real mu200   =(fE+fW+fNE+fSW+fSE+fNW+fTE+fBW+fBE+fTW+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         real mu020   =(fN+fS+fNE+fSW+fSE+fNW+fTN+fBS+fBN+fTS+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         real mu002   =(fT+fB+fTE+fBW+fBE+fTW+fTN+fBS+fBN+fTS+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         real mu110   =(fNE-fSE+fSW-fNW+fTNE-fTSE+fBNE-fBSE+fTSW-fTNW+fBSW-fBNW) * OORho;
         real mu101   =(fTE+fBW-fBE-fTW+fTNE-fBNE+fTSE-fBSE-fTNW+fBNW-fTSW+fBSW) * OORho;
         real mu011   =(fTN+fBS-fBN-fTS+fTNE-fBNE-fTSE+fBSE+fTNW-fBNW-fTSW+fBSW) * OORho;
         real mu210   =(fNE-fSW-fSE+fNW+fTNE+fBNE-fTSE-fBSE+fTNW+fBNW-fTSW-fBSW) * OORho;
         real mu120   =(fNE-fSW+fSE-fNW+fTNE-fTNW+fBNE-fBNW+fTSE-fTSW+fBSE-fBSW) * OORho;
         real mu102   =(fTE-fBW+fBE-fTW+fTNE-fTNW+fBNE-fBNW+fTSE-fTSW+fBSE-fBSW) * OORho;
         real mu111   =(fTNE-fBNE-fTSE+fBSE-fTNW+fBNW+fTSW-fBSW) * OORho;
         real mu201   =(fTE-fBW-fBE+fTW+fTNE-fBNE+fTSE-fBSE+fTNW-fBNW+fTSW-fBSW) * OORho;
         real mu021   =(fTN-fBS-fBN+fTS+fTNE-fBNE+fTSE-fBSE+fTNW-fBNW+fTSW-fBSW) * OORho;
         real mu012   =(fTN-fBS+fBN-fTS+fTNE+fBNE-fTSE-fBSE+fTNW+fBNW-fTSW-fBSW) * OORho;
         real mu220   =(fNE+fSW+fSE+fNW+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         real mu121   =(fTNE-fBNE+fTSE-fBSE-fTNW+fBNW-fTSW+fBSW) * OORho;
         real mu202   =(fTE+fBW+fBE+fTW+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         real mu211   =(fTNE-fBNE-fTSE+fBSE+fTNW-fBNW-fTSW+fBSW) * OORho;
         real mu112   =(fTNE+fBNE-fTSE-fBSE-fTNW-fBNW+fTSW+fBSW) * OORho;
         real mu022   =(fTN+fBS+fBN+fTS+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         real mu221   =(fTNE-fBNE+fTSE-fBSE+fTNW-fBNW+fTSW-fBSW) * OORho;
         real mu122   =(fTNE-fTNW+fBNE-fBNW+fTSE-fTSW+fBSE-fBSW) * OORho;
         real mu212   =(fTNE+fBNE-fTSE-fBSE+fTNW+fBNW-fTSW-fBSW) * OORho;
         real mu222   =(fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         real MzXX,MzYY,MzZZ,MzXY,MzXZ,MzYZ,MzXXY,MzXYY,MzXXZ,MzXZZ,MzYYZ,MzYZZ,MzXYZ,MzXXYY,MzXXZZ,MzYYZZ,MzXXYZ,MzXYYZ,MzXYZZ,MzXXYYZ,MzXXYZZ,MzXYYZZ,MzXXYYZZ;
         {
            //2.
            MzXX      =   mu200-vx2;
            MzXY      =   mu110-vxy;
            MzXZ      =   mu101-vxz;
            MzYY      =   mu020-vy2;
            MzYZ      =   mu011-vyz;
            MzZZ      =   mu002-vz2;

            real pimpmu200 = mu200 + c1o3 * OORho;
            real pimpmu020 = mu020 + c1o3 * OORho;
            real pimpmu002 = mu002 + c1o3 * OORho;

            //3.
            MzXXY     =    two*vx2y - two*vx*mu110 - vy*pimpmu200 + mu210; 
            MzXXZ     =    two*vx2z - two*vx*mu101 - vz*pimpmu200 + mu201; 
            MzXYY     =    two*vxy2 - pimpmu020*vx - two*vy*mu110 + mu120; 
            MzXYZ     =    two*vxyz - mu011*vx-vy*mu101 -vz*mu110 + mu111;
            MzXZZ     =    two*vxz2 - pimpmu002*vx - two*vz*mu101 + mu102; 
            MzYYZ     =    two*vy2z - two*vy*mu011 - vz*pimpmu020 + mu021; 
            MzYZZ     =    two*vyz2 - pimpmu002*vy - two*vz*mu011 + mu012; 

            //4.
            MzXXYY    =   -three*vx2y2+pimpmu020*vx2+four*vxy*mu110-two*vx*mu120+vy2*pimpmu200-two*vy*mu210+mu220; 
            MzXXYZ    =   -three*vx2yz+mu011*vx2+two*vxy*mu101+two*vxz*mu110-two*vx*mu111+vyz*pimpmu200-vy*mu201-vz*mu210+mu211; 
            MzXXZZ    =   -three*vx2z2+pimpmu002*vx2+four*vxz*mu101-two*vx*mu102+vz2*pimpmu200-two*vz*mu201+mu202; 
            MzXYYZ    =   -three*vxy2z+two*vxy*mu011+vxz*pimpmu020-mu021*vx+vy2*mu101+two*vyz*mu110-two*vy*mu111-vz*mu120+mu121; 
            MzXYZZ    =   -three*vxyz2+pimpmu002*vxy+two*vxz*mu011-mu012*vx+two*vyz*mu101-vy*mu102+vz2*mu110-two*vz*mu111+mu112; 
            MzYYZZ    =   -three*vy2z2+pimpmu002*vy2+four*vyz*mu011-two*vy*mu012+vz2*pimpmu020-two*vz*mu021+mu022; 
                             
            real pimpmu220 = mu220 + c1o9 * OORho;
            real pimpmu202 = mu202 + c1o9 * OORho;
            real pimpmu022 = mu022 + c1o9 * OORho;

            //5.
            MzXXYYZ   =    four*vx2y2z-four*vxyz*mu110+four*vxy*mu111+ 
                           two*(vxz*mu120-vxy2*mu101-vx2y*mu011+vyz*mu210-vy*mu211-vx*mu121)+
                           vy2*mu201-vx2z*pimpmu020+mu021*vx2-vz*pimpmu220-vy2z*pimpmu200+mu221; 
            MzXXYZZ   =    four*(vx2yz2-vxyz*mu101+vxz*mu111)+
                           two*(vxy*mu102-vxz2*mu110-vx2z*mu011-vx*mu112-vz*mu211+vyz*mu201)+
                           vz2*mu210-vy*pimpmu202-pimpmu002*vx2y+mu012*vx2-vyz2*pimpmu200+mu212;
            MzXYYZZ   =    four*(vxy2z2-vxyz*mu011+vyz*mu111)+
                           two*(vxy*mu012-vyz2*mu110-vy*mu112+vxz*mu021-vz*mu121-vy2z*mu101)+
                           vy2*mu102+vz2*mu120-vxz2*pimpmu020-pimpmu022*vx-pimpmu002*vxy2+mu122;    

            //6.
            MzXXYYZZ  =   mu222 + pimpmu022*vx2 - two*mu122*vx + pimpmu202*vy2 + pimpmu002*vx2*vy2 - two*mu102*vx*vy2 - two*mu212*vy - two*mu012*vx2*vy + four*mu112*vx*vy + pimpmu220*vz2 + pimpmu020*vx2*vz2 - two*mu120*vx*vz2 + 
                          pimpmu200*vy2*vz2 - five*vx2*vy2*vz2 - two*mu210*vy*vz2 + four*mu110*vx*vy*vz2 - two*mu221*vz - two*mu021*vx2*vz + four*mu121*vx*vz - two*mu201*vy2*vz + four*mu101*vx*vy2*vz + 
                          four*mu211*vy*vz + four*mu011*vx2*vy*vz - eight*mu111*vx*vy*vz;
            //MzXXYYZZ  =   -five*vx2y2z2-eight*vxyz*mu111+vy2*pimpmu202+pimpmu002*vx2y2+vy2z2*pimpmu200+vx2z2*pimpmu020+vz2*pimpmu220+mu222+pimpmu022*vx2+
            //               four*(vx2yz*mu011+vxy2z*mu101+vxyz2*mu110+vxy*mu112+vxz*mu121+vyz*mu211)- 
            //               two*(vy*mu212-vy2z*mu201-vxz2*mu120-vyz2*mu210-vxy2*mu102-vx2z*mu021-vx2y*mu012-vx*mu122-vz*mu221); 
         }
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         real MXXpMYYpMZZ,MXXmMYY, MXXmMZZ, MXXYpMYZZ,MXXZpMYYZ,MXYYpMXZZ, MXXYmMYZZ,MXXZmMYYZ,MXYYmMXZZ;
         {
            //coll faktoren:
            real w1 =  -omega;
            real w2 = -one;//-omega;
            real w3 = -one;
            real w4 = -one;
            real w5 = -one;
            real w6 = -one;
            real w7 = -one;
            real w8 = -one;
            real w9 = -one;
            real w10= -one;

            ////////lin kombi bilden:
            MXXpMYYpMZZ = MzXX + MzYY + MzZZ;
            MXXmMYY     =  MzXX - MzYY;
            MXXmMZZ     =  MzXX - MzZZ;

            MXXYpMYZZ   =  MzXXY+MzYZZ;
            MXXYmMYZZ   =  MzXXY-MzYZZ;
            MXXZpMYYZ   =  MzXXZ+MzYYZ;
            MXXZmMYYZ   =  MzXXZ-MzYYZ;
            MXYYpMXZZ   =  MzXYY+MzXZZ;
            MXYYmMXZZ   =  MzXYY-MzXZZ;

            real MXXYYppp    = MzXXYY + MzXXZZ + MzYYZZ;
            real MXXYYpm2p   = MzXXYY - two*MzXXZZ + MzYYZZ;
            real MXXYYppm2   = MzXXYY + MzXXZZ - two*MzYYZZ;

            //relaxation:
            MXXpMYYpMZZ -= w2*(one-OORho-MXXpMYYpMZZ);

            MzXZ       *=  one+w1;
            MzYZ       *=  one+w1;
            MzXY       *=  one+w1;

            MzXYZ      *=  one+w5;

            MXXmMYY    *=  one+w1;
            MXXmMZZ    *=  one+w1;
                                       
            MXXYpMYZZ  *=  one+w3;
            MXXYmMYZZ  *=  one+w4;
            MXXZpMYYZ  *=  one+w3;
            MXXZmMYYZ  *=  one+w4;
            MXYYpMXZZ  *=  one+w3;
            MXYYmMXZZ  *=  one+w4;

            ////////von Lin Kombis zurueck:
            MzXX  =  c1o3 * (       MXXmMYY +       MXXmMZZ + MXXpMYYpMZZ + OORho);
            MzYY  =  c1o3 * (-two * MXXmMYY +       MXXmMZZ + MXXpMYYpMZZ + OORho);
            MzZZ  =  c1o3 * (       MXXmMYY - two * MXXmMZZ + MXXpMYYpMZZ + OORho);

            MzXXY = (MXXYmMYZZ + MXXYpMYZZ)*c1o2;
            MzYZZ = (MXXYpMYZZ - MXXYmMYZZ)*c1o2;
            MzXYY = (MXYYmMXZZ + MXYYpMXZZ)*c1o2;
            MzXZZ = (MXYYpMXZZ - MXYYmMXZZ)*c1o2;
            MzXXZ = (MXXZmMYYZ + MXXZpMYYZ)*c1o2;
            MzYYZ = (MXXZpMYYZ - MXXZmMYYZ)*c1o2;

            //faktorisierte atraktoren:
            MXXYYppp  -=  w7*MzXX*MzYY + w7*MzXX*MzZZ     + w7*    MzZZ*MzYY - w7*MXXYYppp - w7*c1o3*OORho;
            MXXYYpm2p -=  w6*MzXX*MzYY - w6*two*MzXX*MzZZ + w6*    MzZZ*MzYY - w6*MXXYYpm2p;
            MXXYYppm2 -=  w6*MzXX*MzYY + w6*    MzXX*MzZZ - w6*two*MzZZ*MzYY - w6*MXXYYppm2;
            MzXXYYZZ  -= w10*MzXX*MzYY*MzZZ - w10*c1o27 - w10*MzXXYYZZ;
            MzXYYZ    -=  w8*MzYY*MzXZ - w8*MzXYYZ;
            MzXYZZ    -=  w8*MzZZ*MzXY - w8*MzXYZZ;
            MzXXYZ    -=  w8*MzXX*MzYZ - w8*MzXXYZ;

            MzXXYYZ *= one+w9;
            MzXXYZZ *= one+w9;
            MzXYYZZ *= one+w9;

            MzXXYY =  c1o3 * (MXXYYpm2p + MXXYYppm2 + MXXYYppp);
            MzXXZZ =  c1o3 * (MXXYYppp - MXXYYpm2p);
            MzYYZZ =  c1o3 * (MXXYYppp-MXXYYppm2);
         }

         //2.
         mu200 = vx2 + c1o3 * (MXXmMYY + MXXmMZZ + MXXpMYYpMZZ);
         mu020 = vy2 + c1o3 * (MXXmMZZ +  MXXpMYYpMZZ - two * MXXmMYY);
         mu002 = vz2 + c1o3 * (MXXmMYY - two * MXXmMZZ + MXXpMYYpMZZ);
         mu110 = vxy + MzXY;
         mu101 = vxz + MzXZ;
         mu011 = vyz + MzYZ;

         //3.
         mu111 = vxyz +     vx*MzYZ +     vy*MzXZ + vz*MzXY + MzXYZ;
         mu210 = vx2y + two*vx*MzXY +     vy*MzXX + MzXXY;
         mu120 = vxy2 +     vx*MzYY + two*vy*MzXY + MzXYY;
         mu102 = vxz2 +     vx*MzZZ + two*vz*MzXZ + MzXZZ;
         mu201 = vx2z + two*vx*MzXZ +     vz*MzXX + MzXXZ;
         mu021 = vy2z + two*vy*MzYZ +     vz*MzYY + MzYYZ;
         mu012 = vyz2 +     vy*MzZZ + two*vz*MzYZ + MzYZZ;

         //4.
         mu211 = vx2yz +     vx2*MzYZ + two*vxy*MzXZ + two*vxz*MzXY + two*vx*MzXYZ +     vyz*MzXX +     vy*MzXXZ +     vz*MzXXY + MzXXYZ;
         mu121 = vxy2z + two*vxy*MzYZ +     vxz*MzYY +     vx*MzYYZ +     vy2*MzXZ + two*vyz*MzXY + two*vy*MzXYZ +     vz*MzXYY + MzXYYZ;
         mu112 = vxyz2 +     vxy*MzZZ + two*vxz*MzYZ +     vx*MzYZZ + two*vyz*MzXZ +     vy*MzXZZ +     vz2*MzXY + two*vz*MzXYZ + MzXYZZ;
         mu220 = vx2y2 +     vx2*MzYY + four*vxy*MzXY + two*vx*MzXYY +     vy2*MzXX + two*vy*MzXXY + MzXXYY ;
         mu202 = vx2z2 +     vx2*MzZZ + four*vxz*MzXZ + two*vx*MzXZZ +     vz2*MzXX + two*vz*MzXXZ + MzXXZZ;
         mu022 = vy2z2 +     vy2*MzZZ + four*vyz*MzYZ + two*vy*MzYZZ +     vz2*MzYY + two*vz*MzYYZ + MzYYZZ;

         MzXXYY += c1o9 * OORho;
         MzXXZZ += c1o9 * OORho;
         MzYYZZ += c1o9 * OORho;

         //5.
         mu221 = vx2y2z+vx2z*MzYY+vx2*MzYYZ+vy2z*MzXX+vy2*MzXXZ+vz*MzXXYY+MzXXYYZ+
                  four*(vxyz*MzXY+vxy*MzXYZ)+
                  two*(vxz*MzXYY+vx*MzXYYZ+vx2y*MzYZ+vyz*MzXXY+vy*MzXXYZ+vxy2*MzXZ);
         mu212 = vx2yz2 + vx2y*MzZZ + vx2*MzYZZ + vy*MzXXZZ + vz2*MzXXY + vyz2*MzXX+MzXXYZZ+
                  four*(vxyz*MzXZ + vxz*MzXYZ)+
                  two*(vxy*MzXZZ + vxz2*MzXY + vx*MzXYZZ + vyz*MzXXZ + vx2z*MzYZ + vz*MzXXYZ);
         mu122 = vxy2z2 + vxy2*MzZZ + MzXYYZZ + vz2*MzXYY + vy2*MzXZZ + vx*MzYYZZ + vxz2*MzYY+
                  four*(vxyz*MzYZ + vyz*MzXYZ)+
                  two*(vxy*MzYZZ + vxz*MzYYZ + vy2z*MzXZ + vyz2*MzXY + vy*MzXYZZ + vz*MzXYYZ);

         //6.
         mu222 = vx2y2z2 + vx2y2*MzZZ + vx2z2*MzYY + vx2*MzYYZZ + vy2z2*MzXX + vy2*MzXXZZ + vz2*MzXXYY + MzXXYYZZ + eight*vxyz*MzXYZ +
                  two*(vx2y*MzYZZ +  vx*MzXYYZZ + vz*MzXXYYZ + vyz2*MzXXY + vx2z*MzYYZ + vxy2*MzXZZ + vxz2*MzXYY + vy*MzXXYZZ + vy2z*MzXXZ)+ 
                  four*(vxy2z*MzXZ + vx2yz*MzYZ + vxyz2*MzXY + vxy*MzXYZZ + vxz*MzXYYZ + vyz*MzXXYZ);

         (D.f[ dirE   ])[k   ] =   c1o2*rho*( mu200  - mu220 + mu222 - mu202 +  mu120 - mu122 + mu102 - vx   );   //ke
         (D.f[ dirW   ])[kw  ] =   c1o2*rho*( mu200  - mu220 + mu222 - mu202 -  mu120 + mu122 - mu102 + vx   );   
         (D.f[ dirN   ])[k   ] =   c1o2*rho*( mu210  - mu220 + mu222 - mu212 +  mu020 - mu022 + mu012 - vy   );   //kn
         (D.f[ dirS   ])[ks  ] =   c1o2*rho*(-mu210  - mu220 + mu222 + mu212 +  mu020 - mu022 - mu012 + vy   );   
         (D.f[ dirT   ])[k   ] =   c1o2*rho*(-mu221  + mu222 + mu201 - mu202 +  mu021 - mu022 + mu002 - vz   );   //kt
         (D.f[ dirB   ])[kb  ] =   c1o2*rho*( mu221  + mu222 - mu201 - mu202 -  mu021 - mu022 + mu002 + vz   );   
         (D.f[ dirNE  ])[k   ] =  c1o4*rho*(-mu210  + mu220 - mu222 + mu212 +  mu110 - mu120 + mu122 - mu112);   //kne
         (D.f[ dirSW  ])[ksw ] =  c1o4*rho*( mu210  + mu220 - mu222 - mu212 +  mu110 + mu120 - mu122 - mu112);   
         (D.f[ dirSE  ])[ks  ] =  c1o4*rho*( mu210  + mu220 - mu222 - mu212 -  mu110 - mu120 + mu122 + mu112);   //kse
         (D.f[ dirNW  ])[kw  ] =  c1o4*rho*(-mu210  + mu220 - mu222 + mu212 -  mu110 + mu120 - mu122 + mu112);   //knw
         (D.f[ dirTE  ])[k   ] =  c1o4*rho*( mu221  - mu222 - mu201 + mu202 -  mu121 + mu122 + mu101 - mu102);   //kte
         (D.f[ dirBW  ])[kbw ] =  c1o4*rho*(-mu221  - mu222 + mu201 + mu202 -  mu121 - mu122 + mu101 + mu102);   
         (D.f[ dirBE  ])[kb  ] =  c1o4*rho*(-mu221  - mu222 + mu201 + mu202 +  mu121 + mu122 - mu101 - mu102);   //kbe
         (D.f[ dirTW  ])[kw  ] =  c1o4*rho*( mu221  - mu222 - mu201 + mu202 +  mu121 - mu122 - mu101 + mu102);   //ktw
         (D.f[ dirTN  ])[k   ] =  c1o4*rho*( mu221  - mu222 - mu211 + mu212 -  mu021 + mu022 + mu011 - mu012);   //ktn
         (D.f[ dirBS  ])[kbs ] =  c1o4*rho*(-mu221  - mu222 - mu211 - mu212 +  mu021 + mu022 + mu011 + mu012);   
         (D.f[ dirBN  ])[kb  ] =  c1o4*rho*(-mu221  - mu222 + mu211 + mu212 +  mu021 + mu022 - mu011 - mu012);   //kbn
         (D.f[ dirTS  ])[ks  ] =  c1o4*rho*( mu221  - mu222 + mu211 - mu212 -  mu021 + mu022 - mu011 + mu012);   //kts
         (D.f[ dirZERO])[k   ] =       rho*(-mu200  + mu220 - mu222 + mu202 -  mu020 + mu022 - mu002        )+rho0;   //kzero
         (D.f[ dirTNE ])[k   ] = c1o8*rho*(-mu221  + mu222 + mu211 - mu212 +  mu121 - mu122 - mu111 + mu112);   //ktne
         (D.f[ dirTSE ])[ks  ] = c1o8*rho*(-mu221  + mu222 - mu211 + mu212 +  mu121 - mu122 + mu111 - mu112);   //ktse
         (D.f[ dirBNE ])[kb  ] = c1o8*rho*( mu221  + mu222 - mu211 - mu212 -  mu121 - mu122 + mu111 + mu112);   //kbne
         (D.f[ dirBSE ])[kbs ] = c1o8*rho*( mu221  + mu222 + mu211 + mu212 -  mu121 - mu122 - mu111 - mu112);   //kbse
         (D.f[ dirTNW ])[kw  ] = c1o8*rho*(-mu221  + mu222 + mu211 - mu212 -  mu121 + mu122 + mu111 - mu112);   //ktnw
         (D.f[ dirTSW ])[ksw ] = c1o8*rho*(-mu221  + mu222 - mu211 + mu212 -  mu121 + mu122 - mu111 + mu112);   //ktsw
         (D.f[ dirBNW ])[kbw ] = c1o8*rho*( mu221  + mu222 - mu211 - mu212 +  mu121 + mu122 - mu111 - mu112);   //kbnw
         (D.f[ dirBSW ])[kbsw] = c1o8*rho*( mu221  + mu222 + mu211 + mu212 +  mu121 + mu122 + mu111 + mu112);   


         //////////////////////////////////////////////////////////////////////////
         //BGK
         //////////////////////////////////////////////////////////////////////////
         //(D.f[ dirE   ])[k   ] = fW    ;
         //(D.f[ dirW   ])[kw  ] = fE    ;
         //(D.f[ dirN   ])[k   ] = fS    ;
         //(D.f[ dirS   ])[ks  ] = fN    ;
         //(D.f[ dirT   ])[k   ] = fB    ;
         //(D.f[ dirB   ])[kb  ] = fT    ;
         //(D.f[ dirNE  ])[k   ] = fSW   ;
         //(D.f[ dirSW  ])[ksw ] = fNE   ;
         //(D.f[ dirSE  ])[ks  ] = fNW   ;
         //(D.f[ dirNW  ])[kw  ] = fSE   ;
         //(D.f[ dirTE  ])[k   ] = fBW   ;
         //(D.f[ dirBW  ])[kbw ] = fTE   ;
         //(D.f[ dirBE  ])[kb  ] = fTW   ;
         //(D.f[ dirTW  ])[kw  ] = fBE   ;
         //(D.f[ dirTN  ])[k   ] = fBS   ;
         //(D.f[ dirBS  ])[kbs ] = fTN   ;
         //(D.f[ dirBN  ])[kb  ] = fTS   ;
         //(D.f[ dirTS  ])[ks  ] = fBN   ;
         //(D.f[ dirZERO])[k   ] = fZERO ;
         //(D.f[ dirTNE ])[k   ] = fBSW  ;
         //(D.f[ dirTSE ])[ks  ] = fBNW  ;
         //(D.f[ dirBNE ])[kb  ] = fTSW  ;
         //(D.f[ dirBSE ])[kbs ] = fTNW  ;
         //(D.f[ dirTNW ])[kw  ] = fBSE  ;
         //(D.f[ dirTSW ])[ksw ] = fBNE  ;
         //(D.f[ dirBNW ])[kbw ] = fTSE  ;
         //(D.f[ dirBSW ])[kbsw] = fTNE  ;
      }                                                                                                                    
   }
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Casc_SP_MS_Diff_27(real omega,
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
      BC        =   bcMatD[k];

      if( (BC != GEO_SOLID) && (BC != GEO_VOID))
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
         real fE    =  (D.f[dirE   ])[k  ];//ke
         real fW    =  (D.f[dirW   ])[kw ];
         real fN    =  (D.f[dirN   ])[k  ];//kn
         real fS    =  (D.f[dirS   ])[ks ];
         real fT    =  (D.f[dirT   ])[k  ];//kt
         real fB    =  (D.f[dirB   ])[kb ];
         real fNE   =  (D.f[dirNE  ])[k  ];//kne
         real fSW   =  (D.f[dirSW  ])[ksw];
         real fSE   =  (D.f[dirSE  ])[ks ];//kse
         real fNW   =  (D.f[dirNW  ])[kw ];//knw
         real fTE   =  (D.f[dirTE  ])[k  ];//kte
         real fBW   =  (D.f[dirBW  ])[kbw];
         real fBE   =  (D.f[dirBE  ])[kb ];//kbe
         real fTW   =  (D.f[dirTW  ])[kw ];//ktw
         real fTN   =  (D.f[dirTN  ])[k  ];//ktn
         real fBS   =  (D.f[dirBS  ])[kbs];
         real fBN   =  (D.f[dirBN  ])[kb ];//kbn
         real fTS   =  (D.f[dirTS  ])[ks ];//kts
         real fZERO =  (D.f[dirZERO])[k  ];//kzero
         real fTNE   = (D.f[dirTNE ])[k  ];//ktne
         real fTSW   = (D.f[dirTSW ])[ksw];//ktsw
         real fTSE   = (D.f[dirTSE ])[ks ];//ktse
         real fTNW   = (D.f[dirTNW ])[kw ];//ktnw
         real fBNE   = (D.f[dirBNE ])[kb ];//kbne
         real fBSW   = (D.f[dirBSW ])[kbsw];
         real fBSE   = (D.f[dirBSE ])[kbs];//kbse
         real fBNW   = (D.f[dirBNW ])[kbw];//kbnw
         ////////////////////////////////////////////////////////////////////////////////
         real rho0   =  fZERO+fE+fW+fN+fS+fT+fB+fNE+fSW+fSE+fNW+fTE+fBW+fBE+fTW+fTN+fBS+fBN+fTS+fTNE+fTSW+fTSE+fTNW+fBNE+fBSW+fBSE+fBNW;
         real rho    =  rho0 + one;
         real OORho  =  one/rho;
         real vx     =  (fE -fW +fNE-fSW+fSE-fNW+fTE-fBW+fBE-fTW+ fTNE-fTSW+fTSE-fTNW+ fBNE-fBSW+fBSE-fBNW) * OORho;
         real vy     =  (fN -fS +fNE-fSW-fSE+fNW+fTN-fBS+fBN-fTS+ fTNE-fTSW-fTSE+fTNW+ fBNE-fBSW-fBSE+fBNW) * OORho;
         real vz     =  (fT -fB +fTE-fBW-fBE+fTW+fTN-fBS-fBN+fTS+ fTNE+fTSW+fTSE+fTNW- fBNE-fBSW-fBSE-fBNW) * OORho;

         //////////////////////////////////////////////////////////////////////////
         //BGK
         //////////////////////////////////////////////////////////////////////////
         //real drho   =  fZERO+fE+fW+fN+fS+fT+fB+fNE+fSW+fSE+fNW+fTE+fBW+fBE+fTW+fTN+fBS+fBN+fTS+fTNE+fTSW+fTSE+fTNW+fBNE+fBSW+fBSE+fBNW;
         //real vx1     =  (fE -fW +fNE-fSW+fSE-fNW+fTE-fBW+fBE-fTW+ fTNE-fTSW+fTSE-fTNW+ fBNE-fBSW+fBSE-fBNW);
         //real vx2     =  (fN -fS +fNE-fSW-fSE+fNW+fTN-fBS+fBN-fTS+ fTNE-fTSW-fTSE+fTNW+ fBNE-fBSW-fBSE+fBNW);
         //real vx3     =  (fT -fB +fTE-fBW-fBE+fTW+fTN-fBS-fBN+fTS+ fTNE+fTSW+fTSE+fTNW- fBNE-fBSW-fBSE-fBNW);


         //real cusq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         //fZERO = fZERO *(one+(-omega))-(-omega)*   c8over27*  (drho-cusq);
         //fE    = fE    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq);
         //fW    = fW    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq);
         //fN    = fN    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq);
         //fS    = fS    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq);
         //fT    = fT    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq);
         //fB    = fB    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq);
         //fNE   = fNE   *(one+(-omega))-(-omega)*   c1over54*  (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq);
         //fSW   = fSW   *(one+(-omega))-(-omega)*   c1over54*  (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq);
         //fSE   = fSE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq);
         //fNW   = fNW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq);
         //fTE   = fTE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq);
         //fBW   = fBW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq);
         //fBE   = fBE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq);
         //fTW   = fTW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq);
         //fTN   = fTN   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq);
         //fBS   = fBS   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq);
         //fBN   = fBN   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq);
         //fTS   = fTS   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq);
         //fTNE  = fTNE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq);
         //fBSW  = fBSW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq);
         //fBNE  = fBNE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq);
         //fTSW  = fTSW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq);
         //fTSE  = fTSE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq);
         //fBNW  = fBNW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq);
         //fBSE  = fBSE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq);
         //fTNW  = fTNW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq);


         real vx2    = vx*vx;
         real vy2    = vy*vy;
         real vz2    = vz*vz;
         real vxy    = vx*vy;
         real vxz    = vx*vz;
         real vyz    = vy*vz;
         real vx2y   = vx*vx*vy;
         real vx2z   = vx*vx*vz;
         real vy2z   = vy*vy*vz;
         real vxy2   = vx*vy*vy;
         real vyz2   = vy*vz*vz;
         real vxz2   = vx*vz*vz;
         real vxyz   = vx*vy*vz;
         real vx2y2  = vx*vx*vy*vy;
         real vx2z2  = vx*vx*vz*vz;
         real vy2z2  = vy*vy*vz*vz;
         real vx2yz  = vx*vx*vy*vz;
         real vxyz2  = vx*vy*vz*vz;
         real vxy2z  = vx*vy*vy*vz;
         real vx2y2z = vx*vx*vy*vy*vz;
         real vx2yz2 = vx*vx*vy*vz*vz;
         real vxy2z2 = vx*vy*vy*vz*vz;
         real vx2y2z2= vx*vx*vy*vy*vz*vz;
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         real mu200   =(fE+fW+fNE+fSW+fSE+fNW+fTE+fBW+fBE+fTW+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         real mu020   =(fN+fS+fNE+fSW+fSE+fNW+fTN+fBS+fBN+fTS+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         real mu002   =(fT+fB+fTE+fBW+fBE+fTW+fTN+fBS+fBN+fTS+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         real mu110   =(fNE-fSE+fSW-fNW+fTNE-fTSE+fBNE-fBSE+fTSW-fTNW+fBSW-fBNW) * OORho;
         real mu101   =(fTE+fBW-fBE-fTW+fTNE-fBNE+fTSE-fBSE-fTNW+fBNW-fTSW+fBSW) * OORho;
         real mu011   =(fTN+fBS-fBN-fTS+fTNE-fBNE-fTSE+fBSE+fTNW-fBNW-fTSW+fBSW) * OORho;
         real mu210   =(fNE-fSW-fSE+fNW+fTNE+fBNE-fTSE-fBSE+fTNW+fBNW-fTSW-fBSW) * OORho;
         real mu120   =(fNE-fSW+fSE-fNW+fTNE-fTNW+fBNE-fBNW+fTSE-fTSW+fBSE-fBSW) * OORho;
         real mu102   =(fTE-fBW+fBE-fTW+fTNE-fTNW+fBNE-fBNW+fTSE-fTSW+fBSE-fBSW) * OORho;
         real mu111   =(fTNE-fBNE-fTSE+fBSE-fTNW+fBNW+fTSW-fBSW) * OORho;
         real mu201   =(fTE-fBW-fBE+fTW+fTNE-fBNE+fTSE-fBSE+fTNW-fBNW+fTSW-fBSW) * OORho;
         real mu021   =(fTN-fBS-fBN+fTS+fTNE-fBNE+fTSE-fBSE+fTNW-fBNW+fTSW-fBSW) * OORho;
         real mu012   =(fTN-fBS+fBN-fTS+fTNE+fBNE-fTSE-fBSE+fTNW+fBNW-fTSW-fBSW) * OORho;
         real mu220   =(fNE+fSW+fSE+fNW+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         real mu121   =(fTNE-fBNE+fTSE-fBSE-fTNW+fBNW-fTSW+fBSW) * OORho;
         real mu202   =(fTE+fBW+fBE+fTW+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         real mu211   =(fTNE-fBNE-fTSE+fBSE+fTNW-fBNW-fTSW+fBSW) * OORho;
         real mu112   =(fTNE+fBNE-fTSE-fBSE-fTNW-fBNW+fTSW+fBSW) * OORho;
         real mu022   =(fTN+fBS+fBN+fTS+fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         real mu221   =(fTNE-fBNE+fTSE-fBSE+fTNW-fBNW+fTSW-fBSW) * OORho;
         real mu122   =(fTNE-fTNW+fBNE-fBNW+fTSE-fTSW+fBSE-fBSW) * OORho;
         real mu212   =(fTNE+fBNE-fTSE-fBSE+fTNW+fBNW-fTSW-fBSW) * OORho;
         real mu222   =(fTNE+fBNE+fTSE+fBSE+fTNW+fBNW+fTSW+fBSW) * OORho;
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         real MzXX,MzYY,MzZZ,MzXY,MzXZ,MzYZ,MzXXY,MzXYY,MzXXZ,MzXZZ,MzYYZ,MzYZZ,MzXYZ,MzXXYY,MzXXZZ,MzYYZZ,MzXXYZ,MzXYYZ,MzXYZZ,MzXXYYZ,MzXXYZZ,MzXYYZZ,MzXXYYZZ;
         {
            //2.
            MzXX      =   mu200-vx2;
            MzXY      =   mu110-vxy;
            MzXZ      =   mu101-vxz;
            MzYY      =   mu020-vy2;
            MzYZ      =   mu011-vyz;
            MzZZ      =   mu002-vz2;

            real pimpmu200 = mu200 + c1o3 * OORho;
            real pimpmu020 = mu020 + c1o3 * OORho;
            real pimpmu002 = mu002 + c1o3 * OORho;

            //3.
            MzXXY     =    two*vx2y - two*vx*mu110 - vy*pimpmu200 + mu210; 
            MzXXZ     =    two*vx2z - two*vx*mu101 - vz*pimpmu200 + mu201; 
            MzXYY     =    two*vxy2 - pimpmu020*vx - two*vy*mu110 + mu120; 
            MzXYZ     =    two*vxyz - mu011*vx-vy*mu101 -vz*mu110 + mu111;
            MzXZZ     =    two*vxz2 - pimpmu002*vx - two*vz*mu101 + mu102; 
            MzYYZ     =    two*vy2z - two*vy*mu011 - vz*pimpmu020 + mu021; 
            MzYZZ     =    two*vyz2 - pimpmu002*vy - two*vz*mu011 + mu012; 

            //4.
            MzXXYY    =   -three*vx2y2+pimpmu020*vx2+four*vxy*mu110-two*vx*mu120+vy2*pimpmu200-two*vy*mu210+mu220; 
            MzXXYZ    =   -three*vx2yz+mu011*vx2+two*vxy*mu101+two*vxz*mu110-two*vx*mu111+vyz*pimpmu200-vy*mu201-vz*mu210+mu211; 
            MzXXZZ    =   -three*vx2z2+pimpmu002*vx2+four*vxz*mu101-two*vx*mu102+vz2*pimpmu200-two*vz*mu201+mu202; 
            MzXYYZ    =   -three*vxy2z+two*vxy*mu011+vxz*pimpmu020-mu021*vx+vy2*mu101+two*vyz*mu110-two*vy*mu111-vz*mu120+mu121; 
            MzXYZZ    =   -three*vxyz2+pimpmu002*vxy+two*vxz*mu011-mu012*vx+two*vyz*mu101-vy*mu102+vz2*mu110-two*vz*mu111+mu112; 
            MzYYZZ    =   -three*vy2z2+pimpmu002*vy2+four*vyz*mu011-two*vy*mu012+vz2*pimpmu020-two*vz*mu021+mu022; 

            real pimpmu220 = mu220 + c1o9 * OORho;
            real pimpmu202 = mu202 + c1o9 * OORho;
            real pimpmu022 = mu022 + c1o9 * OORho;

            //5.
            MzXXYYZ   =    four*vx2y2z-four*vxyz*mu110+four*vxy*mu111+ 
               two*(vxz*mu120-vxy2*mu101-vx2y*mu011+vyz*mu210-vy*mu211-vx*mu121)+
               vy2*mu201-vx2z*pimpmu020+mu021*vx2-vz*pimpmu220-vy2z*pimpmu200+mu221; 
            MzXXYZZ   =    four*(vx2yz2-vxyz*mu101+vxz*mu111)+
               two*(vxy*mu102-vxz2*mu110-vx2z*mu011-vx*mu112-vz*mu211+vyz*mu201)+
               vz2*mu210-vy*pimpmu202-pimpmu002*vx2y+mu012*vx2-vyz2*pimpmu200+mu212;
            MzXYYZZ   =    four*(vxy2z2-vxyz*mu011+vyz*mu111)+
               two*(vxy*mu012-vyz2*mu110-vy*mu112+vxz*mu021-vz*mu121-vy2z*mu101)+
               vy2*mu102+vz2*mu120-vxz2*pimpmu020-pimpmu022*vx-pimpmu002*vxy2+mu122;    

            //6.
            MzXXYYZZ  =   -five*vx2y2z2-eight*vxyz*mu111+vy2*pimpmu202+pimpmu002*vx2y2+vy2z2*pimpmu200+vx2z2*pimpmu020+vz2*pimpmu220+mu222+pimpmu022*vx2+
               four*(vx2yz*mu011+vxy2z*mu101+vxyz2*mu110+vxy*mu112+vxz*mu121+vyz*mu211)- 
               two*(vy*mu212-vy2z*mu201-vxz2*mu120-vyz2*mu210-vxy2*mu102-vx2z*mu021-vx2y*mu012-vx*mu122-vz*mu221); 
         }
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         real MXXpMYYpMZZ,MXXmMYY, MXXmMZZ, MXXYpMYZZ,MXXZpMYYZ,MXYYpMXZZ, MXXYmMYZZ,MXXZmMYYZ,MXYYmMXZZ;
         {
            //coll faktoren:
            real w1 =  -omega;
            real w2 = -one;//-omega;
            real w3 = -one;
            real w4 = -one;
            real w5 = -one;
            real w6 = -one;
            real w7 = -one;
            real w8 = -one;
            real w9 = -one;
            real w10= -one;

            ////////lin kombi bilden:
            MXXpMYYpMZZ = MzXX + MzYY + MzZZ;
            MXXmMYY     =  MzXX - MzYY;
            MXXmMZZ     =  MzXX - MzZZ;

            MXXYpMYZZ   =  MzXXY+MzYZZ;
            MXXYmMYZZ   =  MzXXY-MzYZZ;
            MXXZpMYYZ   =  MzXXZ+MzYYZ;
            MXXZmMYYZ   =  MzXXZ-MzYYZ;
            MXYYpMXZZ   =  MzXYY+MzXZZ;
            MXYYmMXZZ   =  MzXYY-MzXZZ;

            real MXXYYppp    = MzXXYY + MzXXZZ + MzYYZZ;
            real MXXYYpm2p   = MzXXYY - two*MzXXZZ + MzYYZZ;
            real MXXYYppm2   = MzXXYY + MzXXZZ - two*MzYYZZ;

            //relaxation:
            MXXpMYYpMZZ -= w2*(one-OORho-MXXpMYYpMZZ);

            MzXZ       *=  one+w1;
            MzYZ       *=  one+w1;
            MzXY       *=  one+w1;

            MzXYZ      *=  one+w5;

            MXXmMYY    *=  one+w1;
            MXXmMZZ    *=  one+w1;

            MXXYpMYZZ  *=  one+w3;
            MXXYmMYZZ  *=  one+w4;
            MXXZpMYYZ  *=  one+w3;
            MXXZmMYYZ  *=  one+w4;
            MXYYpMXZZ  *=  one+w3;
            MXYYmMXZZ  *=  one+w4;

            ////////von Lin Kombis zurueck:
            MzXX  =  c1o3 * (       MXXmMYY +       MXXmMZZ + MXXpMYYpMZZ + OORho);
            MzYY  =  c1o3 * (-two * MXXmMYY +       MXXmMZZ + MXXpMYYpMZZ + OORho);
            MzZZ  =  c1o3 * (       MXXmMYY - two * MXXmMZZ + MXXpMYYpMZZ + OORho);

            MzXXY = (MXXYmMYZZ + MXXYpMYZZ)*c1o2;
            MzYZZ = (MXXYpMYZZ - MXXYmMYZZ)*c1o2;
            MzXYY = (MXYYmMXZZ + MXYYpMXZZ)*c1o2;
            MzXZZ = (MXYYpMXZZ - MXYYmMXZZ)*c1o2;
            MzXXZ = (MXXZmMYYZ + MXXZpMYYZ)*c1o2;
            MzYYZ = (MXXZpMYYZ - MXXZmMYYZ)*c1o2;

            //faktorisierte atraktoren:
            MXXYYppp  -=  w7*MzXX*MzYY + w7*MzXX*MzZZ     + w7*    MzZZ*MzYY - w7*MXXYYppp - w7*c1o3;
            MXXYYpm2p -=  w6*MzXX*MzYY - w6*two*MzXX*MzZZ + w6*    MzZZ*MzYY - w6*MXXYYpm2p;
            MXXYYppm2 -=  w6*MzXX*MzYY + w6*    MzXX*MzZZ - w6*two*MzZZ*MzYY - w6*MXXYYppm2;
            MzXXYYZZ  -= w10*MzXX*MzYY*MzZZ - w10*c1o27 - w10*MzXXYYZZ;
            MzXYYZ    -=  w8*MzYY*MzXZ - w8*MzXYYZ;
            MzXYZZ    -=  w8*MzZZ*MzXY - w8*MzXYZZ;
            MzXXYZ    -=  w8*MzXX*MzYZ - w8*MzXXYZ;

            MzXXYYZ *= one+w9;
            MzXXYZZ *= one+w9;
            MzXYYZZ *= one+w9;

            MzXXYY =  c1o3 * (MXXYYpm2p + MXXYYppm2 + MXXYYppp);
            MzXXZZ =  c1o3 * (MXXYYppp - MXXYYpm2p);
            MzYYZZ =  c1o3 * (MXXYYppp-MXXYYppm2);
         }

         //2.
         mu200 -= vx2 + c1o3 * (MXXmMYY + MXXmMZZ + MXXpMYYpMZZ);
         mu020 -= vy2 + c1o3 * (MXXmMZZ +  MXXpMYYpMZZ - two * MXXmMYY);
         mu002 -= vz2 + c1o3 * (MXXmMYY - two * MXXmMZZ + MXXpMYYpMZZ);
         mu110 -= vxy + MzXY;
         mu101 -= vxz + MzXZ;
         mu011 -= vyz + MzYZ;

         //3.
         mu111 -= vxyz +     vx*MzYZ +     vy*MzXZ + vz*MzXY + MzXYZ;
         mu210 -= vx2y + two*vx*MzXY +     vy*MzXX + MzXXY;
         mu120 -= vxy2 +     vx*MzYY + two*vy*MzXY + MzXYY;
         mu102 -= vxz2 +     vx*MzZZ + two*vz*MzXZ + MzXZZ;
         mu201 -= vx2z + two*vx*MzXZ +     vz*MzXX + MzXXZ;
         mu021 -= vy2z + two*vy*MzYZ +     vz*MzYY + MzYYZ;
         mu012 -= vyz2 +     vy*MzZZ + two*vz*MzYZ + MzYZZ;

         //4.
         mu211 -= vx2yz +     vx2*MzYZ + two*vxy*MzXZ + two*vxz*MzXY + two*vx*MzXYZ +     vyz*MzXX +     vy*MzXXZ +     vz*MzXXY + MzXXYZ;
         mu121 -= vxy2z + two*vxy*MzYZ +     vxz*MzYY +     vx*MzYYZ +     vy2*MzXZ + two*vyz*MzXY + two*vy*MzXYZ +     vz*MzXYY + MzXYYZ;
         mu112 -= vxyz2 +     vxy*MzZZ + two*vxz*MzYZ +     vx*MzYZZ + two*vyz*MzXZ +     vy*MzXZZ +     vz2*MzXY + two*vz*MzXYZ + MzXYZZ;
         mu220 -= vx2y2 +     vx2*MzYY + four*vxy*MzXY + two*vx*MzXYY +     vy2*MzXX + two*vy*MzXXY + MzXXYY ;
         mu202 -= vx2z2 +     vx2*MzZZ + four*vxz*MzXZ + two*vx*MzXZZ +     vz2*MzXX + two*vz*MzXXZ + MzXXZZ;
         mu022 -= vy2z2 +     vy2*MzZZ + four*vyz*MzYZ + two*vy*MzYZZ +     vz2*MzYY + two*vz*MzYYZ + MzYYZZ;

         MzXXYY += c1o9 * OORho;
         MzXXZZ += c1o9 * OORho;
         MzYYZZ += c1o9 * OORho;

         //5.
         mu221 -= vx2y2z+vx2z*MzYY+vx2*MzYYZ+vy2z*MzXX+vy2*MzXXZ+vz*MzXXYY+MzXXYYZ+
            four*(vxyz*MzXY+vxy*MzXYZ)+
            two*(vxz*MzXYY+vx*MzXYYZ+vx2y*MzYZ+vyz*MzXXY+vy*MzXXYZ+vxy2*MzXZ);
         mu212 -= vx2yz2 + vx2y*MzZZ + vx2*MzYZZ + vy*MzXXZZ + vz2*MzXXY + vyz2*MzXX+MzXXYZZ+
            four*(vxyz*MzXZ + vxz*MzXYZ)+
            two*(vxy*MzXZZ + vxz2*MzXY + vx*MzXYZZ + vyz*MzXXZ + vx2z*MzYZ + vz*MzXXYZ);
         mu122 -= vxy2z2 + vxy2*MzZZ + MzXYYZZ + vz2*MzXYY + vy2*MzXZZ + vx*MzYYZZ + vxz2*MzYY+
            four*(vxyz*MzYZ + vyz*MzXYZ)+
            two*(vxy*MzYZZ + vxz*MzYYZ + vy2z*MzXZ + vyz2*MzXY + vy*MzXYZZ + vz*MzXYYZ);

         //6.
         mu222 -= vx2y2z2 + vx2y2*MzZZ + vx2z2*MzYY + vx2*MzYYZZ + vy2z2*MzXX + vy2*MzXXZZ + vz2*MzXXYY + MzXXYYZZ + eight*vxyz*MzXYZ +
            two*(vx2y*MzYZZ +  vx*MzXYYZZ + vz*MzXXYYZ + vyz2*MzXXY + vx2z*MzYYZ + vxy2*MzXZZ + vxz2*MzXYY + vy*MzXXYZZ + vy2z*MzXXZ)+ 
            four*(vxy2z*MzXZ + vx2yz*MzYZ + vxyz2*MzXY + vxy*MzXYZZ + vxz*MzXYYZ + vyz*MzXXYZ);

         (D.f[ dirE   ])[k   ] = fW    -   c1o2*rho*( mu200  - mu220 + mu222 - mu202 +  mu120 - mu122 + mu102        );   //ke
         (D.f[ dirW   ])[kw  ] = fE    -   c1o2*rho*( mu200  - mu220 + mu222 - mu202 -  mu120 + mu122 - mu102        );   
         (D.f[ dirN   ])[k   ] = fS    -   c1o2*rho*( mu210  - mu220 + mu222 - mu212 +  mu020 - mu022 + mu012        );   //kn
         (D.f[ dirS   ])[ks  ] = fN    -   c1o2*rho*(-mu210  - mu220 + mu222 + mu212 +  mu020 - mu022 - mu012        );   
         (D.f[ dirT   ])[k   ] = fB    -   c1o2*rho*(-mu221  + mu222 + mu201 - mu202 +  mu021 - mu022 + mu002        );   //kt
         (D.f[ dirB   ])[kb  ] = fT    -   c1o2*rho*( mu221  + mu222 - mu201 - mu202 -  mu021 - mu022 + mu002        );   
         (D.f[ dirNE  ])[k   ] = fSW   -  c1o4*rho*(-mu210  + mu220 - mu222 + mu212 +  mu110 - mu120 + mu122 - mu112);   //kne
         (D.f[ dirSW  ])[ksw ] = fNE   -  c1o4*rho*( mu210  + mu220 - mu222 - mu212 +  mu110 + mu120 - mu122 - mu112);   
         (D.f[ dirSE  ])[ks  ] = fNW   -  c1o4*rho*( mu210  + mu220 - mu222 - mu212 -  mu110 - mu120 + mu122 + mu112);   //kse
         (D.f[ dirNW  ])[kw  ] = fSE   -  c1o4*rho*(-mu210  + mu220 - mu222 + mu212 -  mu110 + mu120 - mu122 + mu112);   //knw
         (D.f[ dirTE  ])[k   ] = fBW   -  c1o4*rho*( mu221  - mu222 - mu201 + mu202 -  mu121 + mu122 + mu101 - mu102);   //kte
         (D.f[ dirBW  ])[kbw ] = fTE   -  c1o4*rho*(-mu221  - mu222 + mu201 + mu202 -  mu121 - mu122 + mu101 + mu102);   
         (D.f[ dirBE  ])[kb  ] = fTW   -  c1o4*rho*(-mu221  - mu222 + mu201 + mu202 +  mu121 + mu122 - mu101 - mu102);   //kbe
         (D.f[ dirTW  ])[kw  ] = fBE   -  c1o4*rho*( mu221  - mu222 - mu201 + mu202 +  mu121 - mu122 - mu101 + mu102);   //ktw
         (D.f[ dirTN  ])[k   ] = fBS   -  c1o4*rho*( mu221  - mu222 - mu211 + mu212 -  mu021 + mu022 + mu011 - mu012);   //ktn
         (D.f[ dirBS  ])[kbs ] = fTN   -  c1o4*rho*(-mu221  - mu222 - mu211 - mu212 +  mu021 + mu022 + mu011 + mu012);   
         (D.f[ dirBN  ])[kb  ] = fTS   -  c1o4*rho*(-mu221  - mu222 + mu211 + mu212 +  mu021 + mu022 - mu011 - mu012);   //kbn
         (D.f[ dirTS  ])[ks  ] = fBN   -  c1o4*rho*( mu221  - mu222 + mu211 - mu212 -  mu021 + mu022 - mu011 + mu012);   //kts
         (D.f[ dirZERO])[k   ] = fZERO -       rho*(-mu200  + mu220 - mu222 + mu202 -  mu020 + mu022 - mu002        );   //kzero
         (D.f[ dirTNE ])[k   ] = fBSW  - c1o8*rho*(-mu221  + mu222 + mu211 - mu212 +  mu121 - mu122 - mu111 + mu112);   //ktne
         (D.f[ dirTSE ])[ks  ] = fBNW  - c1o8*rho*(-mu221  + mu222 - mu211 + mu212 +  mu121 - mu122 + mu111 - mu112);   //ktse
         (D.f[ dirBNE ])[kb  ] = fTSW  - c1o8*rho*( mu221  + mu222 - mu211 - mu212 -  mu121 - mu122 + mu111 + mu112);   //kbne
         (D.f[ dirBSE ])[kbs ] = fTNW  - c1o8*rho*( mu221  + mu222 + mu211 + mu212 -  mu121 - mu122 - mu111 - mu112);   //kbse
         (D.f[ dirTNW ])[kw  ] = fBSE  - c1o8*rho*(-mu221  + mu222 + mu211 - mu212 -  mu121 + mu122 + mu111 - mu112);   //ktnw
         (D.f[ dirTSW ])[ksw ] = fBNE  - c1o8*rho*(-mu221  + mu222 - mu211 + mu212 -  mu121 + mu122 - mu111 + mu112);   //ktsw
         (D.f[ dirBNW ])[kbw ] = fTSE  - c1o8*rho*( mu221  + mu222 - mu211 - mu212 +  mu121 + mu122 - mu111 - mu112);   //kbnw
         (D.f[ dirBSW ])[kbsw] = fTNE  - c1o8*rho*( mu221  + mu222 + mu211 + mu212 +  mu121 + mu122 + mu111 + mu112);   


         //////////////////////////////////////////////////////////////////////////
         //BGK
         //////////////////////////////////////////////////////////////////////////
         //(D.f[ dirE   ])[k   ] = fW    ;
         //(D.f[ dirW   ])[kw  ] = fE    ;
         //(D.f[ dirN   ])[k   ] = fS    ;
         //(D.f[ dirS   ])[ks  ] = fN    ;
         //(D.f[ dirT   ])[k   ] = fB    ;
         //(D.f[ dirB   ])[kb  ] = fT    ;
         //(D.f[ dirNE  ])[k   ] = fSW   ;
         //(D.f[ dirSW  ])[ksw ] = fNE   ;
         //(D.f[ dirSE  ])[ks  ] = fNW   ;
         //(D.f[ dirNW  ])[kw  ] = fSE   ;
         //(D.f[ dirTE  ])[k   ] = fBW   ;
         //(D.f[ dirBW  ])[kbw ] = fTE   ;
         //(D.f[ dirBE  ])[kb  ] = fTW   ;
         //(D.f[ dirTW  ])[kw  ] = fBE   ;
         //(D.f[ dirTN  ])[k   ] = fBS   ;
         //(D.f[ dirBS  ])[kbs ] = fTN   ;
         //(D.f[ dirBN  ])[kb  ] = fTS   ;
         //(D.f[ dirTS  ])[ks  ] = fBN   ;
         //(D.f[ dirZERO])[k   ] = fZERO ;
         //(D.f[ dirTNE ])[k   ] = fBSW  ;
         //(D.f[ dirTSE ])[ks  ] = fBNW  ;
         //(D.f[ dirBNE ])[kb  ] = fTSW  ;
         //(D.f[ dirBSE ])[kbs ] = fTNW  ;
         //(D.f[ dirTNW ])[kw  ] = fBSE  ;
         //(D.f[ dirTSW ])[ksw ] = fBNE  ;
         //(D.f[ dirBNW ])[kbw ] = fTSE  ;
         //(D.f[ dirBSW ])[kbsw] = fTNE  ;
      }                                                                                                                    
   }
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Casc_SP_27(  real omega,
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
      BC        =   bcMatD[k];

      if( (BC != GEO_SOLID) && (BC != GEO_VOID))
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
       //unsigned int nxny = nx*ny;
       //unsigned int kzero= k;
       //unsigned int ke   = k;
       //unsigned int kw   = k + 1;
       //unsigned int kn   = k;
       //unsigned int ks   = k + nx;
       //unsigned int kt   = k;
       //unsigned int kb   = k + nxny;
       //unsigned int ksw  = k + nx + 1;
       //unsigned int kne  = k;
       //unsigned int kse  = k + nx;
       //unsigned int knw  = k + 1;
       //unsigned int kbw  = k + nxny + 1;
       //unsigned int kte  = k;
       //unsigned int kbe  = k + nxny;
       //unsigned int ktw  = k + 1;
       //unsigned int kbs  = k + nxny + nx;
       //unsigned int ktn  = k;
       //unsigned int kbn  = k + nxny;
       //unsigned int kts  = k + nx;
       //unsigned int ktse = k + nx;
       //unsigned int kbnw = k + nxny + 1;
       //unsigned int ktnw = k + 1;
       //unsigned int kbse = k + nxny + nx;
       //unsigned int ktsw = k + nx + 1;
       //unsigned int kbne = k + nxny;
       //unsigned int ktne = k;
       //unsigned int kbsw = k + nxny + nx + 1;
       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       real f_E,f_W,f_N,f_S,f_T,f_B,f_NE,f_SW,f_SE,f_NW,f_TE,f_BW,f_BE,f_TW,f_TN,f_BS,f_BN,f_TS,f_ZERO, f_TNE,f_TNW,f_TSE,f_TSW, f_BNE,f_BNW,f_BSE,f_BSW;
       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       f_E    =  (D.f[dirE   ])[ke]+c2over27;
       f_W    =  (D.f[dirW   ])[kw]+c2over27;
       f_N    =  (D.f[dirN   ])[kn]+c2over27;
       f_S    =  (D.f[dirS   ])[ks]+c2over27;
       f_T    =  (D.f[dirT   ])[kt]+c2over27;
       f_B    =  (D.f[dirB   ])[kb]+c2over27;
       f_NE   =  (D.f[dirNE  ])[kne]+c1over54;
       f_SW   =  (D.f[dirSW  ])[ksw]+c1over54;
       f_SE   =  (D.f[dirSE  ])[kse]+c1over54;
       f_NW   =  (D.f[dirNW  ])[knw]+c1over54;
       f_TE   =  (D.f[dirTE  ])[kte]+c1over54;
       f_BW   =  (D.f[dirBW  ])[kbw]+c1over54;
       f_BE   =  (D.f[dirBE  ])[kbe]+c1over54;
       f_TW   =  (D.f[dirTW  ])[ktw]+c1over54;
       f_TN   =  (D.f[dirTN  ])[ktn]+c1over54;
       f_BS   =  (D.f[dirBS  ])[kbs]+c1over54;
       f_BN   =  (D.f[dirBN  ])[kbn]+c1over54;
       f_TS   =  (D.f[dirTS  ])[kts]+c1over54;
       f_ZERO =  (D.f[dirZERO])[kzero]+c8over27;
       f_TNE   = (D.f[dirTNE ])[ktne]+c1over216;
       f_TSW   = (D.f[dirTSW ])[ktsw]+c1over216;
       f_TSE   = (D.f[dirTSE ])[ktse]+c1over216;
       f_TNW   = (D.f[dirTNW ])[ktnw]+c1over216;
       f_BNE   = (D.f[dirBNE ])[kbne]+c1over216;
       f_BSW   = (D.f[dirBSW ])[kbsw]+c1over216;
       f_BSE   = (D.f[dirBSE ])[kbse]+c1over216;
       f_BNW   = (D.f[dirBNW ])[kbnw]+c1over216;
       ////////////////////////////////////////////////////////////////////////////////

       if( BC == GEO_FLUID || BC == GEO_VELO)
       {
          {
         //    /////////////////////////version 13.04.
         //   real drho     = f_ZERO+f_E +f_W +f_N +f_S +f_T +f_B +f_NE+f_SW+f_SE+f_NW+f_TE+f_BW+f_BE+f_TW+f_TN+f_BS+f_BN+f_TS+ f_TNE+f_TSW+f_TSE+f_TNW+ f_BNE+f_BSW+f_BSE+f_BNW;
         //   real vx1     =  f_E -f_W +f_NE-f_SW+f_SE-f_NW+f_TE-f_BW+f_BE-f_TW+ f_TNE-f_TSW+f_TSE-f_TNW+ f_BNE-f_BSW+f_BSE-f_BNW;
         //   real vx2     =  f_N -f_S +f_NE-f_SW-f_SE+f_NW+f_TN-f_BS+f_BN-f_TS+ f_TNE-f_TSW-f_TSE+f_TNW+ f_BNE-f_BSW-f_BSE+f_BNW;
         //   real vx3     =  f_T -f_B +f_TE-f_BW-f_BE+f_TW+f_TN-f_BS-f_BN+f_TS+ f_TNE+f_TSW+f_TSE+f_TNW- f_BNE-f_BSW-f_BSE-f_BNW;
         ////   vx1=vx1/(drho); vx2=vx2/(drho); vx3=vx3/(drho);
         //   //////////////////////////////////////////////////////////////////////////
         //   //BGK
         //   //////////////////////////////////////////////////////////////////////////


         //   real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         //   f_ZERO = f_ZERO *(one+(-omega))-(-omega)*   c8over27*  (drho-cu_sq);
         //   f_E    = f_E    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
         //   f_W    = f_W    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
         //   f_N    = f_N    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
         //   f_S    = f_S    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
         //   f_T    = f_T    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
         //   f_B    = f_B    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
         //   f_NE   = f_NE   *(one+(-omega))-(-omega)*   c1over54*  (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         //   f_SW   = f_SW   *(one+(-omega))-(-omega)*   c1over54*  (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         //   f_SE   = f_SE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         //   f_NW   = f_NW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         //   f_TE   = f_TE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         //   f_BW   = f_BW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         //   f_BE   = f_BE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         //   f_TW   = f_TW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         //   f_TN   = f_TN   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         //   f_BS   = f_BS   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         //   f_BN   = f_BN   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         //   f_TS   = f_TS   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         //   f_TNE  = f_TNE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         //   f_BSW  = f_BSW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         //   f_BNE  = f_BNE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         //   f_TSW  = f_TSW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         //   f_TSE  = f_TSE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         //   f_BNW  = f_BNW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         //   f_BSE  = f_BSE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         //   f_TNW  = f_TNW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);

            ///////////////////////version 13.04.
            real rho   =  f_ZERO+f_E +f_W +f_N +f_S +f_T +f_B +f_NE+f_SW+f_SE+f_NW+f_TE+f_BW+f_BE+f_TW+f_TN+f_BS+f_BN+f_TS+ f_TNE+f_TSW+f_TSE+f_TNW+ f_BNE+f_BSW+f_BSE+f_BNW;
            //real rho    =  rho0 + one;
            real vx     =  f_E -f_W +f_NE-f_SW+f_SE-f_NW+f_TE-f_BW+f_BE-f_TW+ f_TNE-f_TSW+f_TSE-f_TNW+ f_BNE-f_BSW+f_BSE-f_BNW;
            real vy     =  f_N -f_S +f_NE-f_SW-f_SE+f_NW+f_TN-f_BS+f_BN-f_TS+ f_TNE-f_TSW-f_TSE+f_TNW+ f_BNE-f_BSW-f_BSE+f_BNW;
            real vz     =  f_T -f_B +f_TE-f_BW-f_BE+f_TW+f_TN-f_BS-f_BN+f_TS+ f_TNE+f_TSW+f_TSE+f_TNW- f_BNE-f_BSW-f_BSE-f_BNW;
            vx=vx/(rho); vy=vy/(rho); vz=vz/(rho);

            real mu[3][3][3];

            mu[2][0][0]=(f_E+f_W+f_NE+f_SW+f_SE+f_NW+f_TE+f_BW+f_BE+f_TW+f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o3);
            mu[0][2][0]=(f_N+f_S+f_NE+f_SW+f_SE+f_NW+f_TN+f_BS+f_BN+f_TS+f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o3);
            mu[0][0][2]=(f_T+f_B+f_TE+f_BW+f_BE+f_TW+f_TN+f_BS+f_BN+f_TS+f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o3);
            mu[1][1][0]=(f_NE-f_SE+f_SW-f_NW+f_TNE-f_TSE+f_BNE-f_BSE+f_TSW-f_TNW+f_BSW-f_BNW)/(rho);
            mu[1][0][1]=(f_TE+f_BW-f_BE-f_TW+f_TNE-f_BNE+f_TSE-f_BSE-f_TNW+f_BNW-f_TSW+f_BSW)/(rho);
            mu[0][1][1]=(f_TN+f_BS-f_BN-f_TS+f_TNE-f_BNE-f_TSE+f_BSE+f_TNW-f_BNW-f_TSW+f_BSW)/(rho);
            mu[2][1][0]=(f_NE-f_SW-f_SE+f_NW+f_TNE+f_BNE-f_TSE-f_BSE+f_TNW+f_BNW-f_TSW-f_BSW)/(rho);
            mu[1][2][0]=(f_NE-f_SW+f_SE-f_NW+f_TNE-f_TNW+f_BNE-f_BNW+f_TSE-f_TSW+f_BSE-f_BSW)/(rho);
            mu[1][0][2]=(f_TE-f_BW+f_BE-f_TW+f_TNE-f_TNW+f_BNE-f_BNW+f_TSE-f_TSW+f_BSE-f_BSW)/(rho);
            mu[1][1][1]=(f_TNE-f_BNE-f_TSE+f_BSE-f_TNW+f_BNW+f_TSW-f_BSW)/(rho);
            mu[2][0][1]=(f_TE-f_BW-f_BE+f_TW+f_TNE-f_BNE+f_TSE-f_BSE+f_TNW-f_BNW+f_TSW-f_BSW)/(rho);
            mu[0][2][1]=(f_TN-f_BS-f_BN+f_TS+f_TNE-f_BNE+f_TSE-f_BSE+f_TNW-f_BNW+f_TSW-f_BSW)/(rho);
            mu[0][1][2]=(f_TN-f_BS+f_BN-f_TS+f_TNE+f_BNE-f_TSE-f_BSE+f_TNW+f_BNW-f_TSW-f_BSW)/(rho);
            mu[2][2][0]=(f_NE+f_SW+f_SE+f_NW+f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o9);
            mu[1][2][1]=(f_TNE-f_BNE+f_TSE-f_BSE-f_TNW+f_BNW-f_TSW+f_BSW)/(rho);
            mu[2][0][2]=(f_TE+f_BW+f_BE+f_TW+f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o9);
            mu[2][1][1]=(f_TNE-f_BNE-f_TSE+f_BSE+f_TNW-f_BNW-f_TSW+f_BSW)/(rho);
            mu[1][1][2]=(f_TNE+f_BNE-f_TSE-f_BSE-f_TNW-f_BNW+f_TSW+f_BSW)/(rho);
            mu[0][2][2]=(f_TN+f_BS+f_BN+f_TS+f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o9);
            mu[2][2][1]=(f_TNE-f_BNE+f_TSE-f_BSE+f_TNW-f_BNW+f_TSW-f_BSW)/(rho);
            mu[1][2][2]=(f_TNE-f_TNW+f_BNE-f_BNW+f_TSE-f_TSW+f_BSE-f_BSW)/(rho);
            mu[2][1][2]=(f_TNE+f_BNE-f_TSE-f_BSE+f_TNW+f_BNW-f_TSW-f_BSW)/(rho);
            mu[2][2][2]=(f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o27);
            mu[0][0][0]=one;
            mu[1][0][0]=vx;
            mu[0][1][0]=vy;
            mu[0][0][1]=vz;

            real M_zXX=zero;      real M_zYY=zero;      real M_zZZ=zero;      real M_zXY=zero;     real M_zXZ=zero;   real M_zYZ=zero;
            real M_zXXY=zero;     real M_zXYY=zero;     real M_zXXZ=zero;     real M_zXZZ=zero;    real M_zYYZ=zero;  real M_zYZZ=zero;  real M_zXYZ=zero;
            real M_zXXYY=zero;    real M_zXXZZ=zero;    real M_zYYZZ=zero;    real M_zXXYZ=zero;   real M_zXYYZ=zero; real M_zXYZZ=zero;
            real M_zXXYYZ=zero;   real M_zXXYZZ=zero;   real M_zXYYZZ=zero;   real M_zXXYYZZ=zero;
            real M_zXYeq=zero;    real M_zXZeq=zero;    real M_zYZeq=zero;    real M_zXYZeq=zero;
            real M_zXXYZeq=zero;  real M_zXYYZeq=zero;  real M_zXYZZeq=zero;
            real M_zXXYYZeq=zero; real M_zXXYZZeq=zero; real M_zXYYZZeq=zero; real M_zXXYYZZeq=zero;

            M_zXXYYZZ= -five*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] * 
                             mu[1][0][0] +     mu[0][0][2] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] * 
                             mu[1][0][0] + four*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][1] *     mu[1][0][0] * 
                             mu[1][0][0] - two*mu[0][1][0] *     mu[0][1][2] *     mu[1][0][0] *     mu[1][0][0] + 
                             mu[0][0][1] *     mu[0][0][1] *     mu[0][2][0] *     mu[1][0][0] *     mu[1][0][0] -
                         two*mu[0][0][1] *     mu[0][2][1] *     mu[1][0][0] *     mu[1][0][0] +     mu[0][2][2] *
                             mu[1][0][0] *     mu[1][0][0] + four*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *
                             mu[1][0][0] *     mu[1][0][1] - two*mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] *
                             mu[1][0][2] + four*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] *
                             mu[1][1][0] - eight*mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][1][1] +
                        four*mu[0][1][0] *     mu[1][0][0] *     mu[1][1][2] - two*mu[0][0][1] *     mu[0][0][1] *
                             mu[1][0][0] *     mu[1][2][0] + four*mu[0][0][1] *     mu[1][0][0] *     mu[1][2][1] -
                         two*mu[1][0][0] *     mu[1][2][2] +     mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *
                             mu[0][1][0] *     mu[2][0][0] - two*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *
                             mu[2][0][1] +     mu[0][1][0] *     mu[0][1][0] *     mu[2][0][2] - two*mu[0][0][1] *
                             mu[0][0][1] *     mu[0][1][0] *     mu[2][1][0] + four*mu[0][0][1] *     mu[0][1][0] *
                             mu[2][1][1] - two*mu[0][1][0] *     mu[2][1][2] +     mu[0][0][1] *     mu[0][0][1] *
                             mu[2][2][0] - two*mu[0][0][1] *     mu[2][2][1] +     mu[2][2][2];  //5
            M_zXXYYZ=   four*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][0] - 
                         two*mu[0][1][0] *     mu[0][1][1] *     mu[1][0][0] *     mu[1][0][0] -     mu[0][0][1] *
                             mu[0][2][0] *     mu[1][0][0] *     mu[1][0][0] +     mu[0][2][1] *     mu[1][0][0] *
                             mu[1][0][0] - two*mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][1] -
                        four*mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][1][0] + four*mu[0][1][0] *
                             mu[1][0][0] *     mu[1][1][1] + two*mu[0][0][1] *     mu[1][0][0] *     mu[1][2][0] - 
                         two*mu[1][0][0] *     mu[1][2][1] -     mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *
                             mu[2][0][0] +     mu[0][1][0] *     mu[0][1][0] *     mu[2][0][1] + two*mu[0][0][1] *
                             mu[0][1][0] *     mu[2][1][0] - two*mu[0][1][0] *     mu[2][1][1] -     mu[0][0][1] *
                             mu[2][2][0] +     mu[2][2][1]; //6
            M_zXXYY=  -three*mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][0] +     mu[0][2][0] *
                             mu[1][0][0] *     mu[1][0][0] + four*mu[0][1][0] *     mu[1][0][0] *     mu[1][1][0] -
                         two*mu[1][0][0] *     mu[1][2][0] +     mu[0][1][0] *     mu[0][1][0] *     mu[2][0][0] -
                         two*mu[0][1][0] *     mu[2][1][0] +     mu[2][2][0]; //7
            M_zXXYZZ=   four*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][0] -
                             mu[0][0][2] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][0] - two*mu[0][0][1] *
                             mu[0][1][1] *     mu[1][0][0] *     mu[1][0][0] +     mu[0][1][2] *     mu[1][0][0] *
                             mu[1][0][0] - four*mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][1] +
                         two*mu[0][1][0] *     mu[1][0][0] *     mu[1][0][2] - two*mu[0][0][1] *     mu[0][0][1] *
                             mu[1][0][0] *     mu[1][1][0] + four*mu[0][0][1] *     mu[1][0][0] *     mu[1][1][1] -
                         two*mu[1][0][0] *     mu[1][1][2] -     mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *
                             mu[2][0][0] + two*mu[0][0][1] *     mu[0][1][0] *     mu[2][0][1] -     mu[0][1][0] *
                             mu[2][0][2] +     mu[0][0][1] *     mu[0][0][1] *     mu[2][1][0] - two*mu[0][0][1] *
                             mu[2][1][1] +     mu[2][1][2]; //8
            M_zXXYZ=  -three*mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][0] +     mu[0][1][1] *
                             mu[1][0][0] *     mu[1][0][0] + two*mu[0][1][0] *     mu[1][0][0] *     mu[1][0][1] +
                         two*mu[0][0][1] *     mu[1][0][0] *     mu[1][1][0] - two*mu[1][0][0] *     mu[1][1][1] +
                             mu[0][0][1] *     mu[0][1][0] *     mu[2][0][0] -     mu[0][1][0] *     mu[2][0][1] - 
                             mu[0][0][1] *     mu[2][1][0] +     mu[2][1][1]; //9
            M_zXXY=      two*mu[0][1][0] *     mu[1][0][0] *     mu[1][0][0] - two*mu[1][0][0] *     mu[1][1][0] -
                             mu[0][1][0] *     mu[2][0][0] +     mu[2][1][0]; //10
            M_zXXZZ=  -three*mu[0][0][1] *     mu[0][0][1] *     mu[1][0][0] *     mu[1][0][0] +     mu[0][0][2] *
                             mu[1][0][0] *     mu[1][0][0] + four*mu[0][0][1] *     mu[1][0][0] *     mu[1][0][1] -
                         two*mu[1][0][0] *     mu[1][0][2] +     mu[0][0][1] *     mu[0][0][1] *     mu[2][0][0] -
                         two*mu[0][0][1] *     mu[2][0][1] +     mu[2][0][2]; //11
            M_zXXZ=      two*mu[0][0][1] *     mu[1][0][0] *     mu[1][0][0] - two*mu[1][0][0] *     mu[1][0][1] -
                             mu[0][0][1] *     mu[2][0][0] +     mu[2][0][1]; //12
            M_zXX=-          mu[1][0][0] *     mu[1][0][0] +     mu[2][0][0]; //13
            M_zXYYZZ=   four*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] -
                             mu[0][0][2] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] - four*mu[0][0][1] * 
                             mu[0][1][0] *     mu[0][1][1] *     mu[1][0][0] + two*mu[0][1][0] *     mu[0][1][2] *
                             mu[1][0][0] -     mu[0][0][1] *     mu[0][0][1] *     mu[0][2][0] *     mu[1][0][0] +
                         two*mu[0][0][1] *     mu[0][2][1] *     mu[1][0][0] -     mu[0][2][2] *     mu[1][0][0] -
                         two*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][1] +     mu[0][1][0] *
                             mu[0][1][0] *     mu[1][0][2] - two*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *
                             mu[1][1][0] + four*mu[0][0][1] *     mu[0][1][0] *     mu[1][1][1] - two*mu[0][1][0] *
                             mu[1][1][2] +     mu[0][0][1] *     mu[0][0][1] *     mu[1][2][0] - two*mu[0][0][1] *
                             mu[1][2][1] +     mu[1][2][2];    //14
            M_zXYYZ=  -three*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] + two*mu[0][1][0] *
                             mu[0][1][1] *     mu[1][0][0] +     mu[0][0][1] *     mu[0][2][0] *     mu[1][0][0] - 
                             mu[0][2][1] *     mu[1][0][0] +     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][1] + 
                         two*mu[0][0][1] *     mu[0][1][0] *     mu[1][1][0] - two*mu[0][1][0] *     mu[1][1][1] - 
                             mu[0][0][1] *     mu[1][2][0] +     mu[1][2][1]; //15
            M_zXYY=      two*mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] -     mu[0][2][0] *     mu[1][0][0] -
                         two*mu[0][1][0] *     mu[1][1][0] +     mu[1][2][0]; //16
            M_zXYZZ=  -three*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] +     mu[0][0][2] *
                             mu[0][1][0] *     mu[1][0][0] + two*mu[0][0][1] *     mu[0][1][1] *     mu[1][0][0] -     
                             mu[0][1][2] *     mu[1][0][0] + two*mu[0][0][1] *     mu[0][1][0] *     mu[1][0][1] -     
                             mu[0][1][0] *     mu[1][0][2] +     mu[0][0][1] *     mu[0][0][1] *     mu[1][1][0] - 
                         two*mu[0][0][1] *     mu[1][1][1] +     mu[1][1][2]; //17
            M_zXYZ=      two*mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] -     mu[0][1][1] *     mu[1][0][0] -
                             mu[0][1][0] *     mu[1][0][1] -     mu[0][0][1] *     mu[1][1][0] +     mu[1][ 1][ 1]; //18
            M_zXY =         -mu[0][1][0] *     mu[1][0][0] +     mu[1][1][0]; //19
            M_zXZZ=      two*mu[0][0][1] *     mu[0][0][1] *     mu[1][0][0] -     mu[0][0][2] *     mu[1][0][0] -
                         two*mu[0][0][1] *     mu[1][0][1] +     mu[1][0][2]; //20
            M_zXZ =         -mu[0][0][1] *     mu[1][0][0] +     mu[1][0][1]; //21
            M_zYYZZ=  -three*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] +     mu[0][0][2] *
                             mu[0][1][0] *     mu[0][1][0] + four*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][1] -
                         two*mu[0][1][0] *     mu[0][1][2] +     mu[0][0][1] *     mu[0][0][1] *     mu[0][2][0] -
                         two*mu[0][0][1] *     mu[0][2][1] +     mu[0][2][2]; //22
            M_zYYZ=      two*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] - two*mu[0][1][0] *     mu[0][1][1] -
                             mu[0][0][1] *     mu[0][2][0] +     mu[0][2][1]; //23
            M_zYY  =        -mu[0][1][0] *     mu[0][1][0] +     mu[0][2][0]; //24
            M_zYZZ =     two*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] -     mu[0][0][2] *     mu[0][1][0] -  
                         two*mu[0][0][1] *     mu[0][1][1] +     mu[0][1][2]; //25
            M_zYZ  =        -mu[0][0][1] *     mu[0][1][0] +     mu[0][1][1]; //26
            M_zZZ  =        -mu[0][0][1] *     mu[0][0][1] +     mu[0][0][2]; //27

            M_zXXYYZZeq=c1o27;
            M_zXXYYZeq=zero;
            M_zXXYZZeq=zero;
            M_zXXYZeq=zero;
            M_zXYYZZeq=zero;
            M_zXYYZeq=zero;
            M_zXYZZeq=zero;

            real MXXpMYYpMZZ,MXXmMYY, MXXmMZZ, MXXYpMYZZ,MXXZpMYYZ,MXYYpMXZZ, MXXYmMYZZ,MXXZmMYYZ,MXYYmMXZZ, MXXYYppp,MXXYYpm2p, MXXYYppm2;
            real MXXpMYYpMZZeq,MXXmMYYeq, MXXmMZZeq, MXXYpMYZZeq,MXXZpMYYZeq,MXYYpMXZZeq, MXXYmMYZZeq,MXXZmMYYZeq,MXYYmMXZZeq, MXXYYpppeq,MXXYYpm2peq, MXXYYppm2eq;

            //coll faktoren:
            real w1,w2,w3,w4,w5,w6,w7,w8,w9,w10;
            w2=-one;//-omega;
            w7=-one;
            w9=-one;
            w10=-one;
            w1=-omega;
            w3=-one;
            w4=-one;
            w5=-one;
            w6=-one;
            w8=-one;
            ////////lin kombi bilden:
            MXXpMYYpMZZ = M_zXX + M_zYY + M_zZZ;
            MXXmMYY = M_zXX - M_zYY;
            MXXmMZZ = M_zXX - M_zZZ;

            MXXYpMYZZ  =  M_zXXY+M_zYZZ;
            MXXYmMYZZ  =  M_zXXY-M_zYZZ;
            MXXZpMYYZ  =  M_zXXZ+M_zYYZ;
            MXXZmMYYZ  =  M_zXXZ-M_zYYZ;
            MXYYpMXZZ  =  M_zXYY+M_zXZZ;
            MXYYmMXZZ  =  M_zXYY-M_zXZZ;

            MXXYYppp  = M_zXXYY + M_zXXZZ + M_zYYZZ;
            MXXYYpm2p = M_zXXYY - 2.0*M_zXXZZ + M_zYYZZ;
            MXXYYppm2 = M_zXXYY + M_zXXZZ - 2.0*M_zYYZZ;
            bool casc=false;
            bool factorCasc=true;
            if(casc)
            {
              MXXpMYYpMZZeq = one;
              MXXmMYYeq     = zero;
              MXXmMZZeq     = zero;

              MXXYpMYZZeq  =  zero;
              MXXYmMYZZeq  =  zero;
              MXXZpMYYZeq  =  zero;
              MXXZmMYYZeq  =  zero;
              MXYYpMXZZeq  =  zero;
              MXYYmMXZZeq  =  zero;

              MXXYYpppeq  = c1o3;
              MXXYYpm2peq = zero;
              MXXYYppm2eq = zero;

              //relaxation:
              MXXpMYYpMZZ -= w2*(MXXpMYYpMZZeq-MXXpMYYpMZZ);
              MXXmMYY -=     w1*(MXXmMYYeq    -MXXmMYY);
              MXXmMZZ -=     w1*(MXXmMZZeq    -MXXmMZZ);

              MXXYpMYZZ  -=  w3*(MXXYpMYZZeq-MXXYpMYZZ);
              MXXYmMYZZ  -=  w4*(MXXYmMYZZeq-MXXYmMYZZ);
              MXXZpMYYZ  -=  w3*(MXXZpMYYZeq-MXXZpMYYZ);
              MXXZmMYYZ  -=  w4*(MXXZmMYYZeq-MXXZmMYYZ);
              MXYYpMXZZ  -=  w3*(MXYYpMXZZeq-MXYYpMXZZ);
              MXYYmMXZZ  -=  w4*(MXYYmMXZZeq-MXYYmMXZZ);

              MXXYYppp  -= w7*(MXXYYpppeq -MXXYYppp );
              MXXYYpm2p -= w6*(MXXYYpm2peq-MXXYYpm2p);
              MXXYYppm2 -= w6*(MXXYYppm2eq-MXXYYppm2);

              M_zXXYYZZ=M_zXXYYZZ-w10*(M_zXXYYZZeq-M_zXXYYZZ);
              M_zXZ=    M_zXZ-    w1*(M_zXZeq   - M_zXZ    );
              M_zYZ=    M_zYZ-    w1*(M_zYZeq   - M_zYZ    );
              M_zXY=    M_zXY-    w1*(M_zXYeq   - M_zXY    );

              M_zXYZ=   M_zXYZ-   w5*(M_zXYZeq  - M_zXYZ   );

              M_zXYYZ=  M_zXYYZ-  w8*(M_zXYYZeq - M_zXYYZ  );
              M_zXYZZ=  M_zXYZZ-  w8*(M_zXYZZeq - M_zXYZZ  );
              M_zXXYZ=  M_zXXYZ-  w8*(M_zXXYZeq - M_zXXYZ  );

              M_zXXYYZ= M_zXXYYZ- w9*(M_zXXYYZeq- M_zXXYYZ );
              M_zXXYZZ= M_zXXYZZ- w9*(M_zXXYZZeq- M_zXXYZZ );
              M_zXYYZZ= M_zXYYZZ- w9*(M_zXYYZZeq- M_zXYYZZ );

              ////////von Lin Kombis zurueck:
              M_zXX = c1o3 * MXXmMYY + c1o3 * MXXmMZZ +  c1o3 * MXXpMYYpMZZ;
              M_zYY =-c2o3 * MXXmMYY + c1o3 * MXXmMZZ +  c1o3 * MXXpMYYpMZZ;
              M_zZZ = c1o3 * MXXmMYY - c2o3 * MXXmMZZ +  c1o3 * MXXpMYYpMZZ;

              M_zXXY = (MXXYmMYZZ + MXXYpMYZZ)*c1o2;
              M_zYZZ = c1o2*(-MXXYmMYZZ + MXXYpMYZZ);
              M_zXYY =(MXYYmMXZZ + MXYYpMXZZ)*c1o2;
              M_zXZZ = c1o2*(-MXYYmMXZZ + MXYYpMXZZ);
              M_zXXZ = (MXXZmMYYZ + MXXZpMYYZ)*c1o2;
              M_zYYZ = c1o2*(-MXXZmMYYZ + MXXZpMYYZ);

			  M_zXXYY =  c1o3* MXXYYpm2p + c1o3* MXXYYppm2 + c1o3*MXXYYppp;
			  M_zXXZZ = -c1o3* MXXYYpm2p + c1o3* MXXYYppp;
			  M_zYYZZ = -c1o3* MXXYYppm2 + c1o3* MXXYYppp;
            }
            else if (factorCasc)
            {
              MXXpMYYpMZZeq = one;
              MXXmMYYeq     = zero;
              MXXmMZZeq     = zero;

              MXXYpMYZZeq  =  zero;
              MXXYmMYZZeq  =  zero;
              MXXZpMYYZeq  =  zero;
              MXXZmMYYZeq  =  zero;
              MXYYpMXZZeq  =  zero;
              MXYYmMXZZeq  =  zero;

              MXXYYpppeq  = c1o3;
              MXXYYpm2peq = zero;
              MXXYYppm2eq = zero;

              //relaxation:
              MXXpMYYpMZZ -= w2*(MXXpMYYpMZZeq-MXXpMYYpMZZ);
              MXXmMYY -=     w1*(MXXmMYYeq    -MXXmMYY);
              MXXmMZZ -=     w1*(MXXmMZZeq    -MXXmMZZ);

              MXXYpMYZZ  -=  w3*(MXXYpMYZZeq-MXXYpMYZZ);
              MXXYmMYZZ  -=  w4*(MXXYmMYZZeq-MXXYmMYZZ);
              MXXZpMYYZ  -=  w3*(MXXZpMYYZeq-MXXZpMYYZ);
              MXXZmMYYZ  -=  w4*(MXXZmMYYZeq-MXXZmMYYZ);
              MXYYpMXZZ  -=  w3*(MXYYpMXZZeq-MXYYpMXZZ);
              MXYYmMXZZ  -=  w4*(MXYYmMXZZeq-MXYYmMXZZ);

              M_zXZ=    M_zXZ-    w1*(M_zXZeq   - M_zXZ    );
              M_zYZ=    M_zYZ-    w1*(M_zYZeq   - M_zYZ    );
              M_zXY=    M_zXY-    w1*(M_zXYeq   - M_zXY    );
              M_zXYZ=   M_zXYZ-   w5*(M_zXYZeq  - M_zXYZ   );

              ////////von Lin Kombis zurueck:
              M_zXX = c1o3 * MXXmMYY + c1o3 * MXXmMZZ +  c1o3 * MXXpMYYpMZZ;
              M_zYY =-c2o3 * MXXmMYY + c1o3 * MXXmMZZ +  c1o3 * MXXpMYYpMZZ;
              M_zZZ = c1o3 * MXXmMYY - c2o3 * MXXmMZZ +  c1o3 * MXXpMYYpMZZ;

              M_zXXY = (MXXYmMYZZ + MXXYpMYZZ)*c1o2;
              M_zYZZ = c1o2*(-MXXYmMYZZ + MXXYpMYZZ);
              M_zXYY =(MXYYmMXZZ + MXYYpMXZZ)*c1o2;
              M_zXZZ = c1o2*(-MXYYmMXZZ + MXYYpMXZZ);
              M_zXXZ = (MXXZmMYYZ + MXXZpMYYZ)*c1o2;
              M_zYYZ = c1o2*(-MXXZmMYYZ + MXXZpMYYZ);

              //faktorisierte atraktoren:
              MXXYYpppeq=M_zXX*M_zYY+M_zXX*M_zZZ+M_zZZ*M_zYY;    
              MXXYYpm2peq=M_zXX*M_zYY-two*M_zXX*M_zZZ+M_zZZ*M_zYY; 
              MXXYYppm2eq=M_zXX*M_zYY+M_zXX*M_zZZ-two*M_zZZ*M_zYY; 
              M_zXYYZeq=M_zYY*M_zXZ;
              M_zXYZZeq=M_zZZ*M_zXY;
              M_zXXYZeq=M_zXX*M_zYZ;

              M_zXXYYZeq=zero;
              M_zXXYZZeq=zero;
              M_zXYYZZeq=zero;

              M_zXXYYZZeq=M_zXX*M_zYY*M_zZZ;

              MXXYYppp  -= w7*(MXXYYpppeq -MXXYYppp );
              MXXYYpm2p -= w6*(MXXYYpm2peq-MXXYYpm2p);
              MXXYYppm2 -= w6*(MXXYYppm2eq-MXXYYppm2);
              M_zXXYYZZ=M_zXXYYZZ-w10*(M_zXXYYZZeq-M_zXXYYZZ);
              M_zXYYZ=  M_zXYYZ-  w8*(M_zXYYZeq - M_zXYYZ  );
              M_zXYZZ=  M_zXYZZ-  w8*(M_zXYZZeq - M_zXYZZ  );
              M_zXXYZ=  M_zXXYZ-  w8*(M_zXXYZeq - M_zXXYZ  );

              M_zXXYYZ= M_zXXYYZ- w9*(M_zXXYYZeq- M_zXXYYZ );
              M_zXXYZZ= M_zXXYZZ- w9*(M_zXXYZZeq- M_zXXYZZ );
              M_zXYYZZ= M_zXYYZZ- w9*(M_zXYYZZeq- M_zXYYZZ );

              M_zXXYY =  c1o3* MXXYYpm2p + c1o3* MXXYYppm2 + c1o3*MXXYYppp;
              M_zXXZZ = -c1o3* MXXYYpm2p + c1o3* MXXYYppp;
              M_zYYZZ = -c1o3* MXXYYppm2 + c1o3* MXXYYppp;
            }
            //real forcingX1=1.11986e-8;//1.41342e-8;//5.48288e-9;//<-96//1.03093e-8;//1.59178e-8;//1.6256e-8; 

            //mu[1][0][0]+=forcingX1;

            mu[2][ 2][ 2]= vx*vx*vy*vy*vz*vz + vx*vx*vy*vy*M_zZZ + four*vx*vx*vy*vz*M_zYZ + two* vx*vx*vy*M_zYZZ +
                           vx*vx*vz*vz*M_zYY + two* vx*vx*vz*M_zYYZ +vx*vx*M_zYYZZ + four*vx*vy*vy*vz*M_zXZ   +
                           two* vx*vy*vy*M_zXZZ + four* vx*vy*vz*vz*M_zXY+   eight* vx*vy*vz*M_zXYZ + four* vx*vy*M_zXYZZ +
                           two* vx*vz*vz*M_zXYY + four* vx*vz*M_zXYYZ + two*vx*M_zXYYZZ +   vy*vy*vz*vz*M_zXX + two* vy*vy*vz*M_zXXZ  +
                           vy*vy*M_zXXZZ     + two* vy*vz*vz*M_zXXY  + four*vy*vz*M_zXXYZ +  two* vy*M_zXXYZZ + vz*vz*M_zXXYY +two* vz*M_zXXYYZ + M_zXXYYZZ;//-(c1o27);
            mu[2][ 2][ 1]= vx*vx*vy*vy*vz + two*vx*vx*vy*M_zYZ + vx*vx*vz*M_zYY +  vx*vx*M_zYYZ + two*vx*vy*vy*M_zXZ +four*vx*vy*vz*M_zXY +four*vx*vy*M_zXYZ+
                           two*vx*vz*M_zXYY + two*vx*M_zXYYZ + vy*vy*vz*M_zXX +  vy*vy*M_zXXZ + two* vy*vz*M_zXXY + two*vy*M_zXXYZ +  vz*M_zXXYY + M_zXXYYZ;
            mu[2][ 2][ 0]= vx*vx*vy*vy + vx*vx*M_zYY + four*vx*vy*M_zXY +   two*vx*M_zXYY + vy*vy*M_zXX + two*vy*M_zXXY +  M_zXXYY ;//-(c1o9);
            mu[2][ 1][ 2]= vx*vx*vy*vz*vz + vx*vx*vy*M_zZZ + two*vx*vx*vz*M_zYZ +   vx*vx*M_zYZZ + four*vx*vy*vz*M_zXZ+ //vy->vz
                           two*vx*vy*M_zXZZ + two*vx*vz*vz*M_zXY +   four*vx*vz*M_zXYZ + two*vx*M_zXYZZ + vy*vz*vz*M_zXX +
                           two*vy*vz*M_zXXZ + vy*M_zXXZZ + vz*vz*M_zXXY +  two*vz*M_zXXYZ + M_zXXYZZ;
            mu[2][ 1][ 1]= vx*vx*vy*vz + vx*vx*M_zYZ + two*vx*vy*M_zXZ + two*vx*vz*M_zXY + two*vx*M_zXYZ + vy*vz*M_zXX+  vy*M_zXXZ + vz*M_zXXY + M_zXXYZ;
            mu[2][ 1][ 0]= vx*vx*vy + two*vx*M_zXY + vy*M_zXX + M_zXXY;
            mu[2][ 0][ 2]= vx*vx*vz*vz + vx*vx*M_zZZ + four*vx*vz*M_zXZ+ two*vx*M_zXZZ + vz*vz*M_zXX+ two*vz*M_zXXZ + M_zXXZZ;//-(c1o9);
            mu[2][ 0][ 1]= vx*vx*vz + two*vx*M_zXZ + vz*M_zXX + M_zXXZ;
            mu[2][ 0][ 0]= vx*vx + M_zXX ;//-(c1o3);
            mu[1][ 2][ 2]= vx*vy*vy*vz*vz + vx*vy*vy*M_zZZ + four*vx*vy*vz*M_zYZ +   two*vx*vy*M_zYZZ + vx*vz*vz*M_zYY + two*vx*vz*M_zYYZ +vx*M_zYYZZ + two*vy*vy*vz*M_zXZ + vy*vy*M_zXZZ +
                           two*vy*vz*vz*M_zXY + four*vy*vz*M_zXYZ + two*vy*M_zXYZZ +vz*vz*M_zXYY + two*vz*M_zXYYZ + M_zXYYZZ;
            mu[1][ 2][ 1]= vx*vy*vy*vz + two*vx*vy*M_zYZ + vx*vz*M_zYY +   vx*M_zYYZ + vy*vy*M_zXZ + two*vy*vz*M_zXY +   two*vy*M_zXYZ + vz*M_zXYY + M_zXYYZ;
            mu[1][ 2][ 0]= vx*vy*vy + vx*M_zYY + two*vy*M_zXY + M_zXYY;
            mu[1][ 1][ 2]= vx*vy*vz*vz+vx*vy*M_zZZ+two*vx*vz*M_zYZ+vx*M_zYZZ+two*vy*vz*M_zXZ+vy*M_zXZZ+vz*vz*M_zXY+two*vz*M_zXYZ+M_zXYZZ;
            mu[1][ 1][ 1]= vx*vy*vz + vx*M_zYZ + vy*M_zXZ + vz*M_zXY+M_zXYZ;
            mu[1][ 1][ 0]= vx*vy+ M_zXY;
            mu[1][ 0][ 2]= vx*vz*vz + vx*M_zZZ + two*vz*M_zXZ + M_zXZZ;
            mu[1][ 0][ 1]= vx*vz + M_zXZ;
            //mu[1][ 0][ 0]= vx;
            mu[0][ 2][ 2]= vy*vy*vz*vz + vy*vy*M_zZZ + four*vy*vz*M_zYZ +two*vy*M_zYZZ + vz*vz*M_zYY + two*vz*M_zYYZ +  M_zYYZZ;//-(c1o9);
            mu[0][ 2][ 1]= vy*vy*vz + two*vy*M_zYZ + vz*M_zYY + M_zYYZ;
            mu[0][ 2][ 0]= vy*vy + M_zYY ;//-(c1o3);
            mu[0][ 1][ 2]= vy*vz*vz+ vy*M_zZZ + two*vz*M_zYZ + M_zYZZ;
            mu[0][ 1][ 1]= vy*vz + M_zYZ;
            //mu[0][ 1][ 0]= vy ;                     //richtig?
            mu[0][ 0][ 2]= vz*vz + M_zZZ;// -(c1o3);
            //mu[0][ 0][ 1]= vz;
            //mu[0][ 0][ 0]=rho0 ;

            f_E =c1o2* (mu[2][0][0] -  mu[2][2][0] + mu[2][2][2] - mu[2][0][2] -   mu[1][2][0] + mu[1][2][2] - mu[1][0][2] +mu[1][0][0]) *(rho);
            f_W =c1o2* (mu[2][0][0] - mu[2][2][0] + mu[2][2][2] - mu[2][0][2] +    mu[1][2][0] - mu[1][2][2] + mu[1][0][2] -mu[1][0][0])*  (rho);
            f_N =c1o2* (-mu[2][1][0] - mu[2][2][0] + mu[2][2][2] + mu[2][1][2] +   mu[0][2][0] - mu[0][2][2] - mu[0][1][2] +mu[0][1][0])* (rho);
            f_S =c1o2* (mu[2][1][0] -  mu[2][2][0] + mu[2][2][2] - mu[2][1][2] +   mu[0][2][0] - mu[0][2][2] + mu[0][1][2] -mu[0][1][0])* (rho);
            f_T =c1o2* (mu[2][2][1] +  mu[2][2][2] - mu[2][0][1] - mu[2][0][2] -   mu[0][2][1] - mu[0][2][2] + mu[0][0][2] +mu[0][0][1])* (rho);
            f_B =c1o2* (-mu[2][2][1] + mu[2][2][2] + mu[2][0][1]  - mu[2][0][2] +  mu[0][2][1] - mu[0][2][2] + mu[0][0][2]-mu[0][0][1])* (rho);
            f_NE =c1o4*( mu[2][1][0]  + mu[2][2][0]- mu[2][2][2] - mu[2][1][2] +  mu[1][1][0]+ mu[1][2][0]- mu[1][2][2] -mu[1][1][2])*  (rho);
            f_SW =c1o4*(-mu[2][1][0] + mu[2][2 ][0]- mu[2][2][2] + mu[2][1][2] +  mu[1][1][0]- mu[1][2][0]+ mu[1][2][2] -mu[1][1][2])*  (rho);
            f_SE =c1o4*(-mu[2][1][0] + mu[2][2 ][0]- mu[2][2][2] + mu[2][1][2] -  mu[1][1][0]+ mu[1][2][0]- mu[1][2][2] +mu[1][1][2])*  (rho);
            f_NW =c1o4*( mu[2][1][0]  + mu[2][2][0]- mu[2][2][2] - mu[2][1][2] -  mu[1][1][0]- mu[1][2][0]+ mu[1][2][2] + mu[1][1][2])* (rho);
            f_TE=c1o4*(-mu[2][2][1] - mu[2][2][2] + mu[2][0][1] + mu[2][0][2] -   mu[1][2][1] - mu[1][2][2] + mu[1][0][1] + mu[1][0][2])*(rho);
            f_BW=c1o4*( mu[2][2][1]  -mu[2][2][2] - mu[2][0][1] + mu[2][0][2] -   mu[1][2][1] + mu[1][2][2] + mu[1][0][1] - mu[1][0][2])*(rho);
            f_BE=c1o4*(mu[2][2][1] - mu[2][2][2] - mu[2][0][1] + mu[2][0][2] +    mu[1][2][1] - mu[1][2][2] - mu[1][0][1] +mu[1][0][2])*  (rho);
            f_TW=c1o4*(-mu[2][2][1] - mu[2][2][2] + mu[2][0][1] + mu[2][0][2] +   mu[1][2][1] + mu[1][2][2] - mu[1][0][1] -mu[1][0][2])* (rho);
            f_TN=c1o4*(-mu[2][2][1] - mu[2][2][2] - mu[2][1][1] - mu[2][1][2] +   mu[0][2][1] + mu[0][2][2] + mu[0][1][1]+mu[0][1][2])*  (rho);
            f_BS=c1o4*( mu[2][2][1] - mu[2][2][2] - mu[2][1][1] + mu[2][1][2] -   mu[0][2][1] + mu[0][2][2] + mu[0][1][1] - mu[0][1][2])*(rho);
            f_BN=c1o4*( mu[2][2][1] - mu[2][2][2] + mu[2][1][1] - mu[2][1][2] -   mu[0][2][1] + mu[0][2][2] - mu[0][1][1] + mu[0][1][2])*(rho);
            f_TS=c1o4*(-mu[2][2][1] - mu[2][2][2] + mu[2][1][1] + mu[2][1][2] +   mu[0][2][1] + mu[0][2][2] - mu[0][1][1] -mu[0][1][2])* (rho);
            f_ZERO=     (-mu[2][0][0] + mu[2][2][0] - mu[2][2][2] + mu[2][0][2] - mu[0][2][0] + mu[0][2][2] - mu[0][0][2] + mu[0][0][0])*(rho);
            f_TNE=c1o8*( mu[2][2][1] + mu[2][2][2] + mu[2][1][1] + mu[2][1][2] + mu[1][2][1] + mu[1][2][2] + mu[1][1][1] + mu[1][1][2])*  (rho);
            f_BNE=c1o8*(-mu[2][2][1] + mu[2][2][2] -mu[2][1][1] + mu[2][1][2]  - mu[1][2][1] + mu[1][2][2] -mu[1][1][1] + mu[1][1][2])*     (rho);
            f_TSE=c1o8*( mu[2][2][1] + mu[2][2][2] - mu[2][1][1] - mu[2][1][2] + mu[1][2][1] + mu[1][2][2] - mu[1][1][1] - mu[1][1][2])*  (rho);
            f_BSE=c1o8*(-mu[2][2][1] + mu[2][2][2] +mu[2][1][1] - mu[2][1][2]  - mu[1][2][1] + mu[1][2][2] +mu[1][1][1] - mu[1][1][2])*     (rho);
            f_TNW=c1o8*( mu[2][2][1] + mu[2][2][2] + mu[2][1][1] + mu[2][1][2] - mu[1][2][1] - mu[1][2][2] - mu[1][1][1] - mu[1][1][2])*  (rho);
            f_BNW=c1o8*(-mu[2][2][1] + mu[2][2][2] -mu[2][1][1] + mu[2][1][2] +  mu[1][2][1] - mu[1][2][2] +mu[1][1][1] - mu[1][1][2])*     (rho);
            f_TSW=c1o8*( mu[2][2][1] + mu[2][2][2] - mu[2][1][1] - mu[2][1][2] - mu[1][2][1] - mu[1][2][2] + mu[1][1][1] + mu[1][1][2])*  (rho);
            f_BSW=c1o8*(-mu[2][2][1] + mu[2][2][2]+mu[2][1][1] - mu[2][1][2]+    mu[1][2][1] - mu[1][2][2]-mu[1][1][1] + mu[1][1][2])*(rho);
          }
        }

       (D.f[ dirE  ])[ke ] = f_W-c2over27;
       (D.f[ dirW  ])[kw ] = f_E-c2over27;

       (D.f[ dirN  ])[kn ] = f_S-c2over27;
       (D.f[ dirS  ])[ks ] = f_N-c2over27;
       (D.f[ dirT  ])[kt ] = f_B-c2over27;
       (D.f[ dirB  ])[kb ] = f_T-c2over27;

       (D.f[ dirNE ])[kne] = f_SW-c1over54;
       (D.f[ dirSW ])[ksw] = f_NE-c1over54;
       (D.f[ dirSE ])[kse] = f_NW-c1over54;
       (D.f[ dirNW ])[knw] = f_SE-c1over54;
       (D.f[ dirTE ])[kte] = f_BW-c1over54;
       (D.f[ dirBW ])[kbw] = f_TE-c1over54;
       (D.f[ dirBE ])[kbe] = f_TW-c1over54;
       (D.f[ dirTW ])[ktw] = f_BE-c1over54;

       (D.f[ dirTN ])[ktn] = f_BS-c1over54;
       (D.f[ dirBS ])[kbs] = f_TN-c1over54;
       (D.f[ dirBN ])[kbn] = f_TS-c1over54;
       (D.f[ dirTS ])[kts] = f_BN-c1over54;

       (D.f[ dirZERO])[k] = f_ZERO-c8over27;

       (D.f[ dirTNE ])[ktne] = f_BSW-c1over216;
       (D.f[ dirTSE ])[ktse] = f_BNW-c1over216;
       (D.f[ dirBNE ])[kbne] = f_TSW-c1over216;
       (D.f[ dirBSE ])[kbse] = f_TNW-c1over216;
       (D.f[ dirTNW ])[ktnw] = f_BSE-c1over216;
       (D.f[ dirTSW ])[ktsw] = f_BNE-c1over216;
       (D.f[ dirBNW ])[kbnw] = f_TSE-c1over216;
       (D.f[ dirBSW ])[kbsw] = f_TNE-c1over216;
      }
     __syncthreads();
     }
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_Casc27(real omega,
                                            unsigned int* bcMatD,
                                            unsigned int* neighborX,
                                            unsigned int* neighborY,
                                            unsigned int* neighborZ,
                                            real* DDStart,
                                            int size_Mat,
                                            bool EvenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   int  k;                   // Zugriff auf arrays im device
   //
   int tx = threadIdx.x;     // Thread index = lokaler i index
   int by = blockIdx.x;      // Block index x
   int bz = blockIdx.y;      // Block index y
   int  x = tx + STARTOFFX;  // Globaler x-Index
   int  y = by + STARTOFFY;  // Globaler y-Index
   int  z = bz + STARTOFFZ;  // Globaler z-Index
   //vom 10.06.:
   const unsigned sizeX = blockDim.x;
   const unsigned sizeY = gridDim.x;

   const unsigned nx = sizeX + 2 * STARTOFFX;
   const unsigned ny = sizeY + 2 * STARTOFFY;

   k = nx*(ny*z + y) + x;
   ////////////////////////////////////////////////////////////////////////////////
   unsigned int BC;
   BC        =   bcMatD[k];

   if( (BC != GEO_SOLID) && (BC != GEO_VOID))
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
      //unsigned int nxny = nx*ny;
      //unsigned int kzero= k;
      //unsigned int ke   = k;
      //unsigned int kw   = k + 1;
      //unsigned int kn   = k;
      //unsigned int ks   = k + nx;
      //unsigned int kt   = k;
      //unsigned int kb   = k + nxny;
      //unsigned int ksw  = k + nx + 1;
      //unsigned int kne  = k;
      //unsigned int kse  = k + nx;
      //unsigned int knw  = k + 1;
      //unsigned int kbw  = k + nxny + 1;
      //unsigned int kte  = k;
      //unsigned int kbe  = k + nxny;
      //unsigned int ktw  = k + 1;
      //unsigned int kbs  = k + nxny + nx;
      //unsigned int ktn  = k;
      //unsigned int kbn  = k + nxny;
      //unsigned int kts  = k + nx;
      //unsigned int ktse = k + nx;
      //unsigned int kbnw = k + nxny + 1;
      //unsigned int ktnw = k + 1;
      //unsigned int kbse = k + nxny + nx;
      //unsigned int ktsw = k + nx + 1;
      //unsigned int kbne = k + nxny;
      //unsigned int ktne = k;
      //unsigned int kbsw = k + nxny + nx + 1;
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f_E,f_W,f_N,f_S,f_T,f_B,f_NE,f_SW,f_SE,f_NW,f_TE,f_BW,f_BE,f_TW,f_TN,f_BS,f_BN,f_TS,f_ZERO, f_TNE,f_TNW,f_TSE,f_TSW, f_BNE,f_BNW,f_BSE,f_BSW;
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      f_E    =  (D.f[dirE   ])[ke]+c2over27;
      f_W    =  (D.f[dirW   ])[kw]+c2over27;
      f_N    =  (D.f[dirN   ])[kn]+c2over27;
      f_S    =  (D.f[dirS   ])[ks]+c2over27;
      f_T    =  (D.f[dirT   ])[kt]+c2over27;
      f_B    =  (D.f[dirB   ])[kb]+c2over27;
      f_NE   =  (D.f[dirNE  ])[kne]+c1over54;
      f_SW   =  (D.f[dirSW  ])[ksw]+c1over54;
      f_SE   =  (D.f[dirSE  ])[kse]+c1over54;
      f_NW   =  (D.f[dirNW  ])[knw]+c1over54;
      f_TE   =  (D.f[dirTE  ])[kte]+c1over54;
      f_BW   =  (D.f[dirBW  ])[kbw]+c1over54;
      f_BE   =  (D.f[dirBE  ])[kbe]+c1over54;
      f_TW   =  (D.f[dirTW  ])[ktw]+c1over54;
      f_TN   =  (D.f[dirTN  ])[ktn]+c1over54;
      f_BS   =  (D.f[dirBS  ])[kbs]+c1over54;
      f_BN   =  (D.f[dirBN  ])[kbn]+c1over54;
      f_TS   =  (D.f[dirTS  ])[kts]+c1over54;
      f_ZERO =  (D.f[dirZERO])[kzero]+c8over27;
      f_TNE   = (D.f[dirTNE ])[ktne]+c1over216;
      f_TSW   = (D.f[dirTSW ])[ktsw]+c1over216;
      f_TSE   = (D.f[dirTSE ])[ktse]+c1over216;
      f_TNW   = (D.f[dirTNW ])[ktnw]+c1over216;
      f_BNE   = (D.f[dirBNE ])[kbne]+c1over216;
      f_BSW   = (D.f[dirBSW ])[kbsw]+c1over216;
      f_BSE   = (D.f[dirBSE ])[kbse]+c1over216;
      f_BNW   = (D.f[dirBNW ])[kbnw]+c1over216;
      ////////////////////////////////////////////////////////////////////////////////

      if( BC == GEO_FLUID || BC == GEO_VELO)
      {
         {
            //    /////////////////////////version 13.04.
            //   real drho     = f_ZERO+f_E +f_W +f_N +f_S +f_T +f_B +f_NE+f_SW+f_SE+f_NW+f_TE+f_BW+f_BE+f_TW+f_TN+f_BS+f_BN+f_TS+ f_TNE+f_TSW+f_TSE+f_TNW+ f_BNE+f_BSW+f_BSE+f_BNW;
            //   real vx1     =  f_E -f_W +f_NE-f_SW+f_SE-f_NW+f_TE-f_BW+f_BE-f_TW+ f_TNE-f_TSW+f_TSE-f_TNW+ f_BNE-f_BSW+f_BSE-f_BNW;
            //   real vx2     =  f_N -f_S +f_NE-f_SW-f_SE+f_NW+f_TN-f_BS+f_BN-f_TS+ f_TNE-f_TSW-f_TSE+f_TNW+ f_BNE-f_BSW-f_BSE+f_BNW;
            //   real vx3     =  f_T -f_B +f_TE-f_BW-f_BE+f_TW+f_TN-f_BS-f_BN+f_TS+ f_TNE+f_TSW+f_TSE+f_TNW- f_BNE-f_BSW-f_BSE-f_BNW;
            ////   vx1=vx1/(drho); vx2=vx2/(drho); vx3=vx3/(drho);
            //   //////////////////////////////////////////////////////////////////////////
            //   //BGK
            //   //////////////////////////////////////////////////////////////////////////


            //   real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

            //   f_ZERO = f_ZERO *(one+(-omega))-(-omega)*   c8over27*  (drho-cu_sq);
            //   f_E    = f_E    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
            //   f_W    = f_W    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
            //   f_N    = f_N    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
            //   f_S    = f_S    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
            //   f_T    = f_T    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
            //   f_B    = f_B    *(one+(-omega))-(-omega)*   c2over27*  (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
            //   f_NE   = f_NE   *(one+(-omega))-(-omega)*   c1over54*  (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
            //   f_SW   = f_SW   *(one+(-omega))-(-omega)*   c1over54*  (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
            //   f_SE   = f_SE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
            //   f_NW   = f_NW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
            //   f_TE   = f_TE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
            //   f_BW   = f_BW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
            //   f_BE   = f_BE   *(one+(-omega))-(-omega)*    c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
            //   f_TW   = f_TW   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
            //   f_TN   = f_TN   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
            //   f_BS   = f_BS   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
            //   f_BN   = f_BN   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
            //   f_TS   = f_TS   *(one+(-omega))-(-omega)*    c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
            //   f_TNE  = f_TNE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
            //   f_BSW  = f_BSW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
            //   f_BNE  = f_BNE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
            //   f_TSW  = f_TSW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
            //   f_TSE  = f_TSE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
            //   f_BNW  = f_BNW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
            //   f_BSE  = f_BSE  *(one+(-omega))-(-omega)*    c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
            //   f_TNW  = f_TNW  *(one+(-omega))-(-omega)*    c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);

            ///////////////////////version 13.04.
            real rho   =  f_ZERO+f_E +f_W +f_N +f_S +f_T +f_B +f_NE+f_SW+f_SE+f_NW+f_TE+f_BW+f_BE+f_TW+f_TN+f_BS+f_BN+f_TS+ f_TNE+f_TSW+f_TSE+f_TNW+ f_BNE+f_BSW+f_BSE+f_BNW;
            //real rho    =  rho0 + one;
            real vx     =  f_E -f_W +f_NE-f_SW+f_SE-f_NW+f_TE-f_BW+f_BE-f_TW+ f_TNE-f_TSW+f_TSE-f_TNW+ f_BNE-f_BSW+f_BSE-f_BNW;
            real vy     =  f_N -f_S +f_NE-f_SW-f_SE+f_NW+f_TN-f_BS+f_BN-f_TS+ f_TNE-f_TSW-f_TSE+f_TNW+ f_BNE-f_BSW-f_BSE+f_BNW;
            real vz     =  f_T -f_B +f_TE-f_BW-f_BE+f_TW+f_TN-f_BS-f_BN+f_TS+ f_TNE+f_TSW+f_TSE+f_TNW- f_BNE-f_BSW-f_BSE-f_BNW;
            vx=vx/(rho); vy=vy/(rho); vz=vz/(rho);

            real mu[3][3][3];

            mu[2][0][0]=(f_E+f_W+f_NE+f_SW+f_SE+f_NW+f_TE+f_BW+f_BE+f_TW+f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o3);
            mu[0][2][0]=(f_N+f_S+f_NE+f_SW+f_SE+f_NW+f_TN+f_BS+f_BN+f_TS+f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o3);
            mu[0][0][2]=(f_T+f_B+f_TE+f_BW+f_BE+f_TW+f_TN+f_BS+f_BN+f_TS+f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o3);
            mu[1][1][0]=(f_NE-f_SE+f_SW-f_NW+f_TNE-f_TSE+f_BNE-f_BSE+f_TSW-f_TNW+f_BSW-f_BNW)/(rho);
            mu[1][0][1]=(f_TE+f_BW-f_BE-f_TW+f_TNE-f_BNE+f_TSE-f_BSE-f_TNW+f_BNW-f_TSW+f_BSW)/(rho);
            mu[0][1][1]=(f_TN+f_BS-f_BN-f_TS+f_TNE-f_BNE-f_TSE+f_BSE+f_TNW-f_BNW-f_TSW+f_BSW)/(rho);
            mu[2][1][0]=(f_NE-f_SW-f_SE+f_NW+f_TNE+f_BNE-f_TSE-f_BSE+f_TNW+f_BNW-f_TSW-f_BSW)/(rho);
            mu[1][2][0]=(f_NE-f_SW+f_SE-f_NW+f_TNE-f_TNW+f_BNE-f_BNW+f_TSE-f_TSW+f_BSE-f_BSW)/(rho);
            mu[1][0][2]=(f_TE-f_BW+f_BE-f_TW+f_TNE-f_TNW+f_BNE-f_BNW+f_TSE-f_TSW+f_BSE-f_BSW)/(rho);
            mu[1][1][1]=(f_TNE-f_BNE-f_TSE+f_BSE-f_TNW+f_BNW+f_TSW-f_BSW)/(rho);
            mu[2][0][1]=(f_TE-f_BW-f_BE+f_TW+f_TNE-f_BNE+f_TSE-f_BSE+f_TNW-f_BNW+f_TSW-f_BSW)/(rho);
            mu[0][2][1]=(f_TN-f_BS-f_BN+f_TS+f_TNE-f_BNE+f_TSE-f_BSE+f_TNW-f_BNW+f_TSW-f_BSW)/(rho);
            mu[0][1][2]=(f_TN-f_BS+f_BN-f_TS+f_TNE+f_BNE-f_TSE-f_BSE+f_TNW+f_BNW-f_TSW-f_BSW)/(rho);
            mu[2][2][0]=(f_NE+f_SW+f_SE+f_NW+f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o9);
            mu[1][2][1]=(f_TNE-f_BNE+f_TSE-f_BSE-f_TNW+f_BNW-f_TSW+f_BSW)/(rho);
            mu[2][0][2]=(f_TE+f_BW+f_BE+f_TW+f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o9);
            mu[2][1][1]=(f_TNE-f_BNE-f_TSE+f_BSE+f_TNW-f_BNW-f_TSW+f_BSW)/(rho);
            mu[1][1][2]=(f_TNE+f_BNE-f_TSE-f_BSE-f_TNW-f_BNW+f_TSW+f_BSW)/(rho);
            mu[0][2][2]=(f_TN+f_BS+f_BN+f_TS+f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o9);
            mu[2][2][1]=(f_TNE-f_BNE+f_TSE-f_BSE+f_TNW-f_BNW+f_TSW-f_BSW)/(rho);
            mu[1][2][2]=(f_TNE-f_TNW+f_BNE-f_BNW+f_TSE-f_TSW+f_BSE-f_BSW)/(rho);
            mu[2][1][2]=(f_TNE+f_BNE-f_TSE-f_BSE+f_TNW+f_BNW-f_TSW-f_BSW)/(rho);
            mu[2][2][2]=(f_TNE+f_BNE+f_TSE+f_BSE+f_TNW+f_BNW+f_TSW+f_BSW)/(rho);//+(c1o27);
            mu[0][0][0]=one;
            mu[1][0][0]=vx;
            mu[0][1][0]=vy;
            mu[0][0][1]=vz;

            real M_zXX=zero;      real M_zYY=zero;      real M_zZZ=zero;      real M_zXY=zero;     real M_zXZ=zero;   real M_zYZ=zero;
            real M_zXXY=zero;     real M_zXYY=zero;     real M_zXXZ=zero;     real M_zXZZ=zero;    real M_zYYZ=zero;  real M_zYZZ=zero;  real M_zXYZ=zero;
            real M_zXXYY=zero;    real M_zXXZZ=zero;    real M_zYYZZ=zero;    real M_zXXYZ=zero;   real M_zXYYZ=zero; real M_zXYZZ=zero;
            real M_zXXYYZ=zero;   real M_zXXYZZ=zero;   real M_zXYYZZ=zero;   real M_zXXYYZZ=zero;
            real M_zXYeq=zero;    real M_zXZeq=zero;    real M_zYZeq=zero;    real M_zXYZeq=zero;
            real M_zXXYZeq=zero;  real M_zXYYZeq=zero;  real M_zXYZZeq=zero;
            real M_zXXYYZeq=zero; real M_zXXYZZeq=zero; real M_zXYYZZeq=zero; real M_zXXYYZZeq=zero;

            M_zXXYYZZ=  -five*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] * 
               mu[1][0][0] +     mu[0][0][2] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] * 
               mu[1][0][0] + four*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][1] *     mu[1][0][0] * 
               mu[1][0][0] - two*mu[0][1][0] *     mu[0][1][2] *     mu[1][0][0] *     mu[1][0][0] + 
               mu[0][0][1] *     mu[0][0][1] *     mu[0][2][0] *     mu[1][0][0] *     mu[1][0][0] -
               two*mu[0][0][1] *     mu[0][2][1] *     mu[1][0][0] *     mu[1][0][0] +     mu[0][2][2] *
               mu[1][0][0] *     mu[1][0][0] + four*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *
               mu[1][0][0] *     mu[1][0][1] - two*mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] *
               mu[1][0][2] + four*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] *
               mu[1][1][0] - eight*mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][1][1] +
               four*mu[0][1][0] *     mu[1][0][0] *     mu[1][1][2] - two*mu[0][0][1] *     mu[0][0][1] *
               mu[1][0][0] *     mu[1][2][0] + four*mu[0][0][1] *     mu[1][0][0] *     mu[1][2][1] -
               two*mu[1][0][0] *     mu[1][2][2] +     mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *
               mu[0][1][0] *     mu[2][0][0] - two*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *
               mu[2][0][1] +     mu[0][1][0] *     mu[0][1][0] *     mu[2][0][2] - two*mu[0][0][1] *
               mu[0][0][1] *     mu[0][1][0] *     mu[2][1][0] + four*mu[0][0][1] *     mu[0][1][0] *
               mu[2][1][1] - two*mu[0][1][0] *     mu[2][1][2] +     mu[0][0][1] *     mu[0][0][1] *
               mu[2][2][0] - two*mu[0][0][1] *     mu[2][2][1] +     mu[2][2][2];  //5
            M_zXXYYZ=    four*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][0] - 
               two*mu[0][1][0] *     mu[0][1][1] *     mu[1][0][0] *     mu[1][0][0] -     mu[0][0][1] *
               mu[0][2][0] *     mu[1][0][0] *     mu[1][0][0] +     mu[0][2][1] *     mu[1][0][0] *
               mu[1][0][0] - two*mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][1] -
               four*mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][1][0] + four*mu[0][1][0] *
               mu[1][0][0] *     mu[1][1][1] + two*mu[0][0][1] *     mu[1][0][0] *     mu[1][2][0] - 
               two*mu[1][0][0] *     mu[1][2][1] -     mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *
               mu[2][0][0] +     mu[0][1][0] *     mu[0][1][0] *     mu[2][0][1] + two*mu[0][0][1] *
               mu[0][1][0] *     mu[2][1][0] - two*mu[0][1][0] *     mu[2][1][1] -     mu[0][0][1] *
               mu[2][2][0] +     mu[2][2][1]; //6
            M_zXXYY=    -three*mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][0] +     mu[0][2][0] *
               mu[1][0][0] *     mu[1][0][0] + four*mu[0][1][0] *     mu[1][0][0] *     mu[1][1][0] -
               two*mu[1][0][0] *     mu[1][2][0] +     mu[0][1][0] *     mu[0][1][0] *     mu[2][0][0] -
               two*mu[0][1][0] *     mu[2][1][0] +     mu[2][2][0]; //7
            M_zXXYZZ=    four*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][0] -
               mu[0][0][2] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][0] - two*mu[0][0][1] *
               mu[0][1][1] *     mu[1][0][0] *     mu[1][0][0] +     mu[0][1][2] *     mu[1][0][0] *
               mu[1][0][0] - four*mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][1] +
               two*mu[0][1][0] *     mu[1][0][0] *     mu[1][0][2] - two*mu[0][0][1] *     mu[0][0][1] *
               mu[1][0][0] *     mu[1][1][0] + four*mu[0][0][1] *     mu[1][0][0] *     mu[1][1][1] -
               two*mu[1][0][0] *     mu[1][1][2] -     mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *
               mu[2][0][0] + two*mu[0][0][1] *     mu[0][1][0] *     mu[2][0][1] -     mu[0][1][0] *
               mu[2][0][2] +     mu[0][0][1] *     mu[0][0][1] *     mu[2][1][0] - two*mu[0][0][1] *
               mu[2][1][1] +     mu[2][1][2]; //8
            M_zXXYZ=    -three*mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] *     mu[1][0][0] +     mu[0][1][1] *
               mu[1][0][0] *     mu[1][0][0] + two*mu[0][1][0] *     mu[1][0][0] *     mu[1][0][1] +
               two*mu[0][0][1] *     mu[1][0][0] *     mu[1][1][0] - two*mu[1][0][0] *     mu[1][1][1] +
               mu[0][0][1] *     mu[0][1][0] *     mu[2][0][0] -     mu[0][1][0] *     mu[2][0][1] - 
               mu[0][0][1] *     mu[2][1][0] +     mu[2][1][1]; //9
            M_zXXY=      two*mu[0][1][0] *     mu[1][0][0] *     mu[1][0][0] - two*mu[1][0][0] *     mu[1][1][0] -
               mu[0][1][0] *     mu[2][0][0] +     mu[2][1][0]; //10
            M_zXXZZ=    -three*mu[0][0][1] *     mu[0][0][1] *     mu[1][0][0] *     mu[1][0][0] +     mu[0][0][2] *
               mu[1][0][0] *     mu[1][0][0] + four*mu[0][0][1] *     mu[1][0][0] *     mu[1][0][1] -
               two*mu[1][0][0] *     mu[1][0][2] +     mu[0][0][1] *     mu[0][0][1] *     mu[2][0][0] -
               two*mu[0][0][1] *     mu[2][0][1] +     mu[2][0][2]; //11
            M_zXXZ=      two*mu[0][0][1] *     mu[1][0][0] *     mu[1][0][0] - two*mu[1][0][0] *     mu[1][0][1] -
               mu[0][0][1] *     mu[2][0][0] +     mu[2][0][1]; //12
            M_zXX=-          mu[1][0][0] *     mu[1][0][0] +     mu[2][0][0]; //13
            M_zXYYZZ=    four*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] -
               mu[0][0][2] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] - four*mu[0][0][1] * 
               mu[0][1][0] *     mu[0][1][1] *     mu[1][0][0] + two*mu[0][1][0] *     mu[0][1][2] *
               mu[1][0][0] -     mu[0][0][1] *     mu[0][0][1] *     mu[0][2][0] *     mu[1][0][0] +
               two*mu[0][0][1] *     mu[0][2][1] *     mu[1][0][0] -     mu[0][2][2] *     mu[1][0][0] -
               two*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][1] +     mu[0][1][0] *
               mu[0][1][0] *     mu[1][0][2] - two*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *
               mu[1][1][0] + four*mu[0][0][1] *     mu[0][1][0] *     mu[1][1][1] - two*mu[0][1][0] *
               mu[1][1][2] +     mu[0][0][1] *     mu[0][0][1] *     mu[1][2][0] - two*mu[0][0][1] *
               mu[1][2][1] +     mu[1][2][2];    //14
            M_zXYYZ=    -three*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] + two*mu[0][1][0] *
               mu[0][1][1] *     mu[1][0][0] +     mu[0][0][1] *     mu[0][2][0] *     mu[1][0][0] - 
               mu[0][2][1] *     mu[1][0][0] +     mu[0][1][0] *     mu[0][1][0] *     mu[1][0][1] + 
               two*mu[0][0][1] *     mu[0][1][0] *     mu[1][1][0] - two*mu[0][1][0] *     mu[1][1][1] - 
               mu[0][0][1] *     mu[1][2][0] +     mu[1][2][1]; //15
            M_zXYY=      two*mu[0][1][0] *     mu[0][1][0] *     mu[1][0][0] -     mu[0][2][0] *     mu[1][0][0] -
               two*mu[0][1][0] *     mu[1][1][0] +     mu[1][2][0]; //16
            M_zXYZZ=    -three*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] +     mu[0][0][2] *
               mu[0][1][0] *     mu[1][0][0] + two*mu[0][0][1] *     mu[0][1][1] *     mu[1][0][0] -     
               mu[0][1][2] *     mu[1][0][0] + two*mu[0][0][1] *     mu[0][1][0] *     mu[1][0][1] -     
               mu[0][1][0] *     mu[1][0][2] +     mu[0][0][1] *     mu[0][0][1] *     mu[1][1][0] - 
               two*mu[0][0][1] *     mu[1][1][1] +     mu[1][1][2]; //17
            M_zXYZ=      two*mu[0][0][1] *     mu[0][1][0] *     mu[1][0][0] -     mu[0][1][1] *     mu[1][0][0] -
               mu[0][1][0] *     mu[1][0][1] -     mu[0][0][1] *     mu[1][1][0] +     mu[1][ 1][ 1]; //18
            M_zXY =         -mu[0][1][0] *     mu[1][0][0] +     mu[1][1][0]; //19
            M_zXZZ=      two*mu[0][0][1] *     mu[0][0][1] *     mu[1][0][0] -     mu[0][0][2] *     mu[1][0][0] -
               two*mu[0][0][1] *     mu[1][0][1] +     mu[1][0][2]; //20
            M_zXZ =         -mu[0][0][1] *     mu[1][0][0] +     mu[1][0][1]; //21
            M_zYYZZ=    -three*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] +     mu[0][0][2] *
               mu[0][1][0] *     mu[0][1][0] + four*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][1] -
               two*mu[0][1][0] *     mu[0][1][2] +     mu[0][0][1] *     mu[0][0][1] *     mu[0][2][0] -
               two*mu[0][0][1] *     mu[0][2][1] +     mu[0][2][2]; //22
            M_zYYZ=      two*mu[0][0][1] *     mu[0][1][0] *     mu[0][1][0] - two*mu[0][1][0] *     mu[0][1][1] -
               mu[0][0][1] *     mu[0][2][0] +     mu[0][2][1]; //23
            M_zYY  =        -mu[0][1][0] *     mu[0][1][0] +     mu[0][2][0]; //24
            M_zYZZ =     two*mu[0][0][1] *     mu[0][0][1] *     mu[0][1][0] -     mu[0][0][2] *     mu[0][1][0] -  
               two*mu[0][0][1] *     mu[0][1][1] +     mu[0][1][2]; //25
            M_zYZ  =        -mu[0][0][1] *     mu[0][1][0] +     mu[0][1][1]; //26
            M_zZZ  =        -mu[0][0][1] *     mu[0][0][1] +     mu[0][0][2]; //27

            M_zXXYYZZeq=c1o27;
            M_zXXYYZeq=zero;
            M_zXXYZZeq=zero;
            M_zXXYZeq=zero;
            M_zXYYZZeq=zero;
            M_zXYYZeq=zero;
            M_zXYZZeq=zero;

            real MXXpMYYpMZZ,MXXmMYY, MXXmMZZ, MXXYpMYZZ,MXXZpMYYZ,MXYYpMXZZ, MXXYmMYZZ,MXXZmMYYZ,MXYYmMXZZ, MXXYYppp,MXXYYpm2p, MXXYYppm2;
            real MXXpMYYpMZZeq,MXXmMYYeq, MXXmMZZeq, MXXYpMYZZeq,MXXZpMYYZeq,MXYYpMXZZeq, MXXYmMYZZeq,MXXZmMYYZeq,MXYYmMXZZeq, MXXYYpppeq,MXXYYpm2peq, MXXYYppm2eq;

            //coll faktoren:
            real w1,w2,w3,w4,w5,w6,w7,w8,w9,w10;
            w2=-one;//-omega;
            w7=-one;
            w9=-one;
            w10=-one;
            w1=-omega;
            w3=-one;
            w4=-one;
            w5=-one;
            w6=-one;
            w8=-one;
            ////////lin kombi bilden:
            MXXpMYYpMZZ = M_zXX + M_zYY + M_zZZ;
            MXXmMYY = M_zXX - M_zYY;
            MXXmMZZ = M_zXX - M_zZZ;

            MXXYpMYZZ  =  M_zXXY+M_zYZZ;
            MXXYmMYZZ  =  M_zXXY-M_zYZZ;
            MXXZpMYYZ  =  M_zXXZ+M_zYYZ;
            MXXZmMYYZ  =  M_zXXZ-M_zYYZ;
            MXYYpMXZZ  =  M_zXYY+M_zXZZ;
            MXYYmMXZZ  =  M_zXYY-M_zXZZ;

            MXXYYppp  = M_zXXYY + M_zXXZZ + M_zYYZZ;
            MXXYYpm2p = M_zXXYY - 2.0*M_zXXZZ + M_zYYZZ;
            MXXYYppm2 = M_zXXYY + M_zXXZZ - 2.0*M_zYYZZ;
            bool casc=false;
            bool factorCasc=true;
            if(casc)
            {
               MXXpMYYpMZZeq = one;
               MXXmMYYeq     = two;
               MXXmMZZeq     = two;

               MXXYpMYZZeq  =  two;
               MXXYmMYZZeq  =  two;
               MXXZpMYYZeq  =  two;
               MXXZmMYYZeq  =  two;
               MXYYpMXZZeq  =  two;
               MXYYmMXZZeq  =  two;

               MXXYYpppeq  = c1o3;
               MXXYYpm2peq = two;
               MXXYYppm2eq = two;

               //relaxation:
               MXXpMYYpMZZ -= w2*(MXXpMYYpMZZeq-MXXpMYYpMZZ);
               MXXmMYY -=     w1*(MXXmMYYeq    -MXXmMYY);
               MXXmMZZ -=     w1*(MXXmMZZeq    -MXXmMZZ);

               MXXYpMYZZ  -=  w3*(MXXYpMYZZeq-MXXYpMYZZ);
               MXXYmMYZZ  -=  w4*(MXXYmMYZZeq-MXXYmMYZZ);
               MXXZpMYYZ  -=  w3*(MXXZpMYYZeq-MXXZpMYYZ);
               MXXZmMYYZ  -=  w4*(MXXZmMYYZeq-MXXZmMYYZ);
               MXYYpMXZZ  -=  w3*(MXYYpMXZZeq-MXYYpMXZZ);
               MXYYmMXZZ  -=  w4*(MXYYmMXZZeq-MXYYmMXZZ);

               MXXYYppp  -= w7*(MXXYYpppeq -MXXYYppp );
               MXXYYpm2p -= w6*(MXXYYpm2peq-MXXYYpm2p);
               MXXYYppm2 -= w6*(MXXYYppm2eq-MXXYYppm2);

               M_zXXYYZZ=M_zXXYYZZ-w10*(M_zXXYYZZeq-M_zXXYYZZ);
               M_zXZ=    M_zXZ-    w1*(M_zXZeq   - M_zXZ    );
               M_zYZ=    M_zYZ-    w1*(M_zYZeq   - M_zYZ    );
               M_zXY=    M_zXY-    w1*(M_zXYeq   - M_zXY    );

               M_zXYZ=   M_zXYZ-   w5*(M_zXYZeq  - M_zXYZ   );

               M_zXYYZ=  M_zXYYZ-  w8*(M_zXYYZeq - M_zXYYZ  );
               M_zXYZZ=  M_zXYZZ-  w8*(M_zXYZZeq - M_zXYZZ  );
               M_zXXYZ=  M_zXXYZ-  w8*(M_zXXYZeq - M_zXXYZ  );

               M_zXXYYZ= M_zXXYYZ- w9*(M_zXXYYZeq- M_zXXYYZ );
               M_zXXYZZ= M_zXXYZZ- w9*(M_zXXYZZeq- M_zXXYZZ );
               M_zXYYZZ= M_zXYYZZ- w9*(M_zXYYZZeq- M_zXYYZZ );

               ////////von Lin Kombis zurueck:
               M_zXX =  c1o3 * MXXmMYY + c1o3 * MXXmMZZ +  c1o3 * MXXpMYYpMZZ;
               M_zYY = -c2o3 * MXXmMYY + c1o3 * MXXmMZZ +  c1o3 * MXXpMYYpMZZ;
               M_zZZ =  c1o3 * MXXmMYY - c2o3 * MXXmMZZ +  c1o3 * MXXpMYYpMZZ;

               M_zXXY = (MXXYmMYZZ + MXXYpMYZZ)*c1o2;
               M_zYZZ = c1o2*(-MXXYmMYZZ + MXXYpMYZZ);
               M_zXYY =(MXYYmMXZZ + MXYYpMXZZ)*c1o2;
               M_zXZZ = c1o2*(-MXYYmMXZZ + MXYYpMXZZ);
               M_zXXZ = (MXXZmMYYZ + MXXZpMYYZ)*c1o2;
               M_zYYZ = c1o2*(-MXXZmMYYZ + MXXZpMYYZ);

               M_zXXYY =  c1o3* MXXYYpm2p + c1o3* MXXYYppm2 + c1o3*MXXYYppp;
               M_zXXZZ = -c1o3* MXXYYpm2p + c1o3* MXXYYppp;
               M_zYYZZ = -c1o3* MXXYYppm2 + c1o3* MXXYYppp;
            }
            else if (factorCasc)
            {
               MXXpMYYpMZZeq = one;
               MXXmMYYeq     = two;
               MXXmMZZeq     = two;

               MXXYpMYZZeq  =  two;
               MXXYmMYZZeq  =  two;
               MXXZpMYYZeq  =  two;
               MXXZmMYYZeq  =  two;
               MXYYpMXZZeq  =  two;
               MXYYmMXZZeq  =  two;

               MXXYYpppeq  =c1o3;
               MXXYYpm2peq = two;
               MXXYYppm2eq = two;

               //relaxation:
               MXXpMYYpMZZ -= w2*(MXXpMYYpMZZeq-MXXpMYYpMZZ);
               MXXmMYY -=     w1*(MXXmMYYeq    -MXXmMYY);
               MXXmMZZ -=     w1*(MXXmMZZeq    -MXXmMZZ);

               MXXYpMYZZ  -=  w3*(MXXYpMYZZeq-MXXYpMYZZ);
               MXXYmMYZZ  -=  w4*(MXXYmMYZZeq-MXXYmMYZZ);
               MXXZpMYYZ  -=  w3*(MXXZpMYYZeq-MXXZpMYYZ);
               MXXZmMYYZ  -=  w4*(MXXZmMYYZeq-MXXZmMYYZ);
               MXYYpMXZZ  -=  w3*(MXYYpMXZZeq-MXYYpMXZZ);
               MXYYmMXZZ  -=  w4*(MXYYmMXZZeq-MXYYmMXZZ);

               M_zXZ=    M_zXZ-    w1*(M_zXZeq   - M_zXZ    );
               M_zYZ=    M_zYZ-    w1*(M_zYZeq   - M_zYZ    );
               M_zXY=    M_zXY-    w1*(M_zXYeq   - M_zXY    );
               M_zXYZ=   M_zXYZ-   w5*(M_zXYZeq  - M_zXYZ   );

               ////////von Lin Kombis zurueck:
               M_zXX =  c1o3 * MXXmMYY + c1o3 * MXXmMZZ +  c1o3 * MXXpMYYpMZZ;
               M_zYY = -c2o3 * MXXmMYY + c1o3 * MXXmMZZ +  c1o3 * MXXpMYYpMZZ;
               M_zZZ =  c1o3 * MXXmMYY - c2o3 * MXXmMZZ +  c1o3 * MXXpMYYpMZZ;

               M_zXXY = (MXXYmMYZZ + MXXYpMYZZ)*c1o2;
               M_zYZZ = c1o2*(-MXXYmMYZZ + MXXYpMYZZ);
               M_zXYY =(MXYYmMXZZ + MXYYpMXZZ)*c1o2;
               M_zXZZ = c1o2*(-MXYYmMXZZ + MXYYpMXZZ);
               M_zXXZ = (MXXZmMYYZ + MXXZpMYYZ)*c1o2;
               M_zYYZ = c1o2*(-MXXZmMYYZ + MXXZpMYYZ);

               //faktorisierte atraktoren:
               MXXYYpppeq=M_zXX*M_zYY+M_zXX*M_zZZ+M_zZZ*M_zYY;    
               MXXYYpm2peq=M_zXX*M_zYY-2.0*M_zXX*M_zZZ+M_zZZ*M_zYY; 
               MXXYYppm2eq=M_zXX*M_zYY+M_zXX*M_zZZ-2.0*M_zZZ*M_zYY; 
               M_zXYYZeq=M_zYY*M_zXZ;
               M_zXYZZeq=M_zZZ*M_zXY;
               M_zXXYZeq=M_zXX*M_zYZ;

               M_zXXYYZeq=zero;
               M_zXXYZZeq=zero;
               M_zXYYZZeq=zero;

               M_zXXYYZZeq=M_zXX*M_zYY*M_zZZ;

               MXXYYppp  -= w7*(MXXYYpppeq -MXXYYppp );
               MXXYYpm2p -= w6*(MXXYYpm2peq-MXXYYpm2p);
               MXXYYppm2 -= w6*(MXXYYppm2eq-MXXYYppm2);
               M_zXXYYZZ=M_zXXYYZZ-w10*(M_zXXYYZZeq-M_zXXYYZZ);
               M_zXYYZ=  M_zXYYZ-  w8*(M_zXYYZeq - M_zXYYZ  );
               M_zXYZZ=  M_zXYZZ-  w8*(M_zXYZZeq - M_zXYZZ  );
               M_zXXYZ=  M_zXXYZ-  w8*(M_zXXYZeq - M_zXXYZ  );

               M_zXXYYZ= M_zXXYYZ- w9*(M_zXXYYZeq- M_zXXYYZ );
               M_zXXYZZ= M_zXXYZZ- w9*(M_zXXYZZeq- M_zXXYZZ );
               M_zXYYZZ= M_zXYYZZ- w9*(M_zXYYZZeq- M_zXYYZZ );

               M_zXXYY =  c1o3* MXXYYpm2p + c1o3* MXXYYppm2 + c1o3*MXXYYppp;
               M_zXXZZ = -c1o3* MXXYYpm2p + c1o3* MXXYYppp;
               M_zYYZZ = -c1o3* MXXYYppm2 + c1o3* MXXYYppp;
            }
            //real forcingX1=1.11986e-8;//1.41342e-8;//5.48288e-9;//<-96//1.03093e-8;//1.59178e-8;//1.6256e-8; 

            //mu[1][0][0]+=forcingX1;

            mu[2][ 2][ 2]= vx*vx*vy*vy*vz*vz + vx*vx*vy*vy*M_zZZ + four*vx*vx*vy*vz*M_zYZ + two* vx*vx*vy*M_zYZZ +
               vx*vx*vz*vz*M_zYY + two* vx*vx*vz*M_zYYZ +vx*vx*M_zYYZZ + four*vx*vy*vy*vz*M_zXZ   +
               two* vx*vy*vy*M_zXZZ + four* vx*vy*vz*vz*M_zXY+   eight* vx*vy*vz*M_zXYZ + four* vx*vy*M_zXYZZ +
               two* vx*vz*vz*M_zXYY + four* vx*vz*M_zXYYZ + two*vx*M_zXYYZZ +   vy*vy*vz*vz*M_zXX + two* vy*vy*vz*M_zXXZ  +
               vy*vy*M_zXXZZ     + two* vy*vz*vz*M_zXXY  + four*vy*vz*M_zXXYZ +  two* vy*M_zXXYZZ + vz*vz*M_zXXYY +two* vz*M_zXXYYZ + M_zXXYYZZ;//-(c1o27);
            mu[2][ 2][ 1]= vx*vx*vy*vy*vz + two*vx*vx*vy*M_zYZ + vx*vx*vz*M_zYY +  vx*vx*M_zYYZ + two*vx*vy*vy*M_zXZ +four*vx*vy*vz*M_zXY +four*vx*vy*M_zXYZ+
               two*vx*vz*M_zXYY + two*vx*M_zXYYZ + vy*vy*vz*M_zXX +  vy*vy*M_zXXZ + two* vy*vz*M_zXXY + two*vy*M_zXXYZ +  vz*M_zXXYY + M_zXXYYZ;
            mu[2][ 2][ 0]= vx*vx*vy*vy + vx*vx*M_zYY + four*vx*vy*M_zXY +   two*vx*M_zXYY + vy*vy*M_zXX + two*vy*M_zXXY +  M_zXXYY ;//-(c1o9);
            mu[2][ 1][ 2]= vx*vx*vy*vz*vz + vx*vx*vy*M_zZZ + two*vx*vx*vz*M_zYZ +   vx*vx*M_zYZZ + four*vx*vy*vz*M_zXZ+ //vy->vz
               two*vx*vy*M_zXZZ + two*vx*vz*vz*M_zXY +   four*vx*vz*M_zXYZ + two*vx*M_zXYZZ + vy*vz*vz*M_zXX +
               two*vy*vz*M_zXXZ + vy*M_zXXZZ + vz*vz*M_zXXY +  two*vz*M_zXXYZ + M_zXXYZZ;
            mu[2][ 1][ 1]= vx*vx*vy*vz + vx*vx*M_zYZ + two*vx*vy*M_zXZ + two*vx*vz*M_zXY + two*vx*M_zXYZ + vy*vz*M_zXX+  vy*M_zXXZ + vz*M_zXXY + M_zXXYZ;
            mu[2][ 1][ 0]= vx*vx*vy + two*vx*M_zXY + vy*M_zXX + M_zXXY;
            mu[2][ 0][ 2]= vx*vx*vz*vz + vx*vx*M_zZZ + four*vx*vz*M_zXZ+ two*vx*M_zXZZ + vz*vz*M_zXX+ two*vz*M_zXXZ + M_zXXZZ;//-(c1o9);
            mu[2][ 0][ 1]= vx*vx*vz + two*vx*M_zXZ + vz*M_zXX + M_zXXZ;
            mu[2][ 0][ 0]= vx*vx + M_zXX ;//-(c1o3);
            mu[1][ 2][ 2]= vx*vy*vy*vz*vz + vx*vy*vy*M_zZZ + four*vx*vy*vz*M_zYZ +   two*vx*vy*M_zYZZ + vx*vz*vz*M_zYY + two*vx*vz*M_zYYZ +vx*M_zYYZZ + two*vy*vy*vz*M_zXZ + vy*vy*M_zXZZ +
               two*vy*vz*vz*M_zXY + four*vy*vz*M_zXYZ + two*vy*M_zXYZZ +vz*vz*M_zXYY + two*vz*M_zXYYZ + M_zXYYZZ;
            mu[1][ 2][ 1]= vx*vy*vy*vz + two*vx*vy*M_zYZ + vx*vz*M_zYY +   vx*M_zYYZ + vy*vy*M_zXZ + two*vy*vz*M_zXY +   two*vy*M_zXYZ + vz*M_zXYY + M_zXYYZ;
            mu[1][ 2][ 0]= vx*vy*vy + vx*M_zYY + two*vy*M_zXY + M_zXYY;
            mu[1][ 1][ 2]= vx*vy*vz*vz+vx*vy*M_zZZ+two*vx*vz*M_zYZ+vx*M_zYZZ+two*vy*vz*M_zXZ+vy*M_zXZZ+vz*vz*M_zXY+two*vz*M_zXYZ+M_zXYZZ;
            mu[1][ 1][ 1]= vx*vy*vz + vx*M_zYZ + vy*M_zXZ + vz*M_zXY+M_zXYZ;
            mu[1][ 1][ 0]= vx*vy+ M_zXY;
            mu[1][ 0][ 2]= vx*vz*vz + vx*M_zZZ + two*vz*M_zXZ + M_zXZZ;
            mu[1][ 0][ 1]= vx*vz + M_zXZ;
            //mu[1][ 0][ 0]= vx;
            mu[0][ 2][ 2]= vy*vy*vz*vz + vy*vy*M_zZZ + four*vy*vz*M_zYZ +two*vy*M_zYZZ + vz*vz*M_zYY + two*vz*M_zYYZ +  M_zYYZZ;//-(c1o9);
            mu[0][ 2][ 1]= vy*vy*vz + two*vy*M_zYZ + vz*M_zYY + M_zYYZ;
            mu[0][ 2][ 0]= vy*vy + M_zYY ;//-(c1o3);
            mu[0][ 1][ 2]= vy*vz*vz+ vy*M_zZZ + two*vz*M_zYZ + M_zYZZ;
            mu[0][ 1][ 1]= vy*vz + M_zYZ;
            //mu[0][ 1][ 0]= vy ;                     //richtig?
            mu[0][ 0][ 2]= vz*vz + M_zZZ;// -(c1o3);
            //mu[0][ 0][ 1]= vz;
            //mu[0][ 0][ 0]=rho0 ;

            f_E =c1o2* (mu[2][0][0] -  mu[2][2][0] + mu[2][2][2] - mu[2][0][2] -   mu[1][2][0] + mu[1][2][2] - mu[1][0][2] +mu[1][0][0]) *(rho);
            f_W =c1o2* (mu[2][0][0] - mu[2][2][0] + mu[2][2][2] - mu[2][0][2] +    mu[1][2][0] - mu[1][2][2] + mu[1][0][2] -mu[1][0][0])*  (rho);
            f_N =c1o2* (-mu[2][1][0] - mu[2][2][0] + mu[2][2][2] + mu[2][1][2] +   mu[0][2][0] - mu[0][2][2] - mu[0][1][2] +mu[0][1][0])* (rho);
            f_S =c1o2* (mu[2][1][0] -  mu[2][2][0] + mu[2][2][2] - mu[2][1][2] +   mu[0][2][0] - mu[0][2][2] + mu[0][1][2] -mu[0][1][0])* (rho);
            f_T =c1o2* (mu[2][2][1] +  mu[2][2][2] - mu[2][0][1] - mu[2][0][2] -   mu[0][2][1] - mu[0][2][2] + mu[0][0][2] +mu[0][0][1])* (rho);
            f_B =c1o2* (-mu[2][2][1] + mu[2][2][2] + mu[2][0][1]  - mu[2][0][2] +  mu[0][2][1] - mu[0][2][2] + mu[0][0][2]-mu[0][0][1])* (rho);
            f_NE =c1o4*( mu[2][1][0]  + mu[2][2][0]- mu[2][2][2] - mu[2][1][2] +  mu[1][1][0]+ mu[1][2][0]- mu[1][2][2] -mu[1][1][2])*  (rho);
            f_SW =c1o4*(-mu[2][1][0] + mu[2][2 ][0]- mu[2][2][2] + mu[2][1][2] +  mu[1][1][0]- mu[1][2][0]+ mu[1][2][2] -mu[1][1][2])*  (rho);
            f_SE =c1o4*(-mu[2][1][0] + mu[2][2 ][0]- mu[2][2][2] + mu[2][1][2] -  mu[1][1][0]+ mu[1][2][0]- mu[1][2][2] +mu[1][1][2])*  (rho);
            f_NW =c1o4*( mu[2][1][0]  + mu[2][2][0]- mu[2][2][2] - mu[2][1][2] -  mu[1][1][0]- mu[1][2][0]+ mu[1][2][2] + mu[1][1][2])* (rho);
            f_TE=c1o4*(-mu[2][2][1] - mu[2][2][2] + mu[2][0][1] + mu[2][0][2] -   mu[1][2][1] - mu[1][2][2] + mu[1][0][1] + mu[1][0][2])*(rho);
            f_BW=c1o4*( mu[2][2][1]  -mu[2][2][2] - mu[2][0][1] + mu[2][0][2] -   mu[1][2][1] + mu[1][2][2] + mu[1][0][1] - mu[1][0][2])*(rho);
            f_BE=c1o4*(mu[2][2][1] - mu[2][2][2] - mu[2][0][1] + mu[2][0][2] +    mu[1][2][1] - mu[1][2][2] - mu[1][0][1] +mu[1][0][2])*  (rho);
            f_TW=c1o4*(-mu[2][2][1] - mu[2][2][2] + mu[2][0][1] + mu[2][0][2] +   mu[1][2][1] + mu[1][2][2] - mu[1][0][1] -mu[1][0][2])* (rho);
            f_TN=c1o4*(-mu[2][2][1] - mu[2][2][2] - mu[2][1][1] - mu[2][1][2] +   mu[0][2][1] + mu[0][2][2] + mu[0][1][1]+mu[0][1][2])*  (rho);
            f_BS=c1o4*( mu[2][2][1] - mu[2][2][2] - mu[2][1][1] + mu[2][1][2] -   mu[0][2][1] + mu[0][2][2] + mu[0][1][1] - mu[0][1][2])*(rho);
            f_BN=c1o4*( mu[2][2][1] - mu[2][2][2] + mu[2][1][1] - mu[2][1][2] -   mu[0][2][1] + mu[0][2][2] - mu[0][1][1] + mu[0][1][2])*(rho);
            f_TS=c1o4*(-mu[2][2][1] - mu[2][2][2] + mu[2][1][1] + mu[2][1][2] +   mu[0][2][1] + mu[0][2][2] - mu[0][1][1] -mu[0][1][2])* (rho);
            f_ZERO=     (-mu[2][0][0] + mu[2][2][0] - mu[2][2][2] + mu[2][0][2] - mu[0][2][0] + mu[0][2][2] - mu[0][0][2] + mu[0][0][0])*(rho);
            f_TNE=c1o8*( mu[2][2][1] + mu[2][2][2] + mu[2][1][1] + mu[2][1][2] + mu[1][2][1] + mu[1][2][2] + mu[1][1][1] + mu[1][1][2])*  (rho);
            f_BNE=c1o8*(-mu[2][2][1] + mu[2][2][2] -mu[2][1][1] + mu[2][1][2]  - mu[1][2][1] + mu[1][2][2] -mu[1][1][1] + mu[1][1][2])*     (rho);
            f_TSE=c1o8*( mu[2][2][1] + mu[2][2][2] - mu[2][1][1] - mu[2][1][2] + mu[1][2][1] + mu[1][2][2] - mu[1][1][1] - mu[1][1][2])*  (rho);
            f_BSE=c1o8*(-mu[2][2][1] + mu[2][2][2] +mu[2][1][1] - mu[2][1][2]  - mu[1][2][1] + mu[1][2][2] +mu[1][1][1] - mu[1][1][2])*     (rho);
            f_TNW=c1o8*( mu[2][2][1] + mu[2][2][2] + mu[2][1][1] + mu[2][1][2] - mu[1][2][1] - mu[1][2][2] - mu[1][1][1] - mu[1][1][2])*  (rho);
            f_BNW=c1o8*(-mu[2][2][1] + mu[2][2][2] -mu[2][1][1] + mu[2][1][2] +  mu[1][2][1] - mu[1][2][2] +mu[1][1][1] - mu[1][1][2])*     (rho);
            f_TSW=c1o8*( mu[2][2][1] + mu[2][2][2] - mu[2][1][1] - mu[2][1][2] - mu[1][2][1] - mu[1][2][2] + mu[1][1][1] + mu[1][1][2])*  (rho);
            f_BSW=c1o8*(-mu[2][2][1] + mu[2][2][2]+mu[2][1][1] - mu[2][1][2]+    mu[1][2][1] - mu[1][2][2]-mu[1][1][1] + mu[1][1][2])*(rho);
         }
      }

      (D.f[ dirE  ])[ke ] = f_W-c2over27;
      (D.f[ dirW  ])[kw ] = f_E-c2over27;

      (D.f[ dirN  ])[kn ] = f_S-c2over27;
      (D.f[ dirS  ])[ks ] = f_N-c2over27;
      (D.f[ dirT  ])[kt ] = f_B-c2over27;
      (D.f[ dirB  ])[kb ] = f_T-c2over27;

      (D.f[ dirNE ])[kne] = f_SW-c1over54;
      (D.f[ dirSW ])[ksw] = f_NE-c1over54;
      (D.f[ dirSE ])[kse] = f_NW-c1over54;
      (D.f[ dirNW ])[knw] = f_SE-c1over54;
      (D.f[ dirTE ])[kte] = f_BW-c1over54;
      (D.f[ dirBW ])[kbw] = f_TE-c1over54;
      (D.f[ dirBE ])[kbe] = f_TW-c1over54;
      (D.f[ dirTW ])[ktw] = f_BE-c1over54;

      (D.f[ dirTN ])[ktn] = f_BS-c1over54;
      (D.f[ dirBS ])[kbs] = f_TN-c1over54;
      (D.f[ dirBN ])[kbn] = f_TS-c1over54;
      (D.f[ dirTS ])[kts] = f_BN-c1over54;

      (D.f[ dirZERO])[k] = f_ZERO-c8over27;

      (D.f[ dirTNE ])[ktne] = f_BSW-c1over216;
      (D.f[ dirTSE ])[ktse] = f_BNW-c1over216;
      (D.f[ dirBNE ])[kbne] = f_TSW-c1over216;
      (D.f[ dirBSE ])[kbse] = f_TNW-c1over216;
      (D.f[ dirTNW ])[ktnw] = f_BSE-c1over216;
      (D.f[ dirTSW ])[ktsw] = f_BNE-c1over216;
      (D.f[ dirBNW ])[kbnw] = f_TSE-c1over216;
      (D.f[ dirBSW ])[kbsw] = f_TNE-c1over216;
   }
   __syncthreads();
}
