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
			// or much sloooowaaaa ... it dep�ndssssss on sadaku
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
			// or much sloooowaaaa ... it dep�ndssssss on sadaku
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
			// or much sloooowaaaa ... it dep�ndssssss on sadaku
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
extern "C" __global__ void LB_Kernel_Kum_New_Comp_SRT_SP_27(
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

		if (BC >= GEO_FLUID)
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

			real rho = one + drho;
			////////////////////////////////////////////////////////////////////////////////////
			//slow
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
			real fx = forces[0] / (pow(two, level)); 
			real fy = forces[1] / (pow(two, level)); 
			real fz = forces[2] / (pow(two, level)); 
			vvx += fx*c1o2;
			vvy += fy*c1o2;
			vvz += fz*c1o2;
			////////////////////////////////////////////////////////////////////////////////////
			//real omega = omega_in;
			////////////////////////////////////////////////////////////////////////////////////
			//fast
			real oMdrho = one; // comp special
			real m0, m1, m2;
			real vx2;
			real vy2;
			real vz2;
			vx2 = vvx*vvx;
			vy2 = vvy*vvy;
			vz2 = vvz*vvz;
			//////////////////////////////////////////////////////////////////////////////////////
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


			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = omega;

			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			real OxyyPxzz = omega; 
			real OxyyMxzz = omega; 
			real Oxyz = omega; 
			////////////////////////////////////////////////////////////
			//4.
			//////////////////////////////
			real O4 = omega;
			////////////////////////////////////////////////////////////
			//5.
			//////////////////////////////
			real O5 = omega;
			////////////////////////////////////////////////////////////
			//6.
			//////////////////////////////
			real O6 = omega;
			////////////////////////////////////////////////////////////


			//central moments to cumulants
			//4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + two * mfbba * mfbab) / rho;	//ab 15.05.2015 verwendet
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + two * mfbba * mfabb) / rho; //ab 15.05.2015 verwendet
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + two * mfbab * mfabb) / rho; //ab 15.05.2015 verwendet

			real CUMcca = mfcca - (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
			real CUMcac = mfcac - (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
			real CUMacc = mfacc - (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

			//5.
			real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)) / rho;
			real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)) / rho;
			real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)) / rho;

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
			//+ c1o27*(one -three/rho +two/(rho*rho)));


			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy = mfcaa - mfaca;
			real mxxMzz = mfcaa - mfaac;

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
				mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz) - three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
				mxxMyy += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
				mxxMzz += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

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
			mfabb += omega * (-mfabb);
			mfbab += omega * (-mfbab);
			mfbba += omega * (-mfbba);

			//////////////////////////////////////////////////////////////////////////
			//limiter-Scheise Teil 3
			//mfabb     += oyz * (-mfabb);
			//mfbab     += oxz * (-mfbab);
			//mfbba     += oxy * (-mfbba);
			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-two*  mxxMyy + mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (mxxMyy - two* mxxMzz + mxxPyyPzz);

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
			//wadjust    = Oxyz+(one-Oxyz)*abs(mfbbb)/(abs(mfbbb)+qudricLimitD);
			//mfbbb     += wadjust * (-mfbbb);
			//wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimitP);
			//mxxyPyzz  += wadjust * (-mxxyPyzz);
			//wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimitM);
			//mxxyMyzz  += wadjust * (-mxxyMyzz);
			//wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimitP);
			//mxxzPyyz  += wadjust * (-mxxzPyyz);
			//wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimitM);
			//mxxzMyyz  += wadjust * (-mxxzMyyz);
			//wadjust    = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimitP);
			//mxyyPxzz  += wadjust * (-mxyyPxzz);
			//wadjust    = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimitM);
			//mxyyMxzz  += wadjust * (-mxyyMxzz);
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
			mfbbb += OxyyMxzz * (-mfbbb);
			mxxyPyzz += OxyyPxzz * (-mxxyPyzz);
			mxxyMyzz += OxyyMxzz * (-mxxyMyzz);
			mxxzPyyz += OxyyPxzz * (-mxxzPyyz);
			mxxzMyyz += OxyyMxzz * (-mxxzMyyz);
			mxyyPxzz += OxyyPxzz * (-mxyyPxzz);
			mxyyMxzz += OxyyMxzz * (-mxyyMxzz);
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

			mfcca = CUMcca + (((mfcaa * mfaca + two * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
			mfcac = CUMcac + (((mfcaa * mfaac + two * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
			mfacc = CUMacc + (((mfaac * mfaca + two * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

	        //5.
			mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)) / rho;
			mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)) / rho;
			mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)) / rho;

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






