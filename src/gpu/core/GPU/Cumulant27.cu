//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
/* Device code */
#include "LBM/LB.h" 
#include "basics/constants/NumericConstants.h"
#include "lbm/constants/D3Q27.h"


using namespace vf::basics::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////
__global__ void LB_Kernel_Kum_AA2016_Comp_Bulk_SP_27(real omega,
																unsigned int* bcMatD,
																unsigned int* neighborX,
																unsigned int* neighborY,
																unsigned int* neighborZ,
																real* DDStart,
																unsigned long long numberOfLBnodes,
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

	if(k<numberOfLBnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			real mfcbb = (D.f[dP00])[k  ];//[ke   ];// +  c2over27 ;(D.f[dP00])[k  ];//ke
			real mfabb = (D.f[dM00])[kw ];//[kw   ];// +  c2over27 ;(D.f[dM00])[kw ];
			real mfbcb = (D.f[d0P0])[k  ];//[kn   ];// +  c2over27 ;(D.f[d0P0])[k  ];//kn
			real mfbab = (D.f[d0M0])[ks ];//[ks   ];// +  c2over27 ;(D.f[d0M0])[ks ];
			real mfbbc = (D.f[d00P])[k  ];//[kt   ];// +  c2over27 ;(D.f[d00P])[k  ];//kt
			real mfbba = (D.f[d00M])[kb ];//[kb   ];// +  c2over27 ;(D.f[d00M])[kb ];
			real mfccb = (D.f[dPP0])[k  ];//[kne  ];// +  c1over54 ;(D.f[dPP0])[k  ];//kne
			real mfaab = (D.f[dMM0])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dMM0])[ksw];
			real mfcab = (D.f[dPM0])[ks ];//[kse  ];// +  c1over54 ;(D.f[dPM0])[ks ];//kse
			real mfacb = (D.f[dMP0])[kw ];//[knw  ];// +  c1over54 ;(D.f[dMP0])[kw ];//knw
			real mfcbc = (D.f[dP0P])[k  ];//[kte  ];// +  c1over54 ;(D.f[dP0P])[k  ];//kte
			real mfaba = (D.f[dM0M])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dM0M])[kbw];
			real mfcba = (D.f[dP0M])[kb ];//[kbe  ];// +  c1over54 ;(D.f[dP0M])[kb ];//kbe
			real mfabc = (D.f[dM0P])[kw ];//[ktw  ];// +  c1over54 ;(D.f[dM0P])[kw ];//ktw
			real mfbcc = (D.f[d0PP])[k  ];//[ktn  ];// +  c1over54 ;(D.f[d0PP])[k  ];//ktn
			real mfbaa = (D.f[d0MM])[kbs];//[kbs  ];// +  c1over54 ;(D.f[d0MM])[kbs];
			real mfbca = (D.f[d0PM])[kb ];//[kbn  ];// +  c1over54 ;(D.f[d0PM])[kb ];//kbn
			real mfbac = (D.f[d0MP])[ks ];//[kts  ];// +  c1over54 ;(D.f[d0MP])[ks ];//kts
			real mfbbb = (D.f[d000])[k  ];//[kzero];// +  c8over27 ;(D.f[d000])[k  ];//kzero
			real mfccc = (D.f[dPPP])[k  ];//[ktne ];// +  c1over216;(D.f[dPPP])[k  ];//ktne
			real mfaac = (D.f[dMMP])[ksw];//[ktsw ];// +  c1over216;(D.f[dMMP])[ksw];//ktsw
			real mfcac = (D.f[dPMP])[ks ];//[ktse ];// +  c1over216;(D.f[dPMP])[ks ];//ktse
			real mfacc = (D.f[dMPP])[kw ];//[ktnw ];// +  c1over216;(D.f[dMPP])[kw ];//ktnw
			real mfcca = (D.f[dPPM])[kb ];//[kbne ];// +  c1over216;(D.f[dPPM])[kb ];//kbne
			real mfaaa = (D.f[dMMM])[kbsw];//[kbsw ];// +  c1over216;(D.f[dMMM])[kbsw];
			real mfcaa = (D.f[dPMM])[kbs];//[kbse ];// +  c1over216;(D.f[dPMM])[kbs];//kbse
			real mfaca = (D.f[dMPM])[kbw];//[kbnw ];// +  c1over216;(D.f[dMPM])[kbw];//kbnw
			////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
							(((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
							((mfabb+mfcbb) + (mfbab+mfbcb) + (mfbba+mfbbc))) + mfbbb;

			real rho = c1o1+drho;
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
			real fx = forces[0] / (pow((double)c2o1, (double)level)); // zero;//0.0032653/(pow(two,level)); //0.000000005;//(two/1600000.0) / 120.0; //
            real fy = forces[1] / (pow((double)c2o1, (double)level)); // zero;
            real fz = forces[2] / (pow((double)c2o1, (double)level)); // zero;
			vvx += fx;
			vvy += fy;
			vvz += fz;
			////////////////////////////////////////////////////////////////////////////////////
			//real omega = omega_in;
			////////////////////////////////////////////////////////////////////////////////////
			real oMdrho = c1o1; // comp special
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
			mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - c2o1*	m1 * vvz + vz2 * m0;
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
			mfaca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += c1o1* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = c4o1 * omega / (c4o1 + c3o1 * gamma * (c2o1 - omega));//one;	//set the bulk viscosity one is high / two is very low and zero is (too) high ... (also called omega 2)

			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			//no bulk viscosity
			//real OxyyPxzz  = eight*(-two+omega)*(one+two*omega)/(-eight-fourteen*omega+seven*omega*omega);//one;
			//real OxyyMxzz  = eight*(-two+omega)*(-seven+four*omega)/(fiftysix-fifty*omega+nine*omega*omega);//one;
			//real Oxyz      = twentyfour*(-two+omega)*(-two-seven*omega+three*omega*omega)/(fourtyeight+c152*omega-c130*omega*omega+twentynine*omega*omega*omega);//one;
			//with bulk viscosity
			real OxyyPxzz  = c8o1*(-c2o1+omega)*(OxxPyyPzz * (-c1o1+c3o1*omega) - c5o1 * omega)/(c8o1 * (c5o1 - c2o1 * omega) * omega + OxxPyyPzz * ( c8o1 + omega * ( c9o1 * omega - c26o1)));//one;
			real OxyyMxzz  = c8o1*(-c2o1+omega)*(omega + OxxPyyPzz*(c3o1 * omega - c7o1))/(OxxPyyPzz * (c56o1 - c42o1 * omega + c9o1 * omega * omega) - c8o1 * omega);//one;
			real Oxyz      = c24o1*(-c2o1+omega)*(c4o1*omega*omega + omega * OxxPyyPzz * (c18o1 - c13o1 * omega) + OxxPyyPzz * OxxPyyPzz * (c2o1 + omega * ( c6o1 * omega - c11o1)))/
								(c16o1 * omega * omega * ( omega - c6o1 ) - c2o1 * omega * OxxPyyPzz *( c216o1 + c5o1  * omega * ( c9o1 * omega - c46o1 )) + 
								OxxPyyPzz * OxxPyyPzz * ( omega * ( c3o1 * omega - c10o1 ) * ( c15o1 * omega - c28o1 ) - c48o1 ));//one;
			////////////////////////////////////////////////////////////
			//4.
			//////////////////////////////
			real O4        = c1o1;
			//////////////////////////////
			//real O4        = omega;//TRT
			////////////////////////////////////////////////////////////
			//5.
			//////////////////////////////
			real O5        = c1o1;
			////////////////////////////////////////////////////////////
			//6.
			//////////////////////////////
			real O6        = c1o1;
			////////////////////////////////////////////////////////////


			//central moments to cumulants
			//4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho;
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho;
				  	 		
			real CUMcca = mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			real CUMcac = mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			real CUMacc = mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			real CUMccc = mfccc + ((-c4o1 *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     c2o1 * (mfcaa * mfaca * mfaac)
							+ c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));




			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;
			
			////////////////////////////////////////////////////////////////////////////
            real Dxy =-c3o1*omega*mfbba;
            real Dxz =-c3o1*omega*mfbab;
            real Dyz =-c3o1*omega*mfabb;

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
				wadjust    = OxxPyyPzz+(c1o1-OxxPyyPzz)*abs(mfaaa  - mxxPyyPzz- c3o1 * (c1o1/OxxPyyPzz - c1o2 ) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz))/
					         (abs(mfaaa  - mxxPyyPzz- c3o1 * (c1o1/OxxPyyPzz - c1o2 ) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz))+qudricLimitOmega2);
 				mxxPyyPzz += wadjust*(mfaaa  - mxxPyyPzz)- c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 				//without limiter (no bulk)
				//mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 				mxxMyy    += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
 				mxxMzz    += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
 
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
			mfaca = c1o3 * (-c2o1*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - c2o1* mxxMzz + mxxPyyPzz);


			//relax
			//////////////////////////////////////////////////////////////////////////
			//das ist der limiter
 			wadjust    = Oxyz+(c1o1-Oxyz)*abs(mfbbb)/(abs(mfbbb)+qudricLimitD);
 			mfbbb     += wadjust * (-mfbbb);
 			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimitP);
 			mxxyPyzz  += wadjust * (-mxxyPyzz);
 			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimitM);
 			mxxyMyzz  += wadjust * (-mxxyMyzz);
 			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimitP);
 			mxxzPyyz  += wadjust * (-mxxzPyyz);
 			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimitM);
 			mxxzMyyz  += wadjust * (-mxxzMyyz);
 			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimitP);
 			mxyyPxzz  += wadjust * (-mxyyPxzz);
 			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimitM);
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
			real factorA = ( c4o1 * omega * omega + c2o1 * omega * OxxPyyPzz * ( omega - c6o1 ) + OxxPyyPzz * OxxPyyPzz * ( omega * ( c10o1 - c3o1 * omega ) - c4o1 )) /
				        ( ( omega - OxxPyyPzz ) * ( OxxPyyPzz * ( c2o1 + c3o1 * omega ) - c8o1 * omega ) );
			real factorB = ( c4o1 * omega * OxxPyyPzz * ( c9o1 * omega - c16o1 ) - c4o1 * omega * omega - c2o1 * OxxPyyPzz * OxxPyyPzz * ( c2o1 + c9o1 * omega * ( omega - c2o1 ))) /
				        ( c3o1 * ( omega - OxxPyyPzz ) * ( OxxPyyPzz * ( c2o1 + c3o1 * omega ) - c8o1 * omega ) );
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
			CUMacc = -O4 * (c1o1 / omega - c1o2) * (dyuy + dzuz) * c2o3 * factorA + (c1o1 - O4) * (CUMacc);
			CUMcac = -O4 * (c1o1 / omega - c1o2) * (dxux + dzuz) * c2o3 * factorA + (c1o1 - O4) * (CUMcac);
			CUMcca = -O4 * (c1o1 / omega - c1o2) * (dyuy + dxux) * c2o3 * factorA + (c1o1 - O4) * (CUMcca);
			CUMbbc = -O4 * (c1o1 / omega - c1o2) * Dxy           * c1o3 * factorB + (c1o1 - O4) * (CUMbbc);
			CUMbcb = -O4 * (c1o1 / omega - c1o2) * Dxz           * c1o3 * factorB + (c1o1 - O4) * (CUMbcb);
			CUMcbb = -O4 * (c1o1 / omega - c1o2) * Dyz           * c1o3 * factorB + (c1o1 - O4) * (CUMcbb);


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
						   
			mfcca = CUMcca + (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			mfcac = CUMcac + (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			mfacc = CUMacc + (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			mfccc = CUMccc - ((-c4o1 *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     c2o1 * (mfcaa * mfaca * mfaac)
							+ c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
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
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + c1o1* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - c2o1* mfaab *  vvz         +  mfaaa                * (c1o1- vz2)              - c1o1* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + c1o1* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - c2o1* mfabb *  vvz         + mfaba * (c1o1- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - c2o1* mfacb *  vvz         +  mfaca                  * (c1o1- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - c2o1* mfbab *  vvz         + mfbaa * (c1o1- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - c2o1* mfbbb *  vvz         + mfbba * (c1o1- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - c2o1* mfbcb *  vvz         + mfbca * (c1o1- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - c2o1* mfcab *  vvz         +  mfcaa                  * (c1o1- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - c2o1* mfcbb *  vvz         + mfcba * (c1o1- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - c2o1* mfccb *  vvz         +  mfcca                  * (c1o1- vz2)              - c1o9 * oMdrho * vz2; 
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
			m1 = -mfaca        - c2o1* mfaba *  vvy         +  mfaaa                  * (c1o1- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - c2o1* mfabb *  vvy         +  mfaab                  * (c1o1- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - c2o1* mfabc *  vvy         +  mfaac                  * (c1o1- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - c2o1* mfbba *  vvy         + mfbaa * (c1o1- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - c2o1* mfbbb *  vvy         + mfbab * (c1o1- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - c2o1* mfbbc *  vvy         + mfbac * (c1o1- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - c2o1* mfcba *  vvy         +  mfcaa                   * (c1o1- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - c2o1* mfcbb *  vvy         +  mfcab                  * (c1o1- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - c2o1* mfcbc *  vvy         +  mfcac                   * (c1o1- vy2)              - c1o18 * oMdrho * vy2; 
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
			m1 = -mfcaa        - c2o1* mfbaa *  vvx         +  mfaaa                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - c2o1* mfbba *  vvx         +  mfaba                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - c2o1* mfbca *  vvx         +  mfaca                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - c2o1* mfbab *  vvx         +  mfaab                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - c2o1* mfbbb *  vvx         +  mfabb                  * (c1o1- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - c2o1* mfbcb *  vvx         +  mfacb                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - c2o1* mfbac *  vvx         +  mfaac                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - c2o1* mfbbc *  vvx         +  mfabc                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - c2o1* mfbcc *  vvx         +  mfacc                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dP00   ])[k   ] = mfabb;//(D.f[ dP00   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dP00   ])[k   ]                                                                     
			(D.f[ dM00   ])[kw  ] = mfcbb;//(D.f[ dM00   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dM00   ])[kw  ]                                                                   
			(D.f[ d0P0   ])[k   ] = mfbab;//(D.f[ d0P0   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ d0P0   ])[k   ]
			(D.f[ d0M0   ])[ks  ] = mfbcb;//(D.f[ d0M0   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ d0M0   ])[ks  ]
			(D.f[ d00P   ])[k   ] = mfbba;//(D.f[ d00P   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ d00P   ])[k   ]
			(D.f[ d00M   ])[kb  ] = mfbbc;//(D.f[ d00M   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ d00M   ])[kb  ]
			(D.f[ dPP0  ])[k   ] = mfaab;//(D.f[ dPP0  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dPP0  ])[k   ]
			(D.f[ dMM0  ])[ksw ] = mfccb;//(D.f[ dMM0  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dMM0  ])[ksw ]
			(D.f[ dPM0  ])[ks  ] = mfacb;//(D.f[ dPM0  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dPM0  ])[ks  ]
			(D.f[ dMP0  ])[kw  ] = mfcab;//(D.f[ dMP0  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dMP0  ])[kw  ]
			(D.f[ dP0P  ])[k   ] = mfaba;//(D.f[ dP0P  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dP0P  ])[k   ]
			(D.f[ dM0M  ])[kbw ] = mfcbc;//(D.f[ dM0M  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dM0M  ])[kbw ]
			(D.f[ dP0M  ])[kb  ] = mfabc;//(D.f[ dP0M  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dP0M  ])[kb  ]
			(D.f[ dM0P  ])[kw  ] = mfcba;//(D.f[ dM0P  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dM0P  ])[kw  ]
			(D.f[ d0PP  ])[k   ] = mfbaa;//(D.f[ d0PP  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ d0PP  ])[k   ]
			(D.f[ d0MM  ])[kbs ] = mfbcc;//(D.f[ d0MM  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ d0MM  ])[kbs ]
			(D.f[ d0PM  ])[kb  ] = mfbac;//(D.f[ d0PM  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ d0PM  ])[kb  ]
			(D.f[ d0MP  ])[ks  ] = mfbca;//(D.f[ d0MP  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ d0MP  ])[ks  ]
			(D.f[ d000])[k   ] = mfbbb;//(D.f[ d000])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ d000])[k   ]
			(D.f[ dPPP ])[k   ] = mfaaa;//(D.f[ dPPP ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dPPP ])[k   ]
			(D.f[ dPMP ])[ks  ] = mfaca;//(D.f[ dPMP ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dPMP ])[ks  ]
			(D.f[ dPPM ])[kb  ] = mfaac;//(D.f[ dPPM ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dPPM ])[kb  ]
			(D.f[ dPMM ])[kbs ] = mfacc;//(D.f[ dPMM ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dPMM ])[kbs ]
			(D.f[ dMPP ])[kw  ] = mfcaa;//(D.f[ dMPP ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dMPP ])[kw  ]
			(D.f[ dMMP ])[ksw ] = mfcca;//(D.f[ dMMP ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dMMP ])[ksw ]
			(D.f[ dMPM ])[kbw ] = mfcac;//(D.f[ dMPM ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dMPM ])[kbw ]
			(D.f[ dMMM ])[kbsw] = mfccc;//(D.f[ dMMM ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dMMM ])[kbsw]
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
__global__ void LB_Kernel_Kum_IsoTest_SP_27( real omega,
														unsigned int* bcMatD,
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														real* DDStart,
														real* dxxUx,
														real* dyyUy,
														real* dzzUz,
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

	if(k<numberOfLBnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			real mfcbb = (D.f[dP00])[k  ];//[ke   ];// +  c2over27 ;(D.f[dP00])[k  ];//ke
			real mfabb = (D.f[dM00])[kw ];//[kw   ];// +  c2over27 ;(D.f[dM00])[kw ];
			real mfbcb = (D.f[d0P0])[k  ];//[kn   ];// +  c2over27 ;(D.f[d0P0])[k  ];//kn
			real mfbab = (D.f[d0M0])[ks ];//[ks   ];// +  c2over27 ;(D.f[d0M0])[ks ];
			real mfbbc = (D.f[d00P])[k  ];//[kt   ];// +  c2over27 ;(D.f[d00P])[k  ];//kt
			real mfbba = (D.f[d00M])[kb ];//[kb   ];// +  c2over27 ;(D.f[d00M])[kb ];
			real mfccb = (D.f[dPP0])[k  ];//[kne  ];// +  c1over54 ;(D.f[dPP0])[k  ];//kne
			real mfaab = (D.f[dMM0])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dMM0])[ksw];
			real mfcab = (D.f[dPM0])[ks ];//[kse  ];// +  c1over54 ;(D.f[dPM0])[ks ];//kse
			real mfacb = (D.f[dMP0])[kw ];//[knw  ];// +  c1over54 ;(D.f[dMP0])[kw ];//knw
			real mfcbc = (D.f[dP0P])[k  ];//[kte  ];// +  c1over54 ;(D.f[dP0P])[k  ];//kte
			real mfaba = (D.f[dM0M])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dM0M])[kbw];
			real mfcba = (D.f[dP0M])[kb ];//[kbe  ];// +  c1over54 ;(D.f[dP0M])[kb ];//kbe
			real mfabc = (D.f[dM0P])[kw ];//[ktw  ];// +  c1over54 ;(D.f[dM0P])[kw ];//ktw
			real mfbcc = (D.f[d0PP])[k  ];//[ktn  ];// +  c1over54 ;(D.f[d0PP])[k  ];//ktn
			real mfbaa = (D.f[d0MM])[kbs];//[kbs  ];// +  c1over54 ;(D.f[d0MM])[kbs];
			real mfbca = (D.f[d0PM])[kb ];//[kbn  ];// +  c1over54 ;(D.f[d0PM])[kb ];//kbn
			real mfbac = (D.f[d0MP])[ks ];//[kts  ];// +  c1over54 ;(D.f[d0MP])[ks ];//kts
			real mfbbb = (D.f[d000])[k  ];//[kzero];// +  c8over27 ;(D.f[d000])[k  ];//kzero
			real mfccc = (D.f[dPPP])[k  ];//[ktne ];// +  c1over216;(D.f[dPPP])[k  ];//ktne
			real mfaac = (D.f[dMMP])[ksw];//[ktsw ];// +  c1over216;(D.f[dMMP])[ksw];//ktsw
			real mfcac = (D.f[dPMP])[ks ];//[ktse ];// +  c1over216;(D.f[dPMP])[ks ];//ktse
			real mfacc = (D.f[dMPP])[kw ];//[ktnw ];// +  c1over216;(D.f[dMPP])[kw ];//ktnw
			real mfcca = (D.f[dPPM])[kb ];//[kbne ];// +  c1over216;(D.f[dPPM])[kb ];//kbne
			real mfaaa = (D.f[dMMM])[kbsw];//[kbsw ];// +  c1over216;(D.f[dMMM])[kbsw];
			real mfcaa = (D.f[dPMM])[kbs];//[kbse ];// +  c1over216;(D.f[dPMM])[kbs];//kbse
			real mfaca = (D.f[dMPM])[kbw];//[kbnw ];// +  c1over216;(D.f[dMPM])[kbw];//kbnw
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
			real oMdrho = c1o1 - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
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
			mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - c2o1*	m1 * vvz + vz2 * m0;
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
			mfaca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += c1o1* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - c2o1*	m1 * vvx + vx2 * m0;
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
			real OxxPyyPzz = c1o1;
			real OxyyPxzz  = c1o1;//two-omega;//eight*(two-omega)/(eight -omega);//one;//omega;//two-omega;//
			real OxyyMxzz  = c1o1;//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
			real O4        = c1o1;
			real O5        = c1o1;
			real O6        = c1o1;

			//Cum 4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + c2o1 * mfbba * mfbab);
			real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + c2o1 * mfbba * mfabb);
			real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + c2o1 * mfbab * mfabb);

			real CUMcca = mfcca - (mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;
			real CUMcac = mfcac - (mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;
			real CUMacc = mfacc - (mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;

			//Cum 5.
			real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
			real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
			real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

			//Cum 6.
			real CUMccc = mfccc  +((-c4o1 *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
							+(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     c2o1 * (mfcaa * mfaca * mfaac)
							+ c16o1 *  mfbba * mfbab * mfabb)
							-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
							-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(c1o1-c2o1* oMdrho)- c1o27* oMdrho * oMdrho*(-c2o1* oMdrho)
							+(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) +c1o27*oMdrho;





			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Iso Test Part 1
			real dxuydyux = -c3o1 * omega * mfbba;
			real dxuzdzux = -c3o1 * omega * mfbab;
			real dyuzdzuy = -c3o1 * omega * mfabb;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			{
				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
				real dyuy = dxux + omega * c3o2 * mxxMyy;
				real dzuz = dxux + omega * c3o2 * mxxMzz;

				//relax
				mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
				mxxMyy    += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
				mxxMzz    += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

			}

			//mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);
			//mxxMyy    += -(-omega) * (-mxxMyy);
			//mxxMzz    += -(-omega) * (-mxxMzz);
			mfabb     += omega * (-mfabb);
			mfbab     += omega * (-mfbab);
			mfbba     += omega * (-mfbba);

			// linear combinations back
			mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-c2o1*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - c2o1* mxxMzz + mxxPyyPzz);

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
			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
			mfbbb     += wadjust * (-mfbbb);
			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
			mxxyPyzz  += wadjust * (-mxxyPyzz);
			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
			mxxyMyzz  += wadjust * (-mxxyMyzz);
			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
			mxxzPyyz  += wadjust * (-mxxzPyyz);
			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
			mxxzMyyz  += wadjust * (-mxxzMyyz);
			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
			mxyyPxzz  += wadjust * (-mxyyPxzz);
			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
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
			CUMbcc = c2o3 * (c1o1 - c1o2 * O5) * (vvy * dxuydyux + vvz * dxuzdzux) * c0o1 + (c1o1 - O5) * CUMbcc;
			CUMcbc = c2o3 * (c1o1 - c1o2 * O5) * (vvx * dxuydyux + vvz * dyuzdzuy) * c0o1 + (c1o1 - O5) * CUMcbc;
			CUMccb = c2o3 * (c1o1 - c1o2 * O5) * (vvx * dxuzdzux + vvy * dyuzdzuy) * c0o1 + (c1o1 - O5) * CUMccb;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			//6.
			//CUMccc += O6 * (-CUMccc);
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Iso Test Part 4
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//new calculation of 6. moment
			CUMccc = O6 * c1o3 * (vx2 + vy2 + vz2) * c0o1 + (c1o1 - O6) * CUMccc;
			// second derivation of ux, uy, uz
			real dxxux = c3o1 * (c3o2 * (CUMbcc - preCUMbcc) - c0o1 * (mxyyPxzz - premxyyPxzz));
			real dyyuy = c3o1 * (c3o2 * (CUMcbc - preCUMcbc) - c0o1 * (mxxyPyzz - premxxyPyzz));
			real dzzuz = c3o1 * (c3o2 * (CUMccb - preCUMccb) - c0o1 * (mxxzPyyz - premxxzPyyz));
			// copy local values to global arrays for paraview files 
			dxxUx[k] = dxxux;
			dyyUy[k] = dyyuy;
			dzzUz[k] = dzzuz;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			//back cumulants to central moments
			//4.
			mfcbb = CUMcbb + ((mfcaa + c1o3 * oMdrho) * mfabb + c2o1 * mfbba * mfbab);
			mfbcb = CUMbcb + ((mfaca + c1o3 * oMdrho) * mfbab + c2o1 * mfbba * mfabb);
			mfbbc = CUMbbc + ((mfaac + c1o3 * oMdrho) * mfbba + c2o1 * mfbab * mfabb); 
						   
			mfcca = CUMcca + (mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;
			mfcac = CUMcac + (mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;
			mfacc = CUMacc + (mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;

			//5.
			mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) * oMdrho;
			mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) * oMdrho;
			mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) * oMdrho;
			
			//6.
			mfccc = CUMccc  -(( -c4o1 *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
							+(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     c2o1 * (mfcaa * mfaca * mfaac)
							+ c16o1 *  mfbba * mfbab * mfabb)
							-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
							-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(c1o1-c2o1* oMdrho)- c1o27* oMdrho * oMdrho*(-c2o1* oMdrho)
							+(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) -c1o27*oMdrho;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			//back
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + c1o1* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - c2o1* mfaab *  vvz         +  mfaaa                * (c1o1- vz2)              - c1o1* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + c1o1* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - c2o1* mfabb *  vvz         + mfaba * (c1o1- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - c2o1* mfacb *  vvz         +  mfaca                  * (c1o1- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - c2o1* mfbab *  vvz         + mfbaa * (c1o1- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - c2o1* mfbbb *  vvz         + mfbba * (c1o1- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - c2o1* mfbcb *  vvz         + mfbca * (c1o1- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - c2o1* mfcab *  vvz         +  mfcaa                  * (c1o1- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - c2o1* mfcbb *  vvz         + mfcba * (c1o1- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - c2o1* mfccb *  vvz         +  mfcca                  * (c1o1- vz2)              - c1o9 * oMdrho * vz2; 
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
			m1 = -mfaca        - c2o1* mfaba *  vvy         +  mfaaa                  * (c1o1- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - c2o1* mfabb *  vvy         +  mfaab                  * (c1o1- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - c2o1* mfabc *  vvy         +  mfaac                  * (c1o1- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - c2o1* mfbba *  vvy         + mfbaa * (c1o1- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - c2o1* mfbbb *  vvy         + mfbab * (c1o1- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - c2o1* mfbbc *  vvy         + mfbac * (c1o1- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - c2o1* mfcba *  vvy         +  mfcaa                   * (c1o1- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - c2o1* mfcbb *  vvy         +  mfcab                  * (c1o1- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - c2o1* mfcbc *  vvy         +  mfcac                   * (c1o1- vy2)              - c1o18 * oMdrho * vy2; 
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
			m1 = -mfcaa        - c2o1* mfbaa *  vvx         +  mfaaa                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - c2o1* mfbba *  vvx         +  mfaba                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - c2o1* mfbca *  vvx         +  mfaca                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - c2o1* mfbab *  vvx         +  mfaab                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - c2o1* mfbbb *  vvx         +  mfabb                  * (c1o1- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - c2o1* mfbcb *  vvx         +  mfacb                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - c2o1* mfbac *  vvx         +  mfaac                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - c2o1* mfbbc *  vvx         +  mfabc                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - c2o1* mfbcc *  vvx         +  mfacc                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dP00   ])[k   ] = mfabb;//(D.f[ dP00   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dP00   ])[k   ]                                                                     
			(D.f[ dM00   ])[kw  ] = mfcbb;//(D.f[ dM00   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dM00   ])[kw  ]                                                                   
			(D.f[ d0P0   ])[k   ] = mfbab;//(D.f[ d0P0   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ d0P0   ])[k   ]
			(D.f[ d0M0   ])[ks  ] = mfbcb;//(D.f[ d0M0   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ d0M0   ])[ks  ]
			(D.f[ d00P   ])[k   ] = mfbba;//(D.f[ d00P   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ d00P   ])[k   ]
			(D.f[ d00M   ])[kb  ] = mfbbc;//(D.f[ d00M   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ d00M   ])[kb  ]
			(D.f[ dPP0  ])[k   ] = mfaab;//(D.f[ dPP0  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dPP0  ])[k   ]
			(D.f[ dMM0  ])[ksw ] = mfccb;//(D.f[ dMM0  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dMM0  ])[ksw ]
			(D.f[ dPM0  ])[ks  ] = mfacb;//(D.f[ dPM0  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dPM0  ])[ks  ]
			(D.f[ dMP0  ])[kw  ] = mfcab;//(D.f[ dMP0  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dMP0  ])[kw  ]
			(D.f[ dP0P  ])[k   ] = mfaba;//(D.f[ dP0P  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dP0P  ])[k   ]
			(D.f[ dM0M  ])[kbw ] = mfcbc;//(D.f[ dM0M  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dM0M  ])[kbw ]
			(D.f[ dP0M  ])[kb  ] = mfabc;//(D.f[ dP0M  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dP0M  ])[kb  ]
			(D.f[ dM0P  ])[kw  ] = mfcba;//(D.f[ dM0P  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dM0P  ])[kw  ]
			(D.f[ d0PP  ])[k   ] = mfbaa;//(D.f[ d0PP  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ d0PP  ])[k   ]
			(D.f[ d0MM  ])[kbs ] = mfbcc;//(D.f[ d0MM  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ d0MM  ])[kbs ]
			(D.f[ d0PM  ])[kb  ] = mfbac;//(D.f[ d0PM  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ d0PM  ])[kb  ]
			(D.f[ d0MP  ])[ks  ] = mfbca;//(D.f[ d0MP  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ d0MP  ])[ks  ]
			(D.f[ d000])[k   ] = mfbbb;//(D.f[ d000])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ d000])[k   ]
			(D.f[ dPPP ])[k   ] = mfaaa;//(D.f[ dPPP ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dPPP ])[k   ]
			(D.f[ dPMP ])[ks  ] = mfaca;//(D.f[ dPMP ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dPMP ])[ks  ]
			(D.f[ dPPM ])[kb  ] = mfaac;//(D.f[ dPPM ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dPPM ])[kb  ]
			(D.f[ dPMM ])[kbs ] = mfacc;//(D.f[ dPMM ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dPMM ])[kbs ]
			(D.f[ dMPP ])[kw  ] = mfcaa;//(D.f[ dMPP ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dMPP ])[kw  ]
			(D.f[ dMMP ])[ksw ] = mfcca;//(D.f[ dMMP ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dMMP ])[ksw ]
			(D.f[ dMPM ])[kbw ] = mfcac;//(D.f[ dMPM ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dMPM ])[kbw ]
			(D.f[ dMMM ])[kbsw] = mfccc;//(D.f[ dMMM ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dMMM ])[kbsw]
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
__global__ void LB_Kernel_Kum_1h_SP_27(  real omega,
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

	if(k<numberOfLBnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			real mfcbb = (D.f[dP00])[k  ];//[ke   ];// +  c2over27 ;(D.f[dP00])[k  ];//ke
			real mfabb = (D.f[dM00])[kw ];//[kw   ];// +  c2over27 ;(D.f[dM00])[kw ];
			real mfbcb = (D.f[d0P0])[k  ];//[kn   ];// +  c2over27 ;(D.f[d0P0])[k  ];//kn
			real mfbab = (D.f[d0M0])[ks ];//[ks   ];// +  c2over27 ;(D.f[d0M0])[ks ];
			real mfbbc = (D.f[d00P])[k  ];//[kt   ];// +  c2over27 ;(D.f[d00P])[k  ];//kt
			real mfbba = (D.f[d00M])[kb ];//[kb   ];// +  c2over27 ;(D.f[d00M])[kb ];
			real mfccb = (D.f[dPP0])[k  ];//[kne  ];// +  c1over54 ;(D.f[dPP0])[k  ];//kne
			real mfaab = (D.f[dMM0])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dMM0])[ksw];
			real mfcab = (D.f[dPM0])[ks ];//[kse  ];// +  c1over54 ;(D.f[dPM0])[ks ];//kse
			real mfacb = (D.f[dMP0])[kw ];//[knw  ];// +  c1over54 ;(D.f[dMP0])[kw ];//knw
			real mfcbc = (D.f[dP0P])[k  ];//[kte  ];// +  c1over54 ;(D.f[dP0P])[k  ];//kte
			real mfaba = (D.f[dM0M])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dM0M])[kbw];
			real mfcba = (D.f[dP0M])[kb ];//[kbe  ];// +  c1over54 ;(D.f[dP0M])[kb ];//kbe
			real mfabc = (D.f[dM0P])[kw ];//[ktw  ];// +  c1over54 ;(D.f[dM0P])[kw ];//ktw
			real mfbcc = (D.f[d0PP])[k  ];//[ktn  ];// +  c1over54 ;(D.f[d0PP])[k  ];//ktn
			real mfbaa = (D.f[d0MM])[kbs];//[kbs  ];// +  c1over54 ;(D.f[d0MM])[kbs];
			real mfbca = (D.f[d0PM])[kb ];//[kbn  ];// +  c1over54 ;(D.f[d0PM])[kb ];//kbn
			real mfbac = (D.f[d0MP])[ks ];//[kts  ];// +  c1over54 ;(D.f[d0MP])[ks ];//kts
			real mfbbb = (D.f[d000])[k  ];//[kzero];// +  c8over27 ;(D.f[d000])[k  ];//kzero
			real mfccc = (D.f[dPPP])[k  ];//[ktne ];// +  c1over216;(D.f[dPPP])[k  ];//ktne
			real mfaac = (D.f[dMMP])[ksw];//[ktsw ];// +  c1over216;(D.f[dMMP])[ksw];//ktsw
			real mfcac = (D.f[dPMP])[ks ];//[ktse ];// +  c1over216;(D.f[dPMP])[ks ];//ktse
			real mfacc = (D.f[dMPP])[kw ];//[ktnw ];// +  c1over216;(D.f[dMPP])[kw ];//ktnw
			real mfcca = (D.f[dPPM])[kb ];//[kbne ];// +  c1over216;(D.f[dPPM])[kb ];//kbne
			real mfaaa = (D.f[dMMM])[kbsw];//[kbsw ];// +  c1over216;(D.f[dMMM])[kbsw];
			real mfcaa = (D.f[dPMM])[kbs];//[kbse ];// +  c1over216;(D.f[dPMM])[kbs];//kbse
			real mfaca = (D.f[dMPM])[kbw];//[kbnw ];// +  c1over216;(D.f[dMPM])[kbw];//kbnw
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
			real oMdrho = c1o1 - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
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
			mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - c2o1*	m1 * vvz + vz2 * m0;
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
			mfaca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += c1o1* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - c2o1*	m1 * vvx + vx2 * m0;
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
			real OxxPyyPzz = c1o1;
			real OxyyPxzz  = c1o1;//omega;//two-omega;//
			real OxyyMxzz  = c2o1-omega;//one;// 
			real O4        = c1o1;
			real O5        = c1o1;
			real O6        = c1o1;

			////Cum 4.
			//real CUMcbb;	real CUMbcb;	real CUMbbc;
			//real CUMcca;	real CUMcac;	real CUMacc;
			////Cum 5.
			//real CUMbcc;	real CUMcbc;	real CUMccb;
			////Cum 6.
			//real CUMccc;

			//Cum four
			real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + c2o1 * mfbba * mfbab);
			real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + c2o1 * mfbba * mfabb);
			real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + c2o1 * mfbab * mfabb);

			real CUMcca = mfcca - (mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;
			real CUMcac = mfcac - (mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;
			real CUMacc = mfacc - (mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;

			//Cum 5.
			//real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) //O(eps^5) 
			//				- c1o3 * (mfbca + mfbac); //O(eps^3)
			//real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) //O(eps^5)
			//				- c1o3 * (mfcba + mfabc); //O(eps^3)
			//real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) //O(eps^5)
			//				- c1o3 * (mfacb + mfcab);//O(eps^3)

			real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
			real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
			real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

			//Cum 6.
			real CUMccc = mfccc  +((-c4o1 *  mfbbb * mfbbb  
				-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				-    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				-     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
				+(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
				+     c2o1 * (mfcaa * mfaca * mfaac)
				+ c16o1 *  mfbba * mfbab * mfabb)
				-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
				-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(c1o1-c2o1* oMdrho)- c1o27* oMdrho * oMdrho*(-c2o1* oMdrho)
				+(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
				+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) +c1o27*oMdrho;





			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;

			{
				real dxux = c0o1;//c1o2 * (-omega) *(mxxMyy + mxxMzz) +  OxxPyyPzz * (mfaaa - mxxPyyPzz);
				real dyuy = c0o1;//dxux + omega * c3o2 * mxxMyy;
				real dzuz = c0o1;//dxux + omega * c3o2 * mxxMzz;

				//relax
				mxxPyyPzz += (OxxPyyPzz)*(mfaaa-c2o1*(ux*vvx+uy*vvy)-ux*ux-uy*uy  - mxxPyyPzz)- c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
				mxxMyy    += omega * (-c2o1*(ux*vvx-uy*vvy)-(ux*ux-uy*uy)-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux + vy2 * dyuy);
				mxxMzz    += omega * (-c2o1*ux*vvx -ux*ux                -mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux + vz2 * dzuz);
			}
			mfabb     += omega * (-uy*vvz-mfabb);
			mfbab     += omega * (-ux*vvz-mfbab);
			mfbba     += omega * (-ux*vvy-uy*vvx-ux*uy-mfbba);

			// linear combinations back
			mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-c2o1*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - c2o1* mxxMzz + mxxPyyPzz);

			//3.
			// linear combinations

			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

			//relax
			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
			mfbbb     += wadjust * (-mfbbb);
			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
			mxxyPyzz  += wadjust * (-mxxyPyzz);
			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
			mxxyMyzz  += wadjust * (-mxxyMyzz);
			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
			mxxzPyyz  += wadjust * (-mxxzPyyz);
			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
			mxxzMyyz  += wadjust * (-mxxzMyyz);
			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
			mxyyPxzz  += wadjust * (-mxyyPxzz);
			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
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
			mfcbb = CUMcbb + ((mfcaa + c1o3 * oMdrho) * mfabb + c2o1 * mfbba * mfbab);
			mfbcb = CUMbcb + ((mfaca + c1o3 * oMdrho) * mfbab + c2o1 * mfbba * mfabb);
			mfbbc = CUMbbc + ((mfaac + c1o3 * oMdrho) * mfbba + c2o1 * mfbab * mfabb); 

			mfcca = CUMcca + (mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;
			mfcac = CUMcac + (mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;
			mfacc = CUMacc + (mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;

			//5.
			mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) * oMdrho;
			mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) * oMdrho;
			mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) * oMdrho;

			//6.
			mfccc = CUMccc  -(( -c4o1 *  mfbbb * mfbbb  
				-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				-    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				-     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
				+(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
				+     c2o1 * (mfcaa * mfaca * mfaac)
				+ c16o1 *  mfbba * mfbab * mfabb)
				-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
				-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(c1o1-c2o1* oMdrho)- c1o27* oMdrho * oMdrho*(-c2o1* oMdrho)
				+(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
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
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + c1o1* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - c2o1* mfaab *  vvz         +  mfaaa                * (c1o1- vz2)              - c1o1* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + c1o1* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - c2o1* mfabb *  vvz         + mfaba * (c1o1- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - c2o1* mfacb *  vvz         +  mfaca                  * (c1o1- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - c2o1* mfbab *  vvz         + mfbaa * (c1o1- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - c2o1* mfbbb *  vvz         + mfbba * (c1o1- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - c2o1* mfbcb *  vvz         + mfbca * (c1o1- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - c2o1* mfcab *  vvz         +  mfcaa                  * (c1o1- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - c2o1* mfcbb *  vvz         + mfcba * (c1o1- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - c2o1* mfccb *  vvz         +  mfcca                  * (c1o1- vz2)              - c1o9 * oMdrho * vz2; 
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
			m1 = -mfaca        - c2o1* mfaba *  vvy         +  mfaaa                  * (c1o1- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - c2o1* mfabb *  vvy         +  mfaab                  * (c1o1- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - c2o1* mfabc *  vvy         +  mfaac                  * (c1o1- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - c2o1* mfbba *  vvy         + mfbaa * (c1o1- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - c2o1* mfbbb *  vvy         + mfbab * (c1o1- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - c2o1* mfbbc *  vvy         + mfbac * (c1o1- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - c2o1* mfcba *  vvy         +  mfcaa                   * (c1o1- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - c2o1* mfcbb *  vvy         +  mfcab                  * (c1o1- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - c2o1* mfcbc *  vvy         +  mfcac                   * (c1o1- vy2)              - c1o18 * oMdrho * vy2; 
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
			m1 = -mfcaa        - c2o1* mfbaa *  vvx         +  mfaaa                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - c2o1* mfbba *  vvx         +  mfaba                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - c2o1* mfbca *  vvx         +  mfaca                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - c2o1* mfbab *  vvx         +  mfaab                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - c2o1* mfbbb *  vvx         +  mfabb                  * (c1o1- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - c2o1* mfbcb *  vvx         +  mfacb                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - c2o1* mfbac *  vvx         +  mfaac                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - c2o1* mfbbc *  vvx         +  mfabc                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - c2o1* mfbcc *  vvx         +  mfacc                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dP00   ])[k   ] = mfabb;//(D.f[ dP00   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dP00   ])[k   ]                                                                     
			(D.f[ dM00   ])[kw  ] = mfcbb;//(D.f[ dM00   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dM00   ])[kw  ]                                                                   
			(D.f[ d0P0   ])[k   ] = mfbab;//(D.f[ d0P0   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ d0P0   ])[k   ]
			(D.f[ d0M0   ])[ks  ] = mfbcb;//(D.f[ d0M0   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ d0M0   ])[ks  ]
			(D.f[ d00P   ])[k   ] = mfbba;//(D.f[ d00P   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ d00P   ])[k   ]
			(D.f[ d00M   ])[kb  ] = mfbbc;//(D.f[ d00M   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ d00M   ])[kb  ]
			(D.f[ dPP0  ])[k   ] = mfaab;//(D.f[ dPP0  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dPP0  ])[k   ]
			(D.f[ dMM0  ])[ksw ] = mfccb;//(D.f[ dMM0  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dMM0  ])[ksw ]
			(D.f[ dPM0  ])[ks  ] = mfacb;//(D.f[ dPM0  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dPM0  ])[ks  ]
			(D.f[ dMP0  ])[kw  ] = mfcab;//(D.f[ dMP0  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dMP0  ])[kw  ]
			(D.f[ dP0P  ])[k   ] = mfaba;//(D.f[ dP0P  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dP0P  ])[k   ]
			(D.f[ dM0M  ])[kbw ] = mfcbc;//(D.f[ dM0M  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dM0M  ])[kbw ]
			(D.f[ dP0M  ])[kb  ] = mfabc;//(D.f[ dP0M  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dP0M  ])[kb  ]
			(D.f[ dM0P  ])[kw  ] = mfcba;//(D.f[ dM0P  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dM0P  ])[kw  ]
			(D.f[ d0PP  ])[k   ] = mfbaa;//(D.f[ d0PP  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ d0PP  ])[k   ]
			(D.f[ d0MM  ])[kbs ] = mfbcc;//(D.f[ d0MM  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ d0MM  ])[kbs ]
			(D.f[ d0PM  ])[kb  ] = mfbac;//(D.f[ d0PM  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ d0PM  ])[kb  ]
			(D.f[ d0MP  ])[ks  ] = mfbca;//(D.f[ d0MP  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ d0MP  ])[ks  ]
			(D.f[ d000])[k   ] = mfbbb;//(D.f[ d000])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ d000])[k   ]
			(D.f[ dPPP ])[k   ] = mfaaa;//(D.f[ dPPP ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dPPP ])[k   ]
			(D.f[ dPMP ])[ks  ] = mfaca;//(D.f[ dPMP ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dPMP ])[ks  ]
			(D.f[ dPPM ])[kb  ] = mfaac;//(D.f[ dPPM ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dPPM ])[kb  ]
			(D.f[ dPMM ])[kbs ] = mfacc;//(D.f[ dPMM ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dPMM ])[kbs ]
			(D.f[ dMPP ])[kw  ] = mfcaa;//(D.f[ dMPP ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dMPP ])[kw  ]
			(D.f[ dMMP ])[ksw ] = mfcca;//(D.f[ dMMP ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dMMP ])[ksw ]
			(D.f[ dMPM ])[kbw ] = mfcac;//(D.f[ dMPM ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dMPM ])[kbw ]
			(D.f[ dMMM ])[kbsw] = mfccc;//(D.f[ dMMM ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dMMM ])[kbsw]
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
__global__ void LB_Kernel_Kum_New_SP_27(     real omega,
														unsigned int* bcMatD,
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														real* DDStart,
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

	if(k<numberOfLBnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			real mfcbb = (D.f[dP00])[k  ];//[ke   ];// +  c2over27 ;(D.f[dP00])[k  ];//ke
			real mfabb = (D.f[dM00])[kw ];//[kw   ];// +  c2over27 ;(D.f[dM00])[kw ];
			real mfbcb = (D.f[d0P0])[k  ];//[kn   ];// +  c2over27 ;(D.f[d0P0])[k  ];//kn
			real mfbab = (D.f[d0M0])[ks ];//[ks   ];// +  c2over27 ;(D.f[d0M0])[ks ];
			real mfbbc = (D.f[d00P])[k  ];//[kt   ];// +  c2over27 ;(D.f[d00P])[k  ];//kt
			real mfbba = (D.f[d00M])[kb ];//[kb   ];// +  c2over27 ;(D.f[d00M])[kb ];
			real mfccb = (D.f[dPP0])[k  ];//[kne  ];// +  c1over54 ;(D.f[dPP0])[k  ];//kne
			real mfaab = (D.f[dMM0])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dMM0])[ksw];
			real mfcab = (D.f[dPM0])[ks ];//[kse  ];// +  c1over54 ;(D.f[dPM0])[ks ];//kse
			real mfacb = (D.f[dMP0])[kw ];//[knw  ];// +  c1over54 ;(D.f[dMP0])[kw ];//knw
			real mfcbc = (D.f[dP0P])[k  ];//[kte  ];// +  c1over54 ;(D.f[dP0P])[k  ];//kte
			real mfaba = (D.f[dM0M])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dM0M])[kbw];
			real mfcba = (D.f[dP0M])[kb ];//[kbe  ];// +  c1over54 ;(D.f[dP0M])[kb ];//kbe
			real mfabc = (D.f[dM0P])[kw ];//[ktw  ];// +  c1over54 ;(D.f[dM0P])[kw ];//ktw
			real mfbcc = (D.f[d0PP])[k  ];//[ktn  ];// +  c1over54 ;(D.f[d0PP])[k  ];//ktn
			real mfbaa = (D.f[d0MM])[kbs];//[kbs  ];// +  c1over54 ;(D.f[d0MM])[kbs];
			real mfbca = (D.f[d0PM])[kb ];//[kbn  ];// +  c1over54 ;(D.f[d0PM])[kb ];//kbn
			real mfbac = (D.f[d0MP])[ks ];//[kts  ];// +  c1over54 ;(D.f[d0MP])[ks ];//kts
			real mfbbb = (D.f[d000])[k  ];//[kzero];// +  c8over27 ;(D.f[d000])[k  ];//kzero
			real mfccc = (D.f[dPPP])[k  ];//[ktne ];// +  c1over216;(D.f[dPPP])[k  ];//ktne
			real mfaac = (D.f[dMMP])[ksw];//[ktsw ];// +  c1over216;(D.f[dMMP])[ksw];//ktsw
			real mfcac = (D.f[dPMP])[ks ];//[ktse ];// +  c1over216;(D.f[dPMP])[ks ];//ktse
			real mfacc = (D.f[dMPP])[kw ];//[ktnw ];// +  c1over216;(D.f[dMPP])[kw ];//ktnw
			real mfcca = (D.f[dPPM])[kb ];//[kbne ];// +  c1over216;(D.f[dPPM])[kb ];//kbne
			real mfaaa = (D.f[dMMM])[kbsw];//[kbsw ];// +  c1over216;(D.f[dMMM])[kbsw];
			real mfcaa = (D.f[dPMM])[kbs];//[kbse ];// +  c1over216;(D.f[dPMM])[kbs];//kbse
			real mfaca = (D.f[dMPM])[kbw];//[kbnw ];// +  c1over216;(D.f[dMPM])[kbw];//kbnw
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
			real oMdrho = c1o1 - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
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
			//real wadjust;
			//real qudricLimit = 0.01f;
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
			mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - c2o1*	m1 * vvz + vz2 * m0;
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
			mfaca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += c1o1* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - c2o1*	m1 * vvx + vx2 * m0;
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
			real OxxPyyPzz = c1o1;
			real OxyyPxzz  = c8o1*(c2o1-omega)/(c8o1 -omega);//one;//two-omega;//eight*(two-omega)/(eight -omega);//one;//omega;//two-omega;//
			real OxyyMxzz  = c8o1*(c2o1-omega)/(c8o1 -omega);//one;//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
			real O4        = c1o1;
			real O5        = c1o1;
			real O6        = c1o1;

			////Cum 4.
			//real CUMcbb;	real CUMbcb;	real CUMbbc;
			//real CUMcca;	real CUMcac;	real CUMacc;
			////Cum 5.
			//real CUMbcc;	real CUMcbc;	real CUMccb;
			////Cum 6.
			//real CUMccc;

			//Cum 4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + c2o1 * mfbba * mfbab); // /rho
			real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + c2o1 * mfbba * mfabb); // /rho
			real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + c2o1 * mfbab * mfabb); // /rho

			real CUMcca = mfcca - ((mfcaa * mfaca + c2o1 * mfbba * mfbba) /* /rho*/ + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho);
			real CUMcac = mfcac - ((mfcaa * mfaac + c2o1 * mfbab * mfbab) /* /rho*/ + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho);
			real CUMacc = mfacc - ((mfaac * mfaca + c2o1 * mfabb * mfabb) /* /rho*/ + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho);

			//Cum 5.
			//real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) //O(eps^5) 
			//				- c1o3 * (mfbca + mfbac); //O(eps^3)
			//real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) //O(eps^5)
			//				- c1o3 * (mfcba + mfabc); //O(eps^3)
			//real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) //O(eps^5)
			//				- c1o3 * (mfacb + mfcab);//O(eps^3)

			real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
			real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
			real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

			//Cum 6.
			real CUMccc = mfccc  +((-c4o1 *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
							+(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     c2o1 * (mfcaa * mfaca * mfaac)
							+ c16o1 *  mfbba * mfbab * mfabb)
							-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
							-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(c1o1-c2o1* oMdrho)- c1o27* oMdrho * oMdrho*(-c2o1* oMdrho)
							+(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) +c1o27*oMdrho;





			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;
			
			//////////////////////////////////////////////////////////////////////////
			//real magicBulk=(CUMacc+CUMcac+CUMcca)*(c1o1/OxxPyyPzz-c1o2)*c3o2*8.;

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
			mfaca = c1o3 * (-c2o1*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - c2o1* mxxMzz + mxxPyyPzz);

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
			mfcbb = CUMcbb + ((mfcaa + c1o3 * oMdrho) * mfabb + c2o1 * mfbba * mfbab);
			mfbcb = CUMbcb + ((mfaca + c1o3 * oMdrho) * mfbab + c2o1 * mfbba * mfabb);
			mfbbc = CUMbbc + ((mfaac + c1o3 * oMdrho) * mfbba + c2o1 * mfbab * mfabb); 
						   
			mfcca = CUMcca + (mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;
			mfcac = CUMcac + (mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;
			mfacc = CUMacc + (mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho;

			//5.
			mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) * oMdrho;
			mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) * oMdrho;
			mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) * oMdrho;
			
			//6.
			mfccc = CUMccc  -(( -c4o1 *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
							+(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     c2o1 * (mfcaa * mfaca * mfaac)
							+ c16o1 *  mfbba * mfbab * mfabb)
							-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
							-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(c1o1-c2o1* oMdrho)- c1o27* oMdrho * oMdrho*(-c2o1* oMdrho)
							+(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) -c1o27*oMdrho;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			//back
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + c1o1* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - c2o1* mfaab *  vvz         +  mfaaa                * (c1o1- vz2)              - c1o1* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + c1o1* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - c2o1* mfabb *  vvz         + mfaba * (c1o1- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - c2o1* mfacb *  vvz         +  mfaca                  * (c1o1- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - c2o1* mfbab *  vvz         + mfbaa * (c1o1- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - c2o1* mfbbb *  vvz         + mfbba * (c1o1- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - c2o1* mfbcb *  vvz         + mfbca * (c1o1- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - c2o1* mfcab *  vvz         +  mfcaa                  * (c1o1- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - c2o1* mfcbb *  vvz         + mfcba * (c1o1- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - c2o1* mfccb *  vvz         +  mfcca                  * (c1o1- vz2)              - c1o9 * oMdrho * vz2; 
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
			m1 = -mfaca        - c2o1* mfaba *  vvy         +  mfaaa                  * (c1o1- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - c2o1* mfabb *  vvy         +  mfaab                  * (c1o1- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - c2o1* mfabc *  vvy         +  mfaac                  * (c1o1- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - c2o1* mfbba *  vvy         + mfbaa * (c1o1- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - c2o1* mfbbb *  vvy         + mfbab * (c1o1- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - c2o1* mfbbc *  vvy         + mfbac * (c1o1- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - c2o1* mfcba *  vvy         +  mfcaa                   * (c1o1- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - c2o1* mfcbb *  vvy         +  mfcab                  * (c1o1- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - c2o1* mfcbc *  vvy         +  mfcac                   * (c1o1- vy2)              - c1o18 * oMdrho * vy2; 
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
			m1 = -mfcaa        - c2o1* mfbaa *  vvx         +  mfaaa                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - c2o1* mfbba *  vvx         +  mfaba                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - c2o1* mfbca *  vvx         +  mfaca                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - c2o1* mfbab *  vvx         +  mfaab                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - c2o1* mfbbb *  vvx         +  mfabb                  * (c1o1- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - c2o1* mfbcb *  vvx         +  mfacb                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - c2o1* mfbac *  vvx         +  mfaac                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - c2o1* mfbbc *  vvx         +  mfabc                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - c2o1* mfbcc *  vvx         +  mfacc                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dP00   ])[k   ] = mfabb;//(D.f[ dP00   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dP00   ])[k   ]                                                                     
			(D.f[ dM00   ])[kw  ] = mfcbb;//(D.f[ dM00   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dM00   ])[kw  ]                                                                   
			(D.f[ d0P0   ])[k   ] = mfbab;//(D.f[ d0P0   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ d0P0   ])[k   ]
			(D.f[ d0M0   ])[ks  ] = mfbcb;//(D.f[ d0M0   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ d0M0   ])[ks  ]
			(D.f[ d00P   ])[k   ] = mfbba;//(D.f[ d00P   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ d00P   ])[k   ]
			(D.f[ d00M   ])[kb  ] = mfbbc;//(D.f[ d00M   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ d00M   ])[kb  ]
			(D.f[ dPP0  ])[k   ] = mfaab;//(D.f[ dPP0  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ dPP0  ])[k   ]
			(D.f[ dMM0  ])[ksw ] = mfccb;//(D.f[ dMM0  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ dMM0  ])[ksw ]
			(D.f[ dPM0  ])[ks  ] = mfacb;//(D.f[ dPM0  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ dPM0  ])[ks  ]
			(D.f[ dMP0  ])[kw  ] = mfcab;//(D.f[ dMP0  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ dMP0  ])[kw  ]
			(D.f[ dP0P  ])[k   ] = mfaba;//(D.f[ dP0P  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ dP0P  ])[k   ]
			(D.f[ dM0M  ])[kbw ] = mfcbc;//(D.f[ dM0M  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ dM0M  ])[kbw ]
			(D.f[ dP0M  ])[kb  ] = mfabc;//(D.f[ dP0M  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ dP0M  ])[kb  ]
			(D.f[ dM0P  ])[kw  ] = mfcba;//(D.f[ dM0P  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ dM0P  ])[kw  ]
			(D.f[ d0PP  ])[k   ] = mfbaa;//(D.f[ d0PP  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ d0PP  ])[k   ]
			(D.f[ d0MM  ])[kbs ] = mfbcc;//(D.f[ d0MM  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ d0MM  ])[kbs ]
			(D.f[ d0PM  ])[kb  ] = mfbac;//(D.f[ d0PM  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ d0PM  ])[kb  ]
			(D.f[ d0MP  ])[ks  ] = mfbca;//(D.f[ d0MP  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ d0MP  ])[ks  ]
			(D.f[ d000])[k   ] = mfbbb;//(D.f[ d000])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ d000])[k   ]
			(D.f[ dPPP ])[k   ] = mfaaa;//(D.f[ dPPP ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ dPPP ])[k   ]
			(D.f[ dPMP ])[ks  ] = mfaca;//(D.f[ dPMP ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ dPMP ])[ks  ]
			(D.f[ dPPM ])[kb  ] = mfaac;//(D.f[ dPPM ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ dPPM ])[kb  ]
			(D.f[ dPMM ])[kbs ] = mfacc;//(D.f[ dPMM ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ dPMM ])[kbs ]
			(D.f[ dMPP ])[kw  ] = mfcaa;//(D.f[ dMPP ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ dMPP ])[kw  ]
			(D.f[ dMMP ])[ksw ] = mfcca;//(D.f[ dMMP ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ dMMP ])[ksw ]
			(D.f[ dMPM ])[kbw ] = mfcac;//(D.f[ dMPM ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ dMPM ])[kbw ]
			(D.f[ dMMM ])[kbsw] = mfccc;//(D.f[ dMMM ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ dMMM ])[kbsw]
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
__global__ void LB_Kernel_Kum_Comp_SP_27(    real omega,
														unsigned int* bcMatD,
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														real* DDStart,
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

	if(k<numberOfLBnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
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
			real f_E     = (D.f[dP00])[ke   ];// +  c2over27 ;
			real f_W     = (D.f[dM00])[kw   ];// +  c2over27 ;
			real f_N     = (D.f[d0P0])[kn   ];// +  c2over27 ;
			real f_S     = (D.f[d0M0])[ks   ];// +  c2over27 ;
			real f_T     = (D.f[d00P])[kt   ];// +  c2over27 ;
			real f_B     = (D.f[d00M])[kb   ];// +  c2over27 ;
			real f_NE    = (D.f[dPP0])[kne  ];// +  c1over54 ;
			real f_SW    = (D.f[dMM0])[ksw  ];// +  c1over54 ;
			real f_SE    = (D.f[dPM0])[kse  ];// +  c1over54 ;
			real f_NW    = (D.f[dMP0])[knw  ];// +  c1over54 ;
			real f_TE    = (D.f[dP0P])[kte  ];// +  c1over54 ;
			real f_BW    = (D.f[dM0M])[kbw  ];// +  c1over54 ;
			real f_BE    = (D.f[dP0M])[kbe  ];// +  c1over54 ;
			real f_TW    = (D.f[dM0P])[ktw  ];// +  c1over54 ;
			real f_TN    = (D.f[d0PP])[ktn  ];// +  c1over54 ;
			real f_BS    = (D.f[d0MM])[kbs  ];// +  c1over54 ;
			real f_BN    = (D.f[d0PM])[kbn  ];// +  c1over54 ;
			real f_TS    = (D.f[d0MP])[kts  ];// +  c1over54 ;
			real f_R     = (D.f[d000])[kzero];// +  c8over27 ;
			real f_TNE   = (D.f[dPPP])[ktne ];// +  c1over216;
			real f_TSW   = (D.f[dMMP])[ktsw ];// +  c1over216;
			real f_TSE   = (D.f[dPMP])[ktse ];// +  c1over216;
			real f_TNW   = (D.f[dMPP])[ktnw ];// +  c1over216;
			real f_BNE   = (D.f[dPPM])[kbne ];// +  c1over216;
			real f_BSW   = (D.f[dMMM])[kbsw ];// +  c1over216;
			real f_BSE   = (D.f[dPMM])[kbse ];// +  c1over216;
			real f_BNW   = (D.f[dMPM])[kbnw ];// +  c1over216;
			////////////////////////////////////////////////////////////////////////////////////
			real fx = c0o1;
			real fy = c0o1;
			real fz = c0o1;
			////////////////////////////////////////////////////////////////////////////////////
			real rho=f_NW+f_W+f_SW+f_S+f_SE+f_E+f_NE+f_N+f_R+f_TN+f_BN+f_TS+f_BS+f_TE+f_BE+f_TW+f_BW+f_TNW+f_BNW+f_TNE+f_BNE+f_TSW+f_BSW+f_TSE+f_BSE+f_T+f_B+c1o1;// ACHTUNG ne EINS !!!!!!!!
			real pix=(f_NE+f_E+f_SE+f_TE+f_BE-f_NW-f_W-f_SW-f_TW-f_BW+f_TNE+f_BNE+f_TSE+f_BSE-f_TNW-f_BNW-f_TSW-f_BSW);
			real piy=(f_NE+f_N+f_NW+f_TN+f_BN-f_SE-f_S-f_SW-f_TS-f_BS+f_TNE+f_BNE+f_TNW+f_BNW-f_TSE-f_BSE-f_TSW-f_BSW);
			real piz=(f_TN+f_TS+f_TW+f_TE+f_T-f_BN-f_BS-f_BW-f_BE-f_B+f_TNE+f_TNW+f_TSE+f_TSW-f_BNE-f_BNW-f_BSE-f_BSW);
			real vvx=pix/rho + fx;
			real vvy=piy/rho + fy;
			real vvz=piz/rho + fz;
			real vx2=vvx*vvx;
			real vy2=vvy*vvy;
			real vz2=vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			real mfaaa = f_BSW;
			real mfaab = f_SW;
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
			real mfcba = f_BE;
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
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - c2o1*	m1 * vvz + vz2 * m0;
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
			mfaca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += c1o1;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9;
			mfbcc = m1 -		m0 * vvx;
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
			real OxyyPxzz  = c1o1;//two+(-omega);//one;
			real OxyyMxzz  = c1o1;//two+(-omega);//one;
			real O4        = c1o1;
			real O5        = c1o1;
			real O6        = c1o1;

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
			real CUMccc = mfccc  +(  -c4o1*  mfbbb * mfbbb  
									-          (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
									-    c4o1* (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
									-     c2o1* (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
									+(   c4o1* (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
									+     c2o1* (mfcaa * mfaca * mfaac)
									+ c16o1*  mfbba * mfbab * mfabb) / (rho * rho)
									-    c1o3* (mfacc + mfcac + mfcca)
									+    c1o9* (mfcaa + mfaca + mfaac)
									+(    c2o1* (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
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
 				mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 				mxxMyy    += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
 				mxxMzz    += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
 
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
			mfaca = c1o3 * (-c2o1 * mxxMyy +       mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - c2o1 * mxxMzz + mxxPyyPzz);

			//3.
			// linear combinations
			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

			//relax
			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
			mfbbb     += wadjust * (-mfbbb);
			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
			mxxyPyzz  += wadjust * (-mxxyPyzz);
			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
			mxxyMyzz  += wadjust * (-mxxyMyzz);
			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
			mxxzPyyz  += wadjust * (-mxxzPyyz);
			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
			mxxzMyyz  += wadjust * (-mxxzMyyz);
			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
			mxyyPxzz  += wadjust * (-mxyyPxzz);
			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
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
			mfccc = CUMccc  -(( -c4o1*  mfbbb * mfbbb  
							-          (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    c4o1* (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     c2o1* (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   c4o1* (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     c2o1* (mfcaa * mfaca * mfaac)
							+ c16o1*  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3* (mfacc + mfcac + mfcca)
							+    c1o9* (mfcaa + mfaca + mfaac)
							+(    c2o1* (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
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
			m1 = -mfaac        - c2o1* mfaab *  vvz         +  mfaaa       * (c1o1- vz2)              - c1o1* vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + 1.) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - c2o1* mfabb *  vvz         + mfaba * (c1o1- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - c2o1* mfacb *  vvz         +  mfaca         * (c1o1- vz2)              - c1o3 * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - c2o1* mfbab *  vvz         + mfbaa * (c1o1- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - c2o1* mfbbb *  vvz         + mfbba * (c1o1- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - c2o1* mfbcb *  vvz         + mfbca * (c1o1- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - c2o1* mfcab *  vvz         +  mfcaa         * (c1o1- vz2)              - c1o3 * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - c2o1* mfcbb *  vvz         + mfcba * (c1o1- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - c2o1* mfccb *  vvz         +  mfcca         * (c1o1- vz2)              - c1o9 * vz2; 
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
			m1 = -mfaca        - c2o1* mfaba *  vvy         +  mfaaa         * (c1o1- vy2)              - c1o6 * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - c2o1* mfabb *  vvy         +  mfaab         * (c1o1- vy2)              - c2o3 * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - c2o1* mfabc *  vvy         +  mfaac         * (c1o1- vy2)              - c1o6 * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - c2o1* mfbba *  vvy         + mfbaa * (c1o1- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - c2o1* mfbbb *  vvy         + mfbab * (c1o1- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - c2o1* mfbbc *  vvy         + mfbac * (c1o1- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - c2o1* mfcba *  vvy         +  mfcaa          * (c1o1- vy2)              - c1o18 * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - c2o1* mfcbb *  vvy         +  mfcab         * (c1o1- vy2)              - c2o9 * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - c2o1* mfcbc *  vvy         +  mfcac          * (c1o1- vy2)              - c1o18 * vy2; 
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
			m1 = -mfcaa        - c2o1* mfbaa *  vvx         +  mfaaa          * (c1o1- vx2)              - c1o36 * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - c2o1* mfbba *  vvx         +  mfaba         * (c1o1- vx2)              - c1o9 * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - c2o1* mfbca *  vvx         +  mfaca          * (c1o1- vx2)              - c1o36 * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - c2o1* mfbab *  vvx         +  mfaab         * (c1o1- vx2)              - c1o9 * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - c2o1* mfbbb *  vvx         +  mfabb         * (c1o1- vx2)              - c4o9 * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - c2o1* mfbcb *  vvx         +  mfacb         * (c1o1- vx2)              - c1o9 * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - c2o1* mfbac *  vvx         +  mfaac          * (c1o1- vx2)              - c1o36 * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - c2o1* mfbbc *  vvx         +  mfabc         * (c1o1- vx2)              - c1o9 * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - c2o1* mfbcc *  vvx         +  mfacc          * (c1o1- vx2)              - c1o36 * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dP00   ])[ke   ] = mfabb;// -  c2over27 ;//                                                                     
			(D.f[ dM00   ])[kw   ] = mfcbb;// -  c2over27 ;                                                                     
			(D.f[ d0P0   ])[kn   ] = mfbab;// -  c2over27 ;
			(D.f[ d0M0   ])[ks   ] = mfbcb;// -  c2over27 ;
			(D.f[ d00P   ])[kt   ] = mfbba;// -  c2over27 ;
			(D.f[ d00M   ])[kb   ] = mfbbc;// -  c2over27 ;
			(D.f[ dPP0  ])[kne  ] = mfaab;// -  c1over54 ;
			(D.f[ dMM0  ])[ksw  ] = mfccb;// -  c1over54 ;
			(D.f[ dPM0  ])[kse  ] = mfacb;// -  c1over54 ;
			(D.f[ dMP0  ])[knw  ] = mfcab;// -  c1over54 ;
			(D.f[ dP0P  ])[kte  ] = mfaba;// -  c1over54 ;
			(D.f[ dM0M  ])[kbw  ] = mfcbc;// -  c1over54 ;
			(D.f[ dP0M  ])[kbe  ] = mfabc;// -  c1over54 ;
			(D.f[ dM0P  ])[ktw  ] = mfcba;// -  c1over54 ;
			(D.f[ d0PP  ])[ktn  ] = mfbaa;// -  c1over54 ;
			(D.f[ d0MM  ])[kbs  ] = mfbcc;// -  c1over54 ;
			(D.f[ d0PM  ])[kbn  ] = mfbac;// -  c1over54 ;
			(D.f[ d0MP  ])[kts  ] = mfbca;// -  c1over54 ;
			(D.f[ d000])[kzero] = mfbbb;// -  c8over27 ;
			(D.f[ dPPP ])[ktne ] = mfaaa;// -  c1over216;
			(D.f[ dPMP ])[ktse ] = mfaca;// -  c1over216;
			(D.f[ dPPM ])[kbne ] = mfaac;// -  c1over216;
			(D.f[ dPMM ])[kbse ] = mfacc;// -  c1over216;
			(D.f[ dMPP ])[ktnw ] = mfcaa;// -  c1over216;
			(D.f[ dMMP ])[ktsw ] = mfcca;// -  c1over216;
			(D.f[ dMPM ])[kbnw ] = mfcac;// -  c1over216;
			(D.f[ dMMM ])[kbsw ] = mfccc;// -  c1over216;
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
__global__ void LB_Kernel_Kum_New_Comp_SRT_SP_27(
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	unsigned long long numberOfLBnodes,
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

	if (k<numberOfLBnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if (BC >= GEO_FLUID)
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
			real mfcbb = (D.f[dP00])[k   ];
			real mfabb = (D.f[dM00])[kw  ];
			real mfbcb = (D.f[d0P0])[k   ];
			real mfbab = (D.f[d0M0])[ks  ];
			real mfbbc = (D.f[d00P])[k   ];
			real mfbba = (D.f[d00M])[kb  ];
			real mfccb = (D.f[dPP0])[k   ];
			real mfaab = (D.f[dMM0])[ksw ];
			real mfcab = (D.f[dPM0])[ks  ];
			real mfacb = (D.f[dMP0])[kw  ];
			real mfcbc = (D.f[dP0P])[k   ];
			real mfaba = (D.f[dM0M])[kbw ];
			real mfcba = (D.f[dP0M])[kb  ];
			real mfabc = (D.f[dM0P])[kw  ];
			real mfbcc = (D.f[d0PP])[k   ];
			real mfbaa = (D.f[d0MM])[kbs ];
			real mfbca = (D.f[d0PM])[kb  ];
			real mfbac = (D.f[d0MP])[ks  ];
			real mfbbb = (D.f[d000])[k   ];
			real mfccc = (D.f[dPPP])[k   ];
			real mfaac = (D.f[dMMP])[ksw ];
			real mfcac = (D.f[dPMP])[ks  ];
			real mfacc = (D.f[dMPP])[kw  ];
			real mfcca = (D.f[dPPM])[kb  ];
			real mfaaa = (D.f[dMMM])[kbsw];
			real mfcaa = (D.f[dPMM])[kbs ];
			real mfaca = (D.f[dMPM])[kbw ];
			////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
				((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

			real rho = c1o1 + drho;
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
            real fx = forces[0] / (pow((double)c2o1, (double)level));
            real fy = forces[1] / (pow((double)c2o1, (double)level));
            real fz = forces[2] / (pow((double)c2o1, (double)level));
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
			//////////////////////////////////////////////////////////////////////////////////////
			//real wadjust;
			//real qudricLimitP = c1o100;// * 0.0001f;
			//real qudricLimitM = c1o100;// * 0.0001f;
			//real qudricLimitD = c1o100;// * 0.001f;
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


			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = omega;

			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			real OxyyPxzz = omega; 
			real OxyyMxzz = omega; 
			//real Oxyz = omega; 
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
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;	//ab 15.05.2015 verwendet
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho; //ab 15.05.2015 verwendet
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho; //ab 15.05.2015 verwendet

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

			//////////////////////////////////////////////////////////////////////////////////////
			real drhoPost =
				((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
					((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
			mfbbb += drho - drhoPost;
			////////////////////////////////////////////////////////////////////////////////////
			(D.f[dP00])[k   ] = mfabb;                                                                   
			(D.f[dM00])[kw  ] = mfcbb;                                                                 
			(D.f[d0P0])[k   ] = mfbab;
			(D.f[d0M0])[ks  ] = mfbcb;
			(D.f[d00P])[k   ] = mfbba;
			(D.f[d00M])[kb  ] = mfbbc;
			(D.f[dPP0])[k   ] = mfaab;
			(D.f[dMM0])[ksw ] = mfccb;
			(D.f[dPM0])[ks  ] = mfacb;
			(D.f[dMP0])[kw  ] = mfcab;
			(D.f[dP0P])[k   ] = mfaba;
			(D.f[dM0M])[kbw ] = mfcbc;
			(D.f[dP0M])[kb  ] = mfabc;
			(D.f[dM0P])[kw  ] = mfcba;
			(D.f[d0PP])[k   ] = mfbaa;
			(D.f[d0MM])[kbs ] = mfbcc;
			(D.f[d0PM])[kb  ] = mfbac;
			(D.f[d0MP])[ks  ] = mfbca;
			(D.f[d000])[k   ] = mfbbb;
			(D.f[dPPP])[k   ] = mfaaa;
			(D.f[dPMP])[ks  ] = mfaca;
			(D.f[dPPM])[kb  ] = mfaac;
			(D.f[dPMM])[kbs ] = mfacc;
			(D.f[dMPP])[kw  ] = mfcaa;
			(D.f[dMMP])[ksw ] = mfcca;
			(D.f[dMPM])[kbw ] = mfcac;
			(D.f[dMMM])[kbsw] = mfccc;
		}
	}
}
////////////////////////////////////////////////////////////////////////////////





