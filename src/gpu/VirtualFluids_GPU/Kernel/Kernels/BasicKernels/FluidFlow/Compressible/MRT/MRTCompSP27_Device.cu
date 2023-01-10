#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void LB_Kernel_MRT_Comp_SP_27(real omega,
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
				D.f[DIR_P00] = &DDStart[DIR_P00 * size_Mat];
				D.f[DIR_M00] = &DDStart[DIR_M00 * size_Mat];
				D.f[DIR_0P0] = &DDStart[DIR_0P0 * size_Mat];
				D.f[DIR_0M0] = &DDStart[DIR_0M0 * size_Mat];
				D.f[DIR_00P] = &DDStart[DIR_00P * size_Mat];
				D.f[DIR_00M] = &DDStart[DIR_00M * size_Mat];
				D.f[DIR_PP0] = &DDStart[DIR_PP0 * size_Mat];
				D.f[DIR_MM0] = &DDStart[DIR_MM0 * size_Mat];
				D.f[DIR_PM0] = &DDStart[DIR_PM0 * size_Mat];
				D.f[DIR_MP0] = &DDStart[DIR_MP0 * size_Mat];
				D.f[DIR_P0P] = &DDStart[DIR_P0P * size_Mat];
				D.f[DIR_M0M] = &DDStart[DIR_M0M * size_Mat];
				D.f[DIR_P0M] = &DDStart[DIR_P0M * size_Mat];
				D.f[DIR_M0P] = &DDStart[DIR_M0P * size_Mat];
				D.f[DIR_0PP] = &DDStart[DIR_0PP * size_Mat];
				D.f[DIR_0MM] = &DDStart[DIR_0MM * size_Mat];
				D.f[DIR_0PM] = &DDStart[DIR_0PM * size_Mat];
				D.f[DIR_0MP] = &DDStart[DIR_0MP * size_Mat];
				D.f[DIR_000] = &DDStart[DIR_000 * size_Mat];
				D.f[DIR_PPP] = &DDStart[DIR_PPP * size_Mat];
				D.f[DIR_MMP] = &DDStart[DIR_MMP * size_Mat];
				D.f[DIR_PMP] = &DDStart[DIR_PMP * size_Mat];
				D.f[DIR_MPP] = &DDStart[DIR_MPP * size_Mat];
				D.f[DIR_PPM] = &DDStart[DIR_PPM * size_Mat];
				D.f[DIR_MMM] = &DDStart[DIR_MMM * size_Mat];
				D.f[DIR_PMM] = &DDStart[DIR_PMM * size_Mat];
				D.f[DIR_MPM] = &DDStart[DIR_MPM * size_Mat];
			}
			else
			{
				D.f[DIR_M00] = &DDStart[DIR_P00 * size_Mat];
				D.f[DIR_P00] = &DDStart[DIR_M00 * size_Mat];
				D.f[DIR_0M0] = &DDStart[DIR_0P0 * size_Mat];
				D.f[DIR_0P0] = &DDStart[DIR_0M0 * size_Mat];
				D.f[DIR_00M] = &DDStart[DIR_00P * size_Mat];
				D.f[DIR_00P] = &DDStart[DIR_00M * size_Mat];
				D.f[DIR_MM0] = &DDStart[DIR_PP0 * size_Mat];
				D.f[DIR_PP0] = &DDStart[DIR_MM0 * size_Mat];
				D.f[DIR_MP0] = &DDStart[DIR_PM0 * size_Mat];
				D.f[DIR_PM0] = &DDStart[DIR_MP0 * size_Mat];
				D.f[DIR_M0M] = &DDStart[DIR_P0P * size_Mat];
				D.f[DIR_P0P] = &DDStart[DIR_M0M * size_Mat];
				D.f[DIR_M0P] = &DDStart[DIR_P0M * size_Mat];
				D.f[DIR_P0M] = &DDStart[DIR_M0P * size_Mat];
				D.f[DIR_0MM] = &DDStart[DIR_0PP * size_Mat];
				D.f[DIR_0PP] = &DDStart[DIR_0MM * size_Mat];
				D.f[DIR_0MP] = &DDStart[DIR_0PM * size_Mat];
				D.f[DIR_0PM] = &DDStart[DIR_0MP * size_Mat];
				D.f[DIR_000] = &DDStart[DIR_000 * size_Mat];
				D.f[DIR_MMM] = &DDStart[DIR_PPP * size_Mat];
				D.f[DIR_PPM] = &DDStart[DIR_MMP * size_Mat];
				D.f[DIR_MPM] = &DDStart[DIR_PMP * size_Mat];
				D.f[DIR_PMM] = &DDStart[DIR_MPP * size_Mat];
				D.f[DIR_MMP] = &DDStart[DIR_PPM * size_Mat];
				D.f[DIR_PPP] = &DDStart[DIR_MMM * size_Mat];
				D.f[DIR_MPP] = &DDStart[DIR_PMM * size_Mat];
				D.f[DIR_PMP] = &DDStart[DIR_MPM * size_Mat];
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
			real mfcbb = (D.f[DIR_P00])[k];//[ke   ];// +  c2over27 ;(D.f[DIR_P00])[k  ];//ke
			real mfabb = (D.f[DIR_M00])[kw];//[kw   ];// +  c2over27 ;(D.f[DIR_M00])[kw ];
			real mfbcb = (D.f[DIR_0P0])[k];//[kn   ];// +  c2over27 ;(D.f[DIR_0P0])[k  ];//kn
			real mfbab = (D.f[DIR_0M0])[ks];//[ks   ];// +  c2over27 ;(D.f[DIR_0M0])[ks ];
			real mfbbc = (D.f[DIR_00P])[k];//[kt   ];// +  c2over27 ;(D.f[DIR_00P])[k  ];//kt
			real mfbba = (D.f[DIR_00M])[kb];//[kb   ];// +  c2over27 ;(D.f[DIR_00M])[kb ];
			real mfccb = (D.f[DIR_PP0])[k];//[kne  ];// +  c1over54 ;(D.f[DIR_PP0])[k  ];//kne
			real mfaab = (D.f[DIR_MM0])[ksw];//[ksw  ];// +  c1over54 ;(D.f[DIR_MM0])[ksw];
			real mfcab = (D.f[DIR_PM0])[ks];//[kse  ];// +  c1over54 ;(D.f[DIR_PM0])[ks ];//kse
			real mfacb = (D.f[DIR_MP0])[kw];//[knw  ];// +  c1over54 ;(D.f[DIR_MP0])[kw ];//knw
			real mfcbc = (D.f[DIR_P0P])[k];//[kte  ];// +  c1over54 ;(D.f[DIR_P0P])[k  ];//kte
			real mfaba = (D.f[DIR_M0M])[kbw];//[kbw  ];// +  c1over54 ;(D.f[DIR_M0M])[kbw];
			real mfcba = (D.f[DIR_P0M])[kb];//[kbe  ];// +  c1over54 ;(D.f[DIR_P0M])[kb ];//kbe
			real mfabc = (D.f[DIR_M0P])[kw];//[ktw  ];// +  c1over54 ;(D.f[DIR_M0P])[kw ];//ktw
			real mfbcc = (D.f[DIR_0PP])[k];//[ktn  ];// +  c1over54 ;(D.f[DIR_0PP])[k  ];//ktn
			real mfbaa = (D.f[DIR_0MM])[kbs];//[kbs  ];// +  c1over54 ;(D.f[DIR_0MM])[kbs];
			real mfbca = (D.f[DIR_0PM])[kb];//[kbn  ];// +  c1over54 ;(D.f[DIR_0PM])[kb ];//kbn
			real mfbac = (D.f[DIR_0MP])[ks];//[kts  ];// +  c1over54 ;(D.f[DIR_0MP])[ks ];//kts
			real mfbbb = (D.f[DIR_000])[k];//[kzero];// +  c8over27 ;(D.f[DIR_000])[k  ];//kzero
			real mfccc = (D.f[DIR_PPP])[k];//[ktne ];// +  c1over216;(D.f[DIR_PPP])[k  ];//ktne
			real mfaac = (D.f[DIR_MMP])[ksw];//[ktsw ];// +  c1over216;(D.f[DIR_MMP])[ksw];//ktsw
			real mfcac = (D.f[DIR_PMP])[ks];//[ktse ];// +  c1over216;(D.f[DIR_PMP])[ks ];//ktse
			real mfacc = (D.f[DIR_MPP])[kw];//[ktnw ];// +  c1over216;(D.f[DIR_MPP])[kw ];//ktnw
			real mfcca = (D.f[DIR_PPM])[kb];//[kbne ];// +  c1over216;(D.f[DIR_PPM])[kb ];//kbne
			real mfaaa = (D.f[DIR_MMM])[kbsw];//[kbsw ];// +  c1over216;(D.f[DIR_MMM])[kbsw];
			real mfcaa = (D.f[DIR_PMM])[kbs];//[kbse ];// +  c1over216;(D.f[DIR_PMM])[kbs];//kbse
			real mfaca = (D.f[DIR_MPM])[kbw];//[kbnw ];// +  c1over216;(D.f[DIR_MPM])[kbw];//kbnw
											////////////////////////////////////////////////////////////////////////////////////
			real rho = (mfccc + mfaaa + mfaca + mfcac + mfacc + mfcaa + mfaac + mfcca +
				mfbac + mfbca + mfbaa + mfbcc + mfabc + mfcba + mfaba + mfcbc + mfacb + mfcab + mfaab + mfccb +
				mfabb + mfcbb + mfbab + mfbcb + mfbba + mfbbc + mfbbb) + c1o1;//  !!!Achtung + one
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
			//fast
			real oMdrho = c1o1; // comp special
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
			real m0, m1, m2;
			real vx2;
			real vy2;
			real vz2;
			vx2 = vvx*vvx;
			vy2 = vvy*vvy;
			vz2 = vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			//real wadjust;
			//real qudricLimit = 0.01f;
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
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaba + mfabc;
			m1 = mfabc - mfaba;
			m0 = m2 + mfabb;
			mfaba = m0;
			m0 += c1o9 * oMdrho;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaca + mfacc;
			m1 = mfacc - mfaca;
			m0 = m2 + mfacb;
			mfaca = m0;
			m0 += c1o36 * oMdrho;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbaa + mfbac;
			m1 = mfbac - mfbaa;
			m0 = m2 + mfbab;
			mfbaa = m0;
			m0 += c1o9 * oMdrho;
			mfbab = m1;
			mfbac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbba + mfbbc;
			m1 = mfbbc - mfbba;
			m0 = m2 + mfbbb;
			mfbba = m0;
			m0 += c4o9 * oMdrho;
			mfbbb = m1;
			mfbbc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbca + mfbcc;
			m1 = mfbcc - mfbca;
			m0 = m2 + mfbcb;
			mfbca = m0;
			m0 += c1o9 * oMdrho;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcaa + mfcac;
			m1 = mfcac - mfcaa;
			m0 = m2 + mfcab;
			mfcaa = m0;
			m0 += c1o36 * oMdrho;
			mfcab = m1;
			mfcac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcba + mfcbc;
			m1 = mfcbc - mfcba;
			m0 = m2 + mfcbb;
			mfcba = m0;
			m0 += c1o9 * oMdrho;
			mfcbb = m1;
			mfcbc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcca + mfccc;
			m1 = mfccc - mfcca;
			m0 = m2 + mfccb;
			mfcca = m0;
			m0 += c1o36 * oMdrho;
			mfccb = m1;
			mfccc = m2;
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
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaab + mfacb;
			m1 = mfacb - mfaab;
			m0 = m2 + mfabb;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaac + mfacc;
			m1 = mfacc - mfaac;
			m0 = m2 + mfabc;
			mfaac = m0;
			m0 += c1o18 * oMdrho;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbaa + mfbca;
			m1 = mfbca - mfbaa;
			m0 = m2 + mfbba;
			mfbaa = m0;
			m0 += c2o3 * oMdrho;
			mfbba = m1;
			mfbca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbab + mfbcb;
			m1 = mfbcb - mfbab;
			m0 = m2 + mfbbb;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfbac + mfbcc;
			m1 = mfbcc - mfbac;
			m0 = m2 + mfbbc;
			mfbac = m0;
			m0 += c2o9 * oMdrho;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcaa + mfcca;
			m1 = mfcca - mfcaa;
			m0 = m2 + mfcba;
			mfcaa = m0;
			m0 += c1o6 * oMdrho;
			mfcba = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcab + mfccb;
			m1 = mfccb - mfcab;
			m0 = m2 + mfcbb;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfcac + mfccc;
			m1 = mfccc - mfcac;
			m0 = m2 + mfcbc;
			mfcac = m0;
			m0 += c1o18 * oMdrho;
			mfcbc = m1;
			mfccc = m2;
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
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaba + mfcba;
			m1 = mfcba - mfaba;
			m0 = m2 + mfbba;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaca + mfcca;
			m1 = mfcca - mfaca;
			m0 = m2 + mfbca;
			mfaca = m0;
			m0 += c1o3 * oMdrho;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaab + mfcab;
			m1 = mfcab - mfaab;
			m0 = m2 + mfbab;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfabb + mfcbb;
			m1 = mfcbb - mfabb;
			m0 = m2 + mfbbb;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfacb + mfccb;
			m1 = mfccb - mfacb;
			m0 = m2 + mfbcb;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfaac + mfcac;
			m1 = mfcac - mfaac;
			m0 = m2 + mfbac;
			mfaac = m0;
			m0 += c1o3 * oMdrho;
			mfbac = m1;
			mfcac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfabc + mfcbc;
			m1 = mfcbc - mfabc;
			m0 = m2 + mfbbc;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m2 = mfacc + mfccc;
			m1 = mfccc - mfacc;
			m0 = m2 + mfbcc;
			mfacc = m0;
			m0 += c1o9 * oMdrho;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// MRT
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = c1o1;
			real OxyyPxzz = c1o1;//two-omega;//eight*(two-omega)/(eight -omega);//one;//omega;//two-omega;//
			real OxyyMxzz = c1o1;//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
			real O4 = c1o1;
			real O5 = c1o1;
			real O6 = c1o1;

			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy = mfcaa - mfaca;
			real mxxMzz = mfcaa - mfaac;

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//incl. correction
			{
				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz + (-c2o1*vx2 + vy2 + vz2)*rho) + c1o2 * OxxPyyPzz * (mfaaa + (vx2 + vy2 + vz2)*rho - mxxPyyPzz);
				real dyuy = dxux + omega * c3o2 * (mxxMyy + (-vx2 + vy2)*rho);
				real dzuz = dxux + omega * c3o2 * (mxxMzz + (-vx2 + vz2)*rho);

				//relax
				mxxPyyPzz += OxxPyyPzz*(mfaaa + (vx2 + vy2 + vz2)*rho - mxxPyyPzz) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
				mxxMyy += omega * ((vx2 - vy2)*rho - mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
				mxxMzz += omega * ((vx2 - vz2)*rho - mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			// 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// 			//no correction
			// 			mxxPyyPzz += OxxPyyPzz*(mfaaa+(vx2+vy2+vz2)*rho-mxxPyyPzz);
			// 			mxxMyy    += -(-omega) * ((vx2-vy2)*rho-mxxMyy);
			// 			mxxMzz    += -(-omega) * ((vx2-vz2)*rho-mxxMzz);
			// 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb += omega * ((vvy*vvz)*rho - mfabb);
			mfbab += omega * ((vvx*vvz)*rho - mfbab);
			mfbba += omega * ((vvx*vvy)*rho - mfbba);

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

			mxxyMyzz += OxyyMxzz*((vx2 - vz2)*vvy*rho - mxxyMyzz);
			mxxzMyyz += OxyyMxzz*((vx2 - vy2)*vvz*rho - mxxzMyyz);
			mxyyMxzz += OxyyMxzz*((vy2 - vz2)*vvx*rho - mxyyMxzz);

			mxxyPyzz += OxyyPxzz*((c2o3 + vx2 + vz2)*vvy*rho - mxxyPyzz);
			mxxzPyyz += OxyyPxzz*((c2o3 + vx2 + vy2)*vvz*rho - mxxzPyyz);
			mxyyPxzz += OxyyPxzz*((c2o3 + vy2 + vz2)*vvx*rho - mxyyPxzz);

			mfbbb += OxyyMxzz * (vvx*vvy*vvz*rho - mfbbb);

			mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
			mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

			//4.
			//mfacc += O4*((c1o3+vy2)*(c1o3+vz2)*rho+c1o9*(mfaaa-one)-mfacc);
			//mfcac += O4*((c1o3+vx2)*(c1o3+vz2)*rho+c1o9*(mfaaa-one)-mfcac);
			//mfcca += O4*((c1o3+vx2)*(c1o3+vy2)*rho+c1o9*(mfaaa-one)-mfcca);
			mfacc += O4*((c1o3 + vy2)*(c1o3 + vz2)*rho - c1o9 - mfacc);
			mfcac += O4*((c1o3 + vx2)*(c1o3 + vz2)*rho - c1o9 - mfcac);
			mfcca += O4*((c1o3 + vx2)*(c1o3 + vy2)*rho - c1o9 - mfcca);

			mfcbb += O4*((c1o3 + vx2)*vvy*vvz*rho - mfcbb);
			mfbcb += O4*((c1o3 + vy2)*vvx*vvz*rho - mfbcb);
			mfbbc += O4*((c1o3 + vz2)*vvx*vvy*rho - mfbbc);

			//5.
			mfbcc += O5*((c1o3 + vy2)*(c1o3 + vz2)*vvx*rho - mfbcc);
			mfcbc += O5*((c1o3 + vx2)*(c1o3 + vz2)*vvy*rho - mfcbc);
			mfccb += O5*((c1o3 + vx2)*(c1o3 + vy2)*vvz*rho - mfccb);

			//6.
			mfccc += O6*((c1o3 + vx2)*(c1o3 + vy2)*(c1o3 + vz2)*rho - c1o27 - mfccc);


			//bad fix
			vvx = c0o1;
			vvy = c0o1;
			vvz = c0o1;
			vx2 = c0o1;
			vy2 = c0o1;
			vz2 = c0o1;
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
			(D.f[DIR_P00])[k] = mfabb;//(D.f[ DIR_P00   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ DIR_P00   ])[k   ]                                                                     
			(D.f[DIR_M00])[kw] = mfcbb;//(D.f[ DIR_M00   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ DIR_M00   ])[kw  ]                                                                   
			(D.f[DIR_0P0])[k] = mfbab;//(D.f[ DIR_0P0   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ DIR_0P0   ])[k   ]
			(D.f[DIR_0M0])[ks] = mfbcb;//(D.f[ DIR_0M0   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ DIR_0M0   ])[ks  ]
			(D.f[DIR_00P])[k] = mfbba;//(D.f[ DIR_00P   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ DIR_00P   ])[k   ]
			(D.f[DIR_00M])[kb] = mfbbc;//(D.f[ DIR_00M   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ DIR_00M   ])[kb  ]
			(D.f[DIR_PP0])[k] = mfaab;//(D.f[ DIR_PP0  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ DIR_PP0  ])[k   ]
			(D.f[DIR_MM0])[ksw] = mfccb;//(D.f[ DIR_MM0  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ DIR_MM0  ])[ksw ]
			(D.f[DIR_PM0])[ks] = mfacb;//(D.f[ DIR_PM0  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ DIR_PM0  ])[ks  ]
			(D.f[DIR_MP0])[kw] = mfcab;//(D.f[ DIR_MP0  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ DIR_MP0  ])[kw  ]
			(D.f[DIR_P0P])[k] = mfaba;//(D.f[ DIR_P0P  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ DIR_P0P  ])[k   ]
			(D.f[DIR_M0M])[kbw] = mfcbc;//(D.f[ DIR_M0M  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ DIR_M0M  ])[kbw ]
			(D.f[DIR_P0M])[kb] = mfabc;//(D.f[ DIR_P0M  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ DIR_P0M  ])[kb  ]
			(D.f[DIR_M0P])[kw] = mfcba;//(D.f[ DIR_M0P  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ DIR_M0P  ])[kw  ]
			(D.f[DIR_0PP])[k] = mfbaa;//(D.f[ DIR_0PP  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ DIR_0PP  ])[k   ]
			(D.f[DIR_0MM])[kbs] = mfbcc;//(D.f[ DIR_0MM  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ DIR_0MM  ])[kbs ]
			(D.f[DIR_0PM])[kb] = mfbac;//(D.f[ DIR_0PM  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ DIR_0PM  ])[kb  ]
			(D.f[DIR_0MP])[ks] = mfbca;//(D.f[ DIR_0MP  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ DIR_0MP  ])[ks  ]
			(D.f[DIR_000])[k] = mfbbb;//(D.f[ DIR_000])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ DIR_000])[k   ]
			(D.f[DIR_PPP])[k] = mfaaa;//(D.f[ DIR_PPP ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ DIR_PPP ])[k   ]
			(D.f[DIR_PMP])[ks] = mfaca;//(D.f[ DIR_PMP ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ DIR_PMP ])[ks  ]
			(D.f[DIR_PPM])[kb] = mfaac;//(D.f[ DIR_PPM ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ DIR_PPM ])[kb  ]
			(D.f[DIR_PMM])[kbs] = mfacc;//(D.f[ DIR_PMM ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ DIR_PMM ])[kbs ]
			(D.f[DIR_MPP])[kw] = mfcaa;//(D.f[ DIR_MPP ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ DIR_MPP ])[kw  ]
			(D.f[DIR_MMP])[ksw] = mfcca;//(D.f[ DIR_MMP ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ DIR_MMP ])[ksw ]
			(D.f[DIR_MPM])[kbw] = mfcac;//(D.f[ DIR_MPM ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ DIR_MPM ])[kbw ]
			(D.f[DIR_MMM])[kbsw] = mfccc;//(D.f[ DIR_MMM ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ DIR_MMM ])[kbsw]
										////////////////////////////////////////////////////////////////////////////////////
		}
	}
}