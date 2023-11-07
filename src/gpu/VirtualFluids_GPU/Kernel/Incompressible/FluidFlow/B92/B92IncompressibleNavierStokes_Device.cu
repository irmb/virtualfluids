#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void B92IncompressibleNavierStokes_Device(real omega,
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
				D.f[dP00] = &DDStart[dP00 * size_Mat];
				D.f[dM00] = &DDStart[dM00 * size_Mat];
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
				D.f[d000] = &DDStart[d000 * size_Mat];
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
				D.f[dM00] = &DDStart[dP00 * size_Mat];
				D.f[dP00] = &DDStart[dM00 * size_Mat];
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
				D.f[d000] = &DDStart[d000 * size_Mat];
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
			real fE = (D.f[dP00])[k];//ke
			real fW = (D.f[dM00])[kw];
			real fN = (D.f[DIR_0P0])[k];//kn
			real fS = (D.f[DIR_0M0])[ks];
			real fT = (D.f[DIR_00P])[k];//kt
			real fB = (D.f[DIR_00M])[kb];
			real fNE = (D.f[DIR_PP0])[k];//kne
			real fSW = (D.f[DIR_MM0])[ksw];
			real fSE = (D.f[DIR_PM0])[ks];//kse
			real fNW = (D.f[DIR_MP0])[kw];//knw
			real fTE = (D.f[DIR_P0P])[k];//kte
			real fBW = (D.f[DIR_M0M])[kbw];
			real fBE = (D.f[DIR_P0M])[kb];//kbe
			real fTW = (D.f[DIR_M0P])[kw];//ktw
			real fTN = (D.f[DIR_0PP])[k];//ktn
			real fBS = (D.f[DIR_0MM])[kbs];
			real fBN = (D.f[DIR_0PM])[kb];//kbn
			real fTS = (D.f[DIR_0MP])[ks];//kts
			real fZERO = (D.f[d000])[k];//kzero
			real fTNE = (D.f[DIR_PPP])[k];//ktne
			real fTSW = (D.f[DIR_MMP])[ksw];//ktsw
			real fTSE = (D.f[DIR_PMP])[ks];//ktse
			real fTNW = (D.f[DIR_MPP])[kw];//ktnw
			real fBNE = (D.f[DIR_PPM])[kb];//kbne
			real fBSW = (D.f[DIR_MMM])[kbsw];
			real fBSE = (D.f[DIR_PMM])[kbs];//kbse
			real fBNW = (D.f[DIR_MPM])[kbw];//kbnw
										   ////////////////////////////////////////////////////////////////////////////////







										   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										   //BGK incomp
										   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real drho = (fTNE + fBSW) + (fTSW + fBNE) + (fTSE + fBNW) + (fTNW + fBSE) + (fNE + fSW) + (fNW + fSE) + (fTE + fBW) + (fBE + fTW) + (fTN + fBS) + (fBN + fTS) + (fE + fW) + (fN + fS) + (fT + fB) + fZERO;
			real vx1 = (fTNE - fBSW) + (fBNE - fTSW) + (fTSE - fBNW) + (fBSE - fTNW) + (fNE - fSW) + (fSE - fNW) + (fTE - fBW) + (fBE - fTW) + (fE - fW);
			real vx2 = (fTNE - fBSW) + (fBNE - fTSW) + (fBNW - fTSE) + (fTNW - fBSE) + (fNE - fSW) + (fNW - fSE) + (fTN - fBS) + (fBN - fTS) + (fN - fS);
			real vx3 = (fTNE - fBSW) + (fTSW - fBNE) + (fTSE - fBNW) + (fTNW - fBSE) + (fTE - fBW) + (fTW - fBE) + (fTN - fBS) + (fTS - fBN) + (fT - fB);
			real cusq = c3o2*(vx1*vx1 + vx2*vx2 + vx3*vx3);
			//////////////////////////////////////////////////////////////////////////                            
			fZERO = fZERO *(c1o1 + (-omega)) - (-omega)*   c8o27*  (drho - cusq);
			fE = fE    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + c3o1*(vx1)+c9o2*(vx1)*(vx1)-cusq);
			fW = fW    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + c3o1*(-vx1) + c9o2*(-vx1)*(-vx1) - cusq);
			fN = fN    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + c3o1*(vx2)+c9o2*(vx2)*(vx2)-cusq);
			fS = fS    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + c3o1*(-vx2) + c9o2*(-vx2)*(-vx2) - cusq);
			fT = fT    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + c3o1*(vx3)+c9o2*(vx3)*(vx3)-cusq);
			fB = fB    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + c3o1*(-vx3) + c9o2*(-vx3)*(-vx3) - cusq);
			fNE = fNE   *(c1o1 + (-omega)) - (-omega)*   c1o54*  (drho + c3o1*(vx1 + vx2) + c9o2*(vx1 + vx2)*(vx1 + vx2) - cusq);
			fSW = fSW   *(c1o1 + (-omega)) - (-omega)*   c1o54*  (drho + c3o1*(-vx1 - vx2) + c9o2*(-vx1 - vx2)*(-vx1 - vx2) - cusq);
			fSE = fSE   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(vx1 - vx2) + c9o2*(vx1 - vx2)*(vx1 - vx2) - cusq);
			fNW = fNW   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(-vx1 + vx2) + c9o2*(-vx1 + vx2)*(-vx1 + vx2) - cusq);
			fTE = fTE   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(vx1 + vx3) + c9o2*(vx1 + vx3)*(vx1 + vx3) - cusq);
			fBW = fBW   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(-vx1 - vx3) + c9o2*(-vx1 - vx3)*(-vx1 - vx3) - cusq);
			fBE = fBE   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(vx1 - vx3) + c9o2*(vx1 - vx3)*(vx1 - vx3) - cusq);
			fTW = fTW   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(-vx1 + vx3) + c9o2*(-vx1 + vx3)*(-vx1 + vx3) - cusq);
			fTN = fTN   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(vx2 + vx3) + c9o2*(vx2 + vx3)*(vx2 + vx3) - cusq);
			fBS = fBS   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(-vx2 - vx3) + c9o2*(-vx2 - vx3)*(-vx2 - vx3) - cusq);
			fBN = fBN   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(vx2 - vx3) + c9o2*(vx2 - vx3)*(vx2 - vx3) - cusq);
			fTS = fTS   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(-vx2 + vx3) + c9o2*(-vx2 + vx3)*(-vx2 + vx3) - cusq);
			fTNE = fTNE  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(vx1 + vx2 + vx3) + c9o2*(vx1 + vx2 + vx3)*(vx1 + vx2 + vx3) - cusq);
			fBSW = fBSW  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(-vx1 - vx2 - vx3) + c9o2*(-vx1 - vx2 - vx3)*(-vx1 - vx2 - vx3) - cusq);
			fBNE = fBNE  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(vx1 + vx2 - vx3) + c9o2*(vx1 + vx2 - vx3)*(vx1 + vx2 - vx3) - cusq);
			fTSW = fTSW  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(-vx1 - vx2 + vx3) + c9o2*(-vx1 - vx2 + vx3)*(-vx1 - vx2 + vx3) - cusq);
			fTSE = fTSE  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(vx1 - vx2 + vx3) + c9o2*(vx1 - vx2 + vx3)*(vx1 - vx2 + vx3) - cusq);
			fBNW = fBNW  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(-vx1 + vx2 - vx3) + c9o2*(-vx1 + vx2 - vx3)*(-vx1 + vx2 - vx3) - cusq);
			fBSE = fBSE  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(vx1 - vx2 - vx3) + c9o2*(vx1 - vx2 - vx3)*(vx1 - vx2 - vx3) - cusq);
			fTNW = fTNW  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(-vx1 + vx2 + vx3) + c9o2*(-vx1 + vx2 + vx3)*(-vx1 + vx2 + vx3) - cusq);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







			//////////////////////////////////////////////////////////////////////////                            
			(D.f[dP00])[k] = fW;
			(D.f[dM00])[kw] = fE;
			(D.f[DIR_0P0])[k] = fS;
			(D.f[DIR_0M0])[ks] = fN;
			(D.f[DIR_00P])[k] = fB;
			(D.f[DIR_00M])[kb] = fT;
			(D.f[DIR_PP0])[k] = fSW;
			(D.f[DIR_MM0])[ksw] = fNE;
			(D.f[DIR_PM0])[ks] = fNW;
			(D.f[DIR_MP0])[kw] = fSE;
			(D.f[DIR_P0P])[k] = fBW;
			(D.f[DIR_M0M])[kbw] = fTE;
			(D.f[DIR_P0M])[kb] = fTW;
			(D.f[DIR_M0P])[kw] = fBE;
			(D.f[DIR_0PP])[k] = fBS;
			(D.f[DIR_0MM])[kbs] = fTN;
			(D.f[DIR_0PM])[kb] = fTS;
			(D.f[DIR_0MP])[ks] = fBN;
			(D.f[d000])[k] = fZERO;
			(D.f[DIR_PPP])[k] = fBSW;
			(D.f[DIR_PMP])[ks] = fBNW;
			(D.f[DIR_PPM])[kb] = fTSW;
			(D.f[DIR_PMM])[kbs] = fTNW;
			(D.f[DIR_MPP])[kw] = fBSE;
			(D.f[DIR_MMP])[ksw] = fBNE;
			(D.f[DIR_MPM])[kbw] = fTSE;
			(D.f[DIR_MMM])[kbsw] = fTNE;
			//////////////////////////////////////////////////////////////////////////                            
		}
	}
}