#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
#include "math.h"


extern "C" __global__ void LB_Kernel_BGK_Comp_SP_27(	real omega,
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
	///////////////////////////////////////////////////////////////////////////////

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
				D.f[dirREST] = &DDStart[dirREST*size_Mat];
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
				D.f[dirREST] = &DDStart[dirREST*size_Mat];
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
			real fE = (D.f[E])[k];//ke
			real fW = (D.f[W])[kw];
			real fN = (D.f[N])[k];//kn
			real fS = (D.f[S])[ks];
			real fT = (D.f[T])[k];//kt
			real fB = (D.f[B])[kb];
			real fNE = (D.f[NE])[k];//kne
			real fSW = (D.f[SW])[ksw];
			real fSE = (D.f[SE])[ks];//kse
			real fNW = (D.f[NW])[kw];//knw
			real fTE = (D.f[TE])[k];//kte
			real fBW = (D.f[BW])[kbw];
			real fBE = (D.f[BE])[kb];//kbe
			real fTW = (D.f[TW])[kw];//ktw
			real fTN = (D.f[TN])[k];//ktn
			real fBS = (D.f[BS])[kbs];
			real fBN = (D.f[BN])[kb];//kbn
			real fTS = (D.f[TS])[ks];//kts
			real fZERO = (D.f[dirREST])[k];//kzero
			real fTNE = (D.f[TNE])[k];//ktne
			real fTSW = (D.f[TSW])[ksw];//ktsw
			real fTSE = (D.f[TSE])[ks];//ktse
			real fTNW = (D.f[TNW])[kw];//ktnw
			real fBNE = (D.f[BNE])[kb];//kbne
			real fBSW = (D.f[BSW])[kbsw];
			real fBSE = (D.f[BSE])[kbs];//kbse
			real fBNW = (D.f[BNW])[kbw];//kbnw
										   ////////////////////////////////////////////////////////////////////////////////
			real drho = (fTNE + fBSW) + (fTSW + fBNE) + (fTSE + fBNW) + (fTNW + fBSE) + (fNE + fSW) + (fNW + fSE) + (fTE + fBW) + (fBE + fTW) + (fTN + fBS) + (fBN + fTS) + (fE + fW) + (fN + fS) + (fT + fB) + fZERO;
			real rho = drho + c1o1;
			real OORho = c1o1 / rho;
			real vx1 = OORho*((fTNE - fBSW) + (fBNE - fTSW) + (fTSE - fBNW) + (fBSE - fTNW) + (fNE - fSW) + (fSE - fNW) + (fTE - fBW) + (fBE - fTW) + (fE - fW));
			real vx2 = OORho*((fTNE - fBSW) + (fBNE - fTSW) + (fBNW - fTSE) + (fTNW - fBSE) + (fNE - fSW) + (fNW - fSE) + (fTN - fBS) + (fBN - fTS) + (fN - fS));
			real vx3 = OORho*((fTNE - fBSW) + (fTSW - fBNE) + (fTSE - fBNW) + (fTNW - fBSE) + (fTE - fBW) + (fTW - fBE) + (fTN - fBS) + (fTS - fBN) + (fT - fB));
			////////////////////////////////////////////////////////////////////////////////







			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//BGK comp
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real cusq = c3o2*(vx1*vx1 + vx2*vx2 + vx3*vx3);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			fZERO = fZERO *(c1o1 + (-omega)) - (-omega)*   c8o27*  (drho - rho * cusq);
			fE = fE    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + rho * (c3o1*(vx1)+c9o2*(vx1)*(vx1)-cusq));
			fW = fW    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + rho * (c3o1*(-vx1) + c9o2*(-vx1)*(-vx1) - cusq));
			fN = fN    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + rho * (c3o1*(vx2)+c9o2*(vx2)*(vx2)-cusq));
			fS = fS    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + rho * (c3o1*(-vx2) + c9o2*(-vx2)*(-vx2) - cusq));
			fT = fT    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + rho * (c3o1*(vx3)+c9o2*(vx3)*(vx3)-cusq));
			fB = fB    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + rho * (c3o1*(-vx3) + c9o2*(-vx3)*(-vx3) - cusq));
			fNE = fNE   *(c1o1 + (-omega)) - (-omega)*   c1o54*  (drho + rho * (c3o1*(vx1 + vx2) + c9o2*(vx1 + vx2)*(vx1 + vx2) - cusq));
			fSW = fSW   *(c1o1 + (-omega)) - (-omega)*   c1o54*  (drho + rho * (c3o1*(-vx1 - vx2) + c9o2*(-vx1 - vx2)*(-vx1 - vx2) - cusq));
			fSE = fSE   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(vx1 - vx2) + c9o2*(vx1 - vx2)*(vx1 - vx2) - cusq));
			fNW = fNW   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(-vx1 + vx2) + c9o2*(-vx1 + vx2)*(-vx1 + vx2) - cusq));
			fTE = fTE   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(vx1 + vx3) + c9o2*(vx1 + vx3)*(vx1 + vx3) - cusq));
			fBW = fBW   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(-vx1 - vx3) + c9o2*(-vx1 - vx3)*(-vx1 - vx3) - cusq));
			fBE = fBE   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(vx1 - vx3) + c9o2*(vx1 - vx3)*(vx1 - vx3) - cusq));
			fTW = fTW   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(-vx1 + vx3) + c9o2*(-vx1 + vx3)*(-vx1 + vx3) - cusq));
			fTN = fTN   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(vx2 + vx3) + c9o2*(vx2 + vx3)*(vx2 + vx3) - cusq));
			fBS = fBS   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(-vx2 - vx3) + c9o2*(-vx2 - vx3)*(-vx2 - vx3) - cusq));
			fBN = fBN   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(vx2 - vx3) + c9o2*(vx2 - vx3)*(vx2 - vx3) - cusq));
			fTS = fTS   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(-vx2 + vx3) + c9o2*(-vx2 + vx3)*(-vx2 + vx3) - cusq));
			fTNE = fTNE  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(vx1 + vx2 + vx3) + c9o2*(vx1 + vx2 + vx3)*(vx1 + vx2 + vx3) - cusq));
			fBSW = fBSW  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(-vx1 - vx2 - vx3) + c9o2*(-vx1 - vx2 - vx3)*(-vx1 - vx2 - vx3) - cusq));
			fBNE = fBNE  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(vx1 + vx2 - vx3) + c9o2*(vx1 + vx2 - vx3)*(vx1 + vx2 - vx3) - cusq));
			fTSW = fTSW  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(-vx1 - vx2 + vx3) + c9o2*(-vx1 - vx2 + vx3)*(-vx1 - vx2 + vx3) - cusq));
			fTSE = fTSE  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(vx1 - vx2 + vx3) + c9o2*(vx1 - vx2 + vx3)*(vx1 - vx2 + vx3) - cusq));
			fBNW = fBNW  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(-vx1 + vx2 - vx3) + c9o2*(-vx1 + vx2 - vx3)*(-vx1 + vx2 - vx3) - cusq));
			fBSE = fBSE  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(vx1 - vx2 - vx3) + c9o2*(vx1 - vx2 - vx3)*(vx1 - vx2 - vx3) - cusq));
			fTNW = fTNW  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(-vx1 + vx2 + vx3) + c9o2*(-vx1 + vx2 + vx3)*(-vx1 + vx2 + vx3) - cusq));
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







			//////////////////////////////////////////////////////////////////////////                            
			(D.f[E])[k] = fW;
			(D.f[W])[kw] = fE;
			(D.f[N])[k] = fS;
			(D.f[S])[ks] = fN;
			(D.f[T])[k] = fB;
			(D.f[B])[kb] = fT;
			(D.f[NE])[k] = fSW;
			(D.f[SW])[ksw] = fNE;
			(D.f[SE])[ks] = fNW;
			(D.f[NW])[kw] = fSE;
			(D.f[TE])[k] = fBW;
			(D.f[BW])[kbw] = fTE;
			(D.f[BE])[kb] = fTW;
			(D.f[TW])[kw] = fBE;
			(D.f[TN])[k] = fBS;
			(D.f[BS])[kbs] = fTN;
			(D.f[BN])[kb] = fTS;
			(D.f[TS])[ks] = fBN;
			(D.f[dirREST])[k] = fZERO;
			(D.f[TNE])[k] = fBSW;
			(D.f[TSE])[ks] = fBNW;
			(D.f[BNE])[kb] = fTSW;
			(D.f[BSE])[kbs] = fTNW;
			(D.f[TNW])[kw] = fBSE;
			(D.f[TSW])[ksw] = fBNE;
			(D.f[BNW])[kbw] = fTSE;
			(D.f[BSW])[kbsw] = fTNE;
			//////////////////////////////////////////////////////////////////////////                            
		}
	}
}