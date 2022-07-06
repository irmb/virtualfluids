#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
#include "math.h"

extern "C" __global__ void LB_Kernel_BGK_Incomp_SP_27(real omega,
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
				D.f[dirREST] = &DDStart[dirREST*size_Mat];
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
				D.f[dirREST] = &DDStart[dirREST*size_Mat];
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
			real fE = (D.f[dirE])[k];//ke
			real fW = (D.f[dirW])[kw];
			real fN = (D.f[dirN])[k];//kn
			real fS = (D.f[dirS])[ks];
			real fT = (D.f[dirT])[k];//kt
			real fB = (D.f[dirB])[kb];
			real fNE = (D.f[dirNE])[k];//kne
			real fSW = (D.f[dirSW])[ksw];
			real fSE = (D.f[dirSE])[ks];//kse
			real fNW = (D.f[dirNW])[kw];//knw
			real fTE = (D.f[dirTE])[k];//kte
			real fBW = (D.f[dirBW])[kbw];
			real fBE = (D.f[dirBE])[kb];//kbe
			real fTW = (D.f[dirTW])[kw];//ktw
			real fTN = (D.f[dirTN])[k];//ktn
			real fBS = (D.f[dirBS])[kbs];
			real fBN = (D.f[dirBN])[kb];//kbn
			real fTS = (D.f[dirTS])[ks];//kts
			real fZERO = (D.f[dirREST])[k];//kzero
			real fTNE = (D.f[dirTNE])[k];//ktne
			real fTSW = (D.f[dirTSW])[ksw];//ktsw
			real fTSE = (D.f[dirTSE])[ks];//ktse
			real fTNW = (D.f[dirTNW])[kw];//ktnw
			real fBNE = (D.f[dirBNE])[kb];//kbne
			real fBSW = (D.f[dirBSW])[kbsw];
			real fBSE = (D.f[dirBSE])[kbs];//kbse
			real fBNW = (D.f[dirBNW])[kbw];//kbnw
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
			(D.f[dirE])[k] = fW;
			(D.f[dirW])[kw] = fE;
			(D.f[dirN])[k] = fS;
			(D.f[dirS])[ks] = fN;
			(D.f[dirT])[k] = fB;
			(D.f[dirB])[kb] = fT;
			(D.f[dirNE])[k] = fSW;
			(D.f[dirSW])[ksw] = fNE;
			(D.f[dirSE])[ks] = fNW;
			(D.f[dirNW])[kw] = fSE;
			(D.f[dirTE])[k] = fBW;
			(D.f[dirBW])[kbw] = fTE;
			(D.f[dirBE])[kb] = fTW;
			(D.f[dirTW])[kw] = fBE;
			(D.f[dirTN])[k] = fBS;
			(D.f[dirBS])[kbs] = fTN;
			(D.f[dirBN])[kb] = fTS;
			(D.f[dirTS])[ks] = fBN;
			(D.f[dirREST])[k] = fZERO;
			(D.f[dirTNE])[k] = fBSW;
			(D.f[dirTSE])[ks] = fBNW;
			(D.f[dirBNE])[kb] = fTSW;
			(D.f[dirBSE])[kbs] = fTNW;
			(D.f[dirTNW])[kw] = fBSE;
			(D.f[dirTSW])[ksw] = fBNE;
			(D.f[dirBNW])[kbw] = fTSE;
			(D.f[dirBSW])[kbsw] = fTNE;
			//////////////////////////////////////////////////////////////////////////                            
		}
	}
}