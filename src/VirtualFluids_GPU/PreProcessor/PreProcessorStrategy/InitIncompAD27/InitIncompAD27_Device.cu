#include "LBM/D3Q27.h"
#include "math.h"
#include "GPU/constant.h"

extern "C" __global__ void LB_Init_Incomp_AD_27(unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int* geoD,
	real* Conc,
	real* ux,
	real* uy,
	real* uz,
	unsigned int size_Mat,
	real* DD27,
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
		BC = geoD[k];

		if (BC != GEO_SOLID && BC != GEO_VOID)
		{
			Distributions27 D27;
			if (EvenOrOdd == true)
			{
				D27.f[dirE] = &DD27[dirE   *size_Mat];
				D27.f[dirW] = &DD27[dirW   *size_Mat];
				D27.f[dirN] = &DD27[dirN   *size_Mat];
				D27.f[dirS] = &DD27[dirS   *size_Mat];
				D27.f[dirT] = &DD27[dirT   *size_Mat];
				D27.f[dirB] = &DD27[dirB   *size_Mat];
				D27.f[dirNE] = &DD27[dirNE  *size_Mat];
				D27.f[dirSW] = &DD27[dirSW  *size_Mat];
				D27.f[dirSE] = &DD27[dirSE  *size_Mat];
				D27.f[dirNW] = &DD27[dirNW  *size_Mat];
				D27.f[dirTE] = &DD27[dirTE  *size_Mat];
				D27.f[dirBW] = &DD27[dirBW  *size_Mat];
				D27.f[dirBE] = &DD27[dirBE  *size_Mat];
				D27.f[dirTW] = &DD27[dirTW  *size_Mat];
				D27.f[dirTN] = &DD27[dirTN  *size_Mat];
				D27.f[dirBS] = &DD27[dirBS  *size_Mat];
				D27.f[dirBN] = &DD27[dirBN  *size_Mat];
				D27.f[dirTS] = &DD27[dirTS  *size_Mat];
				D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
				D27.f[dirTNE] = &DD27[dirTNE *size_Mat];
				D27.f[dirTSW] = &DD27[dirTSW *size_Mat];
				D27.f[dirTSE] = &DD27[dirTSE *size_Mat];
				D27.f[dirTNW] = &DD27[dirTNW *size_Mat];
				D27.f[dirBNE] = &DD27[dirBNE *size_Mat];
				D27.f[dirBSW] = &DD27[dirBSW *size_Mat];
				D27.f[dirBSE] = &DD27[dirBSE *size_Mat];
				D27.f[dirBNW] = &DD27[dirBNW *size_Mat];
			}
			else
			{
				D27.f[dirW] = &DD27[dirE   *size_Mat];
				D27.f[dirE] = &DD27[dirW   *size_Mat];
				D27.f[dirS] = &DD27[dirN   *size_Mat];
				D27.f[dirN] = &DD27[dirS   *size_Mat];
				D27.f[dirB] = &DD27[dirT   *size_Mat];
				D27.f[dirT] = &DD27[dirB   *size_Mat];
				D27.f[dirSW] = &DD27[dirNE  *size_Mat];
				D27.f[dirNE] = &DD27[dirSW  *size_Mat];
				D27.f[dirNW] = &DD27[dirSE  *size_Mat];
				D27.f[dirSE] = &DD27[dirNW  *size_Mat];
				D27.f[dirBW] = &DD27[dirTE  *size_Mat];
				D27.f[dirTE] = &DD27[dirBW  *size_Mat];
				D27.f[dirTW] = &DD27[dirBE  *size_Mat];
				D27.f[dirBE] = &DD27[dirTW  *size_Mat];
				D27.f[dirBS] = &DD27[dirTN  *size_Mat];
				D27.f[dirTN] = &DD27[dirBS  *size_Mat];
				D27.f[dirTS] = &DD27[dirBN  *size_Mat];
				D27.f[dirBN] = &DD27[dirTS  *size_Mat];
				D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
				D27.f[dirBSW] = &DD27[dirTNE *size_Mat];
				D27.f[dirBNE] = &DD27[dirTSW *size_Mat];
				D27.f[dirBNW] = &DD27[dirTSE *size_Mat];
				D27.f[dirBSE] = &DD27[dirTNW *size_Mat];
				D27.f[dirTSW] = &DD27[dirBNE *size_Mat];
				D27.f[dirTNE] = &DD27[dirBSW *size_Mat];
				D27.f[dirTNW] = &DD27[dirBSE *size_Mat];
				D27.f[dirTSE] = &DD27[dirBNW *size_Mat];
			}
			//////////////////////////////////////////////////////////////////////////
			real ConcD = Conc[k];
			real   vx1 = ux[k];
			real   vx2 = uy[k];
			real   vx3 = uz[k];
			//real lambdaD     = -three + sqrt(three);
			//real Diffusivity = c1o20;
			//real Lam         = -(c1o2+one/lambdaD);
			//real nue_d       = Lam/three;
			//real ae          = Diffusivity/nue_d - one;
			//real ux_sq       = vx1 * vx1;
			//real uy_sq       = vx2 * vx2;
			//real uz_sq       = vx3 * vx3;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//D3Q7
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//index
			//unsigned int kzero= k;
			//unsigned int ke   = k;
			//unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			//unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			//unsigned int kb   = neighborZ[k];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//(D7.f[0])[kzero] = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
			//(D7.f[1])[ke   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
			//(D7.f[2])[kw   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
			//(D7.f[3])[kn   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
			//(D7.f[4])[ks   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
			//(D7.f[5])[kt   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
			//(D7.f[6])[kb   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//D3Q27
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real cu_sq = c3o2*(vx1*vx1 + vx2*vx2 + vx3*vx3);

			(D27.f[dirZERO])[kzero] = c8over27* ConcD*(one - cu_sq);
			(D27.f[dirE])[ke] = c2over27* ConcD*(one + three*(vx1)+c9over2*(vx1)*(vx1)-cu_sq);
			(D27.f[dirW])[kw] = c2over27* ConcD*(one + three*(-vx1) + c9over2*(-vx1)*(-vx1) - cu_sq);
			(D27.f[dirN])[kn] = c2over27* ConcD*(one + three*(vx2)+c9over2*(vx2)*(vx2)-cu_sq);
			(D27.f[dirS])[ks] = c2over27* ConcD*(one + three*(-vx2) + c9over2*(-vx2)*(-vx2) - cu_sq);
			(D27.f[dirT])[kt] = c2over27* ConcD*(one + three*(vx3)+c9over2*(vx3)*(vx3)-cu_sq);
			(D27.f[dirB])[kb] = c2over27* ConcD*(one + three*(-vx3) + c9over2*(-vx3)*(-vx3) - cu_sq);
			(D27.f[dirNE])[kne] = c1over54* ConcD*(one + three*(vx1 + vx2) + c9over2*(vx1 + vx2)*(vx1 + vx2) - cu_sq);
			(D27.f[dirSW])[ksw] = c1over54* ConcD*(one + three*(-vx1 - vx2) + c9over2*(-vx1 - vx2)*(-vx1 - vx2) - cu_sq);
			(D27.f[dirSE])[kse] = c1over54* ConcD*(one + three*(vx1 - vx2) + c9over2*(vx1 - vx2)*(vx1 - vx2) - cu_sq);
			(D27.f[dirNW])[knw] = c1over54* ConcD*(one + three*(-vx1 + vx2) + c9over2*(-vx1 + vx2)*(-vx1 + vx2) - cu_sq);
			(D27.f[dirTE])[kte] = c1over54* ConcD*(one + three*(vx1 + vx3) + c9over2*(vx1 + vx3)*(vx1 + vx3) - cu_sq);
			(D27.f[dirBW])[kbw] = c1over54* ConcD*(one + three*(-vx1 - vx3) + c9over2*(-vx1 - vx3)*(-vx1 - vx3) - cu_sq);
			(D27.f[dirBE])[kbe] = c1over54* ConcD*(one + three*(vx1 - vx3) + c9over2*(vx1 - vx3)*(vx1 - vx3) - cu_sq);
			(D27.f[dirTW])[ktw] = c1over54* ConcD*(one + three*(-vx1 + vx3) + c9over2*(-vx1 + vx3)*(-vx1 + vx3) - cu_sq);
			(D27.f[dirTN])[ktn] = c1over54* ConcD*(one + three*(vx2 + vx3) + c9over2*(vx2 + vx3)*(vx2 + vx3) - cu_sq);
			(D27.f[dirBS])[kbs] = c1over54* ConcD*(one + three*(-vx2 - vx3) + c9over2*(-vx2 - vx3)*(-vx2 - vx3) - cu_sq);
			(D27.f[dirBN])[kbn] = c1over54* ConcD*(one + three*(vx2 - vx3) + c9over2*(vx2 - vx3)*(vx2 - vx3) - cu_sq);
			(D27.f[dirTS])[kts] = c1over54* ConcD*(one + three*(-vx2 + vx3) + c9over2*(-vx2 + vx3)*(-vx2 + vx3) - cu_sq);
			(D27.f[dirTNE])[ktne] = c1over216*ConcD*(one + three*(vx1 + vx2 + vx3) + c9over2*(vx1 + vx2 + vx3)*(vx1 + vx2 + vx3) - cu_sq);
			(D27.f[dirBSW])[kbsw] = c1over216*ConcD*(one + three*(-vx1 - vx2 - vx3) + c9over2*(-vx1 - vx2 - vx3)*(-vx1 - vx2 - vx3) - cu_sq);
			(D27.f[dirBNE])[kbne] = c1over216*ConcD*(one + three*(vx1 + vx2 - vx3) + c9over2*(vx1 + vx2 - vx3)*(vx1 + vx2 - vx3) - cu_sq);
			(D27.f[dirTSW])[ktsw] = c1over216*ConcD*(one + three*(-vx1 - vx2 + vx3) + c9over2*(-vx1 - vx2 + vx3)*(-vx1 - vx2 + vx3) - cu_sq);
			(D27.f[dirTSE])[ktse] = c1over216*ConcD*(one + three*(vx1 - vx2 + vx3) + c9over2*(vx1 - vx2 + vx3)*(vx1 - vx2 + vx3) - cu_sq);
			(D27.f[dirBNW])[kbnw] = c1over216*ConcD*(one + three*(-vx1 + vx2 - vx3) + c9over2*(-vx1 + vx2 - vx3)*(-vx1 + vx2 - vx3) - cu_sq);
			(D27.f[dirBSE])[kbse] = c1over216*ConcD*(one + three*(vx1 - vx2 - vx3) + c9over2*(vx1 - vx2 - vx3)*(vx1 - vx2 - vx3) - cu_sq);
			(D27.f[dirTNW])[ktnw] = c1over216*ConcD*(one + three*(-vx1 + vx2 + vx3) + c9over2*(-vx1 + vx2 + vx3)*(-vx1 + vx2 + vx3) - cu_sq);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}
}