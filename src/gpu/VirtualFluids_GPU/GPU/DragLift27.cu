/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////
__global__ void DragLiftPost27(  real* DD, 
											int* k_Q, 
											real* QQ,
											int numberOfBCnodes, 
											double *DragX,
											double *DragY,
											double *DragZ,
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned long long numberOfLBnodes, 
											bool isEvenTimestep)
{
	Distributions27 D;
	if (isEvenTimestep==true)
	{
		D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
		D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
		D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
		D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
		D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
		D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
		D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
		D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
		D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
		D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
		D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
		D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
		D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
		D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
		D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
		D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
		D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
		D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
		D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
		D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
		D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
		D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
		D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
		D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
		D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
		D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
		D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
	} 
	else
	{
		D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
		D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
		D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
		D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
		D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
		D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
		D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
		D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
		D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
		D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
		D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
		D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
		D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
		D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
		D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
		D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
		D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
		D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
		D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
		D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
		D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
		D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
		D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
		D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
		D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
		D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
		D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
	}
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if(k<numberOfBCnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
			*q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
			*q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
			*q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			*q_dirBSE, *q_dirBNW; 
		q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
		q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
		q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
		q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
		q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
		q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
		q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
		q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
		q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
		q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
		q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
		q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
		q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
		q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
		q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
		q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
		q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
		q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
		q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
		q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
		q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
		q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
		q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
		q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
		q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
		q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
		////////////////////////////////////////////////////////////////////////////////
		//index
		unsigned int KQK  = k_Q[k];
		//unsigned int kzero= KQK;
		unsigned int ke   = KQK;
		unsigned int kw   = neighborX[KQK];
		unsigned int kn   = KQK;
		unsigned int ks   = neighborY[KQK];
		unsigned int kt   = KQK;
		unsigned int kb   = neighborZ[KQK];
		unsigned int ksw  = neighborY[kw];
		unsigned int kne  = KQK;
		unsigned int kse  = ks;
		unsigned int knw  = kw;
		unsigned int kbw  = neighborZ[kw];
		unsigned int kte  = KQK;
		unsigned int kbe  = kb;
		unsigned int ktw  = kw;
		unsigned int kbs  = neighborZ[ks];
		unsigned int ktn  = KQK;
		unsigned int kbn  = kb;
		unsigned int kts  = ks;
		unsigned int ktse = ks;
		unsigned int kbnw = kbw;
		unsigned int ktnw = kw;
		unsigned int kbse = kbs;
		unsigned int ktsw = ksw;
		unsigned int kbne = kb;
		unsigned int ktne = KQK;
		unsigned int kbsw = neighborZ[ksw];
		////////////////////////////////////////////////////////////////////////////////
		real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
                f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

		f_W    = (D.f[DIR_P00])[ke   ];
		f_E    = (D.f[DIR_M00])[kw   ];
		f_S    = (D.f[DIR_0P0])[kn   ];
		f_N    = (D.f[DIR_0M0])[ks   ];
		f_B    = (D.f[DIR_00P])[kt   ];
		f_T    = (D.f[DIR_00M])[kb   ];
		f_SW   = (D.f[DIR_PP0])[kne  ];
		f_NE   = (D.f[DIR_MM0])[ksw  ];
		f_NW   = (D.f[DIR_PM0])[kse  ];
		f_SE   = (D.f[DIR_MP0])[knw  ];
		f_BW   = (D.f[DIR_P0P])[kte  ];
		f_TE   = (D.f[DIR_M0M])[kbw  ];
		f_TW   = (D.f[DIR_P0M])[kbe  ];
		f_BE   = (D.f[DIR_M0P])[ktw  ];
		f_BS   = (D.f[DIR_0PP])[ktn  ];
		f_TN   = (D.f[DIR_0MM])[kbs  ];
		f_TS   = (D.f[DIR_0PM])[kbn  ];
		f_BN   = (D.f[DIR_0MP])[kts  ];
		f_BSW  = (D.f[DIR_PPP])[ktne ];
		f_BNE  = (D.f[DIR_MMP])[ktsw ];
		f_BNW  = (D.f[DIR_PMP])[ktse ];
		f_BSE  = (D.f[DIR_MPP])[ktnw ];
		f_TSW  = (D.f[DIR_PPM])[kbne ];
		f_TNE  = (D.f[DIR_MMM])[kbsw ];
		f_TNW  = (D.f[DIR_PMM])[kbse ];
		f_TSE  = (D.f[DIR_MPM])[kbnw ];
		////////////////////////////////////////////////////////////////////////////////
		double	OnE   = c0o1, OnW   = c0o1, OnN   = c0o1, OnS   = c0o1, OnT = c0o1, OnB = c0o1, 
				OnNE  = c0o1, OnSW  = c0o1, OnSE  = c0o1, OnNW  = c0o1, 
				OnTE  = c0o1, OnBW  = c0o1, OnBE  = c0o1, OnTW  = c0o1,
				OnTN  = c0o1, OnBS  = c0o1, OnBN  = c0o1, OnTS  = c0o1, 
				OnTNE = c0o1, OnTSW = c0o1, OnTSE = c0o1, OnTNW = c0o1, 
				OnBNE = c0o1, OnBSW = c0o1, OnBSE = c0o1, OnBNW = c0o1;
		////////////////////////////////////////////////////////////////////////////////
		real q;
		q = q_dirE[k];		if (q>=c0o1 && q<=c1o1) OnE   = c1o1;
		q = q_dirW[k];		if (q>=c0o1 && q<=c1o1) OnW   = c1o1;
		q = q_dirN[k];		if (q>=c0o1 && q<=c1o1) OnN   = c1o1;
		q = q_dirS[k];		if (q>=c0o1 && q<=c1o1) OnS   = c1o1;
		q = q_dirT[k];		if (q>=c0o1 && q<=c1o1) OnT   = c1o1;
		q = q_dirB[k];		if (q>=c0o1 && q<=c1o1) OnB   = c1o1;
		q = q_dirNE[k];		if (q>=c0o1 && q<=c1o1) OnNE  = c1o1;
		q = q_dirSW[k];		if (q>=c0o1 && q<=c1o1) OnSW  = c1o1;
		q = q_dirSE[k];		if (q>=c0o1 && q<=c1o1) OnSE  = c1o1;
		q = q_dirNW[k];		if (q>=c0o1 && q<=c1o1) OnNW  = c1o1;
		q = q_dirTE[k];		if (q>=c0o1 && q<=c1o1) OnTE  = c1o1;
		q = q_dirBW[k];		if (q>=c0o1 && q<=c1o1) OnBW  = c1o1;
		q = q_dirBE[k];		if (q>=c0o1 && q<=c1o1) OnBE  = c1o1;
		q = q_dirTW[k];		if (q>=c0o1 && q<=c1o1) OnTW  = c1o1;
		q = q_dirTN[k];		if (q>=c0o1 && q<=c1o1) OnTN  = c1o1;
		q = q_dirBS[k];		if (q>=c0o1 && q<=c1o1) OnBS  = c1o1;
		q = q_dirBN[k];		if (q>=c0o1 && q<=c1o1) OnBN  = c1o1;
		q = q_dirTS[k];		if (q>=c0o1 && q<=c1o1) OnTS  = c1o1;
		q = q_dirTNE[k];	if (q>=c0o1 && q<=c1o1) OnTNE = c1o1;
		q = q_dirBSW[k];	if (q>=c0o1 && q<=c1o1) OnBSW = c1o1;
		q = q_dirBNE[k];	if (q>=c0o1 && q<=c1o1) OnBNE = c1o1;
		q = q_dirTSW[k];	if (q>=c0o1 && q<=c1o1) OnTSW = c1o1;
		q = q_dirTSE[k];	if (q>=c0o1 && q<=c1o1) OnTSE = c1o1;
		q = q_dirBNW[k];	if (q>=c0o1 && q<=c1o1) OnBNW = c1o1;
		q = q_dirBSE[k];	if (q>=c0o1 && q<=c1o1) OnBSE = c1o1;
		q = q_dirTNW[k];	if (q>=c0o1 && q<=c1o1) OnTNW = c1o1;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double dragX, dragY, dragZ;

		dragX = ((((f_BSE * OnBSE) - (f_TNW * OnTNW))   + 
				  ((f_BNE * OnBNE) - (f_TSW * OnTSW))   + 
				  ((f_TSE * OnTSE) - (f_BNW * OnBNW))   + 
				  ((f_TNE * OnTNE) - (f_BSW * OnBSW)))  + 
				 (((f_BE  * OnBE ) - (f_TW  * OnTW ))   + 
				  ((f_TE  * OnTE ) - (f_BW  * OnBW ))   + 
				  ((f_SE  * OnSE ) - (f_NW  * OnNW ))   + 
				  ((f_NE  * OnNE ) - (f_SW  * OnSW )))) + 
				  ((f_E   * OnE  ) - (f_W   * OnW  ));

		dragY = ((((f_TNW * OnTNW) - (f_BSE * OnBSE))   + 
				  ((f_BNE * OnBNE) - (f_TSW * OnTSW))   + 
				  ((f_BNW * OnBNW) - (f_TSE * OnTSE))   + 
				  ((f_TNE * OnTNE) - (f_BSW * OnBSW)))  + 
				 (((f_BN  * OnBN ) - (f_TS  * OnTS ))   + 
				  ((f_TN  * OnTN ) - (f_BS  * OnBS ))   + 
				  ((f_NW  * OnNW ) - (f_SE  * OnSE ))   + 
				  ((f_NE  * OnNE ) - (f_SW  * OnSW )))) + 
				  ((f_N   * OnN  ) - (f_S   * OnS  ));

		dragZ = ((((f_TNW * OnTNW) - (f_BSE * OnBSE))   + 
				  ((f_TSW * OnTSW) - (f_BNE * OnBNE))   + 
				  ((f_TSE * OnTSE) - (f_BNW * OnBNW))   + 
				  ((f_TNE * OnTNE) - (f_BSW * OnBSW)))  + 
				 (((f_TS  * OnTS ) - (f_BN  * OnBN ))   + 
				  ((f_TN  * OnTN ) - (f_BS  * OnBS ))   + 
				  ((f_TW  * OnTW ) - (f_BE  * OnBE ))   + 
				  ((f_TE  * OnTE ) - (f_BW  * OnBW )))) + 
				  ((f_T   * OnT  ) - (f_B   * OnB  ));
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		DragX[k] = -dragX;
		DragY[k] = -dragY;
		DragZ[k] = -dragZ;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}
}











////////////////////////////////////////////////////////////////////////////////
__global__ void DragLiftPre27(   real* DD, 
											int* k_Q, 
											real* QQ,
											int numberOfBCnodes, 
											double *DragX,
											double *DragY,
											double *DragZ,
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned long long numberOfLBnodes, 
											bool isEvenTimestep)
{
	Distributions27 D;
	if (isEvenTimestep==true)
	{
		D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
		D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
		D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
		D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
		D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
		D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
		D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
		D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
		D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
		D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
		D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
		D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
		D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
		D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
		D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
		D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
		D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
		D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
		D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
		D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
		D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
		D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
		D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
		D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
		D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
		D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
		D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
	} 
	else
	{
		D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
		D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
		D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
		D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
		D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
		D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
		D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
		D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
		D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
		D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
		D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
		D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
		D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
		D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
		D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
		D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
		D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
		D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
		D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
		D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
		D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
		D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
		D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
		D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
		D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
		D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
		D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
	}
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if(k<numberOfBCnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
			*q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
			*q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
			*q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			*q_dirBSE, *q_dirBNW; 
		q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
		q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
		q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
		q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
		q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
		q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
		q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
		q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
		q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
		q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
		q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
		q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
		q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
		q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
		q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
		q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
		q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
		q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
		q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
		q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
		q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
		q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
		q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
		q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
		q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
		q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
		////////////////////////////////////////////////////////////////////////////////
		//index
		unsigned int KQK  = k_Q[k];
		//unsigned int kzero= KQK;
		unsigned int ke   = KQK;
		unsigned int kw   = neighborX[KQK];
		unsigned int kn   = KQK;
		unsigned int ks   = neighborY[KQK];
		unsigned int kt   = KQK;
		unsigned int kb   = neighborZ[KQK];
		unsigned int ksw  = neighborY[kw];
		unsigned int kne  = KQK;
		unsigned int kse  = ks;
		unsigned int knw  = kw;
		unsigned int kbw  = neighborZ[kw];
		unsigned int kte  = KQK;
		unsigned int kbe  = kb;
		unsigned int ktw  = kw;
		unsigned int kbs  = neighborZ[ks];
		unsigned int ktn  = KQK;
		unsigned int kbn  = kb;
		unsigned int kts  = ks;
		unsigned int ktse = ks;
		unsigned int kbnw = kbw;
		unsigned int ktnw = kw;
		unsigned int kbse = kbs;
		unsigned int ktsw = ksw;
		unsigned int kbne = kb;
		unsigned int ktne = KQK;
		unsigned int kbsw = neighborZ[ksw];
		////////////////////////////////////////////////////////////////////////////////
		real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
                f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

		f_E   = (D.f[DIR_P00])[ke   ];
		f_W   = (D.f[DIR_M00])[kw   ];
		f_N   = (D.f[DIR_0P0])[kn   ];
		f_S   = (D.f[DIR_0M0])[ks   ];
		f_T   = (D.f[DIR_00P])[kt   ];
		f_B   = (D.f[DIR_00M])[kb   ];
		f_NE  = (D.f[DIR_PP0])[kne  ];
		f_SW  = (D.f[DIR_MM0])[ksw  ];
		f_SE  = (D.f[DIR_PM0])[kse  ];
		f_NW  = (D.f[DIR_MP0])[knw  ];
		f_TE  = (D.f[DIR_P0P])[kte  ];
		f_BW  = (D.f[DIR_M0M])[kbw  ];
		f_BE  = (D.f[DIR_P0M])[kbe  ];
		f_TW  = (D.f[DIR_M0P])[ktw  ];
		f_TN  = (D.f[DIR_0PP])[ktn  ];
		f_BS  = (D.f[DIR_0MM])[kbs  ];
		f_BN  = (D.f[DIR_0PM])[kbn  ];
		f_TS  = (D.f[DIR_0MP])[kts  ];
		f_TNE = (D.f[DIR_PPP])[ktne ];
		f_TSW = (D.f[DIR_MMP])[ktsw ];
		f_TSE = (D.f[DIR_PMP])[ktse ];
		f_TNW = (D.f[DIR_MPP])[ktnw ];
		f_BNE = (D.f[DIR_PPM])[kbne ];
		f_BSW = (D.f[DIR_MMM])[kbsw ];
		f_BSE = (D.f[DIR_PMM])[kbse ];
		f_BNW = (D.f[DIR_MPM])[kbnw ];
		 ////////////////////////////////////////////////////////////////////////////////
		double	OnE   = c0o1, OnW   = c0o1, OnN   = c0o1, OnS   = c0o1, OnT = c0o1, OnB = c0o1, 
				OnNE  = c0o1, OnSW  = c0o1, OnSE  = c0o1, OnNW  = c0o1, 
				OnTE  = c0o1, OnBW  = c0o1, OnBE  = c0o1, OnTW  = c0o1,
				OnTN  = c0o1, OnBS  = c0o1, OnBN  = c0o1, OnTS  = c0o1, 
				OnTNE = c0o1, OnTSW = c0o1, OnTSE = c0o1, OnTNW = c0o1, 
				OnBNE = c0o1, OnBSW = c0o1, OnBSE = c0o1, OnBNW = c0o1;
		////////////////////////////////////////////////////////////////////////////////
		real q;
		q = q_dirE[k];		if (q>=c0o1 && q<=c1o1) OnW   = c1o1;
		q = q_dirW[k];		if (q>=c0o1 && q<=c1o1) OnE   = c1o1;
		q = q_dirN[k];		if (q>=c0o1 && q<=c1o1) OnS   = c1o1;
		q = q_dirS[k];		if (q>=c0o1 && q<=c1o1) OnN   = c1o1;
		q = q_dirT[k];		if (q>=c0o1 && q<=c1o1) OnB   = c1o1;
		q = q_dirB[k];		if (q>=c0o1 && q<=c1o1) OnT   = c1o1;
		q = q_dirNE[k];		if (q>=c0o1 && q<=c1o1) OnSW  = c1o1;
		q = q_dirSW[k];		if (q>=c0o1 && q<=c1o1) OnNE  = c1o1;
		q = q_dirSE[k];		if (q>=c0o1 && q<=c1o1) OnNW  = c1o1;
		q = q_dirNW[k];		if (q>=c0o1 && q<=c1o1) OnSE  = c1o1;
		q = q_dirTE[k];		if (q>=c0o1 && q<=c1o1) OnBW  = c1o1;
		q = q_dirBW[k];		if (q>=c0o1 && q<=c1o1) OnTE  = c1o1;
		q = q_dirBE[k];		if (q>=c0o1 && q<=c1o1) OnTW  = c1o1;
		q = q_dirTW[k];		if (q>=c0o1 && q<=c1o1) OnBE  = c1o1;
		q = q_dirTN[k];		if (q>=c0o1 && q<=c1o1) OnBS  = c1o1;
		q = q_dirBS[k];		if (q>=c0o1 && q<=c1o1) OnTN  = c1o1;
		q = q_dirBN[k];		if (q>=c0o1 && q<=c1o1) OnTS  = c1o1;
		q = q_dirTS[k];		if (q>=c0o1 && q<=c1o1) OnBN  = c1o1;
		q = q_dirTNE[k];	if (q>=c0o1 && q<=c1o1) OnBSW = c1o1;
		q = q_dirBSW[k];	if (q>=c0o1 && q<=c1o1) OnTNE = c1o1;
		q = q_dirBNE[k];	if (q>=c0o1 && q<=c1o1) OnTSW = c1o1;
		q = q_dirTSW[k];	if (q>=c0o1 && q<=c1o1) OnBNE = c1o1;
		q = q_dirTSE[k];	if (q>=c0o1 && q<=c1o1) OnBNW = c1o1;
		q = q_dirBNW[k];	if (q>=c0o1 && q<=c1o1) OnTSE = c1o1;
		q = q_dirBSE[k];	if (q>=c0o1 && q<=c1o1) OnTNW = c1o1;
		q = q_dirTNW[k];	if (q>=c0o1 && q<=c1o1) OnBSE = c1o1;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double dragX, dragY, dragZ;

		dragX = ((((f_BSE * OnBSE) - (f_TNW * OnTNW))   + 
				  ((f_BNE * OnBNE) - (f_TSW * OnTSW))   + 
				  ((f_TSE * OnTSE) - (f_BNW * OnBNW))   + 
				  ((f_TNE * OnTNE) - (f_BSW * OnBSW)))  + 
				 (((f_BE  * OnBE ) - (f_TW  * OnTW ))   + 
				  ((f_TE  * OnTE ) - (f_BW  * OnBW ))   + 
				  ((f_SE  * OnSE ) - (f_NW  * OnNW ))   + 
				  ((f_NE  * OnNE ) - (f_SW  * OnSW )))) + 
				  ((f_E   * OnE  ) - (f_W   * OnW  ));

		dragY = ((((f_TNW * OnTNW) - (f_BSE * OnBSE))   + 
				  ((f_BNE * OnBNE) - (f_TSW * OnTSW))   + 
				  ((f_BNW * OnBNW) - (f_TSE * OnTSE))   + 
				  ((f_TNE * OnTNE) - (f_BSW * OnBSW)))  + 
				 (((f_BN  * OnBN ) - (f_TS  * OnTS ))   + 
				  ((f_TN  * OnTN ) - (f_BS  * OnBS ))   + 
				  ((f_NW  * OnNW ) - (f_SE  * OnSE ))   + 
				  ((f_NE  * OnNE ) - (f_SW  * OnSW )))) + 
				  ((f_N   * OnN  ) - (f_S   * OnS  ));

		dragZ = ((((f_TNW * OnTNW) - (f_BSE * OnBSE))   + 
				  ((f_TSW * OnTSW) - (f_BNE * OnBNE))   + 
				  ((f_TSE * OnTSE) - (f_BNW * OnBNW))   + 
				  ((f_TNE * OnTNE) - (f_BSW * OnBSW)))  + 
				 (((f_TS  * OnTS ) - (f_BN  * OnBN ))   + 
				  ((f_TN  * OnTN ) - (f_BS  * OnBS ))   + 
				  ((f_TW  * OnTW ) - (f_BE  * OnBE ))   + 
				  ((f_TE  * OnTE ) - (f_BW  * OnBW )))) + 
				  ((f_T   * OnT  ) - (f_B   * OnB  ));
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		DragX[k] = -dragX;
		DragY[k] = -dragY;
		DragZ[k] = -dragZ;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}
}
