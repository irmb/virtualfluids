/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;

////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void DragLiftPost27(  real* DD, 
											int* k_Q, 
											real* QQ,
											int numberOfBCnodes, 
											double *DragX,
											double *DragY,
											double *DragZ,
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int size_Mat, 
											bool isEvenTimestep)
{
	Distributions27 D;
	if (isEvenTimestep==true)
	{
		D.f[E   ] = &DD[E   *size_Mat];
		D.f[W   ] = &DD[W   *size_Mat];
		D.f[N   ] = &DD[N   *size_Mat];
		D.f[S   ] = &DD[S   *size_Mat];
		D.f[T   ] = &DD[T   *size_Mat];
		D.f[B   ] = &DD[B   *size_Mat];
		D.f[NE  ] = &DD[NE  *size_Mat];
		D.f[SW  ] = &DD[SW  *size_Mat];
		D.f[SE  ] = &DD[SE  *size_Mat];
		D.f[NW  ] = &DD[NW  *size_Mat];
		D.f[TE  ] = &DD[TE  *size_Mat];
		D.f[BW  ] = &DD[BW  *size_Mat];
		D.f[BE  ] = &DD[BE  *size_Mat];
		D.f[TW  ] = &DD[TW  *size_Mat];
		D.f[TN  ] = &DD[TN  *size_Mat];
		D.f[BS  ] = &DD[BS  *size_Mat];
		D.f[BN  ] = &DD[BN  *size_Mat];
		D.f[TS  ] = &DD[TS  *size_Mat];
		D.f[REST] = &DD[REST*size_Mat];
		D.f[TNE ] = &DD[TNE *size_Mat];
		D.f[TSW ] = &DD[TSW *size_Mat];
		D.f[TSE ] = &DD[TSE *size_Mat];
		D.f[TNW ] = &DD[TNW *size_Mat];
		D.f[BNE ] = &DD[BNE *size_Mat];
		D.f[BSW ] = &DD[BSW *size_Mat];
		D.f[BSE ] = &DD[BSE *size_Mat];
		D.f[BNW ] = &DD[BNW *size_Mat];
	} 
	else
	{
		D.f[W   ] = &DD[E   *size_Mat];
		D.f[E   ] = &DD[W   *size_Mat];
		D.f[S   ] = &DD[N   *size_Mat];
		D.f[N   ] = &DD[S   *size_Mat];
		D.f[B   ] = &DD[T   *size_Mat];
		D.f[T   ] = &DD[B   *size_Mat];
		D.f[SW  ] = &DD[NE  *size_Mat];
		D.f[NE  ] = &DD[SW  *size_Mat];
		D.f[NW  ] = &DD[SE  *size_Mat];
		D.f[SE  ] = &DD[NW  *size_Mat];
		D.f[BW  ] = &DD[TE  *size_Mat];
		D.f[TE  ] = &DD[BW  *size_Mat];
		D.f[TW  ] = &DD[BE  *size_Mat];
		D.f[BE  ] = &DD[TW  *size_Mat];
		D.f[BS  ] = &DD[TN  *size_Mat];
		D.f[TN  ] = &DD[BS  *size_Mat];
		D.f[TS  ] = &DD[BN  *size_Mat];
		D.f[BN  ] = &DD[TS  *size_Mat];
		D.f[REST] = &DD[REST*size_Mat];
		D.f[TNE ] = &DD[BSW *size_Mat];
		D.f[TSW ] = &DD[BNE *size_Mat];
		D.f[TSE ] = &DD[BNW *size_Mat];
		D.f[TNW ] = &DD[BSE *size_Mat];
		D.f[BNE ] = &DD[TSW *size_Mat];
		D.f[BSW ] = &DD[TNE *size_Mat];
		D.f[BSE ] = &DD[TNW *size_Mat];
		D.f[BNW ] = &DD[TSE *size_Mat];
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
		q_dirE   = &QQ[E   * numberOfBCnodes];
		q_dirW   = &QQ[W   * numberOfBCnodes];
		q_dirN   = &QQ[N   * numberOfBCnodes];
		q_dirS   = &QQ[S   * numberOfBCnodes];
		q_dirT   = &QQ[T   * numberOfBCnodes];
		q_dirB   = &QQ[B   * numberOfBCnodes];
		q_dirNE  = &QQ[NE  * numberOfBCnodes];
		q_dirSW  = &QQ[SW  * numberOfBCnodes];
		q_dirSE  = &QQ[SE  * numberOfBCnodes];
		q_dirNW  = &QQ[NW  * numberOfBCnodes];
		q_dirTE  = &QQ[TE  * numberOfBCnodes];
		q_dirBW  = &QQ[BW  * numberOfBCnodes];
		q_dirBE  = &QQ[BE  * numberOfBCnodes];
		q_dirTW  = &QQ[TW  * numberOfBCnodes];
		q_dirTN  = &QQ[TN  * numberOfBCnodes];
		q_dirBS  = &QQ[BS  * numberOfBCnodes];
		q_dirBN  = &QQ[BN  * numberOfBCnodes];
		q_dirTS  = &QQ[TS  * numberOfBCnodes];
		q_dirTNE = &QQ[TNE * numberOfBCnodes];
		q_dirTSW = &QQ[TSW * numberOfBCnodes];
		q_dirTSE = &QQ[TSE * numberOfBCnodes];
		q_dirTNW = &QQ[TNW * numberOfBCnodes];
		q_dirBNE = &QQ[BNE * numberOfBCnodes];
		q_dirBSW = &QQ[BSW * numberOfBCnodes];
		q_dirBSE = &QQ[BSE * numberOfBCnodes];
		q_dirBNW = &QQ[BNW * numberOfBCnodes];
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

		f_W    = (D.f[E   ])[ke   ];
		f_E    = (D.f[W   ])[kw   ];
		f_S    = (D.f[N   ])[kn   ];
		f_N    = (D.f[S   ])[ks   ];
		f_B    = (D.f[T   ])[kt   ];
		f_T    = (D.f[B   ])[kb   ];
		f_SW   = (D.f[NE  ])[kne  ];
		f_NE   = (D.f[SW  ])[ksw  ];
		f_NW   = (D.f[SE  ])[kse  ];
		f_SE   = (D.f[NW  ])[knw  ];
		f_BW   = (D.f[TE  ])[kte  ];
		f_TE   = (D.f[BW  ])[kbw  ];
		f_TW   = (D.f[BE  ])[kbe  ];
		f_BE   = (D.f[TW  ])[ktw  ];
		f_BS   = (D.f[TN  ])[ktn  ];
		f_TN   = (D.f[BS  ])[kbs  ];
		f_TS   = (D.f[BN  ])[kbn  ];
		f_BN   = (D.f[TS  ])[kts  ];
		f_BSW  = (D.f[TNE ])[ktne ];
		f_BNE  = (D.f[TSW ])[ktsw ];
		f_BNW  = (D.f[TSE ])[ktse ];
		f_BSE  = (D.f[TNW ])[ktnw ];
		f_TSW  = (D.f[BNE ])[kbne ];
		f_TNE  = (D.f[BSW ])[kbsw ];
		f_TNW  = (D.f[BSE ])[kbse ];
		f_TSE  = (D.f[BNW ])[kbnw ];
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
extern "C" __global__ void DragLiftPre27(   real* DD, 
											int* k_Q, 
											real* QQ,
											int numberOfBCnodes, 
											double *DragX,
											double *DragY,
											double *DragZ,
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int size_Mat, 
											bool isEvenTimestep)
{
	Distributions27 D;
	if (isEvenTimestep==true)
	{
		D.f[E   ] = &DD[E   *size_Mat];
		D.f[W   ] = &DD[W   *size_Mat];
		D.f[N   ] = &DD[N   *size_Mat];
		D.f[S   ] = &DD[S   *size_Mat];
		D.f[T   ] = &DD[T   *size_Mat];
		D.f[B   ] = &DD[B   *size_Mat];
		D.f[NE  ] = &DD[NE  *size_Mat];
		D.f[SW  ] = &DD[SW  *size_Mat];
		D.f[SE  ] = &DD[SE  *size_Mat];
		D.f[NW  ] = &DD[NW  *size_Mat];
		D.f[TE  ] = &DD[TE  *size_Mat];
		D.f[BW  ] = &DD[BW  *size_Mat];
		D.f[BE  ] = &DD[BE  *size_Mat];
		D.f[TW  ] = &DD[TW  *size_Mat];
		D.f[TN  ] = &DD[TN  *size_Mat];
		D.f[BS  ] = &DD[BS  *size_Mat];
		D.f[BN  ] = &DD[BN  *size_Mat];
		D.f[TS  ] = &DD[TS  *size_Mat];
		D.f[REST] = &DD[REST*size_Mat];
		D.f[TNE ] = &DD[TNE *size_Mat];
		D.f[TSW ] = &DD[TSW *size_Mat];
		D.f[TSE ] = &DD[TSE *size_Mat];
		D.f[TNW ] = &DD[TNW *size_Mat];
		D.f[BNE ] = &DD[BNE *size_Mat];
		D.f[BSW ] = &DD[BSW *size_Mat];
		D.f[BSE ] = &DD[BSE *size_Mat];
		D.f[BNW ] = &DD[BNW *size_Mat];
	} 
	else
	{
		D.f[W   ] = &DD[E   *size_Mat];
		D.f[E   ] = &DD[W   *size_Mat];
		D.f[S   ] = &DD[N   *size_Mat];
		D.f[N   ] = &DD[S   *size_Mat];
		D.f[B   ] = &DD[T   *size_Mat];
		D.f[T   ] = &DD[B   *size_Mat];
		D.f[SW  ] = &DD[NE  *size_Mat];
		D.f[NE  ] = &DD[SW  *size_Mat];
		D.f[NW  ] = &DD[SE  *size_Mat];
		D.f[SE  ] = &DD[NW  *size_Mat];
		D.f[BW  ] = &DD[TE  *size_Mat];
		D.f[TE  ] = &DD[BW  *size_Mat];
		D.f[TW  ] = &DD[BE  *size_Mat];
		D.f[BE  ] = &DD[TW  *size_Mat];
		D.f[BS  ] = &DD[TN  *size_Mat];
		D.f[TN  ] = &DD[BS  *size_Mat];
		D.f[TS  ] = &DD[BN  *size_Mat];
		D.f[BN  ] = &DD[TS  *size_Mat];
		D.f[REST] = &DD[REST*size_Mat];
		D.f[TNE ] = &DD[BSW *size_Mat];
		D.f[TSW ] = &DD[BNE *size_Mat];
		D.f[TSE ] = &DD[BNW *size_Mat];
		D.f[TNW ] = &DD[BSE *size_Mat];
		D.f[BNE ] = &DD[TSW *size_Mat];
		D.f[BSW ] = &DD[TNE *size_Mat];
		D.f[BSE ] = &DD[TNW *size_Mat];
		D.f[BNW ] = &DD[TSE *size_Mat];
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
		q_dirE   = &QQ[E   * numberOfBCnodes];
		q_dirW   = &QQ[W   * numberOfBCnodes];
		q_dirN   = &QQ[N   * numberOfBCnodes];
		q_dirS   = &QQ[S   * numberOfBCnodes];
		q_dirT   = &QQ[T   * numberOfBCnodes];
		q_dirB   = &QQ[B   * numberOfBCnodes];
		q_dirNE  = &QQ[NE  * numberOfBCnodes];
		q_dirSW  = &QQ[SW  * numberOfBCnodes];
		q_dirSE  = &QQ[SE  * numberOfBCnodes];
		q_dirNW  = &QQ[NW  * numberOfBCnodes];
		q_dirTE  = &QQ[TE  * numberOfBCnodes];
		q_dirBW  = &QQ[BW  * numberOfBCnodes];
		q_dirBE  = &QQ[BE  * numberOfBCnodes];
		q_dirTW  = &QQ[TW  * numberOfBCnodes];
		q_dirTN  = &QQ[TN  * numberOfBCnodes];
		q_dirBS  = &QQ[BS  * numberOfBCnodes];
		q_dirBN  = &QQ[BN  * numberOfBCnodes];
		q_dirTS  = &QQ[TS  * numberOfBCnodes];
		q_dirTNE = &QQ[TNE * numberOfBCnodes];
		q_dirTSW = &QQ[TSW * numberOfBCnodes];
		q_dirTSE = &QQ[TSE * numberOfBCnodes];
		q_dirTNW = &QQ[TNW * numberOfBCnodes];
		q_dirBNE = &QQ[BNE * numberOfBCnodes];
		q_dirBSW = &QQ[BSW * numberOfBCnodes];
		q_dirBSE = &QQ[BSE * numberOfBCnodes];
		q_dirBNW = &QQ[BNW * numberOfBCnodes];
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

		f_E   = (D.f[E   ])[ke   ];
		f_W   = (D.f[W   ])[kw   ];
		f_N   = (D.f[N   ])[kn   ];
		f_S   = (D.f[S   ])[ks   ];
		f_T   = (D.f[T   ])[kt   ];
		f_B   = (D.f[B   ])[kb   ];
		f_NE  = (D.f[NE  ])[kne  ];
		f_SW  = (D.f[SW  ])[ksw  ];
		f_SE  = (D.f[SE  ])[kse  ];
		f_NW  = (D.f[NW  ])[knw  ];
		f_TE  = (D.f[TE  ])[kte  ];
		f_BW  = (D.f[BW  ])[kbw  ];
		f_BE  = (D.f[BE  ])[kbe  ];
		f_TW  = (D.f[TW  ])[ktw  ];
		f_TN  = (D.f[TN  ])[ktn  ];
		f_BS  = (D.f[BS  ])[kbs  ];
		f_BN  = (D.f[BN  ])[kbn  ];
		f_TS  = (D.f[TS  ])[kts  ];
		f_TNE = (D.f[TNE ])[ktne ];
		f_TSW = (D.f[TSW ])[ktsw ];
		f_TSE = (D.f[TSE ])[ktse ];
		f_TNW = (D.f[TNW ])[ktnw ];
		f_BNE = (D.f[BNE ])[kbne ];
		f_BSW = (D.f[BSW ])[kbsw ];
		f_BSE = (D.f[BSE ])[kbse ];
		f_BNW = (D.f[BNW ])[kbnw ];
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
