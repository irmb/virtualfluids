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
											bool evenOrOdd)
{
	Distributions27 D;
	if (evenOrOdd==true)
	{
		D.f[dirE   ] = &DD[dirE   *size_Mat];
		D.f[dirW   ] = &DD[dirW   *size_Mat];
		D.f[dirN   ] = &DD[dirN   *size_Mat];
		D.f[dirS   ] = &DD[dirS   *size_Mat];
		D.f[dirT   ] = &DD[dirT   *size_Mat];
		D.f[dirB   ] = &DD[dirB   *size_Mat];
		D.f[dirNE  ] = &DD[dirNE  *size_Mat];
		D.f[dirSW  ] = &DD[dirSW  *size_Mat];
		D.f[dirSE  ] = &DD[dirSE  *size_Mat];
		D.f[dirNW  ] = &DD[dirNW  *size_Mat];
		D.f[dirTE  ] = &DD[dirTE  *size_Mat];
		D.f[dirBW  ] = &DD[dirBW  *size_Mat];
		D.f[dirBE  ] = &DD[dirBE  *size_Mat];
		D.f[dirTW  ] = &DD[dirTW  *size_Mat];
		D.f[dirTN  ] = &DD[dirTN  *size_Mat];
		D.f[dirBS  ] = &DD[dirBS  *size_Mat];
		D.f[dirBN  ] = &DD[dirBN  *size_Mat];
		D.f[dirTS  ] = &DD[dirTS  *size_Mat];
		D.f[dirZERO] = &DD[dirZERO*size_Mat];
		D.f[dirTNE ] = &DD[dirTNE *size_Mat];
		D.f[dirTSW ] = &DD[dirTSW *size_Mat];
		D.f[dirTSE ] = &DD[dirTSE *size_Mat];
		D.f[dirTNW ] = &DD[dirTNW *size_Mat];
		D.f[dirBNE ] = &DD[dirBNE *size_Mat];
		D.f[dirBSW ] = &DD[dirBSW *size_Mat];
		D.f[dirBSE ] = &DD[dirBSE *size_Mat];
		D.f[dirBNW ] = &DD[dirBNW *size_Mat];
	} 
	else
	{
		D.f[dirW   ] = &DD[dirE   *size_Mat];
		D.f[dirE   ] = &DD[dirW   *size_Mat];
		D.f[dirS   ] = &DD[dirN   *size_Mat];
		D.f[dirN   ] = &DD[dirS   *size_Mat];
		D.f[dirB   ] = &DD[dirT   *size_Mat];
		D.f[dirT   ] = &DD[dirB   *size_Mat];
		D.f[dirSW  ] = &DD[dirNE  *size_Mat];
		D.f[dirNE  ] = &DD[dirSW  *size_Mat];
		D.f[dirNW  ] = &DD[dirSE  *size_Mat];
		D.f[dirSE  ] = &DD[dirNW  *size_Mat];
		D.f[dirBW  ] = &DD[dirTE  *size_Mat];
		D.f[dirTE  ] = &DD[dirBW  *size_Mat];
		D.f[dirTW  ] = &DD[dirBE  *size_Mat];
		D.f[dirBE  ] = &DD[dirTW  *size_Mat];
		D.f[dirBS  ] = &DD[dirTN  *size_Mat];
		D.f[dirTN  ] = &DD[dirBS  *size_Mat];
		D.f[dirTS  ] = &DD[dirBN  *size_Mat];
		D.f[dirBN  ] = &DD[dirTS  *size_Mat];
		D.f[dirZERO] = &DD[dirZERO*size_Mat];
		D.f[dirTNE ] = &DD[dirBSW *size_Mat];
		D.f[dirTSW ] = &DD[dirBNE *size_Mat];
		D.f[dirTSE ] = &DD[dirBNW *size_Mat];
		D.f[dirTNW ] = &DD[dirBSE *size_Mat];
		D.f[dirBNE ] = &DD[dirTSW *size_Mat];
		D.f[dirBSW ] = &DD[dirTNE *size_Mat];
		D.f[dirBSE ] = &DD[dirTNW *size_Mat];
		D.f[dirBNW ] = &DD[dirTSE *size_Mat];
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
		unsigned int sizeQ = numberOfBCnodes;
		////////////////////////////////////////////////////////////////////////////////
		real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
			*q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
			*q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
			*q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			*q_dirBSE, *q_dirBNW; 
		q_dirE   = &QQ[dirE   *sizeQ];
		q_dirW   = &QQ[dirW   *sizeQ];
		q_dirN   = &QQ[dirN   *sizeQ];
		q_dirS   = &QQ[dirS   *sizeQ];
		q_dirT   = &QQ[dirT   *sizeQ];
		q_dirB   = &QQ[dirB   *sizeQ];
		q_dirNE  = &QQ[dirNE  *sizeQ];
		q_dirSW  = &QQ[dirSW  *sizeQ];
		q_dirSE  = &QQ[dirSE  *sizeQ];
		q_dirNW  = &QQ[dirNW  *sizeQ];
		q_dirTE  = &QQ[dirTE  *sizeQ];
		q_dirBW  = &QQ[dirBW  *sizeQ];
		q_dirBE  = &QQ[dirBE  *sizeQ];
		q_dirTW  = &QQ[dirTW  *sizeQ];
		q_dirTN  = &QQ[dirTN  *sizeQ];
		q_dirBS  = &QQ[dirBS  *sizeQ];
		q_dirBN  = &QQ[dirBN  *sizeQ];
		q_dirTS  = &QQ[dirTS  *sizeQ];
		q_dirTNE = &QQ[dirTNE *sizeQ];
		q_dirTSW = &QQ[dirTSW *sizeQ];
		q_dirTSE = &QQ[dirTSE *sizeQ];
		q_dirTNW = &QQ[dirTNW *sizeQ];
		q_dirBNE = &QQ[dirBNE *sizeQ];
		q_dirBSW = &QQ[dirBSW *sizeQ];
		q_dirBSE = &QQ[dirBSE *sizeQ];
		q_dirBNW = &QQ[dirBNW *sizeQ];
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

		f_W    = (D.f[dirE   ])[ke   ];
		f_E    = (D.f[dirW   ])[kw   ];
		f_S    = (D.f[dirN   ])[kn   ];
		f_N    = (D.f[dirS   ])[ks   ];
		f_B    = (D.f[dirT   ])[kt   ];
		f_T    = (D.f[dirB   ])[kb   ];
		f_SW   = (D.f[dirNE  ])[kne  ];
		f_NE   = (D.f[dirSW  ])[ksw  ];
		f_NW   = (D.f[dirSE  ])[kse  ];
		f_SE   = (D.f[dirNW  ])[knw  ];
		f_BW   = (D.f[dirTE  ])[kte  ];
		f_TE   = (D.f[dirBW  ])[kbw  ];
		f_TW   = (D.f[dirBE  ])[kbe  ];
		f_BE   = (D.f[dirTW  ])[ktw  ];
		f_BS   = (D.f[dirTN  ])[ktn  ];
		f_TN   = (D.f[dirBS  ])[kbs  ];
		f_TS   = (D.f[dirBN  ])[kbn  ];
		f_BN   = (D.f[dirTS  ])[kts  ];
		f_BSW  = (D.f[dirTNE ])[ktne ];
		f_BNE  = (D.f[dirTSW ])[ktsw ];
		f_BNW  = (D.f[dirTSE ])[ktse ];
		f_BSE  = (D.f[dirTNW ])[ktnw ];
		f_TSW  = (D.f[dirBNE ])[kbne ];
		f_TNE  = (D.f[dirBSW ])[kbsw ];
		f_TNW  = (D.f[dirBSE ])[kbse ];
		f_TSE  = (D.f[dirBNW ])[kbnw ];
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
											bool evenOrOdd)
{
	Distributions27 D;
	if (evenOrOdd==true)
	{
		D.f[dirE   ] = &DD[dirE   *size_Mat];
		D.f[dirW   ] = &DD[dirW   *size_Mat];
		D.f[dirN   ] = &DD[dirN   *size_Mat];
		D.f[dirS   ] = &DD[dirS   *size_Mat];
		D.f[dirT   ] = &DD[dirT   *size_Mat];
		D.f[dirB   ] = &DD[dirB   *size_Mat];
		D.f[dirNE  ] = &DD[dirNE  *size_Mat];
		D.f[dirSW  ] = &DD[dirSW  *size_Mat];
		D.f[dirSE  ] = &DD[dirSE  *size_Mat];
		D.f[dirNW  ] = &DD[dirNW  *size_Mat];
		D.f[dirTE  ] = &DD[dirTE  *size_Mat];
		D.f[dirBW  ] = &DD[dirBW  *size_Mat];
		D.f[dirBE  ] = &DD[dirBE  *size_Mat];
		D.f[dirTW  ] = &DD[dirTW  *size_Mat];
		D.f[dirTN  ] = &DD[dirTN  *size_Mat];
		D.f[dirBS  ] = &DD[dirBS  *size_Mat];
		D.f[dirBN  ] = &DD[dirBN  *size_Mat];
		D.f[dirTS  ] = &DD[dirTS  *size_Mat];
		D.f[dirZERO] = &DD[dirZERO*size_Mat];
		D.f[dirTNE ] = &DD[dirTNE *size_Mat];
		D.f[dirTSW ] = &DD[dirTSW *size_Mat];
		D.f[dirTSE ] = &DD[dirTSE *size_Mat];
		D.f[dirTNW ] = &DD[dirTNW *size_Mat];
		D.f[dirBNE ] = &DD[dirBNE *size_Mat];
		D.f[dirBSW ] = &DD[dirBSW *size_Mat];
		D.f[dirBSE ] = &DD[dirBSE *size_Mat];
		D.f[dirBNW ] = &DD[dirBNW *size_Mat];
	} 
	else
	{
		D.f[dirW   ] = &DD[dirE   *size_Mat];
		D.f[dirE   ] = &DD[dirW   *size_Mat];
		D.f[dirS   ] = &DD[dirN   *size_Mat];
		D.f[dirN   ] = &DD[dirS   *size_Mat];
		D.f[dirB   ] = &DD[dirT   *size_Mat];
		D.f[dirT   ] = &DD[dirB   *size_Mat];
		D.f[dirSW  ] = &DD[dirNE  *size_Mat];
		D.f[dirNE  ] = &DD[dirSW  *size_Mat];
		D.f[dirNW  ] = &DD[dirSE  *size_Mat];
		D.f[dirSE  ] = &DD[dirNW  *size_Mat];
		D.f[dirBW  ] = &DD[dirTE  *size_Mat];
		D.f[dirTE  ] = &DD[dirBW  *size_Mat];
		D.f[dirTW  ] = &DD[dirBE  *size_Mat];
		D.f[dirBE  ] = &DD[dirTW  *size_Mat];
		D.f[dirBS  ] = &DD[dirTN  *size_Mat];
		D.f[dirTN  ] = &DD[dirBS  *size_Mat];
		D.f[dirTS  ] = &DD[dirBN  *size_Mat];
		D.f[dirBN  ] = &DD[dirTS  *size_Mat];
		D.f[dirZERO] = &DD[dirZERO*size_Mat];
		D.f[dirTNE ] = &DD[dirBSW *size_Mat];
		D.f[dirTSW ] = &DD[dirBNE *size_Mat];
		D.f[dirTSE ] = &DD[dirBNW *size_Mat];
		D.f[dirTNW ] = &DD[dirBSE *size_Mat];
		D.f[dirBNE ] = &DD[dirTSW *size_Mat];
		D.f[dirBSW ] = &DD[dirTNE *size_Mat];
		D.f[dirBSE ] = &DD[dirTNW *size_Mat];
		D.f[dirBNW ] = &DD[dirTSE *size_Mat];
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
		unsigned int sizeQ = numberOfBCnodes;
		////////////////////////////////////////////////////////////////////////////////
		real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
			*q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
			*q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
			*q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			*q_dirBSE, *q_dirBNW; 
		q_dirE   = &QQ[dirE   *sizeQ];
		q_dirW   = &QQ[dirW   *sizeQ];
		q_dirN   = &QQ[dirN   *sizeQ];
		q_dirS   = &QQ[dirS   *sizeQ];
		q_dirT   = &QQ[dirT   *sizeQ];
		q_dirB   = &QQ[dirB   *sizeQ];
		q_dirNE  = &QQ[dirNE  *sizeQ];
		q_dirSW  = &QQ[dirSW  *sizeQ];
		q_dirSE  = &QQ[dirSE  *sizeQ];
		q_dirNW  = &QQ[dirNW  *sizeQ];
		q_dirTE  = &QQ[dirTE  *sizeQ];
		q_dirBW  = &QQ[dirBW  *sizeQ];
		q_dirBE  = &QQ[dirBE  *sizeQ];
		q_dirTW  = &QQ[dirTW  *sizeQ];
		q_dirTN  = &QQ[dirTN  *sizeQ];
		q_dirBS  = &QQ[dirBS  *sizeQ];
		q_dirBN  = &QQ[dirBN  *sizeQ];
		q_dirTS  = &QQ[dirTS  *sizeQ];
		q_dirTNE = &QQ[dirTNE *sizeQ];
		q_dirTSW = &QQ[dirTSW *sizeQ];
		q_dirTSE = &QQ[dirTSE *sizeQ];
		q_dirTNW = &QQ[dirTNW *sizeQ];
		q_dirBNE = &QQ[dirBNE *sizeQ];
		q_dirBSW = &QQ[dirBSW *sizeQ];
		q_dirBSE = &QQ[dirBSE *sizeQ];
		q_dirBNW = &QQ[dirBNW *sizeQ];
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

		f_E   = (D.f[dirE   ])[ke   ];
		f_W   = (D.f[dirW   ])[kw   ];
		f_N   = (D.f[dirN   ])[kn   ];
		f_S   = (D.f[dirS   ])[ks   ];
		f_T   = (D.f[dirT   ])[kt   ];
		f_B   = (D.f[dirB   ])[kb   ];
		f_NE  = (D.f[dirNE  ])[kne  ];
		f_SW  = (D.f[dirSW  ])[ksw  ];
		f_SE  = (D.f[dirSE  ])[kse  ];
		f_NW  = (D.f[dirNW  ])[knw  ];
		f_TE  = (D.f[dirTE  ])[kte  ];
		f_BW  = (D.f[dirBW  ])[kbw  ];
		f_BE  = (D.f[dirBE  ])[kbe  ];
		f_TW  = (D.f[dirTW  ])[ktw  ];
		f_TN  = (D.f[dirTN  ])[ktn  ];
		f_BS  = (D.f[dirBS  ])[kbs  ];
		f_BN  = (D.f[dirBN  ])[kbn  ];
		f_TS  = (D.f[dirTS  ])[kts  ];
		f_TNE = (D.f[dirTNE ])[ktne ];
		f_TSW = (D.f[dirTSW ])[ktsw ];
		f_TSE = (D.f[dirTSE ])[ktse ];
		f_TNW = (D.f[dirTNW ])[ktnw ];
		f_BNE = (D.f[dirBNE ])[kbne ];
		f_BSW = (D.f[dirBSW ])[kbsw ];
		f_BSE = (D.f[dirBSE ])[kbse ];
		f_BNW = (D.f[dirBNW ])[kbnw ];
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
