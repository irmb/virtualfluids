/* Device code */
#include "LBM/D3Q27.h"
#include "GPU/constant.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void CalcCP27(real* DD, 
									int* cpIndex, 
									int nonCp, 
									double *cpPress,
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

	if(k<nonCp)
	{
		////////////////////////////////////////////////////////////////////////////////
		//index
		unsigned int KQK  = cpIndex[k];
		unsigned int kzero= KQK;
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
		double PressCP;

		PressCP  =   (D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
                     (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
                     (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
                     (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
                     (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                     (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
                     (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                     (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
                     (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
                     (D.f[dirZERO])[kzero]+ 
                     (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                     (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                     (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
                     (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw];
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		cpPress[k] = PressCP;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}
}

