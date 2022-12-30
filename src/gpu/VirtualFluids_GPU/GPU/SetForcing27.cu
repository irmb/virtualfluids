/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////
__global__ void GetVeloforForcing27( real* DD, 
												int* bcIndex, 
												int nonAtBC, 
												real* Vx,
												real* Vy,
												real* Vz,
												unsigned int* neighborX,
												unsigned int* neighborY,
												unsigned int* neighborZ,
												unsigned long long numberOfLBnodes, 
												bool isEvenTimestep)
{
	Distributions27 D;
	if (isEvenTimestep==false)
	{
		D.f[DIR_P00   ] = &DD[DIR_P00   *numberOfLBnodes];
		D.f[DIR_M00   ] = &DD[DIR_M00   *numberOfLBnodes];
		D.f[DIR_0P0   ] = &DD[DIR_0P0   *numberOfLBnodes];
		D.f[DIR_0M0   ] = &DD[DIR_0M0   *numberOfLBnodes];
		D.f[DIR_00P   ] = &DD[DIR_00P   *numberOfLBnodes];
		D.f[DIR_00M   ] = &DD[DIR_00M   *numberOfLBnodes];
		D.f[DIR_PP0  ] = &DD[DIR_PP0  *numberOfLBnodes];
		D.f[DIR_MM0  ] = &DD[DIR_MM0  *numberOfLBnodes];
		D.f[DIR_PM0  ] = &DD[DIR_PM0  *numberOfLBnodes];
		D.f[DIR_MP0  ] = &DD[DIR_MP0  *numberOfLBnodes];
		D.f[DIR_P0P  ] = &DD[DIR_P0P  *numberOfLBnodes];
		D.f[DIR_M0M  ] = &DD[DIR_M0M  *numberOfLBnodes];
		D.f[DIR_P0M  ] = &DD[DIR_P0M  *numberOfLBnodes];
		D.f[DIR_M0P  ] = &DD[DIR_M0P  *numberOfLBnodes];
		D.f[DIR_0PP  ] = &DD[DIR_0PP  *numberOfLBnodes];
		D.f[DIR_0MM  ] = &DD[DIR_0MM  *numberOfLBnodes];
		D.f[DIR_0PM  ] = &DD[DIR_0PM  *numberOfLBnodes];
		D.f[DIR_0MP  ] = &DD[DIR_0MP  *numberOfLBnodes];
		D.f[DIR_000] = &DD[DIR_000*numberOfLBnodes];
		D.f[DIR_PPP ] = &DD[DIR_PPP *numberOfLBnodes];
		D.f[DIR_MMP ] = &DD[DIR_MMP *numberOfLBnodes];
		D.f[DIR_PMP ] = &DD[DIR_PMP *numberOfLBnodes];
		D.f[DIR_MPP ] = &DD[DIR_MPP *numberOfLBnodes];
		D.f[DIR_PPM ] = &DD[DIR_PPM *numberOfLBnodes];
		D.f[DIR_MMM ] = &DD[DIR_MMM *numberOfLBnodes];
		D.f[DIR_PMM ] = &DD[DIR_PMM *numberOfLBnodes];
		D.f[DIR_MPM ] = &DD[DIR_MPM *numberOfLBnodes];
	} 
	else
	{
		D.f[DIR_M00   ] = &DD[DIR_P00   *numberOfLBnodes];
		D.f[DIR_P00   ] = &DD[DIR_M00   *numberOfLBnodes];
		D.f[DIR_0M0   ] = &DD[DIR_0P0   *numberOfLBnodes];
		D.f[DIR_0P0   ] = &DD[DIR_0M0   *numberOfLBnodes];
		D.f[DIR_00M   ] = &DD[DIR_00P   *numberOfLBnodes];
		D.f[DIR_00P   ] = &DD[DIR_00M   *numberOfLBnodes];
		D.f[DIR_MM0  ] = &DD[DIR_PP0  *numberOfLBnodes];
		D.f[DIR_PP0  ] = &DD[DIR_MM0  *numberOfLBnodes];
		D.f[DIR_MP0  ] = &DD[DIR_PM0  *numberOfLBnodes];
		D.f[DIR_PM0  ] = &DD[DIR_MP0  *numberOfLBnodes];
		D.f[DIR_M0M  ] = &DD[DIR_P0P  *numberOfLBnodes];
		D.f[DIR_P0P  ] = &DD[DIR_M0M  *numberOfLBnodes];
		D.f[DIR_M0P  ] = &DD[DIR_P0M  *numberOfLBnodes];
		D.f[DIR_P0M  ] = &DD[DIR_M0P  *numberOfLBnodes];
		D.f[DIR_0MM  ] = &DD[DIR_0PP  *numberOfLBnodes];
		D.f[DIR_0PP  ] = &DD[DIR_0MM  *numberOfLBnodes];
		D.f[DIR_0MP  ] = &DD[DIR_0PM  *numberOfLBnodes];
		D.f[DIR_0PM  ] = &DD[DIR_0MP  *numberOfLBnodes];
		D.f[DIR_000] = &DD[DIR_000*numberOfLBnodes];
		D.f[DIR_PPP ] = &DD[DIR_MMM *numberOfLBnodes];
		D.f[DIR_MMP ] = &DD[DIR_PPM *numberOfLBnodes];
		D.f[DIR_PMP ] = &DD[DIR_MPM *numberOfLBnodes];
		D.f[DIR_MPP ] = &DD[DIR_PMM *numberOfLBnodes];
		D.f[DIR_PPM ] = &DD[DIR_MMP *numberOfLBnodes];
		D.f[DIR_MMM ] = &DD[DIR_PPP *numberOfLBnodes];
		D.f[DIR_PMM ] = &DD[DIR_MPP *numberOfLBnodes];
		D.f[DIR_MPM ] = &DD[DIR_PMP *numberOfLBnodes];
	}
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////
	if(k < nonAtBC)
	{
		////////////////////////////////////////////////////////////////////////////////
		//index
		unsigned int KQK  = bcIndex[k];
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
		real mfcbb = (D.f[DIR_P00   ])[ke   ];
		real mfabb = (D.f[DIR_M00   ])[kw   ];
		real mfbcb = (D.f[DIR_0P0   ])[kn   ];
		real mfbab = (D.f[DIR_0M0   ])[ks   ];
		real mfbbc = (D.f[DIR_00P   ])[kt   ];
		real mfbba = (D.f[DIR_00M   ])[kb   ];
		real mfccb = (D.f[DIR_PP0  ])[kne  ];
		real mfaab = (D.f[DIR_MM0  ])[ksw  ];
		real mfcab = (D.f[DIR_PM0  ])[kse  ];
		real mfacb = (D.f[DIR_MP0  ])[knw  ];
		real mfcbc = (D.f[DIR_P0P  ])[kte  ];
		real mfaba = (D.f[DIR_M0M  ])[kbw  ];
		real mfcba = (D.f[DIR_P0M  ])[kbe  ];
		real mfabc = (D.f[DIR_M0P  ])[ktw  ];
		real mfbcc = (D.f[DIR_0PP  ])[ktn  ];
		real mfbaa = (D.f[DIR_0MM  ])[kbs  ];
		real mfbca = (D.f[DIR_0PM  ])[kbn  ];
		real mfbac = (D.f[DIR_0MP  ])[kts  ];
		real mfbbb = (D.f[DIR_000])[kzero];
		real mfccc = (D.f[DIR_PPP ])[ktne ];
		real mfaac = (D.f[DIR_MMP ])[ktsw ];
		real mfcac = (D.f[DIR_PMP ])[ktse ];
		real mfacc = (D.f[DIR_MPP ])[ktnw ];
		real mfcca = (D.f[DIR_PPM ])[kbne ];
		real mfaaa = (D.f[DIR_MMM ])[kbsw ];
		real mfcaa = (D.f[DIR_PMM ])[kbse ];
		real mfaca = (D.f[DIR_MPM ])[kbnw ];
		////////////////////////////////////////////////////////////////////////////////////
		real rho   = (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
					 	 mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
						 mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb + c1o1);//!!!!Achtung + one
		////////////////////////////////////////////////////////////////////////////////////
		real vx =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
			         (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
				       (mfcbb-mfabb))/ rho;
		real vy =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
			         (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
			           (mfbcb-mfbab)) / rho;
		real vz =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
			         (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
			           (mfbbc-mfbba)) / rho;
		////////////////////////////////////////////////////////////////////////////////////
		Vx[k] = vx;
		Vy[k] = vy;
		Vz[k] = vz;
		////////////////////////////////////////////////////////////////////////////////////
	}
}

