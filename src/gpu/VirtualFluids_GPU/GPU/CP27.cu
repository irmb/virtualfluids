/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////
__global__ void CalcCP27(real* DD, 
									int* cpIndex, 
									int nonCp, 
									double *cpPress,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									unsigned long long numberOfLBnodes, 
									bool isEvenTimestep)
{
	Distributions27 D;
	if (isEvenTimestep==true)
	{
		D.f[dP00] = &DD[dP00 * numberOfLBnodes];
		D.f[dM00] = &DD[dM00 * numberOfLBnodes];
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
		D.f[d000] = &DD[d000 * numberOfLBnodes];
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
		D.f[dM00] = &DD[dP00 * numberOfLBnodes];
		D.f[dP00] = &DD[dM00 * numberOfLBnodes];
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
		D.f[d000] = &DD[d000 * numberOfLBnodes];
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

		PressCP  =   (D.f[dP00])[ke  ]+ (D.f[dM00])[kw  ]+ 
                     (D.f[DIR_0P0])[kn  ]+ (D.f[DIR_0M0])[ks  ]+
                     (D.f[DIR_00P])[kt  ]+ (D.f[DIR_00M])[kb  ]+
                     (D.f[DIR_PP0])[kne ]+ (D.f[DIR_MM0])[ksw ]+
                     (D.f[DIR_PM0])[kse ]+ (D.f[DIR_MP0])[knw ]+
                     (D.f[DIR_P0P])[kte ]+ (D.f[DIR_M0M])[kbw ]+
                     (D.f[DIR_P0M])[kbe ]+ (D.f[DIR_M0P])[ktw ]+
                     (D.f[DIR_0PP])[ktn ]+ (D.f[DIR_0MM])[kbs ]+
                     (D.f[DIR_0PM])[kbn ]+ (D.f[DIR_0MP])[kts ]+
                     (D.f[d000])[kzero]+ 
                     (D.f[DIR_PPP])[ktne]+ (D.f[DIR_MMP])[ktsw]+ 
                     (D.f[DIR_PMP])[ktse]+ (D.f[DIR_MPP])[ktnw]+ 
                     (D.f[DIR_PPM])[kbne]+ (D.f[DIR_MMM])[kbsw]+ 
                     (D.f[DIR_PMM])[kbse]+ (D.f[DIR_MPM])[kbnw];
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		cpPress[k] = PressCP;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}
}

