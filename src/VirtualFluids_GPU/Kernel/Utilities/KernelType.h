#ifndef KERNEL_TYPE_H
#define KERNEL_TYPE_H


enum KernelType {
	LB_BGKCompSP27,
	LB_BGKPlusCompSP27,
	LB_CascadeCompSP27,
	LB_CumulantCompSP27,
	LB_CumulantK17Comp,
	LB_CumulantK17BulkComp,
	LB_CumulantAll4CompSP27,
	LB_CumulantK18Comp,
	LB_CumulantK20Comp,
	LB_CumulantK15Comp,
	LB_CumulantK15BulkComp,
	LB_CumulantK15SpongeComp,
	LB_MRTCompSP27,
	
	LB_BGKIncompSP27,
	LB_BGKPlusIncompSP27,
	LB_CascadeIncompSP27,
	LB_Cumulant1hIncompSP27,
	LB_CumulantIsoIncompSP27,
	LB_CumulantK15Incomp,
	LB_MRTIncompSP27,
	
	LB_PMCumulantOneCompSP27,
	
	LB_WaleCumulantAA2016CompSP27,
	LB_WaleCumulantAA2016DebugCompSP27,
	LB_WaleCumulantOneCompSP27,
	LB_WaleBySoniMalavCumulantOneCompSP27
};

enum ADKernelType
{
	LB_ADComp27,
	LB_ADComp7,

	LB_ADIncomp27,
	LB_ADIncomp7
};
#endif