#ifndef KERNEL_TYPE_H
#define KERNEL_TYPE_H


enum KernelType {
	LB_BGKCompSP27,
	LB_BGKPlusCompSP27,
	LB_CascadeCompSP27,
	LB_CumulantCompSP27,
	LB_CumulantAA2016CompSP27,
	LB_CumulantAA2016CompBulkSP27,
	LB_CumulantAll4CompSP27,
	LB_CumulantF3CompSP27,
	LB_CumulantF32018CompSP27,
	LB_CumulantOneCompSP27,
	LB_CumulantOneCompBulkSP27,
	LB_CumulantOneCompSpongeSP27,
	LB_MRTCompSP27,
	
	LB_BGKIncompSP27,
	LB_BGKPlusIncompSP27,
	LB_CascadeIncompSP27,
	LB_Cumulant1hIncompSP27,
	LB_CumulantIsoIncompSP27,
	LB_CumulantOneIncompSP27,
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