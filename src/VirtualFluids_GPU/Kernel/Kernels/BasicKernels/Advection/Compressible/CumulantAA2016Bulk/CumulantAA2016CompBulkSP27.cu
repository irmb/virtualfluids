#include "CumulantAA2016CompBulkSP27.h"

#include "CumulantAA2016CompBulkSP27_Device.cuh"
#include "Parameter\Parameter.h"

std::shared_ptr<CumulantAA2016CompBulkSP27> CumulantAA2016CompBulkSP27::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantAA2016CompBulkSP27>(new CumulantAA2016CompBulkSP27(para, level));
}

void CumulantAA2016CompBulkSP27::run()
{
	int size_Array = para->getParD(level)->size_Array_SP;
	int numberOfThreads = para->getParD(level)->numberofthreads;

	int Grid = size_Array / numberOfThreads;
	dim3 grid(Grid, 1, 1);
	dim3 threads(numberOfThreads, 1, 1);

	LB_Kernel_Cum_AA2016_Comp_Bulk_SP_27 << < grid, threads >> >(	para->getParD(level)->omega,
																	para->getParD(level)->geoSP,
																	para->getParD(level)->neighborX_SP,
																	para->getParD(level)->neighborY_SP,
																	para->getParD(level)->neighborZ_SP,
																	para->getParD(level)->d0SP.f[0],
																	para->getParD(level)->size_Mat_SP,
																	level,
																	para->getForcesDev(),
																	para->getParD(level)->evenOrOdd);
	getLastCudaError("LB_Kernel_Kum_AA2016_Comp_Bulk_SP_27 execution failed");
}

CumulantAA2016CompBulkSP27::CumulantAA2016CompBulkSP27(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}

CumulantAA2016CompBulkSP27::CumulantAA2016CompBulkSP27()
{
}
