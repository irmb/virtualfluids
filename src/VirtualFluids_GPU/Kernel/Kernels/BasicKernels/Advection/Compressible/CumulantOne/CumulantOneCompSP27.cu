#include "CumulantOneCompSP27.h"

#include "CumulantOneCompSP27_Device.cuh"

#include "Parameter/Parameter.h"

std::shared_ptr<CumulantOneCompSP27> CumulantOneCompSP27::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantOneCompSP27>(new CumulantOneCompSP27(para, level));
}

void CumulantOneCompSP27::run()
{
	int numberOfThreads = para->getParD(level)->numberofthreads;
	int size_Mat = para->getParD(level)->size_Mat_SP;

	int Grid = (size_Mat / numberOfThreads) + 1;
	int Grid1, Grid2;
	if (Grid > 512)
	{
		Grid1 = 512;
		Grid2 = (Grid / Grid1) + 1;
	}
	else
	{
		Grid1 = 1;
		Grid2 = Grid;
	}
	dim3 grid(Grid1, Grid2, 1);
	dim3 threads(numberOfThreads, 1, 1);

	LB_Kernel_Kum_One_Comp_SP_27 << < grid, threads >> >(	para->getParD(level)->omega,
															para->getParD(level)->geoSP,
															para->getParD(level)->neighborX_SP,
															para->getParD(level)->neighborY_SP,
															para->getParD(level)->neighborZ_SP,
															para->getParD(level)->d0SP.f[0],
															size_Mat,
															level,
															para->getForcesDev(),
															para->getParD(level)->evenOrOdd);
	getLastCudaError("LB_Kernel_Kum_New_Comp_SP_27 execution failed");
}

CumulantOneCompSP27::CumulantOneCompSP27(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}