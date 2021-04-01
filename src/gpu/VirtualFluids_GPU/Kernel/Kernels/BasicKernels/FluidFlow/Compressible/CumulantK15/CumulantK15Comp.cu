#include "CumulantK15Comp.h"

#include "CumulantK15Comp_Device.cuh"

#include "Parameter/Parameter.h"

std::shared_ptr<CumulantK15Comp> CumulantK15Comp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK15Comp>(new CumulantK15Comp(para, level));
}

void CumulantK15Comp::run()
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

	LB_Kernel_CumulantK15Comp <<< grid, threads >>>(para->getParD(level)->omega,
													para->getParD(level)->geoSP,
													para->getParD(level)->neighborX_SP,
													para->getParD(level)->neighborY_SP,
													para->getParD(level)->neighborZ_SP,
													para->getParD(level)->d0SP.f[0],
													size_Mat,
													level,
													para->getForcesDev(),
													para->getParD(level)->evenOrOdd);
	getLastCudaError("LB_Kernel_CumulantK15Comp execution failed");
}

CumulantK15Comp::CumulantK15Comp(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}