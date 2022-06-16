#include "MRTIncompSP27.h"

#include "MRTIncompSP27_Device.cuh"
#include "Parameter/Parameter.h"

std::shared_ptr<MRTIncompSP27> MRTIncompSP27::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<MRTIncompSP27>(new MRTIncompSP27(para, level));
}

void MRTIncompSP27::run()
{
	int size_Mat = para->getParD(level)->numberOfNodes;
	int numberOfThreads = para->getParD(level)->numberofthreads;

	int Grid = (size_Mat / numberOfThreads) + 1;
	int Grid1, Grid2;
	if (Grid>512)
	{
		Grid1 = 512;
		Grid2 = (Grid / Grid1) + 1;
	}
	else
	{
		Grid1 = 1;
		Grid2 = Grid;
	}
	dim3 grid(Grid1, Grid2);
	dim3 threads(numberOfThreads, 1, 1);

	LB_Kernel_MRT_Incomp_SP_27 << < grid, threads >> >(	para->getParD(level)->omega,
													para->getParD(level)->typeOfGridNode,
													para->getParD(level)->neighborX,
													para->getParD(level)->neighborY,
													para->getParD(level)->neighborZ,
													para->getParD(level)->distributions.f[0],
													para->getParD(level)->numberOfNodes,
													para->getParD(level)->isEvenTimestep);
	getLastCudaError("LB_Kernel_MRT_SP_27 execution failed");
}

MRTIncompSP27::MRTIncompSP27(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitSP27);

	myKernelGroup = BasicKernel;
}

MRTIncompSP27::MRTIncompSP27()
{
}
