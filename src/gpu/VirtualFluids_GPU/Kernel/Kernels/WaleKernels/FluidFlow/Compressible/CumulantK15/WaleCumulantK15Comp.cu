#include "WaleCumulantK15Comp.h"

#include "WaleCumulantK15Comp_Device.cuh"
#include "Parameter/Parameter.h"

std::shared_ptr<WaleCumulantK15Comp> WaleCumulantK15Comp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<WaleCumulantK15Comp>(new WaleCumulantK15Comp(para, level));
}

void WaleCumulantK15Comp::run()
{
	int size_Mat = para->getParD(level)->numberOfNodes;
	int numberOfThreads = para->getParD(level)->numberofthreads;

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

	LB_Kernel_WaleCumulantK15Comp << < grid, threads >> >(	para->getParD(level)->omega,
																para->getParD(level)->typeOfGridNode,
																para->getParD(level)->neighborX,
																para->getParD(level)->neighborY,
																para->getParD(level)->neighborZ,
																para->getParD(level)->neighborInverse,
																para->getParD(level)->velocityX,
																para->getParD(level)->velocityY,
																para->getParD(level)->velocityZ,
																para->getParD(level)->distributions.f[0],
																para->getParD(level)->turbViscosity,
																para->getParD(level)->numberOfNodes,
																level,
																para->getTimestepOfCoarseLevel(),
																para->getForcesDev(),
																para->getParD(level)->isEvenTimestep);
	getLastCudaError("LB_Kernel_WaleCumulantK15Comp execution failed");
}

WaleCumulantK15Comp::WaleCumulantK15Comp(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicWaleKernel;
}

WaleCumulantK15Comp::WaleCumulantK15Comp()
{
}
