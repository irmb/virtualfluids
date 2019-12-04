#include "WaleCumulantK17Comp.h"

#include "WaleCumulantK17Comp_Device.cuh"
#include "Parameter\Parameter.h"

std::shared_ptr<WaleCumulantK17Comp> WaleCumulantK17Comp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<WaleCumulantK17Comp>(new WaleCumulantK17Comp(para, level));
}

void WaleCumulantK17Comp::run()
{
	int size_Mat = para->getParD(level)->size_Mat_SP;
	int numberOfThreads = para->getParD(level)->numberofthreads;

	//int Grid = size_Array / numberOfThreads;
	//dim3 grid(Grid, 1, 1);
	//dim3 threads(numberOfThreads, 1, 1 );

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

	LB_Kernel_WaleCumulantK17Comp <<< grid, threads >>>(para->getParD(level)->omega,
														para->getParD(level)->geoSP,
														para->getParD(level)->neighborX_SP,
														para->getParD(level)->neighborY_SP,
														para->getParD(level)->neighborZ_SP,
														para->getParD(level)->neighborWSB_SP,
														para->getParD(level)->vx_SP,
														para->getParD(level)->vy_SP,
														para->getParD(level)->vz_SP,
														para->getParD(level)->d0SP.f[0],
														para->getParD(level)->turbViscosity,
														para->getParD(level)->size_Mat_SP,
														level,
														para->getTimestepOfCoarseLevel(),
														para->getForcesDev(),
                                                        para->getQuadricLimitersDev(),
														para->getParD(level)->evenOrOdd);
	getLastCudaError("LB_Kernel_WaleCumulantK17Comp execution failed");
}

WaleCumulantK17Comp::WaleCumulantK17Comp(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicWaleKernel;
}

WaleCumulantK17Comp::WaleCumulantK17Comp()
{
}
