#include "WaleCumulantAA2016DebugCompSP27.h"

#include "WaleCumulantAA2016DebugCompSP27_Device.cuh"
#include "Parameter\Parameter.h"

std::shared_ptr<WaleCumulantAA2016DebugCompSP27> WaleCumulantAA2016DebugCompSP27::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<WaleCumulantAA2016DebugCompSP27>(new WaleCumulantAA2016DebugCompSP27(para, level));
}

void WaleCumulantAA2016DebugCompSP27::run()
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

	LB_Kernel_Wale_Cum_AA2016_Debug_Comp_SP_27 << < grid, threads >> >(
																		para->getParD(level)->omega,
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
																		para->getParD(level)->gSij,
																		para->getParD(level)->gSDij,
																		para->getParD(level)->gDxvx,
																		para->getParD(level)->gDyvx,
																		para->getParD(level)->gDzvx,
																		para->getParD(level)->gDxvy,
																		para->getParD(level)->gDyvy,
																		para->getParD(level)->gDzvy,
																		para->getParD(level)->gDxvz,
																		para->getParD(level)->gDyvz,
																		para->getParD(level)->gDzvz,
																		para->getParD(level)->size_Mat_SP,
																		level,
																		para->getForcesDev(),
																		para->getParD(level)->evenOrOdd);
	getLastCudaError("LB_Kernel_Wale_Cum_AA2016_Debug_Comp_SP_27 execution failed");
}

WaleCumulantAA2016DebugCompSP27::WaleCumulantAA2016DebugCompSP27(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicWaleKernel;
}

WaleCumulantAA2016DebugCompSP27::WaleCumulantAA2016DebugCompSP27()
{
}
