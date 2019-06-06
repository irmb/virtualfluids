#include "WaleCumulantOneCompSP27.h"

#include "WaleCumulantOneCompSP27_Device.cuh"
#include "Parameter\Parameter.h"

std::shared_ptr<WaleCumulantOneCompSP27> WaleCumulantOneCompSP27::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<WaleCumulantOneCompSP27>(new WaleCumulantOneCompSP27(para, level));
}

void WaleCumulantOneCompSP27::run()
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

	LB_Kernel_Wale_Cum_One_Comp_SP_27 << < grid, threads >> >(	para->getParD(level)->omega,
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
																para->getForcesDev(),
																para->getParD(level)->evenOrOdd);
	getLastCudaError("LB_Kernel_Wale_Cum_One_Comp_SP_27 execution failed");
}

WaleCumulantOneCompSP27::WaleCumulantOneCompSP27(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicWaleKernel;
}

WaleCumulantOneCompSP27::WaleCumulantOneCompSP27()
{
}
