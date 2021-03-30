#include "CumulantAll4CompSP27.h"

#include "CumulantAll4CompSP27_Device.cuh"

#include "Parameter/Parameter.h"

std::shared_ptr<CumulantAll4CompSP27> CumulantAll4CompSP27::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantAll4CompSP27>(new CumulantAll4CompSP27(para, level));
}

void CumulantAll4CompSP27::run()
{
	int numberOfThreads = para->getParD(level)->numberofthreads;
	int size_Mat = para->getParD(level)->size_Mat_SP;

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

	LB_Kernel_Cumulant_D3Q27All4 << < grid, threads >> >(	para->getParD(level)->omega,
															para->getParD(level)->geoSP,
															para->getParD(level)->neighborX_SP,
															para->getParD(level)->neighborY_SP,
															para->getParD(level)->neighborZ_SP,
															para->getParD(level)->d0SP.f[0],
															size_Mat,
															level,
															para->getForcesDev(),
                                                            para->getQuadricLimitersDev(),
															para->getParD(level)->evenOrOdd);
	getLastCudaError("LB_Kernel_Cumulant_D3Q27All4 execution failed");
}

CumulantAll4CompSP27::CumulantAll4CompSP27(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}