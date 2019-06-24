#include "CumulantF3CompSP27.h"

#include "CumulantF3CompSP27_Device.cuh"

#include "Parameter/Parameter.h"

std::shared_ptr<CumulantF3CompSP27> CumulantF3CompSP27::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantF3CompSP27>(new CumulantF3CompSP27(para, level));
}

void CumulantF3CompSP27::run()
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

	LB_Kernel_Cumulant_D3Q27F3 << < grid, threads >> >(	para->getParD(level)->omega,
														para->getParD(level)->geoSP,
														para->getParD(level)->neighborX_SP,
														para->getParD(level)->neighborY_SP,
														para->getParD(level)->neighborZ_SP,
														para->getParD(level)->d0SP.f[0],
														para->getParD(level)->g6.g[0],
														para->getParD(level)->size_Mat_SP,
														level,
														para->getForcesDev(),
														para->getParD(level)->evenOrOdd);
	getLastCudaError("LB_Kernel_Cumulant_D3Q27F3 execution failed");
}

CumulantF3CompSP27::CumulantF3CompSP27(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);
	myPreProcessorTypes.push_back(InitF3);

	myKernelGroup = F3Kernel;
}