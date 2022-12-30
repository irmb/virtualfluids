#include "CumulantK18Comp.h"

#include "CumulantK18Comp_Device.cuh"

#include "Parameter/Parameter.h"

std::shared_ptr<CumulantK18Comp> CumulantK18Comp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK18Comp>(new CumulantK18Comp(para, level));
}

void CumulantK18Comp::run()
{
	int numberOfThreads = para->getParD(level)->numberofthreads;
	int size_Mat = (int)para->getParD(level)->numberOfNodes;

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

	LB_Kernel_CumulantK18Comp <<< grid, threads >>>(
		para->getParD(level)->omega,
		para->getParD(level)->typeOfGridNode,
		para->getParD(level)->neighborX,
		para->getParD(level)->neighborY,
		para->getParD(level)->neighborZ,
		para->getParD(level)->distributions.f[0],
		para->getParD(level)->g6.g[0],
		para->getParD(level)->numberOfNodes,
		level,
		para->getForcesDev(),
        para->getQuadricLimitersDev(),
		para->getParD(level)->isEvenTimestep);
	getLastCudaError("LB_Kernel_CumulantK18Comp execution failed");
}

CumulantK18Comp::CumulantK18Comp(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);
	myPreProcessorTypes.push_back(InitF3);

	myKernelGroup = F3Kernel;
}