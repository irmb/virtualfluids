#include "ADIncomp7.h"

#include "ADIncomp7_Device.cuh"
#include "Parameter/Parameter.h"

std::shared_ptr<ADIncomp7> ADIncomp7::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<ADIncomp7>(new ADIncomp7(para, level));
}

void ADIncomp7::run()
{
	int size_Mat = (int)para->getParD(level)->numberOfNodes;
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

	LB_Kernel_AD_Incomp_7 <<< grid, threads >>>(
		para->getParD(level)->diffusivity, 
		para->getParD(level)->typeOfGridNode,
		para->getParD(level)->neighborX, 
		para->getParD(level)->neighborY, 
		para->getParD(level)->neighborZ,
		para->getParD(level)->distributions.f[0],
		para->getParD(level)->distributionsAD7.f[0], 
		para->getParD(level)->numberOfNodes,
		para->getParD(level)->isEvenTimestep);
	getLastCudaError("LB_Kernel_AD_Incomp_7 execution failed");
}

ADIncomp7::ADIncomp7(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitIncompAD7);

	myKernelGroup = ADKernel7;
}

ADIncomp7::ADIncomp7()
{
}
