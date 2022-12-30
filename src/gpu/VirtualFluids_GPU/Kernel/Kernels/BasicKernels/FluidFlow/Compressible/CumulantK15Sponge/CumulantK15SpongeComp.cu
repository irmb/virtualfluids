#include "CumulantK15SpongeComp.h"

#include "CumulantK15SpongeComp_Device.cuh"
#include "Parameter/Parameter.h"

std::shared_ptr<CumulantK15SpongeComp> CumulantK15SpongeComp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK15SpongeComp>(new CumulantK15SpongeComp(para, level));
}

void CumulantK15SpongeComp::run()
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

	LB_Kernel_CumulantK15SpongeComp <<< grid, threads >>>(
		para->getParD(level)->omega,
		para->getParD(level)->typeOfGridNode,
		para->getParD(level)->neighborX,
		para->getParD(level)->neighborY,
		para->getParD(level)->neighborZ,
		para->getParD(level)->coordinateX,
		para->getParD(level)->coordinateY,
		para->getParD(level)->coordinateZ,
		para->getParD(level)->distributions.f[0],
		para->getParD(level)->numberOfNodes,
		para->getParD(level)->isEvenTimestep);
	getLastCudaError("LB_Kernel_CumulantK15SpongeComp execution failed");
}

CumulantK15SpongeComp::CumulantK15SpongeComp(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}

CumulantK15SpongeComp::CumulantK15SpongeComp()
{
}
