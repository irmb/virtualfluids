#include "CumulantK15BulkComp.h"

#include "CumulantK15BulkComp_Device.cuh"
#include "Parameter/Parameter.h"

std::shared_ptr<CumulantK15BulkComp> CumulantK15BulkComp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK15BulkComp>(new CumulantK15BulkComp(para, level));
}

void CumulantK15BulkComp::run()
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

	LB_Kernel_CumulantK15BulkComp <<< grid, threads >>>(
		para->getParD(level)->omega,
		para->getParD(level)->typeOfGridNode,
		para->getParD(level)->neighborX,
		para->getParD(level)->neighborY,
		para->getParD(level)->neighborZ,
		para->getParD(level)->distributions.f[0],
		para->getParD(level)->numberOfNodes,
		level,
		para->getForcesDev(),
		para->getParD(level)->isEvenTimestep);
	getLastCudaError("LB_Kernel_CumulantK15BulkComp execution failed");
}

CumulantK15BulkComp::CumulantK15BulkComp(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}

CumulantK15BulkComp::CumulantK15BulkComp()
{
}
