#include "ADComp7.h"

#include "ADComp7_Device.cuh"
#include "Parameter/Parameter.h"

std::shared_ptr<ADComp7> ADComp7::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<ADComp7>(new ADComp7(para, level));
}

void ADComp7::run()
{
	int size_Mat = para->getParD(level)->size_Mat_SP;
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

	LB_Kernel_AD_Comp_7 << < grid, threads >> >(	para->getParD(level)->diffusivity,
											para->getParD(level)->geoSP,
											para->getParD(level)->neighborX_SP,
											para->getParD(level)->neighborY_SP,
											para->getParD(level)->neighborZ_SP,
											para->getParD(level)->d0SP.f[0], 
											para->getParD(level)->d7.f[0], 
											para->getParD(level)->size_Mat_SP,
											para->getParD(level)->evenOrOdd);
	getLastCudaError("LB_Kernel_ThS7 execution failed");
}

ADComp7::ADComp7(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompAD7);

	myKernelGroup = ADKernel7;
}

ADComp7::ADComp7()
{
}
