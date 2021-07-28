#include "CumulantK17CompChimSparse.h"

#include "Parameter/Parameter.h"
#include "CumulantK17CompChimSparse_Device.cuh"

std::shared_ptr<CumulantK17CompChimSparse> CumulantK17CompChimSparse::getNewInstance(std::shared_ptr<Parameter> para,
                                                                               int level)
{
    return std::shared_ptr<CumulantK17CompChimSparse>(new CumulantK17CompChimSparse(para, level));
}

void CumulantK17CompChimSparse::run()
{
	int numberOfThreads = para->getParD(level)->numberofthreads;
    int size_Mat = para->getParD(level)->numberOfFluidNodes;

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

	LB_Kernel_CumulantK17CompChimSparse <<< grid, threads >>>(
		para->getParD(level)->omega,
		para->getParD(level)->neighborX_SP,
		para->getParD(level)->neighborY_SP,
		para->getParD(level)->neighborZ_SP,
		para->getParD(level)->d0SP.f[0],
		para->getParD(level)->size_Mat_SP,
		level,
		para->getForcesDev(),
        para->getQuadricLimitersDev(),
		para->getParD(level)->evenOrOdd,
        para->getParD(level)->fluidNodeIndices,
		para->getParD(level)->numberOfFluidNodes);
	getLastCudaError("LB_Kernel_CumulantK17CompChim execution failed");
}

CumulantK17CompChimSparse::CumulantK17CompChimSparse(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}