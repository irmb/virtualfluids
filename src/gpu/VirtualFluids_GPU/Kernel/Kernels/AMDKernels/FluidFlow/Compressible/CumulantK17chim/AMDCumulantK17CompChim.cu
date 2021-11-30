#include "AMDCumulantK17CompChim.h"

#include "Parameter/Parameter.h"
#include "AMDCumulantK17CompChim_Device.cuh"

std::shared_ptr<AMDCumulantK17CompChim> AMDCumulantK17CompChim::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<AMDCumulantK17CompChim>(new AMDCumulantK17CompChim(para,level));
}

void AMDCumulantK17CompChim::run()
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
	//printf("size mat before %d, after %ld  \n", para->getParD(level)->size_Mat_SP, (unsigned long)para->getParD(level)->size_Mat_SP);
	LB_Kernel_AMDCumulantK17CompChim <<< grid, threads >>>(
		para->getParD(level)->omega,
		para->getSGSConstant(),
		para->getParD(level)->geoSP,
		para->getParD(level)->neighborX_SP,
		para->getParD(level)->neighborY_SP,
		para->getParD(level)->neighborZ_SP,
		para->getParD(level)->neighborWSB_SP,
		para->getParD(level)->d0SP.f[0],
		para->getParD(level)->vx_SP,
		para->getParD(level)->vy_SP,
		para->getParD(level)->vz_SP,
		(unsigned long)para->getParD(level)->size_Mat_SP,
		level,
		para->getIsBodyForce(),
		para->getForcesDev(),
		para->getParD(level)->forceX_SP,
		para->getParD(level)->forceY_SP,
		para->getParD(level)->forceZ_SP,
        para->getQuadricLimitersDev(),
		para->getParD(level)->evenOrOdd);
	getLastCudaError("LB_Kernel_CumulantK17CompChim execution failed");
}

AMDCumulantK17CompChim::AMDCumulantK17CompChim(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}