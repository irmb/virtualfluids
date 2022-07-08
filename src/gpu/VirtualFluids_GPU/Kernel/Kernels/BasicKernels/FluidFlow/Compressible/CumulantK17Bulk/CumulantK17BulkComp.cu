#include "CumulantK17BulkComp.h"

#include "CumulantK17BulkComp_Device.cuh"
#include "Parameter/Parameter.h"

std::shared_ptr<CumulantK17BulkComp> CumulantK17BulkComp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK17BulkComp>(new CumulantK17BulkComp(para, level));
}

void CumulantK17BulkComp::run()
{
	int size_Array = para->getParD(level)->size_Array_SP;
	int numberOfThreads = para->getParD(level)->numberofthreads;

	int Grid = size_Array / numberOfThreads;
	dim3 grid(Grid, 1, 1);
	dim3 threads(numberOfThreads, 1, 1);

	LB_Kernel_CumulantK17BulkComp << < grid, threads >> >(	para->getParD(level)->omega,
																	para->getParD(level)->typeOfGridNode,
																	para->getParD(level)->neighborX,
																	para->getParD(level)->neighborY,
																	para->getParD(level)->neighborZ,
																	para->getParD(level)->distributions.f[0],
																	para->getParD(level)->numberOfNodes,
																	level,
																	para->getForcesDev(),
                                                                    para->getQuadricLimitersDev(),
																	para->getParD(level)->isEvenTimestep);
	getLastCudaError("LB_Kernel_CumulantK17BulkComp execution failed");
}

CumulantK17BulkComp::CumulantK17BulkComp(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}

CumulantK17BulkComp::CumulantK17BulkComp()
{
}
