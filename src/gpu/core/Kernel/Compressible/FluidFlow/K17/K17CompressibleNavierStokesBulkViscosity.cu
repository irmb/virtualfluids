#include "K17CompressibleNavierStokesBulkViscosity.h"

#include "K17CompressibleNavierStokesBulkViscosity_Device.cuh"
#include "Parameter/Parameter.h"

std::shared_ptr<K17CompressibleNavierStokesBulkViscosity> K17CompressibleNavierStokesBulkViscosity::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<K17CompressibleNavierStokesBulkViscosity>(new K17CompressibleNavierStokesBulkViscosity(para, level));
}

void K17CompressibleNavierStokesBulkViscosity::run()
{
	int size_Array = para->getParD(level)->size_Array_SP;
	int numberOfThreads = para->getParD(level)->numberofthreads;

	int Grid = size_Array / numberOfThreads;
	dim3 grid(Grid, 1, 1);
	dim3 threads(numberOfThreads, 1, 1);

	K17CompressibleNavierStokesBulkViscosity_Device << < grid, threads >> >(
		para->getParD(level)->omega,
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
	getLastCudaError("K17CompressibleNavierStokesBulkViscosity_Device execution failed");
}

K17CompressibleNavierStokesBulkViscosity::K17CompressibleNavierStokesBulkViscosity(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	
}

K17CompressibleNavierStokesBulkViscosity::K17CompressibleNavierStokesBulkViscosity()
{
}
