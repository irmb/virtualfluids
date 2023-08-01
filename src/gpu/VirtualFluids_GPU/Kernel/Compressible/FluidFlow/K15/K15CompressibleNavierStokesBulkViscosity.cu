#include "K15CompressibleNavierStokesBulkViscosity.h"

#include "K15CompressibleNavierStokesBulkViscosity_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<K15CompressibleNavierStokesBulkViscosity> K15CompressibleNavierStokesBulkViscosity::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<K15CompressibleNavierStokesBulkViscosity>(new K15CompressibleNavierStokesBulkViscosity(para, level));
}

void K15CompressibleNavierStokesBulkViscosity::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    K15CompressibleNavierStokesBulkViscosity_Device <<<grid.grid, grid.threads>>>(
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

K15CompressibleNavierStokesBulkViscosity::K15CompressibleNavierStokesBulkViscosity(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	
}

K15CompressibleNavierStokesBulkViscosity::K15CompressibleNavierStokesBulkViscosity()
{
}
