#include "K20CompressibleNavierStokes.h"

#include "K20CompressibleNavierStokes_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<K20CompressibleNavierStokes> K20CompressibleNavierStokes::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<K20CompressibleNavierStokes>(new K20CompressibleNavierStokes(para, level));
}

void K20CompressibleNavierStokes::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    K20CompressibleNavierStokes_Device <<< grid.grid, grid.threads >>>(
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
    getLastCudaError("LB_Kernel_CumulantK20Comp execution failed");
}

K20CompressibleNavierStokes::K20CompressibleNavierStokes(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);
	myPreProcessorTypes.push_back(InitF3);

	
}