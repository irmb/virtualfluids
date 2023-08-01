#include "K15IncompressibleNavierStokes.h"

#include "K15IncompressibleNavierStokes_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<K15IncompressibleNavierStokes> K15IncompressibleNavierStokes::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<K15IncompressibleNavierStokes>(new K15IncompressibleNavierStokes(para, level));
}

void K15IncompressibleNavierStokes::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    K15IncompressibleNavierStokes_Device <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_CumulantK15Incomp execution failed");
}

K15IncompressibleNavierStokes::K15IncompressibleNavierStokes(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitSP27);

	
}

K15IncompressibleNavierStokes::K15IncompressibleNavierStokes()
{
}
