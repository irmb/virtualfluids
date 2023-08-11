#include "K08CompressibleNavierStokes.h"

#include "K08CompressibleNavierStokes_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<K08CompressibleNavierStokes> K08CompressibleNavierStokes::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<K08CompressibleNavierStokes>(new K08CompressibleNavierStokes(para, level));
}

void K08CompressibleNavierStokes::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    K08CompressibleNavierStokes_Device <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_Cum_Comp_SP_27 execution failed");
}


K08CompressibleNavierStokes::K08CompressibleNavierStokes(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	
}

K08CompressibleNavierStokes::K08CompressibleNavierStokes()
{
}