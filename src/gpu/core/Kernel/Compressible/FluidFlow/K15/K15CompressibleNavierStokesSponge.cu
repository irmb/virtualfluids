#include "K15CompressibleNavierStokesSponge.h"

#include "K15CompressibleNavierStokesSponge_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<K15CompressibleNavierStokesSponge> K15CompressibleNavierStokesSponge::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<K15CompressibleNavierStokesSponge>(new K15CompressibleNavierStokesSponge(para, level));
}

void K15CompressibleNavierStokesSponge::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    K15CompressibleNavierStokesSponge_Device <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->coordinateX,
        para->getParD(level)->coordinateY,
        para->getParD(level)->coordinateZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_CumulantK15SpongeComp execution failed");
}

K15CompressibleNavierStokesSponge::K15CompressibleNavierStokesSponge(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	
}

K15CompressibleNavierStokesSponge::K15CompressibleNavierStokesSponge()
{
}
