#include "F16CompressibleAdvectionDiffusion.h"

#include "F16CompressibleAdvectionDiffusion_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<F16CompressibleAdvectionDiffusion> F16CompressibleAdvectionDiffusion::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::shared_ptr<F16CompressibleAdvectionDiffusion>(new F16CompressibleAdvectionDiffusion(para, level));
}

void F16CompressibleAdvectionDiffusion::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    F16CompressibleAdvectionDiffusion_Device<<< grid.grid, grid.threads >>>(
        para->getParD(level)->diffusivity,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->distributionsAD.f[0],
        para->getParD(level)->numberOfNodes, 
        para->getParD(level)->forcing,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("F16CompressibleAdvectionDiffusion_Device execution failed");
}

F16CompressibleAdvectionDiffusion::F16CompressibleAdvectionDiffusion(std::shared_ptr<Parameter> para, int level)
{
    this->para = para;
    this->level = level;

    myPreProcessorTypes.push_back(InitCompAD27);

}

F16CompressibleAdvectionDiffusion::F16CompressibleAdvectionDiffusion()
{
}
