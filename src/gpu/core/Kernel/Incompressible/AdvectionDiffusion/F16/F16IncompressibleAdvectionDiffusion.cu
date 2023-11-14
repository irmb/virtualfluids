#include "F16IncompressibleAdvectionDiffusion.h"

#include "F16IncompressibleAdvectionDiffusion_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<F16IncompressibleAdvectionDiffusion> F16IncompressibleAdvectionDiffusion::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::shared_ptr<F16IncompressibleAdvectionDiffusion>(new F16IncompressibleAdvectionDiffusion(para, level));
}

void F16IncompressibleAdvectionDiffusion::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    F16IncompressibleAdvectionDiffusion_Device<<< grid.grid, grid.threads >>>(
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
    getLastCudaError("F16IncompressibleAdvectionDiffusion_Device execution failed");
}

F16IncompressibleAdvectionDiffusion::F16IncompressibleAdvectionDiffusion(std::shared_ptr<Parameter> para, int level)
{
    this->para = para;
    this->level = level;

    myPreProcessorTypes.push_back(InitAdvectionDiffusionIncompressible);

}

F16IncompressibleAdvectionDiffusion::F16IncompressibleAdvectionDiffusion()
{
}
