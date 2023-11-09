#include "B12CompressibleAdvectionDiffusionD3Q7.h"

#include "B12CompressibleAdvectionDiffusionD3Q7_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<B12CompressibleAdvectionDiffusionD3Q7> B12CompressibleAdvectionDiffusionD3Q7::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<B12CompressibleAdvectionDiffusionD3Q7>(new B12CompressibleAdvectionDiffusionD3Q7(para, level));
}

void B12CompressibleAdvectionDiffusionD3Q7::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    B12CompressibleAdvectionDiffusionD3Q7_Device<<< grid.grid, grid.threads >>>(
        para->getParD(level)->diffusivity,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0], 
        para->getParD(level)->distributionsAD7.f[0], 
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("B12CompressibleAdvectionDiffusionD3Q7_Device execution failed");
}

B12CompressibleAdvectionDiffusionD3Q7::B12CompressibleAdvectionDiffusionD3Q7(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompAD7);

}

B12CompressibleAdvectionDiffusionD3Q7::B12CompressibleAdvectionDiffusionD3Q7()
{
}
