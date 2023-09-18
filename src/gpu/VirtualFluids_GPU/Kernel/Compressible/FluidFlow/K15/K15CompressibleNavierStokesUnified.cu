#include "K15CompressibleNavierStokesUnified.h"

#include <stdexcept>

#include "../RunLBMKernel.cuh"

#include "Parameter/Parameter.h"

#include <lbm/CumulantChimera.h>

namespace vf
{
namespace gpu
{

K15CompressibleNavierStokesUnified::K15CompressibleNavierStokesUnified(std::shared_ptr<Parameter> para, int level)
    : KernelImp(para, level)
{
#ifndef BUILD_CUDA_LTO
    throw std::invalid_argument(
        "To use the CumulantK15Unified kernel, pass -DBUILD_CUDA_LTO=ON to cmake. Requires: CUDA 11.2 & cc 5.0");
#endif

    myPreProcessorTypes.push_back(InitCompSP27);

    

    this->cudaGrid = cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);
}

void K15CompressibleNavierStokesUnified::run()
{
    GPUKernelParameter kernelParameter{
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        (int)para->getParD(level)->numberOfNodes,
        para->getParD(level)->forcing,
        para->getParD(level)->isEvenTimestep };

    auto lambda = [] __device__(lbm::KernelParameter parameter) {
        return lbm::cumulantChimera(parameter, lbm::setRelaxationRatesK15);
    };

    vf::gpu::runKernel<<<cudaGrid.grid, cudaGrid.threads>>>(lambda, kernelParameter);

    getLastCudaError("LB_Kernel_CumulantK15Comp execution failed");
}


}
}