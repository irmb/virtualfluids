#include "CumulantK17Unified.h"

#include <stdexcept>

#include "Parameter/Parameter.h"
#include "../RunLBMKernel.cuh"

#include <lbm/CumulantChimera.h>

namespace vf
{
namespace gpu
{


CumulantK17Unified::CumulantK17Unified(std::shared_ptr<Parameter> para, int level)
    : KernelImp(para, level)
{
#ifndef BUILD_CUDA_LTO
    throw std::invalid_argument("To use the CumulantK17Unified kernel, pass -DBUILD_CUDA_LTO=ON to cmake. Requires: CUDA 11.2 & cc 5.0");
#endif

    myPreProcessorTypes.push_back(InitCompSP27);

    myKernelGroup = BasicKernel;

    this->cudaGrid = cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);
}



void CumulantK17Unified::run()
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
        return lbm::cumulantChimera(parameter, lbm::setRelaxationRatesK17);
    };

    runKernel<<<cudaGrid.grid, cudaGrid.threads>>>(lambda, kernelParameter);

    getLastCudaError("LB_Kernel_CumulantK17Unified execution failed");
}


}
}
