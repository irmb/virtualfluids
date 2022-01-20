#include "BGKUnified.h"

#include <stdexcept>

#include "Parameter/Parameter.h"
#include "../RunLBMKernel.cuh"

#include <lbm/BGK.h>


namespace vf
{
namespace gpu
{


BGKUnified::BGKUnified(std::shared_ptr<Parameter> para, int level) 
    : KernelImp(para, level)
{
#ifndef BUILD_CUDA_LTO
    throw std::invalid_argument("To use the BKGUnified kernel, pass -DBUILD_CUDA_LTO=ON to cmake. Requires: CUDA 11.2 & cc 5.0");
#endif

    myPreProcessorTypes.push_back(InitCompSP27);

    myKernelGroup = BasicKernel;

    this->cudaGrid = cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->size_Mat_SP);
}


void BGKUnified::run()
{
    GPUKernelParameter kernelParameter{ para->getParD(level)->omega,
                                                 para->getParD(level)->geoSP,
                                                 para->getParD(level)->neighborX_SP,
                                                 para->getParD(level)->neighborY_SP,
                                                 para->getParD(level)->neighborZ_SP,
                                                 para->getParD(level)->d0SP.f[0],
                                                 (int)para->getParD(level)->size_Mat_SP,
                                                 nullptr, /* forces not used in bgk kernel */
                                                 para->getParD(level)->evenOrOdd };

    auto lambda = [] __device__(lbm::KernelParameter parameter) {
        return lbm::bgk(parameter);
    };

    runKernel<<<cudaGrid.grid, cudaGrid.threads>>>(lambda, kernelParameter);

    getLastCudaError("LB_Kernel_BGKUnified execution failed");
}


}
}
