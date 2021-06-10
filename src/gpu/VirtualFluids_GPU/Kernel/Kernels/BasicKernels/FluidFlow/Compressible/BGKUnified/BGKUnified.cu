#include "BGKUnified.h"

#include "Parameter/Parameter.h"
#include "../CumulantKernel.cuh"
#include "Kernel/Utilities/CudaGrid.h"
#include <stdexcept>

#include <lbm/BGK.h>


std::shared_ptr<BGKUnified> BGKUnified::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::make_shared<BGKUnified>(para, level);
}

void BGKUnified::run()
{
    vf::gpu::LBMKernelParameter kernelParameter{ para->getParD(level)->omega,
                                                 para->getParD(level)->geoSP,
                                                 para->getParD(level)->neighborX_SP,
                                                 para->getParD(level)->neighborY_SP,
                                                 para->getParD(level)->neighborZ_SP,
                                                 para->getParD(level)->d0SP.f[0],
                                                 (int)para->getParD(level)->size_Mat_SP,
                                                 nullptr, /* forces not used in bgk kernel */
                                                 para->getParD(level)->evenOrOdd };

    auto lambda = [] __device__(vf::lbm::CumulantChimeraParameter parameter) {
        return vf::lbm::bgk(parameter);
    };

    vf::gpu::cumulantKernel<<<cudaGrid.grid, cudaGrid.threads>>>(lambda, kernelParameter);

    getLastCudaError("LB_Kernel_BGKUnified execution failed");
}

BGKUnified::BGKUnified(std::shared_ptr<Parameter> para, int level)
{
#ifndef BUILD_CUDA_LTO
    throw std::invalid_argument("To use the BKGUnified kernel, pass -DBUILD_CUDA_LTO=ON to cmake. Requires: CUDA 11.2 & cc 5.0");
#endif

    this->para  = para;
    this->level = level;

    myPreProcessorTypes.push_back(InitCompSP27);

    myKernelGroup = BasicKernel;

    this->cudaGrid = vf::gpu::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->size_Mat_SP);
}
