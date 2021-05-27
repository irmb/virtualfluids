#include "CumulantK15Unified.h"

#include "../CumulantKernel.cuh"

#include "Parameter/Parameter.h"

#include <lbm/Distribution27.h>

#include <lbm/CumulantChimera.h>

std::shared_ptr<CumulantK15Unified> CumulantK15Unified::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::make_shared<CumulantK15Unified>(para, level);
}

void CumulantK15Unified::run()
{
    vf::gpu::LBMKernelParameter kernelParameter{ para->getParD(level)->omega,
                                                 para->getParD(level)->geoSP,
                                                 para->getParD(level)->neighborX_SP,
                                                 para->getParD(level)->neighborY_SP,
                                                 para->getParD(level)->neighborZ_SP,
                                                 para->getParD(level)->d0SP.f[0],
                                                 (int)para->getParD(level)->size_Mat_SP,
                                                 level,
                                                 para->getForcesDev(),
                                                 para->getParD(level)->evenOrOdd };

    auto lambda = [] __device__(vf::lbm::CumulantChimeraParameter parameter) {
        return vf::lbm::cumulantChimera(parameter, vf::lbm::setRelaxationRatesK15);
    };

    vf::gpu::cumulantKernel<<<cudaGrid.grid, cudaGrid.threads>>>(lambda, kernelParameter);

    getLastCudaError("LB_Kernel_CumulantK15Comp execution failed");
}

CumulantK15Unified::CumulantK15Unified(std::shared_ptr<Parameter> para, int level)
{
#ifndef BUILD_CUDA_LTO
    throw std::invalid_argument(
        "To use the CumulantK15Unified kernel, pass -DBUILD_CUDA_LTO=ON to cmake. Requires: CUDA 11.2 & cc 5.0");
#endif

    this->para  = para;
    this->level = level;

    myPreProcessorTypes.push_back(InitCompSP27);

    myKernelGroup = BasicKernel;

    this->cudaGrid = vf::gpu::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->size_Mat_SP);
}
