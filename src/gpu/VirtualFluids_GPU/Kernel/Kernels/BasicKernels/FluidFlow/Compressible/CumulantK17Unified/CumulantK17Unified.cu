#include "CumulantK17Unified.h"

#include "CumulantK17Unified_Device.cuh"
#include "Parameter/Parameter.h"

#include <stdexcept>

std::shared_ptr<CumulantK17Unified> CumulantK17Unified::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::shared_ptr<CumulantK17Unified>(new CumulantK17Unified(para, level));
}

void CumulantK17Unified::run()
{
    int numberOfThreads = para->getParD(level)->numberofthreads;
    int size_Mat        = para->getParD(level)->size_Mat_SP;

    int Grid = (size_Mat / numberOfThreads) + 1;
    int Grid1, Grid2;
    if (Grid > 512) {
        Grid1 = 512;
        Grid2 = (Grid / Grid1) + 1;
    } else {
        Grid1 = 1;
        Grid2 = Grid;
    }
    dim3 grid(Grid1, Grid2);
    dim3 threads(numberOfThreads, 1, 1);

    vf::gpu::LB_Kernel_CumulantK17Unified<<<grid, threads>>>(
        para->getParD(level)->omega,
        para->getParD(level)->geoSP,
        para->getParD(level)->neighborX_SP,
        para->getParD(level)->neighborY_SP,
        para->getParD(level)->neighborZ_SP,
        para->getParD(level)->d0SP.f[0],
        para->getParD(level)->size_Mat_SP,
        level,
        para->getForcesDev(),
        para->getParD(level)->evenOrOdd);

    getLastCudaError("LB_Kernel_CumulantK17Unified execution failed");
}

CumulantK17Unified::CumulantK17Unified(std::shared_ptr<Parameter> para, int level)
{
#ifndef BUILD_CUDA_LTO
    throw std::invalid_argument("To use the CumulantK17Unified kernel, pass -DBUILD_CUDA_LTO=ON to cmake. Requires: CUDA 11.2 & cc 5.0");
#endif

    this->para  = para;
    this->level = level;

    myPreProcessorTypes.push_back(InitCompSP27);

    myKernelGroup = BasicKernel;
}