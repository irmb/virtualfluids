#ifndef GPU_CUMULANT_KERNEL_H
#define GPU_CUMULANT_KERNEL_H


#include <DataTypes.h>
#include <cuda_runtime.h>

#include <lbm/Distribution27.h>
#include <lbm/CumulantChimeraParameter.h>

#include "Kernel/Utilities/DistributionHelper.cuh"

namespace vf
{
namespace gpu
{


struct LBMKernelParameter
{
    real omega;
    unsigned int* typeOfGridNode;
    unsigned int* neighborX;
    unsigned int* neighborY;
    unsigned int* neighborZ;
    real* distributions;
    int size_Mat;
    real* forces;
    bool isEvenTimestep;
};

template<typename KernelFunctor>
__global__ void cumulantKernel(KernelFunctor kernel, LBMKernelParameter kernelParameter)
{
    const uint k = vf::gpu::getNodeIndex();
    const uint nodeType = kernelParameter.typeOfGridNode[k];

    if (!vf::gpu::isValidFluidNode(k, kernelParameter.size_Mat, nodeType))
        return;

    vf::gpu::DistributionWrapper distributionWrapper {
        kernelParameter.distributions,
        kernelParameter.size_Mat,
        kernelParameter.isEvenTimestep,
        k,
        kernelParameter.neighborX,
        kernelParameter.neighborY,
        kernelParameter.neighborZ
    };

    lbm::CumulantChimeraParameter chimeraParameter {distributionWrapper.distribution, kernelParameter.omega, kernelParameter.forces};
    kernel(chimeraParameter);

    distributionWrapper.write();
}

}
}

#endif