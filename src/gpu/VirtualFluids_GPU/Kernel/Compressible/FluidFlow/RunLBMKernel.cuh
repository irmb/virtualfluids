#ifndef GPU_CUMULANT_KERNEL_H
#define GPU_CUMULANT_KERNEL_H


#include <DataTypes.h>
#include <cuda_runtime.h>

#include "lbm/KernelParameter.h"

#include "Kernel/Utilities/DistributionHelper.cuh"

namespace vf
{
namespace gpu
{


struct GPUKernelParameter
{
    real omega;
    unsigned int* typeOfGridNode;
    unsigned int* neighborX;
    unsigned int* neighborY;
    unsigned int* neighborZ;
    real* distributions;
    int numberOfLBnodes;
    real* forces;
    bool isEvenTimestep;
};

template<typename KernelFunctor>
__global__ void runKernel(KernelFunctor kernel, GPUKernelParameter kernelParameter)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    if(nodeIndex >= kernelParameter.numberOfLBnodes)
        return;

    if (!isValidFluidNode(kernelParameter.typeOfGridNode[nodeIndex]))
        return;

    DistributionWrapper distributionWrapper {
        kernelParameter.distributions,
        (unsigned int)kernelParameter.numberOfLBnodes,
        kernelParameter.isEvenTimestep,
        nodeIndex,
        kernelParameter.neighborX,
        kernelParameter.neighborY,
        kernelParameter.neighborZ
    };

    lbm::KernelParameter parameter {distributionWrapper.distribution, kernelParameter.omega, kernelParameter.forces};
    kernel(parameter);

    distributionWrapper.write();
}

}
}

#endif