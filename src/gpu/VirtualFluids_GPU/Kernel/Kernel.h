#ifndef GPU_KERNEL_H
#define GPU_KERNEL_H

#include <vector>

#include "Kernel/Utilities/KernelGroup.h"
#include "PreProcessor/PreProcessorType.h"
#include "Parameter/CudaStreamManager.h"

#include <helper_cuda.h>

class Kernel
{
public:
    virtual ~Kernel()  = default;
    virtual void run() = 0;
    virtual void runOnIndices(const unsigned int *indices, unsigned int size_indices, CudaStreamIndex streamIdx=CudaStreamIndex::Legacy) = 0; //if stream == -1: run on default stream

    virtual bool checkParameter()                                = 0;
    virtual std::vector<PreProcessorType> getPreProcessorTypes() = 0;
    virtual KernelGroup getKernelGroup()                         = 0;
};
#endif
