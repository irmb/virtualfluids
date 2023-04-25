#ifndef GPU_KERNEL_H
#define GPU_KERNEL_H

#include <vector>

#include "LBM/LB.h" 

#include "PreProcessor/PreProcessorType.h"
#include "Parameter/CudaStreamManager.h"

#include <helper_cuda.h>

class Kernel
{
public:
    virtual ~Kernel()  = default;
    virtual void run() = 0;
    virtual void runOnIndices(const unsigned int *indices, unsigned int size_indices, CollisionTemplate collisionTemplate, CudaStreamIndex streamIdx=CudaStreamIndex::Legacy) = 0;

    virtual bool checkParameter()                                = 0;
    virtual std::vector<PreProcessorType> getPreProcessorTypes() = 0;
};
#endif
