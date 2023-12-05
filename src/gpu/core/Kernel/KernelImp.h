#ifndef KERNEL_IMP_H
#define KERNEL_IMP_H

#include "Calculation/Calculation.h" 

#include "Kernel.h"

#include <memory>

#include <cuda_helper/CudaGrid.h>

class Parameter;
class CudaStreamManager; 
class KernelImp : public Kernel
{
public:
    virtual void run() = 0;
    virtual void runOnIndices(const unsigned int *indices, unsigned int size_indices, CollisionTemplate collisionTemplate, CudaStreamIndex streamIndex=CudaStreamIndex::Legacy);

    std::vector<PreProcessorType> getPreProcessorTypes();

    bool getKernelUsesFluidNodeIndices();

protected:
    KernelImp(std::shared_ptr<Parameter> para, int level);
    KernelImp();

    std::shared_ptr<Parameter> para;
    int level;
    std::vector<PreProcessorType> myPreProcessorTypes;
    vf::cuda::CudaGrid cudaGrid;

    bool kernelUsesFluidNodeIndices = false;
};

#endif
