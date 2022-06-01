#ifndef KERNEL_IMP_H
#define KERNEL_IMP_H

#include "Kernel.h"

#include <memory>

#include <cuda/CudaGrid.h>

class CheckParameterStrategy;
class Parameter;

class KernelImp : public Kernel
{
public:
    virtual void run() = 0;
    virtual void runOnIndices(const unsigned int *indices, unsigned int size_indices, int stream = -1);

    bool checkParameter();
    std::vector<PreProcessorType> getPreProcessorTypes();
    KernelGroup getKernelGroup();

    void setCheckParameterStrategy(std::shared_ptr<CheckParameterStrategy> strategy);

protected:
    KernelImp(std::shared_ptr<Parameter> para, int level);
    KernelImp();

    std::shared_ptr<Parameter> para;
    std::shared_ptr<CheckParameterStrategy> checkStrategy;
    int level;
    std::vector<PreProcessorType> myPreProcessorTypes;
    KernelGroup myKernelGroup;

    vf::cuda::CudaGrid cudaGrid;
    
    std::unique_ptr<std::pair<dim3, dim3>> calcGridDimensions(unsigned int size_Mat, int numberOfThreads);
};

#endif
