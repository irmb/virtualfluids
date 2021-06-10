#ifndef GPU_KERNEL_H
#define GPU_KERNEL_H

#include <vector>

#include "Kernel/Utilities/KernelGroup.h"
#include "PreProcessor/PreProcessorType.h"

class Kernel
{
public:
    virtual ~Kernel()  = default;
    virtual void run() = 0;

    virtual bool checkParameter()                                = 0;
    virtual std::vector<PreProcessorType> getPreProcessorTypes() = 0;
    virtual KernelGroup getKernelGroup()                         = 0;
};
#endif
