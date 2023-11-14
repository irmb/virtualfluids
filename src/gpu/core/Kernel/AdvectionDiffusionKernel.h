#ifndef AdvectionDiffusionKernel_H
#define AdvectionDiffusionKernel_H

#include "Kernel/KernelImp.h"

class AdvectionDiffusionKernel : public KernelImp
{
public:
    virtual ~AdvectionDiffusionKernel() = default;
    virtual void run() = 0;
};
#endif
