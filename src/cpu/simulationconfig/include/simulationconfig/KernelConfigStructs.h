#ifndef VIRTUALFLUIDSPYTHONBINDINGS_KERNELCONFIGSTRUCTS_H
#define VIRTUALFLUIDSPYTHONBINDINGS_KERNELCONFIGSTRUCTS_H

#include <string>
#include <LBM/LBMSystem.h>

struct LBMKernelConfig {
    KernelFactory::KernelType kernelType;
    bool useForcing = false;
    LBMReal forcingX1{};
    LBMReal forcingX2{};
    LBMReal forcingX3{};

    explicit LBMKernelConfig(KernelFactory::KernelType kernelType) : kernelType(kernelType)
    {
    }
};

#endif //VIRTUALFLUIDSPYTHONBINDINGS_KERNELCONFIGSTRUCTS_H
