#ifndef VIRTUALFLUIDSPYTHONBINDINGS_KERNELFACTORY_H
#define VIRTUALFLUIDSPYTHONBINDINGS_KERNELFACTORY_H

#include <LBM/LBMKernel.h>
#include "AbstractLBMSystem.h"


class KernelFactory {
public:
    enum KernelType {
        BGK,
        COMPRESSIBLE_CUMULANT_4TH_ORDER_VISCOSITY /*,
        COMPRESSIBLE_CUMULANT,
        CUMULANT_K17,
        INCOMPRESSIBLE_CUMULANT,
        INCOMPRESSIBLE_CUMULANT_WITH_SPONGE_LAYER,
        INIT_DENSITIY,
        ET_D3Q27_BGK
        */
    };

    KernelFactory() = default;

    std::shared_ptr<LBMKernel> makeKernel(KernelType kernelType);

    std::shared_ptr<AbstractLBMSystem> makeLBMSystem(KernelType type);
};


#endif //VIRTUALFLUIDSPYTHONBINDINGS_KERNELFACTORY_H
