#include <LBM/LBMKernel.h>
#include <LBM/K17CompressibleNavierStokes.h>
#include <LBM/B92IncompressibleNavierStokes.h>
#include <simulationconfig/D3Q27LBMSystem.h>
#include "simulationconfig/KernelFactory.h"

std::shared_ptr<LBMKernel> KernelFactory::makeKernel(KernelType kernelType)
{
    switch (kernelType) {
        case BGK:
            return std::shared_ptr<LBMKernel>(new B92IncompressibleNavierStokes());
        case COMPRESSIBLE_CUMULANT_4TH_ORDER_VISCOSITY:
            return std::shared_ptr<LBMKernel>(new K17CompressibleNavierStokes());
        default:
            throw std::logic_error("No such kernel type");
    }
}

std::shared_ptr<AbstractLBMSystem> KernelFactory::makeLBMSystem(KernelType type)
{
    return std::shared_ptr<AbstractLBMSystem>(new D3Q27LBMSystem());
}
