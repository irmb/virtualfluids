#include "KernelFactoryImp.h"

#include <logger/Logger.h>

#include "Parameter/Parameter.h"

#include "Kernel/KernelTypes.h"

//LBM kernel (compressible)
#include "Kernel/Compressible/NavierStokes/B92/B92CompressibleNavierStokes.h"
#include "Kernel/Compressible/NavierStokes/B15/B15CompressibleNavierStokesBGKplus.h"
#include "Kernel/Compressible/NavierStokes/K15/K15CompressibleNavierStokes.h"
#include "Kernel/Compressible/NavierStokes/K15/K15CompressibleNavierStokesBulkViscosity.h"
#include "Kernel/Compressible/NavierStokes/K15/K15CompressibleNavierStokesSponge.h"
#include "Kernel/Compressible/NavierStokes/K17/K17CompressibleNavierStokes.h"
#include "Kernel/Compressible/NavierStokes/K17/K17CompressibleNavierStokesChimeraLegacy.h"
#include "Kernel/Compressible/NavierStokes/K17/K17CompressibleNavierStokesBulkViscosity.h"
#include "Kernel/Compressible/NavierStokes/K17/K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants.h"

//LBM kernel (inkompressible)
#include "Kernel/Incompressible/NavierStokes/B92/B92IncompressibleNavierStokes.h"
#include "Kernel/Incompressible/NavierStokes/B15/B15IncompressibleNavierStokesBGKplus.h"
#include "Kernel/Incompressible/NavierStokes/K15/K15IncompressibleNavierStokesRotatingVelocityField.h"
#include "Kernel/Incompressible/NavierStokes/K15/K15IncompressibleNavierStokesIsoCheck.h"
#include "Kernel/Incompressible/NavierStokes/K15/K15IncompressibleNavierStokes.h"

//advection diffusion kernel (compressible)
#include "Kernel/Compressible/AdvectionDiffusion/D3Q27/F16/F16CompressibleAdvectionDiffusion.h"

//advection diffusion kernel (incompressible)
#include "Kernel/Incompressible/AdvectionDiffusion/D3Q27/F16/F16IncompressibleAdvectionDiffusion.h"

#include <lbm/collision/TurbulentViscosity.h>

using namespace vf;

std::vector<std::shared_ptr<Kernel>> KernelFactoryImp::makeKernels(std::shared_ptr<Parameter> para)
{
    std::vector< std::shared_ptr< Kernel>> kernels;
    for (int level = 0; level <= para->getMaxLevel(); level++)
        kernels.push_back(makeKernel(para, para->getMainKernel(), level));

    if (para->getMaxLevel() > 0)
        if (para->getMultiKernelOn())
            for (std::size_t i = 0; i < para->getMultiKernelLevel().size(); i++)
                setKernelAtLevel(kernels, para, para->getMultiKernel().at(i), para->getMultiKernelLevel().at(i));
    return kernels;
}

std::vector<std::shared_ptr<AdvectionDiffusionKernel>> KernelFactoryImp::makeAdvDifKernels(std::shared_ptr<Parameter> para)
{
    std::vector< std::shared_ptr< AdvectionDiffusionKernel>> aDKernels;
    for (int level = 0; level <= para->getMaxLevel(); level++)
        aDKernels.push_back(makeAdvDifKernel(para, para->getADKernel(), level));
    return aDKernels;
}

void KernelFactoryImp::setPorousMedia(std::vector<std::shared_ptr<PorousMedia>> pm)
{
    this->pm = pm;
}

void KernelFactoryImp::setKernelAtLevel(std::vector<std::shared_ptr<Kernel>> kernels, std::shared_ptr<Parameter> para, std::string kernel, int level)
{
    kernels.at(level) = makeKernel(para, kernel, level);
}

std::shared_ptr<Kernel> KernelFactoryImp::makeKernel(std::shared_ptr<Parameter> para, std::string kernel, int level)
{
    VF_LOG_INFO("Instantiating Kernel: {}", kernel);
    std::shared_ptr<KernelImp> newKernel;

    if (kernel == collisionKernel::compressible::BGK) {
        newKernel     = B92CompressibleNavierStokes::getNewInstance(para, level);               // compressible
    } else if (kernel == collisionKernel::compressible::BGKPlus) {
        newKernel     = B15CompressibleNavierStokesBGKplus::getNewInstance(para, level);
    } else if (kernel == collisionKernel::compressible::K17CompressibleNavierStokes){
        switch(para->getTurbulenceModel())
        {
            case lbm::TurbulenceModel::AMD:
                newKernel = K17CompressibleNavierStokes<lbm::TurbulenceModel::AMD>::getNewInstance(para, level);
                break;
            case lbm::TurbulenceModel::Smagorinsky:
                newKernel = K17CompressibleNavierStokes<lbm::TurbulenceModel::Smagorinsky>::getNewInstance(para, level);
                break;
            case lbm::TurbulenceModel::QR:
                newKernel = K17CompressibleNavierStokes<lbm::TurbulenceModel::QR>::getNewInstance(para, level);
                break;
            case lbm::TurbulenceModel::None:
                newKernel = K17CompressibleNavierStokes<lbm::TurbulenceModel::None>::getNewInstance(para, level);
                break;
            default:
                throw std::runtime_error("Unknown turbulence model!");
            break;
        }
    } else if (kernel == collisionKernel::compressible::K15CompressibleNavierStokes) {
        newKernel     = K15CompressibleNavierStokes::getNewInstance(para, level);
    }                                                                           //===============
    else if (  kernel == collisionKernel::incompressible::BGK) {                // incompressible
        newKernel     = B92IncompressibleNavierStokes::getNewInstance(para, level);             //     ||
    } else if (kernel == collisionKernel::incompressible::BGKPlus) {
        newKernel     = B15IncompressibleNavierStokesBGKplus::getNewInstance(para, level);
    } else if (kernel == collisionKernel::incompressible::CumulantK15) {          //     /\      //
        newKernel     = K15IncompressibleNavierStokes::getNewInstance(para, level);           //     ||
    }                                                                             //===============
    else {
        throw std::runtime_error("KernelFactory does not know the KernelType.");
    }
    para->setKernelNeedsFluidNodeIndicesToRun(newKernel->getKernelUsesFluidNodeIndices());
    return newKernel;
}

std::shared_ptr<AdvectionDiffusionKernel> KernelFactoryImp::makeAdvDifKernel(std::shared_ptr<Parameter> para, std::string kernel, int level)
{
    std::shared_ptr<AdvectionDiffusionKernel> newKernel;

    if (kernel == "ADComp27") {
        newKernel     = F16CompressibleAdvectionDiffusion::getNewInstance(para, level);
    } else if (kernel == "ADIncomp27") {
        newKernel     = F16IncompressibleAdvectionDiffusion::getNewInstance(para, level);
    } else {
        throw std::runtime_error("KernelFactory does not know the KernelType.");
    }

    return newKernel;
}
