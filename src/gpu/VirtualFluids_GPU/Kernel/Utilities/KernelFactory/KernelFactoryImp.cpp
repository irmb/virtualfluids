#include "KernelFactoryImp.h"

#include <logger/Logger.h>

#include "Parameter/Parameter.h"

#include "Kernel/Utilities/KernelTypes.h"

//LBM kernel (compressible)
#include "Kernel/Compressible/FluidFlow/B92/B92CompressibleNavierStokes.h"
#include "Kernel/Compressible/FluidFlow/B15/B15CompressibleNavierStokesBGKplusUnified.h"
#include "Kernel/Compressible/FluidFlow/B15/B15CompressibleNavierStokesBGKplus.h"
#include "Kernel/Compressible/FluidFlow/M02/M02CompressibleNavierStokes.h"
#include "Kernel/Compressible/FluidFlow/C06/C06CompressibleNavierStokes.h"
#include "Kernel/Compressible/FluidFlow/K08/K08CompressibleNavierStokes.h"
#include "Kernel/Compressible/FluidFlow/K15/K15CompressibleNavierStokes.h"
#include "Kernel/Compressible/FluidFlow/K15/K15CompressibleNavierStokesBulkViscosity.h"
#include "Kernel/Compressible/FluidFlow/K15/K15CompressibleNavierStokesSponge.h"
#include "Kernel/Compressible/FluidFlow/K15/K15CompressibleNavierStokesUnified.h"
#include "Kernel/Compressible/FluidFlow/K17/K17CompressibleNavierStokes.h"
#include "Kernel/Compressible/FluidFlow/K17/K17CompressibleNavierStokesChimeraLegacy.h"
#include "Kernel/Compressible/FluidFlow/K17/K17CompressibleNavierStokesBulkViscosity.h"
#include "Kernel/Compressible/FluidFlow/K17/K17CompressibleNavierStokesUnified.h"
#include "Kernel/Compressible/FluidFlow/K17/K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants.h"
#include "Kernel/Compressible/FluidFlow/K18/K18CompressibleNavierStokes.h"
#include "Kernel/Compressible/FluidFlow/K20/K20CompressibleNavierStokes.h"

//LBM kernel (inkompressible)
#include "Kernel/Incompressible/FluidFlow/B92/B92IncompressibleNavierStokes.h"
#include "Kernel/Incompressible/FluidFlow/B15/B15IncompressibleNavierStokesBGKplus.h"
#include "Kernel/Incompressible/FluidFlow/C06/C06IncompressibleNavierStokes.h"
#include "Kernel/Incompressible/FluidFlow/K15/K15IncompressibleNavierStokesRotatingVelocityField.h"
#include "Kernel/Incompressible/FluidFlow/K15/K15IncompressibleNavierStokesIsoCheck.h"
#include "Kernel/Incompressible/FluidFlow/K15/K15IncompressibleNavierStokes.h"
#include "Kernel/Incompressible/FluidFlow/M02/M02IncompressibleNavierStokes.h"

//advection diffusion kernel (compressible)
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Compressible/Mod27/ADComp27/ADComp27.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Compressible/Mod7/ADComp7/ADComp7.h"

//advection diffusion kernel (incompressible)
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Incompressible/Mod27/ADIncomp27/ADIncomp27.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Incompressible/Mod7/ADIncomp7/ADIncomp7.h"

//porous media kernel
#include "Kernel/Kernels/PorousMediaKernels/FluidFlow/Compressible/CumulantOne/PMCumulantOneCompSP27.h"

//wale kernel
#include "Kernel/Kernels/WaleKernels/FluidFlow/Compressible/CumulantK17/WaleCumulantK17Comp.h"
#include "Kernel/Kernels/WaleKernels/FluidFlow/Compressible/CumulantK17Debug/WaleCumulantK17DebugComp.h"
#include "Kernel/Kernels/WaleKernels/FluidFlow/Compressible/CumulantK15/WaleCumulantK15Comp.h"
#include "Kernel/Kernels/WaleKernels/FluidFlow/Compressible/CumulantK15BySoniMalav/WaleBySoniMalavCumulantK15Comp.h"

//strategies
#include "Kernel/Compressible/FluidFlow/FluidFlowCompStrategy.h"
#include "Kernel/Incompressible/FluidFlow/FluidFlowIncompStrategy.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Compressible/Mod27/ADMod27CompStrategy.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Compressible/Mod7/ADMod7CompStrategy.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Incompressible/Mod27/ADMod27IncompStrategy.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Incompressible/Mod7/ADMod7IncompStrategy.h"
#include "Kernel/Kernels/PorousMediaKernels/FluidFlow/Compressible/PMFluidFlowCompStrategy.h"
#include "Kernel/Kernels/WaleKernels/FluidFlow/Compressible/WaleFluidFlowCompStrategy.h"

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

std::vector<std::shared_ptr<ADKernel>> KernelFactoryImp::makeAdvDifKernels(std::shared_ptr<Parameter> para)
{
    std::vector< std::shared_ptr< ADKernel>> aDKernels;
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
    std::shared_ptr<CheckParameterStrategy> checkStrategy;

    if (kernel == collisionKernel::compressible::BGK) {
        newKernel     = B92CompressibleNavierStokes::getNewInstance(para, level);               // compressible
        checkStrategy = FluidFlowCompStrategy::getInstance();                   //      ||
    } else if (kernel == collisionKernel::compressible::BGKUnified) {           //      \/
        newKernel     = std::make_shared<vf::gpu::B15CompressibleNavierStokesBGKplusUnified>(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::BGKPlus) {
        newKernel     = B15CompressibleNavierStokesBGKplus::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::MRT) {
        newKernel     = M02CompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::Cascade) {
        newKernel     = C06CompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::CumulantClassic) {
        newKernel     = K08CompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::CumulantK15Unified) {
        newKernel     = std::make_shared<vf::gpu::K15CompressibleNavierStokesUnified>(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::K17CompressibleNavierStokesUnified) {
        newKernel     = std::make_shared<vf::gpu::K17CompressibleNavierStokesUnified>(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::K17CompressibleNavierStokesBulkViscosity) {
        newKernel     = K17CompressibleNavierStokesBulkViscosity::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::K17CompressibleNavierStokesChimeraLegacy) {
        newKernel     = K17CompressibleNavierStokesChimeraLegacy::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::K17CompressibleNavierStokes){
        switch(para->getTurbulenceModel())
        {
            case TurbulenceModel::AMD:
                newKernel = K17CompressibleNavierStokes<TurbulenceModel::AMD>::getNewInstance(para, level);
                break;
            case TurbulenceModel::Smagorinsky:
                newKernel = K17CompressibleNavierStokes<TurbulenceModel::Smagorinsky>::getNewInstance(para, level);
                break;
            case TurbulenceModel::QR:
                newKernel = K17CompressibleNavierStokes<TurbulenceModel::QR>::getNewInstance(para, level);
                break;
            case TurbulenceModel::None:
                newKernel = K17CompressibleNavierStokes<TurbulenceModel::None>::getNewInstance(para, level);
                break;
            default:
                throw std::runtime_error("Unknown turbulence model!");
            break;
        }
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::CumulantAll4SP27) {
        newKernel     = K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::CumulantK18) {
        newKernel     = K18CompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::CumulantK20) {
        newKernel     = K20CompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::K15CompressibleNavierStokes) {
        newKernel     = K15CompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::K15CompressibleNavierStokesBulk) {
        newKernel     = K15CompressibleNavierStokesBulkViscosity::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::compressible::K15CompressibleNavierStokesSponge) { //     /\      //
        newKernel     = K15CompressibleNavierStokesSponge::getNewInstance(para, level);     //     ||
        checkStrategy = FluidFlowCompStrategy::getInstance();                   // compressible
    }                                                                           //===============
    else if (  kernel == collisionKernel::incompressible::BGK) {                // incompressible
        newKernel     = B92IncompressibleNavierStokes::getNewInstance(para, level);             //     ||
        checkStrategy = FluidFlowIncompStrategy::getInstance();                 //     \/
    } else if (kernel == collisionKernel::incompressible::BGKPlus) {
        newKernel     = B15IncompressibleNavierStokesBGKplus::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == collisionKernel::incompressible::MRT) {
        newKernel     = M02IncompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == collisionKernel::incompressible::Cascade) {
        newKernel     = C06IncompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == collisionKernel::incompressible::Cumulant1h) {
        newKernel     = K15IncompressibleNavierStokesRotatingVelocityField::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    //} else if (kernel == collisionKernel::incompressible::CumulantIsometric) {
    //    newKernel     = K15IncompressibleNavierStokesIsoTest::getNewInstance(para, level);
    //    checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == collisionKernel::incompressible::CumulantK15) {          //     /\      //
        newKernel     = K15IncompressibleNavierStokes::getNewInstance(para, level);           //     ||
        checkStrategy = FluidFlowIncompStrategy::getInstance();                   // incompressible
    }                                                                             //===============
    else if (kernel == collisionKernel::porousMedia::CumulantOne) {               // porous media
        newKernel     = PMCumulantOneCompSP27::getNewInstance(para, pm, level);   //     ||
        checkStrategy = PMFluidFlowCompStrategy::getInstance();                   // porous media
    }                                                                             //===============
    else if (kernel == collisionKernel::wale::CumulantK17) {                      // wale model
        newKernel     = WaleCumulantK17Comp::getNewInstance(para, level);         //     ||
        checkStrategy = WaleFluidFlowCompStrategy::getInstance();                 //     \/
    } else if (kernel == collisionKernel::wale::CumulantK17Debug) {
        newKernel     = WaleCumulantK17DebugComp::getNewInstance(para, level);
        checkStrategy = WaleFluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::wale::CumulantK15) {
        newKernel     = WaleCumulantK15Comp::getNewInstance(para, level);
        checkStrategy = WaleFluidFlowCompStrategy::getInstance();
    } else if (kernel == collisionKernel::wale::CumulantK15SoniMalav) {              //     /\      //
        newKernel     = WaleBySoniMalavCumulantK15Comp::getNewInstance(para, level); //     ||
        checkStrategy = WaleFluidFlowCompStrategy::getInstance();                    // wale model
    }                                                                                //===============
    else {
        throw std::runtime_error("KernelFactory does not know the KernelType.");
    }
    newKernel->setCheckParameterStrategy(checkStrategy);
    para->setKernelNeedsFluidNodeIndicesToRun(newKernel->getKernelUsesFluidNodeIndices());
    return newKernel;
}

std::shared_ptr<ADKernel> KernelFactoryImp::makeAdvDifKernel(std::shared_ptr<Parameter> para, std::string kernel, int level)
{
    std::shared_ptr<ADKernel> newKernel;
    std::shared_ptr<CheckParameterStrategy> checkStrategy;

    if (kernel == "ADComp27") {
        newKernel     = ADComp27::getNewInstance(para, level);
        checkStrategy = ADMod27CompStrategy::getInstance();
    } else if(kernel == "ADComp7") {
        newKernel     = ADComp7::getNewInstance(para, level);
        checkStrategy = ADMod7CompStrategy::getInstance();
    } else if (kernel == "ADIncomp27") {
        newKernel     = ADIncomp27::getNewInstance(para, level);
        checkStrategy = ADMod7IncompStrategy::getInstance();
    } else if (kernel == "ADIncomp7") {
        newKernel     = ADIncomp7::getNewInstance(para, level);
        checkStrategy = ADMod7IncompStrategy::getInstance();
    } else {
        throw std::runtime_error("KernelFactory does not know the KernelType.");
    }

    if (newKernel) {
        newKernel->setCheckParameterStrategy(checkStrategy);
        return newKernel;
    }
    else
        throw  std::runtime_error("KernelFactory does not know the KernelType.");
}
