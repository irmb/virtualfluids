#include "KernelFactoryImp.h"

#include <logger/Logger.h>

#include "Parameter/Parameter.h"

#include "Kernel/Utilities/KernelTypes.h"

//LBM kernel (compressible)
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/BGK/BGKCompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/BGKUnified/BGKUnified.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/BGKPlus/BGKPlusCompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/Cascade/CascadeCompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/Cumulant/CumulantCompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17Unified/K17CompressibleNavierStokesUnified.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17chim/K17CompressibleNavierStokesChimeraLegacy.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17/K17CompressibleNavierStokes.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17Bulk/K17CompressibleNavierStokesBulkViscosity.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantAll4/CumulantAll4CompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK18/CumulantK18Comp.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK20/CumulantK20Comp.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK15/CumulantK15Comp.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK15Unified/CumulantK15Unified.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK15Bulk/CumulantK15BulkComp.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK15Sponge/CumulantK15SpongeComp.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/MRT/MRTCompSP27.h"

//LBM kernel (inkompressible)
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/BGK/BGKIncompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/BGKPlus/BGKPlusIncompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/Cascade/CascadeIncompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/Cumulant1hSP27/Cumulant1hIncompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/CumulantIsoSP27/CumulantIsoIncompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/CumulantK15/CumulantK15Incomp.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/MRT/MRTIncompSP27.h"

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
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/FluidFlowCompStrategy.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/FluidFlowIncompStrategy.h"
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

    if (kernel == CollisionKernel::Compressible::BGK) {
        newKernel     = BGKCompSP27::getNewInstance(para, level);               // compressible
        checkStrategy = FluidFlowCompStrategy::getInstance();                   //      ||
    } else if (kernel == CollisionKernel::Compressible::BGKUnified) {           //      \/
        newKernel     = std::make_shared<vf::gpu::BGKUnified>(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::BGKPlus) {
        newKernel     = BGKPlusCompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::MRT) {
        newKernel     = MRTCompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::Cascade) {
        newKernel     = CascadeCompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantClassic) {
        newKernel     = CumulantCompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK15Unified) {
        newKernel     = std::make_shared<vf::gpu::CumulantK15Unified>(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK17Unified) {
        newKernel     = std::make_shared<vf::gpu::K17CompressibleNavierStokesUnified>(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK17Bulk) {
        newKernel     = K17CompressibleNavierStokesBulkViscosity::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK17Chim) {
        newKernel     = K17CompressibleNavierStokesChimeraLegacy::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK17){
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
    } else if (kernel == CollisionKernel::Compressible::CumulantAll4SP27) {
        newKernel     = CumulantAll4CompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK18) {
        newKernel     = CumulantK18Comp::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK20) {
        newKernel     = CumulantK20Comp::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK15) {
        newKernel     = CumulantK15Comp::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK15Bulk) {
        newKernel     = CumulantK15BulkComp::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK15Sponge) {    //     /\      //
        newKernel     = CumulantK15SpongeComp::getNewInstance(para, level);     //     ||
        checkStrategy = FluidFlowCompStrategy::getInstance();                   // compressible
    }                                                                           //===============
    else if (  kernel == CollisionKernel::Incompressible::BGK) {                // incompressible
        newKernel     = BGKIncompSP27::getNewInstance(para, level);             //     ||
        checkStrategy = FluidFlowIncompStrategy::getInstance();                 //     \/
    } else if (kernel == CollisionKernel::Incompressible::BGKPlus) {
        newKernel     = BGKPlusIncompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Incompressible::MRT) {
        newKernel     = MRTIncompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Incompressible::Cascade) {
        newKernel     = CascadeIncompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Incompressible::Cumulant1h) {
        newKernel     = Cumulant1hIncompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Incompressible::CumulantIsometric) {
        newKernel     = CumulantIsoIncompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Incompressible::CumulantK15) {          //     /\      //
        newKernel     = CumulantK15Incomp::getNewInstance(para, level);           //     ||
        checkStrategy = FluidFlowIncompStrategy::getInstance();                   // incompressible
    }                                                                             //===============
    else if (kernel == CollisionKernel::PorousMedia::CumulantOne) {               // porous media
        newKernel     = PMCumulantOneCompSP27::getNewInstance(para, pm, level);   //     ||
        checkStrategy = PMFluidFlowCompStrategy::getInstance();                   // porous media
    }                                                                             //===============
    else if (kernel == CollisionKernel::Wale::CumulantK17) {                      // wale model
        newKernel     = WaleCumulantK17Comp::getNewInstance(para, level);         //     ||
        checkStrategy = WaleFluidFlowCompStrategy::getInstance();                 //     \/
    } else if (kernel == CollisionKernel::Wale::CumulantK17Debug) {
        newKernel     = WaleCumulantK17DebugComp::getNewInstance(para, level);
        checkStrategy = WaleFluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Wale::CumulantK15) {
        newKernel     = WaleCumulantK15Comp::getNewInstance(para, level);
        checkStrategy = WaleFluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Wale::CumulantK15SoniMalav) {              //     /\      //
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
