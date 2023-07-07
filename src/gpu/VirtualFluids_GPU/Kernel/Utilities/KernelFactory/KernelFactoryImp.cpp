#include "KernelFactoryImp.h"

#include <logger/Logger.h>

#include "Parameter/Parameter.h"

#include "Kernel/Utilities/KernelTypes.h"

//LBM kernel (compressible)
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/B92CompressibleNavierStokes/B92CompressibleNavierStokes.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/B15CompressibleNavierStokesBGKplusUnified/B15CompressibleNavierStokesBGKplusUnified.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/B15CompressibleNavierStokesBGKplus/B15CompressibleNavierStokesBGKplus.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/M02CompressibleNavierStokes/M02CompressibleNavierStokes.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/C06CompressibleNavierStokes/C06CompressibleNavierStokes.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/K08CompressibleNavierStokes/K08CompressibleNavierStokes.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/K15CompressibleNavierStokes/K15CompressibleNavierStokes.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/K15CompressibleNavierStokesBulkViscosity/K15CompressibleNavierStokesBulkViscosity.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/K15CompressibleNavierStokesSponge/K15CompressibleNavierStokesSponge.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/K15CompressibleNavierStokesUnified/K15CompressibleNavierStokesUnified.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/K17CompressibleNavierStokes/K17CompressibleNavierStokes.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/K17CompressibleNavierStokesChimeraLegacy/K17CompressibleNavierStokesChimeraLegacy.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/K17CompressibleNavierStokesBulkViscosity/K17CompressibleNavierStokesBulkViscosity.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/K17CompressibleNavierStokesUnified/K17CompressibleNavierStokesUnified.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants/K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/K18CompressibleNavierStokes/K18CompressibleNavierStokes.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/K20CompressibleNavierStokes/K20CompressibleNavierStokes.h"

//LBM kernel (inkompressible)
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/B92IncompressibleNavierStokes/B92IncompressibleNavierStokes.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/B15IncompressibleNavierStokesBGKplus/B15IncompressibleNavierStokesBGKplus.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/C06IncompressibleNavierStokes/C06IncompressibleNavierStokes.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/K15IncompressibleNavierStokesRotatingVelocityField/K15IncompressibleNavierStokesRotatingVelocityField.h"
//#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/K15IncompressibleNavierStokesIsoTest/K15IncompressibleNavierStokesIsoTest.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/K15IncompressibleNavierStokes/K15IncompressibleNavierStokes.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/M02IncompressibleNavierStokes/M02IncompressibleNavierStokes.h"

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
        newKernel     = B92CompressibleNavierStokes::getNewInstance(para, level);               // compressible
        checkStrategy = FluidFlowCompStrategy::getInstance();                   //      ||
    } else if (kernel == CollisionKernel::Compressible::BGKUnified) {           //      \/
        newKernel     = std::make_shared<vf::gpu::B15CompressibleNavierStokesBGKplusUnified>(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::BGKPlus) {
        newKernel     = B15CompressibleNavierStokesBGKplus::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::MRT) {
        newKernel     = M02CompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::Cascade) {
        newKernel     = C06CompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantClassic) {
        newKernel     = K08CompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK15Unified) {
        newKernel     = std::make_shared<vf::gpu::K15CompressibleNavierStokesUnified>(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::K17CompressibleNavierStokesUnified) {
        newKernel     = std::make_shared<vf::gpu::K17CompressibleNavierStokesUnified>(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::K17CompressibleNavierStokesBulkViscosity) {
        newKernel     = K17CompressibleNavierStokesBulkViscosity::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::K17CompressibleNavierStokesChimeraLegacy) {
        newKernel     = K17CompressibleNavierStokesChimeraLegacy::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::K17CompressibleNavierStokes){
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
        newKernel     = K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK18) {
        newKernel     = K18CompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::CumulantK20) {
        newKernel     = K20CompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::K15CompressibleNavierStokes) {
        newKernel     = K15CompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::K15CompressibleNavierStokesBulk) {
        newKernel     = K15CompressibleNavierStokesBulkViscosity::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Compressible::K15CompressibleNavierStokesSponge) { //     /\      //
        newKernel     = K15CompressibleNavierStokesSponge::getNewInstance(para, level);     //     ||
        checkStrategy = FluidFlowCompStrategy::getInstance();                   // compressible
    }                                                                           //===============
    else if (  kernel == CollisionKernel::Incompressible::BGK) {                // incompressible
        newKernel     = B92IncompressibleNavierStokes::getNewInstance(para, level);             //     ||
        checkStrategy = FluidFlowIncompStrategy::getInstance();                 //     \/
    } else if (kernel == CollisionKernel::Incompressible::BGKPlus) {
        newKernel     = B15IncompressibleNavierStokesBGKplus::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Incompressible::MRT) {
        newKernel     = M02IncompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Incompressible::Cascade) {
        newKernel     = C06IncompressibleNavierStokes::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Incompressible::Cumulant1h) {
        newKernel     = K15IncompressibleNavierStokesRotatingVelocityField::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    //} else if (kernel == CollisionKernel::Incompressible::CumulantIsometric) {
    //    newKernel     = K15IncompressibleNavierStokesIsoTest::getNewInstance(para, level);
    //    checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == CollisionKernel::Incompressible::CumulantK15) {          //     /\      //
        newKernel     = K15IncompressibleNavierStokes::getNewInstance(para, level);           //     ||
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
