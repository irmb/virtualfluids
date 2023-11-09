#include "KernelFactoryImp.h"

#include <logger/Logger.h>

#include "Parameter/Parameter.h"

#include "Kernel/KernelTypes.h"

//LBM kernel (compressible)
#include "Kernel/Compressible/NavierStokes/B92/B92CompressibleNavierStokes.h"
#include "Kernel/Compressible/NavierStokes/B15/B15CompressibleNavierStokesBGKplus.h"
#include "Kernel/Compressible/NavierStokes/M02/M02CompressibleNavierStokes.h"
#include "Kernel/Compressible/NavierStokes/C06/C06CompressibleNavierStokes.h"
#include "Kernel/Compressible/NavierStokes/K08/K08CompressibleNavierStokes.h"
#include "Kernel/Compressible/NavierStokes/K15/K15CompressibleNavierStokes.h"
#include "Kernel/Compressible/NavierStokes/K15/K15CompressibleNavierStokesBulkViscosity.h"
#include "Kernel/Compressible/NavierStokes/K15/K15CompressibleNavierStokesSponge.h"
#include "Kernel/Compressible/NavierStokes/K17/K17CompressibleNavierStokes.h"
#include "Kernel/Compressible/NavierStokes/K17/K17CompressibleNavierStokesChimeraLegacy.h"
#include "Kernel/Compressible/NavierStokes/K17/K17CompressibleNavierStokesBulkViscosity.h"
#include "Kernel/Compressible/NavierStokes/K17/K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants.h"
#include "Kernel/Compressible/NavierStokes/K18/K18CompressibleNavierStokes.h"
#include "Kernel/Compressible/NavierStokes/K20/K20CompressibleNavierStokes.h"

//LBM kernel (inkompressible)
#include "Kernel/Incompressible/NavierStokes/B92/B92IncompressibleNavierStokes.h"
#include "Kernel/Incompressible/NavierStokes/B15/B15IncompressibleNavierStokesBGKplus.h"
#include "Kernel/Incompressible/NavierStokes/C06/C06IncompressibleNavierStokes.h"
#include "Kernel/Incompressible/NavierStokes/K15/K15IncompressibleNavierStokesRotatingVelocityField.h"
#include "Kernel/Incompressible/NavierStokes/K15/K15IncompressibleNavierStokesIsoCheck.h"
#include "Kernel/Incompressible/NavierStokes/K15/K15IncompressibleNavierStokes.h"
#include "Kernel/Incompressible/NavierStokes/M02/M02IncompressibleNavierStokes.h"

//advection diffusion kernel (compressible)
#include "Kernel/Compressible/AdvectionDiffusion/D3Q27/F16/F16CompressibleAdvectionDiffusion.h"
#include "Kernel/Compressible/AdvectionDiffusion/D3Q7/B12/B12CompressibleAdvectionDiffusionD3Q7.h"

//advection diffusion kernel (incompressible)
#include "Kernel/Incompressible/AdvectionDiffusion/D3Q27/F16/F16IncompressibleAdvectionDiffusion.h"
#include "Kernel/Incompressible/AdvectionDiffusion/D3Q7/B12/B12IncompressibleAdvectionDiffusionD3Q7.h"

//porous media kernel
#include "Kernel/Compressible/NavierStokes/K15PorousMedia/K15CompressibleNavierStokesPorousMedia.h"

//wale kernel+
#include "Kernel/Compressible/NavierStokes/K17Wale/K17CompressibleNavierStokesWale.h"
#include "Kernel/Compressible/NavierStokes/K17WaleDebug/K17CompressibleNavierStokesWaleDebug.h"
#include "Kernel/Compressible/NavierStokes/K15Wale/K15CompressibleNavierStokesWale.h"
#include "Kernel/Compressible/NavierStokes/K15WaleSM/K15CompressibleNavierStokesWaleBySoniMalav.h"

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

    if (kernel == collisionKernel::compressible::BGK) {
        newKernel     = B92CompressibleNavierStokes::getNewInstance(para, level);               // compressible
    } else if (kernel == collisionKernel::compressible::BGKPlus) {
        newKernel     = B15CompressibleNavierStokesBGKplus::getNewInstance(para, level);
    } else if (kernel == collisionKernel::compressible::MRT) {
        newKernel     = M02CompressibleNavierStokes::getNewInstance(para, level);
    } else if (kernel == collisionKernel::compressible::Cascade) {
        newKernel     = C06CompressibleNavierStokes::getNewInstance(para, level);
    } else if (kernel == collisionKernel::compressible::CumulantClassic) {
        newKernel     = K08CompressibleNavierStokes::getNewInstance(para, level);
    } else if (kernel == collisionKernel::compressible::K17CompressibleNavierStokesBulkViscosity) {
        newKernel     = K17CompressibleNavierStokesBulkViscosity::getNewInstance(para, level);
    } else if (kernel == collisionKernel::compressible::K17CompressibleNavierStokesChimeraLegacy) {
        newKernel     = K17CompressibleNavierStokesChimeraLegacy::getNewInstance(para, level);
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
    } else if (kernel == collisionKernel::compressible::CumulantAll4SP27) {
        newKernel     = K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants::getNewInstance(para, level);
    } else if (kernel == collisionKernel::compressible::CumulantK18) {
        newKernel     = K18CompressibleNavierStokes::getNewInstance(para, level);
    } else if (kernel == collisionKernel::compressible::CumulantK20) {
        newKernel     = K20CompressibleNavierStokes::getNewInstance(para, level);
    } else if (kernel == collisionKernel::compressible::K15CompressibleNavierStokes) {
        newKernel     = K15CompressibleNavierStokes::getNewInstance(para, level);
    } else if (kernel == collisionKernel::compressible::K15CompressibleNavierStokesBulk) {
        newKernel     = K15CompressibleNavierStokesBulkViscosity::getNewInstance(para, level);
    } else if (kernel == collisionKernel::compressible::K15CompressibleNavierStokesSponge) { //     /\      //
        newKernel     = K15CompressibleNavierStokesSponge::getNewInstance(para, level);     //     ||
    }                                                                           //===============
    else if (  kernel == collisionKernel::incompressible::BGK) {                // incompressible
        newKernel     = B92IncompressibleNavierStokes::getNewInstance(para, level);             //     ||
    } else if (kernel == collisionKernel::incompressible::BGKPlus) {
        newKernel     = B15IncompressibleNavierStokesBGKplus::getNewInstance(para, level);
    } else if (kernel == collisionKernel::incompressible::MRT) {
        newKernel     = M02IncompressibleNavierStokes::getNewInstance(para, level);
    } else if (kernel == collisionKernel::incompressible::Cascade) {
        newKernel     = C06IncompressibleNavierStokes::getNewInstance(para, level);
    } else if (kernel == collisionKernel::incompressible::Cumulant1h) {
        newKernel     = K15IncompressibleNavierStokesRotatingVelocityField::getNewInstance(para, level);
    } else if (kernel == collisionKernel::incompressible::CumulantK15) {          //     /\      //
        newKernel     = K15IncompressibleNavierStokes::getNewInstance(para, level);           //     ||
    }                                                                             //===============
    else if (kernel == collisionKernel::porousMedia::CumulantOne) {               // porous media
        newKernel     = K15CompressibleNavierStokesPorousMedia::getNewInstance(para, pm, level);   //     ||
    }                                                                             //===============
    else if (kernel == collisionKernel::wale::CumulantK17) {                      // wale model
        newKernel     = K17CompressibleNavierStokesWale::getNewInstance(para, level);         //     ||
    } else if (kernel == collisionKernel::wale::CumulantK17Debug) {
        newKernel     = K17CompressibleNavierStokesWaleDebug::getNewInstance(para, level);
    } else if (kernel == collisionKernel::wale::CumulantK15) {
        newKernel     = K15CompressibleNavierStokesWale::getNewInstance(para, level);
    } else if (kernel == collisionKernel::wale::CumulantK15SoniMalav) {              //     /\      //
        newKernel     = K15CompressibleNavierStokesWaleBySoniMalav::getNewInstance(para, level); //     ||
    }                                                                                //===============
    else {
        throw std::runtime_error("KernelFactory does not know the KernelType.");
    }
    para->setKernelNeedsFluidNodeIndicesToRun(newKernel->getKernelUsesFluidNodeIndices());
    return newKernel;
}

std::shared_ptr<ADKernel> KernelFactoryImp::makeAdvDifKernel(std::shared_ptr<Parameter> para, std::string kernel, int level)
{
    std::shared_ptr<ADKernel> newKernel;

    if (kernel == "ADComp27") {
        newKernel     = F16CompressibleAdvectionDiffusion::getNewInstance(para, level);
    } else if(kernel == "ADComp7") {
        newKernel     = B12CompressibleAdvectionDiffusionD3Q7::getNewInstance(para, level);
    } else if (kernel == "ADIncomp27") {
        newKernel     = F16IncompressibleAdvectionDiffusion::getNewInstance(para, level);
    } else if (kernel == "ADIncomp7") {
        newKernel     = B12IncompressibleAdvectionDiffusionD3Q7::getNewInstance(para, level);
    } else {
        throw std::runtime_error("KernelFactory does not know the KernelType.");
    }

    return newKernel;
}
