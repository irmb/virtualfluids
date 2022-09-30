#include "KernelFactoryImp.h"

#include "Parameter/Parameter.h"

//LBM kernel (compressible)
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/BGK/BGKCompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/BGKUnified/BGKUnified.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/BGKPlus/BGKPlusCompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/Cascade/CascadeCompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/Cumulant/CumulantCompSP27.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17/CumulantK17Comp.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17Unified/CumulantK17Unified.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17chim/CumulantK17CompChim.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17chimStream/CumulantK17CompChimStream.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17chimRedesigned/CumulantK17CompChimRedesigned.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17Bulk/CumulantK17BulkComp.h"
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

//turbulent viscosity kernel
#include "Kernel/Kernels/TurbulentViscosityKernels/FluidFlow/Compressible/CumulantK17chim/TurbulentViscosityCumulantK17CompChim.h"

//strategies
#include "Kernel/Kernels/BasicKernels/FluidFlow/Compressible/FluidFlowCompStrategy.h"
#include "Kernel/Kernels/BasicKernels/FluidFlow/Incompressible/FluidFlowIncompStrategy.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Compressible/Mod27/ADMod27CompStrategy.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Compressible/Mod7/ADMod7CompStrategy.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Incompressible/Mod27/ADMod27IncompStrategy.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Incompressible/Mod7/ADMod7IncompStrategy.h"
#include "Kernel/Kernels/PorousMediaKernels/FluidFlow/Compressible/PMFluidFlowCompStrategy.h"
#include "Kernel/Kernels/WaleKernels/FluidFlow/Compressible/WaleFluidFlowCompStrategy.h"

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
    printf("Instantiating Kernel: %s\n", kernel.c_str());
	std::shared_ptr<KernelImp> newKernel;
	std::shared_ptr<CheckParameterStrategy> checkStrategy;

    if (kernel == "BGKCompSP27") {
        newKernel     = BGKCompSP27::getNewInstance(para, level);   // compressible
        checkStrategy = FluidFlowCompStrategy::getInstance();       //      ||
    } else if (kernel == "BGKUnified") {                            //      \/
        newKernel     = std::make_shared<vf::gpu::BGKUnified>(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "BGKPlusCompSP27") {
        newKernel     = BGKPlusCompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "MRTCompSP27") {
        newKernel     = MRTCompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CascadeCompSP27") {
        newKernel     = CascadeCompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantCompSP27") {
        newKernel     = CumulantCompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantK17Comp") {
        newKernel     = CumulantK17Comp::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantK15Unified") {
        newKernel     = std::make_shared<vf::gpu::CumulantK15Unified>(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantK17Unified") {
        newKernel     = std::make_shared<vf::gpu::CumulantK17Unified>(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantK17BulkComp") {
        newKernel     = CumulantK17BulkComp::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantK17CompChim") {
        newKernel     = CumulantK17CompChim::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantK17CompChimStream") {
        newKernel     = CumulantK17CompChimStream::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantK17CompChimRedesigned") {
        newKernel     = CumulantK17CompChimRedesigned::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantAll4CompSP27") {
        newKernel     = CumulantAll4CompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantK18Comp") {
        newKernel     = CumulantK18Comp::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantK20Comp") {
        newKernel     = CumulantK20Comp::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantK15Comp") {
        newKernel     = CumulantK15Comp::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantK15BulkComp") {
        newKernel     = CumulantK15BulkComp::getNewInstance(para, level);
        checkStrategy = FluidFlowCompStrategy::getInstance();
    } else if (kernel == "CumulantK15SpongeComp") {                             //     /\      //
        newKernel     = CumulantK15SpongeComp::getNewInstance(para, level);     //	   ||
        checkStrategy = FluidFlowCompStrategy::getInstance();                   // compressible
    }																			//===============
	else if (  kernel == "BGKIncompSP27") {										// incompressible
        newKernel     = BGKIncompSP27::getNewInstance(para, level);				//	   ||
        checkStrategy = FluidFlowIncompStrategy::getInstance();                 //     \/
    } else if (kernel == "BGKPlusIncompSP27") {
        newKernel     = BGKPlusIncompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == "MRTIncompSP27") {
        newKernel     = MRTIncompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == "CascadeIncompSP27") {
        newKernel     = CascadeIncompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == "Cumulant1hIncompSP27") {
        newKernel     = Cumulant1hIncompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == "CumulantIsoIncompSP27") {
        newKernel     = CumulantIsoIncompSP27::getNewInstance(para, level);
        checkStrategy = FluidFlowIncompStrategy::getInstance();
    } else if (kernel == "CumulantK15Incomp") {									//     /\      //
        newKernel     = CumulantK15Incomp::getNewInstance(para, level);			//	   ||
        checkStrategy = FluidFlowIncompStrategy::getInstance();                 // incompressible
    }																			//===============
	else if (kernel == "PMCumulantOneCompSP27") {								// porous media
        newKernel     = PMCumulantOneCompSP27::getNewInstance(para, pm, level);	//	   ||
        checkStrategy = PMFluidFlowCompStrategy::getInstance();                 // porous media
    }                                                                           //===============
    else if (kernel == "WaleCumulantK17Comp") {                                 // wale model
        newKernel     = WaleCumulantK17Comp::getNewInstance(para, level);       //	   ||
        checkStrategy = WaleFluidFlowCompStrategy::getInstance();               //     \/
    } else if (kernel == "WaleCumulantK17DebugComp") {
        newKernel     = WaleCumulantK17DebugComp::getNewInstance(para, level);
        checkStrategy = WaleFluidFlowCompStrategy::getInstance();
    } else if (kernel == "WaleCumulantK15Comp") {
        newKernel     = WaleCumulantK15Comp::getNewInstance(para, level);
        checkStrategy = WaleFluidFlowCompStrategy::getInstance();
    } else if (kernel == "WaleBySoniMalavCumulantK15Comp") {                    //     /\      //
        newKernel     = WaleBySoniMalavCumulantK15Comp::getNewInstance(para, level);// ||
        checkStrategy = WaleFluidFlowCompStrategy::getInstance();               // wale model
    }                                                                          //===============
    else if (kernel == "TurbulentViscosityCumulantK17CompChim"){               // compressible with turbulent viscosity
        switch(para->getTurbulenceModel())                                     //       ||          
        {                                                                      //       \/      //
            case TurbulenceModel::AMD:
                newKernel = TurbulentViscosityCumulantK17CompChim<TurbulenceModel::AMD>::getNewInstance(para, level);   
                break;
            case TurbulenceModel::Smagorinsky:
                newKernel = TurbulentViscosityCumulantK17CompChim<TurbulenceModel::Smagorinsky>::getNewInstance(para, level);  
                break;
            case TurbulenceModel::QR:
                newKernel = TurbulentViscosityCumulantK17CompChim<TurbulenceModel::QR>::getNewInstance(para, level);  
                break;
            case TurbulenceModel::None:
                newKernel = TurbulentViscosityCumulantK17CompChim<TurbulenceModel::None>::getNewInstance(para, level); 
                break;
            default:
                throw std::runtime_error("Unknown turbulence model!");
            break;                                                              
        }                                                                       
        checkStrategy = FluidFlowCompStrategy::getInstance();
                                                                                //     /\      //
                                                                                //     ||    
                                                                                // compressible with turbulent viscosity  
                                                                                //===============         
    }
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
        checkStrategy = ADMod7CompStrategy::getInstance();
    } else if (kernel == "ADIncomp7") {
        newKernel     = ADIncomp7::getNewInstance(para, level);
        checkStrategy = ADMod7CompStrategy::getInstance();
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
