#include "KernelFactoryImp.h"

#include "Parameter/Parameter.h"

#include "Kernel/Kernels/BasicKernels/Advection/Compressible/BGK/BGKCompSP27.h"
#include "Kernel/Kernels/BasicKernels/Advection/Compressible/BGKPlus/BGKPlusCompSP27.h"
#include "Kernel/Kernels/BasicKernels/Advection/Compressible/Cascade/CascadeCompSP27.h"
#include "Kernel/Kernels/BasicKernels/Advection/Compressible/Cumulant/CumulantCompSP27.h"
#include "Kernel/Kernels/BasicKernels/Advection/Compressible/CumulantK17/CumulantK17Comp.h"
#include "Kernel/Kernels/BasicKernels/Advection/Compressible/CumulantK17Bulk/CumulantK17BulkComp.h"
#include "Kernel/Kernels/BasicKernels/Advection/Compressible/CumulantAll4/CumulantAll4CompSP27.h"
#include "Kernel/Kernels/BasicKernels/Advection/Compressible/CumulantK18/CumulantK18Comp.h"
#include "Kernel/Kernels/BasicKernels/Advection/Compressible/CumulantK20/CumulantK20Comp.h"
#include "Kernel/Kernels/BasicKernels/Advection/Compressible/CumulantK15/CumulantK15Comp.h"
#include "Kernel/Kernels/BasicKernels/Advection/Compressible/CumulantK15Bulk/CumulantK15BulkComp.h"
#include "Kernel/Kernels/BasicKernels/Advection/Compressible/CumulantK15Sponge/CumulantK15SpongeComp.h"
#include "Kernel/Kernels/BasicKernels/Advection/Compressible/MRT/MRTCompSP27.h"

#include "Kernel/Kernels/BasicKernels/Advection/Incompressible/BGK/BGKIncompSP27.h"
#include "Kernel/Kernels/BasicKernels/Advection/Incompressible/BGKPlus/BGKPlusIncompSP27.h"
#include "Kernel/Kernels/BasicKernels/Advection/Incompressible/Cascade/CascadeIncompSP27.h"
#include "Kernel/Kernels/BasicKernels/Advection/Incompressible/Cumulant1hSP27/Cumulant1hIncompSP27.h"
#include "Kernel/Kernels/BasicKernels/Advection/Incompressible/CumulantIsoSP27/CumulantIsoIncompSP27.h"
#include "Kernel/Kernels/BasicKernels/Advection/Incompressible/CumulantK15/CumulantK15Incomp.h"
#include "Kernel/Kernels/BasicKernels/Advection/Incompressible/MRT/MRTIncompSP27.h"

#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Compressible/Mod27/ADComp27/ADComp27.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Compressible/Mod7/ADComp7/ADComp7.h"

#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Incompressible/Mod27/ADIncomp27/ADIncomp27.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Incompressible/Mod7/ADIncomp7/ADIncomp7.h"

#include "Kernel/Kernels/PorousMediaKernels/Advection/Compressible/CumulantOne/PMCumulantOneCompSP27.h"

#include "Kernel/Kernels/WaleKernels/Advection/Compressible/CumulantK17/WaleCumulantK17Comp.h"
#include "Kernel/Kernels/WaleKernels/Advection/Compressible/CumulantK17Debug/WaleCumulantK17DebugComp.h"
#include "Kernel/Kernels/WaleKernels/Advection/Compressible/CumulantK15/WaleCumulantK15Comp.h"
#include "Kernel/Kernels/WaleKernels/Advection/Compressible/CumulantK15BySoniMalav/WaleBySoniMalavCumulantK15Comp.h"

#include "Kernel/Kernels/BasicKernels/Advection/Compressible/AdvecCompStrategy.h"
#include "Kernel/Kernels/BasicKernels/Advection/Incompressible/AdvecIncompStrategy.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Compressible/Mod27/ADMod27CompStrategy.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Compressible/Mod7/ADMod7CompStrategy.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Incompressible/Mod27/ADMod27IncompStrategy.h"
#include "Kernel/Kernels/BasicKernels/AdvectionDiffusion/Incompressible/Mod7/ADMod7IncompStrategy.h"
#include "Kernel/Kernels/PorousMediaKernels/Advection/Compressible/PMAdvecCompStrategy.h"
#include "Kernel/Kernels/WaleKernels/Advection/Compressible/WaleAdvecCompStrategy.h"


std::shared_ptr<KernelFactoryImp> KernelFactoryImp::getInstance()
{
	static std::shared_ptr<KernelFactoryImp> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<KernelFactoryImp>(new KernelFactoryImp());
	return uniqueInstance;
}

std::vector<std::shared_ptr<Kernel>> KernelFactoryImp::makeKernels(std::shared_ptr<Parameter> para)
{
	std::vector< std::shared_ptr< Kernel>> kernels;
	for (int level = 0; level <= para->getMaxLevel(); level++)
		kernels.push_back(makeKernel(para, para->getMainKernel(), level));

	if (para->getMaxLevel() > 0)
		if (para->getMultiKernelOn())
			for (int i = 0; i < para->getMultiKernelLevel().size(); i++)
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

void KernelFactoryImp::setKernelAtLevel(std::vector<std::shared_ptr<Kernel>> kernels, std::shared_ptr<Parameter> para, KernelType kernel, int level)
{
	kernels.at(level) = makeKernel(para, kernel, level);
}

std::shared_ptr<Kernel> KernelFactoryImp::makeKernel(std::shared_ptr<Parameter> para, KernelType kernel, int level)
{
	std::shared_ptr<KernelImp> newKernel;
	std::shared_ptr<CheckParameterStrategy> checkStrategy;

	switch (kernel)
	{
	case LB_BGKCompSP27:
		newKernel = BGKCompSP27::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;
	case LB_BGKPlusCompSP27:
		newKernel = BGKPlusCompSP27::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;
	case LB_CascadeCompSP27:
		newKernel = CascadeCompSP27::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;
	case LB_CumulantCompSP27:
		newKernel = CumulantCompSP27::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;
	case LB_CumulantK17Comp:
		newKernel = CumulantK17Comp::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;
	case LB_CumulantK17BulkComp:
		newKernel = CumulantK17BulkComp::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;
	case LB_CumulantAll4CompSP27:
		newKernel = CumulantAll4CompSP27::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;
	case LB_CumulantK18Comp:
		newKernel = CumulantK18Comp::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;
	case LB_CumulantK20Comp:
		newKernel = CumulantK20Comp::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;
	case LB_CumulantK15Comp:
		newKernel = CumulantK15Comp::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;
	case LB_CumulantK15BulkComp:
		newKernel = CumulantK15Comp::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;
	case LB_CumulantK15SpongeComp:
		newKernel = CumulantK15SpongeComp::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;
	case LB_MRTCompSP27:
		newKernel = MRTCompSP27::getNewInstance(para, level);
		checkStrategy = AdvecCompStrategy::getInstance();
		break;


	case LB_BGKIncompSP27:
		newKernel = BGKIncompSP27::getNewInstance(para, level);
		checkStrategy = AdvecIncompStrategy::getInstance();
		break;
	case LB_BGKPlusIncompSP27:
		newKernel = BGKPlusIncompSP27::getNewInstance(para, level);
		checkStrategy = AdvecIncompStrategy::getInstance();
		break;
	case LB_CascadeIncompSP27:
		newKernel = CascadeIncompSP27::getNewInstance(para, level);
		checkStrategy = AdvecIncompStrategy::getInstance();
		break;
	case LB_Cumulant1hIncompSP27:
		newKernel = Cumulant1hIncompSP27::getNewInstance(para, level);
		checkStrategy = AdvecIncompStrategy::getInstance();
		break;
	case LB_CumulantIsoIncompSP27:
		newKernel = CumulantIsoIncompSP27::getNewInstance(para, level);
		checkStrategy = AdvecIncompStrategy::getInstance();
		break;
	case LB_CumulantK15Incomp:
		newKernel = CumulantK15Incomp::getNewInstance(para, level);
		checkStrategy = AdvecIncompStrategy::getInstance();
		break;
	case LB_MRTIncompSP27:
		newKernel = MRTIncompSP27::getNewInstance(para, level);
		checkStrategy = AdvecIncompStrategy::getInstance();
		break;

	case LB_PMCumulantOneCompSP27:
		newKernel = PMCumulantOneCompSP27::getNewInstance(para, pm, level);
		checkStrategy = PMAdvecCompStrategy::getInstance();
		break;



	case LB_WaleCumulantK17Comp:
		newKernel = WaleCumulantK17Comp::getNewInstance(para, level);
		checkStrategy = WaleAdvecCompStrategy::getInstance();
		break;
	case LB_WaleCumulantK17DebugComp:
		newKernel = WaleCumulantK17DebugComp::getNewInstance(para, level);
		checkStrategy = WaleAdvecCompStrategy::getInstance();
		break;
	case LB_WaleCumulantK15Comp:
		newKernel = WaleCumulantK15Comp::getNewInstance(para, level);
		checkStrategy = WaleAdvecCompStrategy::getInstance();
		break;
	case LB_WaleBySoniMalavCumulantK15Comp:
		newKernel = WaleBySoniMalavCumulantK15Comp::getNewInstance(para, level);
		checkStrategy = WaleAdvecCompStrategy::getInstance();
		break;
	default:
		break;
	}

	if (newKernel) {
		newKernel->setCheckParameterStrategy(checkStrategy);
		return newKernel;
	}
	else
		throw  std::runtime_error("KernelFactory does not know the KernelType.");

	
}

std::shared_ptr<ADKernel> KernelFactoryImp::makeAdvDifKernel(std::shared_ptr<Parameter> para, ADKernelType kernel, int level)
{
	std::shared_ptr<ADKernel> newKernel;
	std::shared_ptr<CheckParameterStrategy> checkStrategy;

	switch (kernel)
	{
	case LB_ADComp27:
		newKernel = ADComp27::getNewInstance(para, level);
		checkStrategy = ADMod27CompStrategy::getInstance();
		break;
	case LB_ADComp7:
		newKernel = ADComp7::getNewInstance(para, level);
		checkStrategy = ADMod7CompStrategy::getInstance();
		break;
	case LB_ADIncomp27:
		newKernel = ADIncomp27::getNewInstance(para, level);
		checkStrategy = ADMod27IncompStrategy::getInstance();
		break;
	case LB_ADIncomp7:
		newKernel = ADIncomp7::getNewInstance(para, level);
		checkStrategy = ADMod7IncompStrategy::getInstance();
		break;
	default:
		break;
	}

	if (newKernel) {
		newKernel->setCheckParameterStrategy(checkStrategy);
		return newKernel;
	}
	else
		throw  std::runtime_error("KernelFactory does not know the KernelType.");
}

KernelFactoryImp::KernelFactoryImp()
{

}
