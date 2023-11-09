#ifndef KERNEL_FACTORY_IMP_H
#define KERNEL_FACTORY_IMP_H

#include "KernelFactory.h"

class PorousMedia;

class KernelFactoryImp : public KernelFactory
{
public:
	std::vector< std::shared_ptr< Kernel>> makeKernels(std::shared_ptr<Parameter> para);
	std::vector< std::shared_ptr< AdvectionDiffusionKernel>> makeAdvDifKernels(std::shared_ptr<Parameter> para);

	void setPorousMedia(std::vector<std::shared_ptr<PorousMedia>> pm);

	std::shared_ptr<Kernel> makeKernel(std::shared_ptr<Parameter> para, std::string kernel, int level);
	std::shared_ptr<AdvectionDiffusionKernel> makeAdvDifKernel(std::shared_ptr<Parameter> para, std::string kernel, int level);
	void setKernelAtLevel(std::vector< std::shared_ptr<Kernel>> kernels, std::shared_ptr<Parameter> para, std::string kernel, int level);

	std::vector<std::shared_ptr<PorousMedia>> pm;
};
#endif