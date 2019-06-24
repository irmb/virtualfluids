#ifndef KERNEL_FACTORY_IMP_H
#define KERNEL_FACTORY_IMP_H

#include "KernelFactory.h"

class PorousMedia;

class VF_PUBLIC KernelFactoryImp : public KernelFactory
{
public:
	static std::shared_ptr< KernelFactoryImp> getInstance();
	std::vector< std::shared_ptr< Kernel>> makeKernels(std::shared_ptr<Parameter> para);
	std::vector< std::shared_ptr< ADKernel>> makeAdvDifKernels(std::shared_ptr<Parameter> para);

	void setPorousMedia(std::vector<std::shared_ptr<PorousMedia>> pm);

private:
	KernelFactoryImp();

	std::shared_ptr< Kernel> makeKernel(std::shared_ptr<Parameter> para, KernelType kernel, int level);
	std::shared_ptr< ADKernel> makeAdvDifKernel(std::shared_ptr<Parameter> para, ADKernelType kernel, int level);
	void setKernelAtLevel(std::vector< std::shared_ptr<Kernel>> kernels, std::shared_ptr<Parameter> para, KernelType kernel, int level);

	std::vector<std::shared_ptr<PorousMedia>> pm;
};
#endif