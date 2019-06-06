#ifndef KERNEL_CONFIGURATION_IMP_H
#define KERNEL_CONFIGURATION_IMP_H

#include "KernelConfiguration.h"

#include <memory>

class KernelConfigurationImp : public KernelConfiguration
{
public:
	KernelType getMainKernel();
	bool getMultiKernelOn();
	std::vector<int> getMultiKernelLevel();
	std::vector<KernelType> getMultiKernel();

	static std::shared_ptr<KernelConfigurationImp> getNewInstance(KernelType kernel);
	static std::shared_ptr<KernelConfigurationImp> getNewInstance(KernelType kernel, std::vector<int> multiKernelLevel, std::vector<KernelType> multiKernel);

private:
	KernelConfigurationImp(KernelType kernel);
	KernelConfigurationImp(KernelType kernel, std::vector<int> multiKernelLevel, std::vector<KernelType> multiKernel);
	KernelConfigurationImp() {};

	KernelType mainKernel;
	bool multiKernelOn;
	std::vector<int> multiKernelLevel;
	std::vector<KernelType> multiKernel;
};
#endif 