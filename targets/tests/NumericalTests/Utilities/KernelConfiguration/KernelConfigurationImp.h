#ifndef KERNEL_CONFIGURATION_IMP_H
#define KERNEL_CONFIGURATION_IMP_H

#include "KernelConfiguration.h"

#include <memory>

class KernelConfigurationImp : public KernelConfiguration
{
public:
	std::string getMainKernel();
	bool getMultiKernelOn();
	std::vector<int> getMultiKernelLevel();
	std::vector<std::string> getMultiKernelName();

	static std::shared_ptr<KernelConfigurationImp> getNewInstance(std::string kernelName);
	static std::shared_ptr<KernelConfigurationImp> getNewInstance(std::string mainKernelName, std::vector<int> multiKernelLevel, std::vector<std::string> multiKernelName);

private:
	KernelConfigurationImp(std::string kernelName);
	KernelConfigurationImp(std::string kernelName, std::vector<int> multiKernelLevel, std::vector<std::string> multiKernelName);
	KernelConfigurationImp() {};

	std::string mainKernelName;
	bool multiKernelOn;
	std::vector<int> multiKernelLevel;
	std::vector<std::string> multiKernelName;
};
#endif 