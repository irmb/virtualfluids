#include "KernelConfigurationImp.h"

std::string KernelConfigurationImp::getMainKernel() {
	return mainKernel;
}

bool KernelConfigurationImp::getMultiKernelOn()
{
	return multiKernelOn;
}

std::vector<int> KernelConfigurationImp::getMultiKernelLevel()
{
	return multiKernelLevel;
}

std::vector<std::string> KernelConfigurationImp::getMultiKernel() {
	return multiKernel;
}

std::shared_ptr<KernelConfigurationImp> KernelConfigurationImp::getNewInstance(std::string kernelName)
{
	return std::shared_ptr<KernelConfigurationImp>(new KernelConfigurationImp(kernelName));
}

std::shared_ptr<KernelConfigurationImp> KernelConfigurationImp::getNewInstance(std::string kernel, std::vector<int> multiKernelLevel, std::vector<std::string> multiKernelName)
{
	return std::shared_ptr<KernelConfigurationImp>(new KernelConfigurationImp(kernel, multiKernelLevel, multiKernelName));
}

KernelConfigurationImp::KernelConfigurationImp(std::string kernel)
{
	this->mainKernel = kernel;
	multiKernelOn = false;
	multiKernelLevel.resize(0);
	multiKernel.resize(0);
}

KernelConfigurationImp::KernelConfigurationImp(std::string mainKernel, std::vector<int> multiKernelLevel, std::vector<std::string> multiKernel)
{
	this->mainKernel = mainKernel;
	multiKernelOn = true;
	this->multiKernelLevel = multiKernelLevel;
	this->multiKernel = multiKernel;
}