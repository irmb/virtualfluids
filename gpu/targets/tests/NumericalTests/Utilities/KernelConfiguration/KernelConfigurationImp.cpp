#include "KernelConfigurationImp.h"

KernelType KernelConfigurationImp::getMainKernel()
{
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

std::vector<KernelType> KernelConfigurationImp::getMultiKernel()
{
	return multiKernel;
}

std::shared_ptr<KernelConfigurationImp> KernelConfigurationImp::getNewInstance(KernelType kernelName)
{
	return std::shared_ptr<KernelConfigurationImp>(new KernelConfigurationImp(kernelName));
}

std::shared_ptr<KernelConfigurationImp> KernelConfigurationImp::getNewInstance(KernelType kernel, std::vector<int> multiKernelLevel, std::vector<KernelType> multiKernelName)
{
	return std::shared_ptr<KernelConfigurationImp>(new KernelConfigurationImp(kernel, multiKernelLevel, multiKernelName));
}

KernelConfigurationImp::KernelConfigurationImp(KernelType kernel)
{
	this->mainKernel = kernel;
	multiKernelOn = false;
	multiKernelLevel.resize(0);
	multiKernel.resize(0);
}

KernelConfigurationImp::KernelConfigurationImp(KernelType mainKernel, std::vector<int> multiKernelLevel, std::vector<KernelType> multiKernel)
{
	this->mainKernel = mainKernel;
	multiKernelOn = true;
	this->multiKernelLevel = multiKernelLevel;
	this->multiKernel = multiKernel;
}