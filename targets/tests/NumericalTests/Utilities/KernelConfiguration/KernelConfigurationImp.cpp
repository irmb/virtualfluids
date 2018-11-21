#include "KernelConfigurationImp.h"

std::string KernelConfigurationImp::getMainKernel()
{
	return mainKernelName;
}

bool KernelConfigurationImp::getMultiKernelOn()
{
	return multiKernelOn;
}

std::vector<int> KernelConfigurationImp::getMultiKernelLevel()
{
	return multiKernelLevel;
}

std::vector<std::string> KernelConfigurationImp::getMultiKernelName()
{
	return multiKernelName;
}

std::shared_ptr<KernelConfigurationImp> KernelConfigurationImp::getNewInstance(std::string kernelName)
{
	return std::shared_ptr<KernelConfigurationImp>(new KernelConfigurationImp(kernelName));
}

std::shared_ptr<KernelConfigurationImp> KernelConfigurationImp::getNewInstance(std::string kernelName, std::vector<int> multiKernelLevel, std::vector<std::string> multiKernelName)
{
	return std::shared_ptr<KernelConfigurationImp>(new KernelConfigurationImp(kernelName, multiKernelLevel, multiKernelName));
}

KernelConfigurationImp::KernelConfigurationImp(std::string kernelName)
{
	this->mainKernelName = kernelName;
	multiKernelOn = false;
	multiKernelLevel.resize(0);
	multiKernelName.resize(0);
}

KernelConfigurationImp::KernelConfigurationImp(std::string mainKernelName, std::vector<int> multiKernelLevel, std::vector<std::string> multiKernelName)
{
	this->mainKernelName = mainKernelName;
	multiKernelOn = true;
	this->multiKernelLevel = multiKernelLevel;
	this->multiKernelName = multiKernelName;
}