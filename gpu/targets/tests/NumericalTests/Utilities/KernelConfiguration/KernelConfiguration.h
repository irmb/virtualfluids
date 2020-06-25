#ifndef KERNEL_CONFIGURATION_H
#define KERNEL_CONFIGURATION_H

#include <vector>

#include "VirtualFluids_GPU/Kernel//Utilities/KernelType.h"

class KernelConfiguration
{
public:
	virtual KernelType getMainKernel() = 0;
	virtual bool getMultiKernelOn() = 0;
	virtual	std::vector<int> getMultiKernelLevel() = 0;
	virtual std::vector<KernelType> getMultiKernel() = 0;

private:	

};
#endif 