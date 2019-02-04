#ifndef KERNEL_CONFIGURATION_H
#define KERNEL_CONFIGURATION_H

#include <vector>

class KernelConfiguration
{
public:
	virtual std::string getMainKernel() = 0;
	virtual bool getMultiKernelOn() = 0;
	virtual	std::vector<int> getMultiKernelLevel() = 0;
	virtual std::vector<std::string> getMultiKernelName() = 0;

private:	

};
#endif 