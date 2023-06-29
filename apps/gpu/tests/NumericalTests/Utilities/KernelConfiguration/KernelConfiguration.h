#ifndef KERNEL_CONFIGURATION_H
#define KERNEL_CONFIGURATION_H

#include <string>
#include <vector>

class KernelConfiguration
{
public:
    virtual ~KernelConfiguration() = default;
    virtual std::string getMainKernel() = 0;
	virtual bool getMultiKernelOn() = 0;
	virtual	std::vector<int> getMultiKernelLevel() = 0;
    virtual std::vector<std::string> getMultiKernel() = 0;

private:	

};
#endif 