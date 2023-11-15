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
    std::vector<std::string> getMultiKernel();

    static std::shared_ptr<KernelConfigurationImp> getNewInstance(std::string kernel);
    static std::shared_ptr<KernelConfigurationImp> getNewInstance(std::string kernel, std::vector<int> multiKernelLevel, std::vector<std::string> multiKernel);

private:
    KernelConfigurationImp(std::string kernel);
    KernelConfigurationImp(std::string kernel, std::vector<int> multiKernelLevel, std::vector<std::string> multiKernel);
    KernelConfigurationImp() {};

    std::string mainKernel;
    bool multiKernelOn;
    std::vector<int> multiKernelLevel;
    std::vector<std::string> multiKernel;
};
#endif 