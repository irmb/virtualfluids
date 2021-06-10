#include "KernelImp.h"

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"

bool KernelImp::checkParameter() 
{ 
    return checkStrategy->checkParameter(para);
}

std::vector<PreProcessorType> KernelImp::getPreProcessorTypes() 
{ 
    return myPreProcessorTypes;
}

KernelGroup KernelImp::getKernelGroup() 
{ 
    return myKernelGroup; 
}

void KernelImp::setCheckParameterStrategy(std::shared_ptr<CheckParameterStrategy> strategy)
{
    this->checkStrategy = strategy;
}

KernelImp::KernelImp(std::shared_ptr<Parameter> para, int level) : para(para), level(level) {}