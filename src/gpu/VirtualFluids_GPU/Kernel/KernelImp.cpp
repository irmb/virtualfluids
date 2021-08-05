#include "KernelImp.h"

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


void KernelImp::runOnIndices(const unsigned int *indices, unsigned int size_indices, int stream)
{
    printf("Method not implemented for this Kernel \n");
}

bool KernelImp::checkParameter() { 
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

KernelImp::KernelImp() {}