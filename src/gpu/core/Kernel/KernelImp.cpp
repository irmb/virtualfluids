#include "KernelImp.h"

#include "Calculation/Calculation.h" 


void KernelImp::runOnIndices(const unsigned int *indices, unsigned int size_indices, CollisionTemplate collisionTemplate, CudaStreamIndex streamIndex)
{
    printf("Method not implemented for this Kernel \n");
}

std::vector<PreProcessorType> KernelImp::getPreProcessorTypes() 
{ 
    return myPreProcessorTypes;
}

bool KernelImp::getKernelUsesFluidNodeIndices(){
    return this->kernelUsesFluidNodeIndices;
}

KernelImp::KernelImp(std::shared_ptr<Parameter> para, int level) : para(para), level(level) {}

KernelImp::KernelImp() {}