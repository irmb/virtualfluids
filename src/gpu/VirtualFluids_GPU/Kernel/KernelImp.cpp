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

std::unique_ptr<std::pair<dim3, dim3>> KernelImp::calcGridDimensions(unsigned int size_Mat, int numberOfThreads)
{
    int Grid = (size_Mat / numberOfThreads) + 1;
    int Grid1, Grid2;
    if (Grid > 512) {
        Grid1 = 512;
        Grid2 = (Grid / Grid1) + 1;
    } else {
        Grid1 = 1;
        Grid2 = Grid;
    }
    dim3 grid(Grid1, Grid2);
    dim3 threads(numberOfThreads, 1, 1);
    std::pair<dim3, dim3> dimensions(grid, threads);
    return std::make_unique<std::pair<dim3, dim3>>(dimensions);
}
