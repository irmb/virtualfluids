#include "GridGpuBuilder.h"

#include <GridGenerator/grid/GridWrapper/GridWrapperGPU/GridWrapperGPU.h>
#include <GridGenerator/geometries/BoundingBox/BoundingBox.cuh>


GridGpuBuilder::GridGpuBuilder() : GridBuilderImp()
{

}

GridGpuBuilder::~GridGpuBuilder()
{

}

void GridGpuBuilder::createGridKernels(std::string distribution)
{
    for (int i = 0; i < rankTasks.size(); i += 2)
    {
        int level = rankTasks[i];
        int index = rankTasks[i + 1];
        this->gridKernels[level][index] = std::shared_ptr<GridWrapperGPU>(new GridWrapperGPU(this->boxes[level][index], distribution));
    }
}
