#include "GridCpuBuilder.h"

#include <GridGenerator/grid/GridWrapper/GridWrapperCPU/GridWrapperCPU.h>
#include <GridGenerator/geometries/BoundingBox/BoundingBox.cuh>


GridCpuBuilder::GridCpuBuilder() : GridBuilderImp()
{

}

GridCpuBuilder::~GridCpuBuilder()
{

}

void GridCpuBuilder::createGridKernels(std::string distribution)
{
    for (int i = 0; i < rankTasks.size(); i += 2)
    {
        int level = rankTasks[i];
        int index = rankTasks[i + 1];
        this->gridKernels[level][index] = std::shared_ptr<GridWrapperCPU>(new GridWrapperCPU(this->boxes[level][index], distribution));
    }
}
