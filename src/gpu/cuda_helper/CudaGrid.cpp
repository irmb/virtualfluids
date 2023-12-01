#include "CudaGrid.h"

#include <logger/Logger.h>

namespace vf::cuda
{

CudaGrid::CudaGrid(unsigned int numberOfThreads, unsigned int numberOfEntities)
    : threads { dim3(numberOfThreads, 1, 1) }, grid { getCudaGrid(numberOfThreads, numberOfEntities) }
{
}
void CudaGrid::print() const
{
    VF_LOG_INFO("blocks: ({},{},{}), threads: ({},{},{})", grid.x, grid.y, grid.z, threads.x, threads.y, threads.z);
}

dim3 getCudaGrid(const unsigned int numberOfThreads, const unsigned int numberOfEntities)
{
    unsigned int Grid = (numberOfEntities / numberOfThreads) + 1;
    unsigned int Grid1, Grid2;
    if (Grid > 512) {
        Grid1 = 512;
        Grid2 = (Grid / Grid1) + 1;
    } else {
        Grid1 = 1;
        Grid2 = Grid;
    }
    dim3 cudaGrid(Grid1, Grid2);
    return cudaGrid;
}

} // namespace vf::cuda
