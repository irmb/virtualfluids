#include "CudaGrid.h"



namespace vf
{
namespace gpu
{

CudaGrid::CudaGrid(unsigned int numberOfThreads, unsigned int size_matrix)
{
    int Grid = (size_matrix / numberOfThreads) + 1;
    int Grid1, Grid2;
    if (Grid > 512) {
        Grid1 = 512;
        Grid2 = (Grid / Grid1) + 1;
    } else {
        Grid1 = 1;
        Grid2 = Grid;
    }
    
    grid = dim3(Grid1, Grid2);
    threads = dim3(numberOfThreads, 1, 1);
}

}
}
