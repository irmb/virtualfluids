#ifndef GPU_CUDA_GRID_H
#define GPU_CUDA_GRID_H


#include <cuda_runtime.h>

namespace vf
{
namespace gpu
{


struct CudaGrid 
{
    dim3 threads;
    dim3 grid;

    CudaGrid(unsigned int numberOfEntities, unsigned int threadsPerBlock);
    CudaGrid() = default;
};

}
}

#endif