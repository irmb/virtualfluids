#ifndef CUDA_GRID_H
#define CUDA_GRID_H


#include <cuda_runtime.h>

namespace vf::cuda
{

struct CudaGrid
{
    dim3 threads;
    dim3 grid;

    CudaGrid(unsigned int numberOfThreads, unsigned int numberOfEntities);
    CudaGrid() = default;

    void print() const;
};

dim3 getCudaGrid(const unsigned int numberOfThreads, const unsigned int numberOfEntities);

} // namespace vf::cuda

#endif
