#ifndef  CudaRunKernel_HPP
#define  CudaRunKernel_HPP

#include <cuda_runtime.h>

#include "CudaUtility/CudaUtility.h"

template<typename Functor, typename... TArgs>
void runKernel(Functor kernel, const CudaUtility::CudaGrid& grid, TArgs... args)
{
    kernel<<< grid.blocks, grid.threads >>>(args...);
}

#endif
