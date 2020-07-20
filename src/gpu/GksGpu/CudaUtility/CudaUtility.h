#ifndef  CudaUtilExtern_H
#define  CudaUtilExtern_H

#include <cuda.h>
#include <cuda_runtime.h>

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"

namespace GksGpu {

class VIRTUALFLUIDS_GPU_EXPORT CudaUtility
{
public:

    struct CudaGrid 
    {
        dim3 threads;
        dim3 blocks;

        uint numberOfEntities;

        cudaStream_t stream;

        CudaGrid( uint numberOfEntities, uint threadsPerBlock, cudaStream_t stream = 0 );
    };

    static cudaStream_t computeStream;
    static cudaStream_t communicationStream;

    static void printCudaMemoryUsage();

    static int getCudaDeviceCount();

    static void setCudaDevice( int device );

    static int getCudaDevice(  );

    static void synchronizeCudaDevice();

    static void synchronizeCudaStream( cudaStream_t stream );
};

} // namespace GksGpu

#endif
