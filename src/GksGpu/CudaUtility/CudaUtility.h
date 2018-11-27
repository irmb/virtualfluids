#ifndef  CudaUtilExtern_H
#define  CudaUtilExtern_H

#include <cuda.h>
#include <cuda_runtime.h>

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"

class VF_PUBLIC CudaUtility
{
public:

    struct CudaGrid 
    {
        dim3 threads;
        dim3 blocks;

        uint numberOfEntities;

        CudaGrid( uint numberOfEntities, uint threadsPerBlock );
    };

    static void printCudaMemoryUsage();

    static void setCudaDevice( int device );
};

#endif
