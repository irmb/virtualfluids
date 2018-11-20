#ifndef  CudaUtilExtern_H
#define  CudaUtilExtern_H

#include <VirtualFluidsDefinitions.h>

class VF_PUBLIC CudaUtility
{
public:

    static struct CudaGrid 
    {
        dim3 threads;
        dim3 blocks;
    };

    static void printCudaMemoryUsage();

    static void setCudaDevice( int device );
};

#endif
