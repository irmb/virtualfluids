#ifndef CUDA_DEFINES_H
#define CUDA_DEFINES_H

#include "cuda_runtime.h"
#include <stdio.h>

#define HOST __host__
#define DEVICE __device__
#define GLOBAL __global__
#define CONSTANT __constant__


#define HOSTDEVICE HOST DEVICE 

static void printCudaInformation(int i) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    printf(" --- General Information for device %d ---\n", i);
    printf("Name: %s\n", prop.name);
    printf("Compute capability: %d.%d\n", prop.major, prop.minor);
    printf("Clock rate: %d\n", prop.clockRate);
    printf("Device copy overlap: ");
    if (prop.deviceOverlap)
        printf("Enabled\n");
    else
        printf("Disabled\n");
    printf("Kernel execition timeout : ");
    if (prop.kernelExecTimeoutEnabled)
        printf("Enabled\n");
    else
        printf("Disabled\n");
    printf(" --- Memory Information for device %d ---\n", i);
    printf("Total global mem: %llu\n", prop.totalGlobalMem);
    printf("Total constant Mem: %zd\n", prop.totalConstMem);
    printf("Max mem pitch: %zd\n", prop.memPitch);
    printf("Texture Alignment: %zd\n", prop.textureAlignment);
    printf("max Texture 1D: %ld\n", prop.maxTexture1D);
    printf("max Texture 2D: %ld, %ld\n", prop.maxTexture2D[0], prop.maxTexture2D[1]);
    printf("max Texture 3D: %ld, %ld, %ld\n", prop.maxTexture3D[0], prop.maxTexture3D[1], prop.maxTexture3D[2]);
    printf(" --- MP Information for device %d ---\n", i);
    printf("Multiprocessor count: %d\n",
        prop.multiProcessorCount);
    printf("Shared mem per mp: %zd\n", prop.sharedMemPerBlock);
    printf("Registers per mp: %d\n", prop.regsPerBlock);
    printf("Threads in warp: %d\n", prop.warpSize);
    printf("Max threads per block: %d\n",
        prop.maxThreadsPerBlock);
    printf("Max thread dimensions: (%d, %d, %d)\n",
        prop.maxThreadsDim[0], prop.maxThreadsDim[1],
        prop.maxThreadsDim[2]);
    printf("Max grid dimensions: (%d, %d, %d)\n",
        prop.maxGridSize[0], prop.maxGridSize[1],
        prop.maxGridSize[2]);
    printf(" --- -------------------------------- ---\n");
    printf("\n");

    cudaSetDevice(i);
    size_t free;
    size_t total;
    cudaMemGetInfo(&free, &total);
    printf("Free: %llu Bytes, Total: %llu Bytes\n", free, total);
    printf("Free: %llu MB, Total: %llu MB\n", free / 1000 / 1000, total / 1000 / 1000);
    //cudaDeviceSetLimit(cudaLimitMallocHeapSize, free);
}

#endif 
