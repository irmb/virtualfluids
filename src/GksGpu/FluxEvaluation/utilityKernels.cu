#include <cstdio>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "DataBase/DataBaseStruct.h"

#define THREADS_PER_BLOCK 128

__global__ void dummyKernel();

__global__ void debugKernel( const DataBaseStruct dataBase );

//////////////////////////////////////////////////////////////////////////

__global__ void dummyKernel(  )
{
    printf("I am thread %d.\n", threadIdx.x);
}

__global__ void debugKernel( const DataBaseStruct dataBase )
{   
    if( threadIdx.x + blockIdx.x == 0 ){
        printf("numberOfCells     : %d\n", dataBase.numberOfCells     );
        printf("numberOfFaces     : %d\n", dataBase.numberOfFaces     );
        printf("\n");
        printf("\n");
        printf("faceToCell: %p\n", dataBase.faceToCell);
        printf("faceCenter: %p\n", dataBase.faceCenter);
        printf("cellCenter: %p\n", dataBase.cellCenter);
        printf("\n");
        printf("data      : %p\n", dataBase.data      );
        printf("dataNew   : %p\n", dataBase.dataUpdate);
    }
}
