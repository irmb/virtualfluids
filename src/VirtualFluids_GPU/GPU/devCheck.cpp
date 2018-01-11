#include "devCheck.h"

#include <stdio.h>
#include <stdlib.h> 

#include <cuda_runtime.h>

  
int devCheck(int gpudevice) 
{ 
   int device_count=0; 
    int device;  // used with  cudaGetDevice() to verify cudaSetDevice() 

   // get the number of non-emulation devices  detected 
   cudaGetDeviceCount( &device_count); 
   if (gpudevice > device_count) 
   { 
       printf("gpudevice >=  device_count ... exiting\n"); 
       exit(1); 
   } 
   cudaError_t cudareturn; 
   cudaDeviceProp deviceProp; 

   // cudaGetDeviceProperties() is also  demonstrated in the deviceQuery/ example
   // of the sdk projects directory 
   cudaGetDeviceProperties(&deviceProp, gpudevice); 
    printf("[deviceProp.major.deviceProp.minor] = [%d.%d]\n", 
      deviceProp.major, deviceProp.minor); 

   if (deviceProp.major > 999) 
   { 
       printf("warning, CUDA Device  Emulation (CPU) detected, exiting\n"); 
               exit(1); 
   } 

   // choose a cuda device for kernel  execution 
   cudareturn=cudaSetDevice(gpudevice); 
   if (cudareturn == cudaErrorInvalidDevice) 
   { 
       perror("cudaSetDevice returned  cudaErrorInvalidDevice"); 
 } 
 else 
 { 
     // double check that device was  properly selected 
     cudaGetDevice(&device); 
     //printf("cudaGetDevice()=%d\n",device); 
     return device;
 } 
   return -1;
} 
