//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file devCheck.cpp
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
#include "devCheck.h"

#include <stdio.h>
#include <stdlib.h> 

#include <cuda_runtime.h>


int devCheck(int gpudevice)
{
	int device_count = 0;
	int device;  // used with  cudaGetDevice() to verify cudaSetDevice() 

   // get the number of non-emulation devices  detected 
	cudaGetDeviceCount(&device_count);
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
	printf("[compute capability] = [%d.%d]\n",
		deviceProp.major, deviceProp.minor);

	if (deviceProp.major > 999)
	{
		printf("warning, CUDA Device  Emulation (CPU) detected, exiting\n");
		exit(1);
	}

	// choose a cuda device for kernel  execution 
	cudareturn = cudaSetDevice(gpudevice);
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
