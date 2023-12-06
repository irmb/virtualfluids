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
//! \author Soeren Peters
//=======================================================================================
#include "DeviceInfo.h"

#include <cstdio>
#include <cstdlib>

#include <cuda_runtime.h>

#include <logger/Logger.h>

namespace vf::cuda
{

void verifyNumberOfDevices(int deviceId)
{
    int device_count = 0;
    cudaError_t errorId = cudaGetDeviceCount(&device_count);
    if(errorId != cudaSuccess) {
        VF_LOG_CRITICAL("Device {}: Error while accessing the device count: {}", deviceId, cudaGetErrorString(errorId));
    }
    if (deviceId > device_count) {
        throw std::runtime_error("chosen gpudevice >=  device_count ... exiting\n");
    }
}

void verifyComputeCapability(int deviceId)
{
    cudaDeviceProp deviceProp;
    cudaError_t errorId = cudaGetDeviceProperties(&deviceProp, deviceId);

    if(errorId != cudaSuccess){
        VF_LOG_CRITICAL("Device {}: Error while accessing the device properties occurs: {}", deviceId, cudaGetErrorString(errorId));
    }

    VF_LOG_INFO("[compute capability] = [{}.{}]\n", deviceProp.major, deviceProp.minor);

    if (deviceProp.major > 999) {
        throw std::runtime_error("Warning, CUDA Device Emulation (CPU) detected, exiting\n");
    }
}

void setCudaDevice(int deviceId)
{
    // choose a cuda device for kernel execution
    cudaError_t errorId = cudaSetDevice(deviceId);
    if (errorId != cudaSuccess) {
        VF_LOG_CRITICAL("Device {}: Error while setting the device to: {}", deviceId, cudaGetErrorString(errorId));
    } else {
        int device;
        // double check that device was properly selected
        errorId = cudaGetDevice(&device);
        if(errorId != cudaSuccess) {
            VF_LOG_CRITICAL("Device {}: Error while getting the device: {}", deviceId, cudaGetErrorString(errorId));
        }
    }
}

void verifyAndSetDevice(int deviceId)
{
    verifyNumberOfDevices(deviceId);
    verifyComputeCapability(deviceId);

    setCudaDevice(deviceId);
}



void printCudaInformation(int deviceId) 
{
    cudaDeviceProp prop;
    cudaError_t errorId = cudaGetDeviceProperties(&prop, deviceId);

    if(errorId != cudaSuccess){
        VF_LOG_CRITICAL("Device {}: Error while accessing the device properties occurs: {}", deviceId, cudaGetErrorString(errorId));
    }

    printf(" --- General Information for device %d ---\n", deviceId);
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
    printf(" --- Memory Information for device %d ---\n", deviceId);
    printf("Total global mem: %zu\n", prop.totalGlobalMem);
    printf("Total constant Mem: %zu\n", prop.totalConstMem);
    printf("Max mem pitch: %zu\n", prop.memPitch);
    printf("Texture Alignment: %zu\n", prop.textureAlignment);
    printf("max Texture 1D: %d\n", prop.maxTexture1D);
    printf("max Texture 2D: %d, %d\n", prop.maxTexture2D[0], prop.maxTexture2D[1]);
    printf("max Texture 3D: %d, %d, %d\n", prop.maxTexture3D[0], prop.maxTexture3D[1], prop.maxTexture3D[2]);
    printf(" --- MP Information for device %d ---\n", deviceId);
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

    cudaSetDevice(deviceId);
    size_t free;
    size_t total;
    cudaMemGetInfo(&free, &total);
    printf("Free: %zu Bytes, Total: %zu Bytes\n", free, total);
    printf("Free: %zu MB, Total: %zu MB\n", free / 1000 / 1000, total / 1000 / 1000);
}

} // namespace vf::cuda
