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
//=======================================================================================
#include "CudaStreamManager.h"
#include <helper_cuda.h>
#include <iostream>

CudaStreamManager::CudaStreamManager() {}

CudaStreamManager::~CudaStreamManager() {}

void CudaStreamManager::launchStreams(uint numberOfStreams)
{
    cudaStreams.resize(numberOfStreams);
    for (cudaStream_t &stream : cudaStreams)
        cudaStreamCreate(&stream);
}

void CudaStreamManager::terminateStreams()
{
    for (cudaStream_t &stream : cudaStreams)
        cudaStreamDestroy(stream);
}

cudaStream_t &CudaStreamManager::getStream(uint streamIndex)
{
    return cudaStreams[streamIndex]; }

void CudaStreamManager::createCudaEvents()
{
    checkCudaErrors(cudaEventCreateWithFlags(&startBulkKernel, cudaEventDisableTiming));
}

void CudaStreamManager::destroyCudaEvents() 
{ 
    checkCudaErrors(cudaEventDestroy(startBulkKernel)); 
}

void CudaStreamManager::triggerStartBulkKernel(int streamIndex)
{
    checkCudaErrors(cudaEventRecord(startBulkKernel, cudaStreams[streamIndex]));
}

void CudaStreamManager::triggerEventByName(std::string eventName, int streamIndex)
{
    if (eventName == "startBulkKernel")
        checkCudaErrors(cudaEventRecord(startBulkKernel, cudaStreams[streamIndex]));
    else
        std::cout << "unknown event name" << std::endl;
}

void CudaStreamManager::waitOnStartBulkKernelEvent(int streamIndex)
{
    checkCudaErrors(cudaStreamWaitEvent(cudaStreams[streamIndex], startBulkKernel));
}
