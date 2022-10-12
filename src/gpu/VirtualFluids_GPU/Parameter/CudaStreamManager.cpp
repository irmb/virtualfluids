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
#include <cuda_runtime.h>
#include <cuda.h>
#include <iostream>

void CudaStreamManager::registerStream(StreamIndex streamIndex)
{
    cudaStreams.emplace(streamIndex, nullptr);
}

void CudaStreamManager::launchStreams()
{
    for (auto &stream : cudaStreams)
        cudaStreamCreate(&stream.second);
}

void CudaStreamManager::terminateStreams()
{
    for (auto &stream : cudaStreams)
        cudaStreamDestroy(stream.second);
}

bool CudaStreamManager::streamIsRegistered(StreamIndex streamIndex)
{
    return cudaStreams.count(streamIndex) == 1;
}

cudaStream_t &CudaStreamManager::getStream(StreamIndex streamIndex)
{ 
    return streamIsRegistered(streamIndex) ? cudaStreams[streamIndex] : legacyStream;
}


void CudaStreamManager::createCudaEvents()
{
    checkCudaErrors(cudaEventCreateWithFlags(&startBulkKernel, cudaEventDisableTiming));
}

void CudaStreamManager::destroyCudaEvents() 
{ 
    checkCudaErrors(cudaEventDestroy(startBulkKernel)); 
}

void CudaStreamManager::triggerStartBulkKernel(StreamIndex streamIndex)
{
    checkCudaErrors(cudaEventRecord(startBulkKernel, cudaStreams[streamIndex]));
}

void CudaStreamManager::waitOnStartBulkKernelEvent(StreamIndex streamIndex)
{
    checkCudaErrors(cudaStreamWaitEvent(cudaStreams[streamIndex], startBulkKernel));
}
