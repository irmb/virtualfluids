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

void CudaStreamManager::registerStream(CudaStreamIndex streamIndex)
{
    if(streamIndex != CudaStreamIndex::Legacy)
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

cudaStream_t &CudaStreamManager::getStream(CudaStreamIndex streamIndex, uint multiStreamIndex)
{
    if(streamIndex == CudaStreamIndex::Legacy)  return legacyStream;
    if(streamIsRegistered(streamIndex))
    {
        auto it = cudaStreams.find(streamIndex);
        for(uint idx=0; idx<multiStreamIndex; idx++) it++;
        return it->second;
    }
    return legacyStream;
}

bool CudaStreamManager::streamIsRegistered(CudaStreamIndex streamIndex)
{
    return cudaStreams.count(streamIndex) > 0;
}

void CudaStreamManager::createCudaEvents()
{
    checkCudaErrors(cudaEventCreateWithFlags(&startBulkKernel, cudaEventDisableTiming));
}

void CudaStreamManager::destroyCudaEvents() 
{ 
    checkCudaErrors(cudaEventDestroy(startBulkKernel)); 
}

void CudaStreamManager::triggerStartBulkKernel(CudaStreamIndex streamIndex, uint multiStreamIndex)
{
    checkCudaErrors(cudaEventRecord(startBulkKernel, getStream(streamIndex, multiStreamIndex)));
}

void CudaStreamManager::waitOnStartBulkKernelEvent(CudaStreamIndex streamIndex, uint multiStreamIndex)
{
    checkCudaErrors(cudaStreamWaitEvent(getStream(streamIndex, multiStreamIndex), startBulkKernel));
}
