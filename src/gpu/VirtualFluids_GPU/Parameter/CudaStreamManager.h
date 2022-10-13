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
#ifndef STREAM_MANAGER_H
#define STREAM_MANAGER_H

#include <map>
#include <cuda.h>
#include <cuda_runtime.h>
enum class CudaStreamIndex
    {
        Legacy,
        Bulk,
        Border,
        Precursor
    };
class CudaStreamManager
{
public:
    
private:
    std::map<CudaStreamIndex, cudaStream_t> cudaStreams;
    cudaEvent_t startBulkKernel = NULL;
    cudaStream_t legacyStream = CU_STREAM_LEGACY;


public:
    void registerStream(CudaStreamIndex streamIndex);
    void launchStreams();
    void terminateStreams();
    cudaStream_t &getStream(CudaStreamIndex streamIndex);

    bool streamIsRegistered(CudaStreamIndex streamIndex);
    // Events
    void createCudaEvents();
    void destroyCudaEvents();
    void triggerStartBulkKernel(CudaStreamIndex streamIndex);
    void waitOnStartBulkKernelEvent(CudaStreamIndex streamIndex);
};

#endif
