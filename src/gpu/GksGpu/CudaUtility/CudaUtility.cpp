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
//! \file CudaUtility.cpp
//! \ingroup CudaUtility
//! \author Stephan Lenz
//=======================================================================================
#include "CudaUtility.h"

#include <iostream>
#include <cuda_runtime.h>
#include <cuda.h>
#include <helper_cuda.h>

#include "Core/DataTypes.h"
#include "Core/Logger/Logger.h"

cudaStream_t CudaUtility::computeStream = nullptr;
cudaStream_t CudaUtility::communicationStream = nullptr;

CudaUtility::CudaGrid::CudaGrid( uint numberOfEntities, uint threadsPerBlock, cudaStream_t stream )
{
    this->numberOfEntities = numberOfEntities;
    this->threads.x = threadsPerBlock;
    this->blocks.x  = ( numberOfEntities + threadsPerBlock - 1 ) / threadsPerBlock;

    this->stream = stream;
}

void CudaUtility::printCudaMemoryUsage()
{
    size_t free_byte ;
    size_t total_byte ;

    checkCudaErrors( cudaMemGetInfo( &free_byte, &total_byte ) );

    double free_db = (double)free_byte ;
    double total_db = (double)total_byte ;
    double used_db = total_db - free_db ;

    *logging::out << logging::Logger::INFO_HIGH << "GPU memory usage:" << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "    used  = " << used_db /1024.0/1024.0/1024.0 << " GB\n";
    *logging::out << logging::Logger::INFO_HIGH << "    free  = " << free_db /1024.0/1024.0/1024.0 << " GB\n";
    *logging::out << logging::Logger::INFO_HIGH << "    total = " << total_db/1024.0/1024.0/1024.0 << " GB\n";
}

int CudaUtility::getCudaDeviceCount()
{    
    int deviceCount = 0;
    checkCudaErrors( cudaGetDeviceCount(&deviceCount) );
    return deviceCount;
}

void CudaUtility::setCudaDevice(int device)
{    
    checkCudaErrors( cudaSetDevice( device ) );
    checkCudaErrors( cudaGetDevice( &device ) );

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device);

    *logging::out << logging::Logger::INFO_HIGH << "Set device " << device << ": " << prop.name << "\n";

    // set communication stream on high priority, such that it can interleave the compute stream
    // the non blocking flag disable implicit synchronization with the default thread '0'
    // based on https://fenix.tecnico.ulisboa.pt/downloadFile/563568428758047/CUDA_StreamsEvents.pdf
    // slide 5
    int priority_high, priority_low;
    cudaDeviceGetStreamPriorityRange(&priority_low , &priority_high ) ;

    // the flag needs to be cudaStreamDefault to ensure synchronization with default stream
    //cudaStreamCreateWithPriority (&communicationStream, cudaStreamDefault, priority_high );
    //cudaStreamCreateWithPriority (&computeStream      , cudaStreamDefault, priority_low  );
    cudaStreamCreateWithPriority (&communicationStream, cudaStreamNonBlocking, priority_high );
    cudaStreamCreateWithPriority (&computeStream      , cudaStreamNonBlocking, priority_low  );
}

int CudaUtility::getCudaDevice()
{
    int device;
    checkCudaErrors( cudaGetDevice( &device ) );

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device);

    *logging::out << logging::Logger::INFO_HIGH << "The current device " << device << ": " << prop.name << "\n";

    return device;
}

void CudaUtility::synchronizeCudaDevice()
{
    checkCudaErrors( cudaDeviceSynchronize() );
}

void CudaUtility::synchronizeCudaStream(cudaStream_t stream)
{
    checkCudaErrors( cudaStreamSynchronize(stream) );
}
