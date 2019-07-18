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
