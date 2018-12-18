#include "CudaUtility.h"

#include <iostream>
#include <cuda_runtime.h>
#include <cuda.h>
#include <helper_cuda.h>

#include "Core/DataTypes.h"
#include "Core/Logger/Logger.h"

CudaUtility::CudaGrid::CudaGrid( uint numberOfEntities, uint threadsPerBlock )
{
    this->numberOfEntities = numberOfEntities;
    this->threads.x = threadsPerBlock;
    this->blocks.x  = ( numberOfEntities + threadsPerBlock - 1 ) / threadsPerBlock;
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

void CudaUtility::setCudaDevice(int device)
{    
    checkCudaErrors( cudaSetDevice( device ) );
    checkCudaErrors( cudaGetDevice( &device ) );

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device);

    *logging::out << logging::Logger::INFO_HIGH << "Set device " << device << ": " << prop.name << "\n";
}
