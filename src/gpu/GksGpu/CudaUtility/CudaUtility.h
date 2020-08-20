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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \file CudaUtility.h
//! \ingroup CudaUtility
//! \author Stephan Lenz
//=======================================================================================
#ifndef  CudaUtilExtern_H
#define  CudaUtilExtern_H

#include <cuda.h>
#include <cuda_runtime.h>

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"

#include "GksGpu_export.h"

//! \brief comprises some utility functions and classes for CUDA
class GKSGPU_EXPORT CudaUtility
{
public:

    //! \brief collects CUDA launch parameters for the runKernel Function
    struct CudaGrid 
    {
        dim3 threads;           //!< number of threads per block
        dim3 blocks;            //!< number of blocks

        uint numberOfEntities;  //!< total number of entities, smaller or equal threads * blocks

        cudaStream_t stream;    //!< CUDA stream

        CudaGrid( uint numberOfEntities, uint threadsPerBlock, cudaStream_t stream = 0 );
    };

    static cudaStream_t computeStream;          //!< Stream for compute tasks
    static cudaStream_t communicationStream;    //!< Stream for communication tasks, if communication hiding is enabled

    //! prints the used device memory to the logger
    static void printCudaMemoryUsage();

    // \return number of devices on the present machine
    static int getCudaDeviceCount();

    // \param device index of device that should be used
    static void setCudaDevice( int device );

    // \return current device index
    static int getCudaDevice(  );

    // encapsulated cudaDeviceSynchronize
    static void synchronizeCudaDevice();

    // encapsulates cudaStreamSynchronize
    // \param stream stream to be synchronized
    static void synchronizeCudaStream( cudaStream_t stream );
};

#endif
