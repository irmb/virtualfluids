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
//! \file CudaRunKernel.hpp
//! \ingroup CudaUtility
//! \author Stephan Lenz
//=======================================================================================
#ifndef  CudaRunKernel_HPP
#define  CudaRunKernel_HPP

#include <string>
#include <device_launch_parameters.h>

#include "CudaUtility/CudaUtility.h"

template<typename KernelFunctor, typename FunctionFunctor, typename... TArgs>
void runKernel(KernelFunctor kernel, FunctionFunctor function, std::string deviceType, const CudaUtility::CudaGrid& grid, TArgs... args)
{
    if( grid.numberOfEntities == 0 ) return;

    if( deviceType == "GPU" )
    {
        kernel<<< grid.blocks, grid.threads, 0, grid.stream >>>( args..., grid.numberOfEntities );
    }
    else
    {
        for( uint index = 0; index < grid.numberOfEntities; index++ )
        {
            function( args..., index );
        }
    }
}

#endif
