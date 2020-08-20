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
//! \file CellProperties.cuh
//! \ingroup CellProperties
//! \author Stephan Lenz
//=======================================================================================
#ifndef CellProperties_H
#define CellProperties_H

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#define __host__
#define __device__
#endif

//////////////////////////////////////////////////////////////////////////

#define CELL_PROPERTIES_DEFAULT        (0u)
#define CELL_PROPERTIES_GHOST          (1u)
#define CELL_PROPERTIES_WALL           (2u)
#define CELL_PROPERTIES_FINE_GHOST     (4u)
#define CELL_PROPERTIES_IS_FLUX_BC     (8u)
#define CELL_PROPERTIES_IS_INSULATED   (16u)
#define CELL_PROPERTIES_5              (32u)
#define CELL_PROPERTIES_6              (64u)
#define CELL_PROPERTIES_7              (128u)

//////////////////////////////////////////////////////////////////////////

//! stores several cell properties as bitmap
typedef unsigned char CellProperties;

//////////////////////////////////////////////////////////////////////////

//! sets a cell property
//! \param left   \ref CellProperties object in which the property should be set
//! \param right  cell property that should be set, usually one of the above preprocessor defines
__host__ __device__ inline void setCellProperties( CellProperties& left, const CellProperties& right )
{
    left |= right;
}

//! unsets a cell property
//! \param left   \ref CellProperties object in which the property should be unset
//! \param right  cell property that should be unset, usually one of the above preprocessor defines
__host__ __device__ inline void unsetCellProperties( CellProperties& left, const CellProperties& right )
{
    left &= ~right;
}

//! queries a specific cell property
//! \param left   \ref CellProperties object in which the property should be queried
//! \param right  cell property that should be queried, usually one of the above preprocessor defines
//! \return true if the cell property is set and false otherwise
__host__ __device__ inline bool isCellProperties( const CellProperties& left, const CellProperties& right )
{
    return (left & right) == right;
}

//////////////////////////////////////////////////////////////////////////

#endif

