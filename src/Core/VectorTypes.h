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
//! \file VectorTypes.h
//! \ingroup Core
//! \author Stephan Lenz
//=======================================================================================
#ifndef VECTORTYPES_H
#define VECTORTYPES_H

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#define __host__
#define __device__
#endif

#include <cmath>

#include "VirtualFluidsDefinitions.h"

#include "DataTypes.h"
#include "RealConstants.h"

struct VF_PUBLIC Vec3 {
    real x, y, z; 

    __host__ __device__ Vec3(real x, real y, real z) : x(x), y(y), z(z) {}
    __host__ __device__ Vec3() : x(c0o1), y(c0o1), z(c0o1) {}

    __host__ __device__ real length() {
        return std::sqrt( x*x + y*y + z*z );
    }
};

VF_PUBLIC Vec3 operator+( Vec3& left, Vec3& right );
VF_PUBLIC Vec3 operator-( Vec3& left, Vec3& right );
VF_PUBLIC Vec3 operator*( real scalar, Vec3& vec );

#endif
