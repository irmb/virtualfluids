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
//! \file Transformation.cuh
//! \ingroup FluxComputation
//! \author Stephan Lenz
//=======================================================================================
#ifndef Transformation_CUH
#define Transformation_CUH

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! \brief Transforms vector from global frame of reference to local frame of reference
//!
//! The transformation is based on the Cartesian directions.
//! 
//! \param[in,out] vector       vector that is transformed
//! \param[in]     direction    char with 'x', 'y' or 'z'
__host__ __device__ inline void transformGlobalToLocal(Vec3& vector, const char direction)
{
    if( direction == 'x' ) return;

    if( direction == 'y' )
    {
        Vec3 tmp = vector;
    
        vector.x = tmp.y;
        vector.y = tmp.z;
        vector.z = tmp.x;

        return;
    }

    if( direction == 'z' )
    {
        Vec3 tmp = vector;
    
        vector.x = tmp.z;
        vector.y = tmp.x;
        vector.z = tmp.y;

        return;
    }
}

//! \brief Transforms vector from local frame of reference to global frame of reference
//!
//! The transformation is based on the Cartesian directions.
//! 
//! \param[in,out] vector       vector that is transformed
//! \param[in]     direction    char with 'x', 'y' or 'z'
__host__ __device__ inline void transformLocalToGlobal(Vec3& vector, const char direction)
{
    if( direction == 'x' ) return;

    if( direction == 'y' )
    {
        Vec3 tmp;
    
        tmp.y = vector.x;
        tmp.z = vector.y;
        tmp.x = vector.z;

        vector = tmp;

        return;
    }

    if( direction == 'z' )
    {
        Vec3 tmp;
    
        tmp.z = vector.x;
        tmp.x = vector.y;
        tmp.y = vector.z;

        vector = tmp;

        return;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! \brief Transforms momentum components of \ref ConservedVariables from global frame of reference to local frame of reference
//!
//! The transformation is based on the Cartesian directions.
//! 
//! \param[in,out] cons         \ref ConservedVariables object that is transformed
//! \param[in]     direction    char with 'x', 'y' or 'z'
__host__ __device__ inline void transformGlobalToLocal(ConservedVariables& cons, const char direction)
{
    Vec3 vector( cons.rhoU, cons.rhoV, cons.rhoW );

    transformGlobalToLocal( vector, direction );

    cons.rhoU = vector.x;
    cons.rhoV = vector.y;
    cons.rhoW = vector.z;
}

//! \brief Transforms velocity components of \ref PrimitiveVariables from global frame of reference to local frame of reference
//!
//! The transformation is based on the Cartesian directions.
//! 
//! \param[in,out] prim         \ref PrimitiveVariables object that is transformed
//! \param[in]     direction    char with 'x', 'y' or 'z'
__host__ __device__ inline void transformGlobalToLocal(PrimitiveVariables& prim, const char direction)
{
    Vec3 vector( prim.U, prim.V, prim.W );

    transformGlobalToLocal( vector, direction );

    prim.U = vector.x;
    prim.V = vector.y;
    prim.W = vector.z;
}

//////////////////////////////////////////////////////////////////////////

//! \brief Transforms momentum components of \ref ConservedVariables from local frame of reference to global frame of reference
//!
//! The transformation is based on the Cartesian directions.
//! 
//! \param[in,out]  cons         \ref ConservedVariables object that is transformed
//! \param[in]      direction    char with 'x', 'y' or 'z'
__host__ __device__ inline void transformLocalToGlobal(ConservedVariables& cons, const char direction)
{
    Vec3 vector( cons.rhoU, cons.rhoV, cons.rhoW );

    transformLocalToGlobal( vector, direction );

    cons.rhoU = vector.x;
    cons.rhoV = vector.y;
    cons.rhoW = vector.z;
}

//! \brief Transforms velocity components of \ref PrimitiveVariables from local frame of reference to global frame of reference
//!
//! The transformation is based on the Cartesian directions.
//! 
//! \param[in,out] prim         \ref PrimitiveVariables object that is transformed
//! \param[in]     direction    char with 'x', 'y' or 'z'
__host__ __device__ inline void transformLocalToGlobal(PrimitiveVariables& prim, const char direction)
{
    Vec3 vector( prim.U, prim.V, prim.W );

    transformLocalToGlobal( vector, direction );

    prim.U = vector.x;
    prim.V = vector.y;
    prim.W = vector.z;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
