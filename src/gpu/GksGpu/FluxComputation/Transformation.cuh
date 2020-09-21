#ifndef Transformation_CUH
#define Transformation_CUH


#include "GksGpu_export.h"

#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

namespace GksGpu {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

__host__ __device__ inline void transformGlobalToLocal(ConservedVariables& cons, const char direction)
{
    Vec3 vector( cons.rhoU, cons.rhoV, cons.rhoW );

    transformGlobalToLocal( vector, direction );

    cons.rhoU = vector.x;
    cons.rhoV = vector.y;
    cons.rhoW = vector.z;
}

__host__ __device__ inline void transformGlobalToLocal(PrimitiveVariables& prim, const char direction)
{
    Vec3 vector( prim.U, prim.V, prim.W );

    transformGlobalToLocal( vector, direction );

    prim.U = vector.x;
    prim.V = vector.y;
    prim.W = vector.z;
}

//////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void transformLocalToGlobal(ConservedVariables& cons, const char direction)
{
    Vec3 vector( cons.rhoU, cons.rhoV, cons.rhoW );

    transformLocalToGlobal( vector, direction );

    cons.rhoU = vector.x;
    cons.rhoV = vector.y;
    cons.rhoW = vector.z;
}

__host__ __device__ inline void transformLocalToGlobal(PrimitiveVariables& prim, const char direction)
{
    Vec3 vector( prim.U, prim.V, prim.W );

    transformLocalToGlobal( vector, direction );

    prim.U = vector.x;
    prim.V = vector.y;
    prim.W = vector.z;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace GksGpu

#endif
