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

    Vec3 operator+( Vec3& right );
    Vec3 operator-( Vec3& right );
};

//VF_PUBLIC Vec3 operator+( Vec3& left, Vec3& right );
//VF_PUBLIC Vec3 operator-( Vec3& left, Vec3& right );
VF_PUBLIC Vec3 operator*( real scalar, Vec3& vec );

#endif