#ifndef VECTORTYPES_H
#define VECTORTYPES_H

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif
#endif

#include <cmath>

#include "basics_export.h"

#include "DataTypes.h"

struct BASICS_EXPORT Vec3 {
    real x{ 0. }, y{ 0. }, z{ 0. };

    __host__ __device__ Vec3(real x, real y, real z) : x(x), y(y), z(z) {}
    Vec3() = default;

    __host__ __device__ real length() { return std::sqrt(x * x + y * y + z * z); }

    Vec3 operator+(Vec3 &right);
    Vec3 operator-(Vec3 &right);
};

// BASICS_EXPORT Vec3 operator+( Vec3& left, Vec3& right );
// BASICS_EXPORT Vec3 operator-( Vec3& left, Vec3& right );
BASICS_EXPORT Vec3 operator*(real scalar, Vec3 &vec);

#endif
