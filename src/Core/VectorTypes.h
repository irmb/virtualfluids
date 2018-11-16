#ifndef VECTORTYPES_H
#define VECTORTYPES_H

#include <cmath>

#include "VirtualFluidsDefinitions.h"

#include "DataTypes.h"
#include "RealConstants.h"

struct VF_PUBLIC Vec3 {
    real x, y, z; 

    Vec3(real x, real y, real z) : x(x), y(y), z(z) {}
    Vec3() : x(zero), y(zero), z(zero) {}

    real length() {
        return std::sqrt( x*x + y*y + z*z );
    }
};

VF_PUBLIC Vec3 operator+( Vec3& left, Vec3& right );
VF_PUBLIC Vec3 operator-( Vec3& left, Vec3& right );
VF_PUBLIC Vec3 operator*( real scalar, Vec3& vec );

#endif
