#include "VectorTypes.h"

Vec3 operator+( Vec3& left, Vec3& right ){
    return Vec3( left.x + right.x, 
                 left.y + right.y, 
                 left.z + right.z );
}

Vec3 operator-( Vec3& left, Vec3& right ){
    return Vec3( left.x - right.x, 
                 left.y - right.y, 
                 left.z - right.z );
}

Vec3 operator*( real scalar, Vec3& vec ){
    return Vec3( scalar * vec.x, 
                 scalar * vec.y, 
                 scalar * vec.z );
}
