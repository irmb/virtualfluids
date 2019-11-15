#include "VectorTypes.h"

//Vec3 Vec3::operator+( Vec3& left, Vec3& right ){
Vec3 Vec3::operator+( Vec3& right ){
    return Vec3( this->x + right.x, 
                 this->y + right.y, 
                 this->z + right.z );
}

Vec3 Vec3::operator-( Vec3& right ){
    return Vec3( this->x - right.x, 
                 this->y - right.y, 
                 this->z - right.z );
}

Vec3 operator*( real scalar, Vec3& vec ){
    return Vec3( scalar * vec.x, 
                 scalar * vec.y, 
                 scalar * vec.z );
}
