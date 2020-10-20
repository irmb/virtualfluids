/*
 *  Author: S. Peters
 *  mail: peters@irmb.tu-bs.de
 */
#ifndef PE_ADAPTER_H
#define PE_ADAPTER_H

#include "Vector3D.h"
#include <pe/basic.h>

class PeConverter
{
public:
    static Vector3D convert(walberla::pe::Vec3 vec3) { return Vector3D(vec3[0], vec3[1], vec3[2]); }

    static walberla::pe::Vec3 convert(const Vector3D &vec3) { return walberla::pe::Vec3(vec3[0], vec3[1], vec3[2]); }
};

#endif
