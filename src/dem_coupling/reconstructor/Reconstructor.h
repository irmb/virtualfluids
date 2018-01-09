/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef RECONSTRCUTOR_H
#define RECONSTRCUTOR_H

#include "boost/shared_ptr.hpp"

#include "VirtualFluidsBasics/basics/utilities/Vector3D.h"

class ILBMKernel;
class PhysicsEngineGeometryAdapter;

class Reconstructor
{
public:
    virtual ~Reconstructor() {}

    virtual void reconstructNode(const int &x1, const int &x2, const int &x3, const Vector3D& worldCoordinates, std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry, std::shared_ptr<ILBMKernel> kernel) const = 0;

};



#endif

