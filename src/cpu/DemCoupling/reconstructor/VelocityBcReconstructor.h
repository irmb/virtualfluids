/*
 *  Author: S. Peters
 *  mail: peters@irmb.tu-bs.de
 */
#ifndef VELOCITY_BC_RECONSTRUCTOR_H
#define VELOCITY_BC_RECONSTRUCTOR_H

#include "UbTuple.h"

#include "Reconstructor.h"

class ILBMKernel;
class PhysicsEngineGeometryAdapter;

class VelocityBcReconstructor : public Reconstructor
{
public:
    virtual ~VelocityBcReconstructor() {}

    void reconstructNode(const int &x1, const int &x2, const int &x3, const Vector3D &worldCoordinates,
                         std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry,
                         std::shared_ptr<ILBMKernel> kernel) const override;
};

#endif
