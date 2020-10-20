/*
 *  Author: S. Peters
 *  mail: peters@irmb.tu-bs.de
 */
#ifndef PHYSICS_ENGINE_GEOMETRY_ADAPTER_H
#define PHYSICS_ENGINE_GEOMETRY_ADAPTER_H

#include "Vector3D.h"

enum class State { PIN, UNPIN };

class PhysicsEngineGeometryAdapter
{
public:
    virtual ~PhysicsEngineGeometryAdapter() {}

    virtual void addForce(const Vector3D &force)   = 0;
    virtual void addTorque(const Vector3D &torque) = 0;

    virtual void setForce(const Vector3D &force)   = 0;
    virtual void setTorque(const Vector3D &torque) = 0;

    virtual void addForceAtPosition(const Vector3D &force, const Vector3D &position) = 0;
    virtual void setLinearVelolocity(const Vector3D &velocity)                       = 0;
    virtual void setAngularVelocity(const Vector3D &velocity)                        = 0;

    virtual void resetForceAndTorque() = 0;

    virtual Vector3D getPosition() const                                   = 0;
    virtual Vector3D getVelocityAtPosition(const Vector3D &position) const = 0;
    virtual Vector3D getLinearVelocity() const                             = 0;
    virtual Vector3D getAngularVelocity() const                            = 0;

    virtual Vector3D getForce() const  = 0;
    virtual Vector3D getTorque() const = 0;

    virtual void changeState(State state) = 0;
};

#endif
