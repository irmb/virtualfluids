#include "DummyPhysicsEngineGeometryAdapter.h"



void DummyPhysicsEngineGeometryAdapter::addForce(const Vector3D& force)
{

}

void DummyPhysicsEngineGeometryAdapter::addTorque(const Vector3D& torque)
{

}

void DummyPhysicsEngineGeometryAdapter::setForce(const Vector3D& force)
{

}

void DummyPhysicsEngineGeometryAdapter::setTorque(const Vector3D& torque)
{

}

void DummyPhysicsEngineGeometryAdapter::addForceAtPosition(const Vector3D& force, const Vector3D& position)
{

}

void DummyPhysicsEngineGeometryAdapter::setLinearVelolocity(const Vector3D& velocity)
{
    this->velocity = velocity;
}

void DummyPhysicsEngineGeometryAdapter::setAngularVelocity(const Vector3D& velocity)
{

}

void DummyPhysicsEngineGeometryAdapter::resetForceAndTorque()
{

}

Vector3D DummyPhysicsEngineGeometryAdapter::getVelocityAtPosition(const Vector3D& position) const
{
    return velocity;
}

Vector3D DummyPhysicsEngineGeometryAdapter::getLinearVelocity() const
{
    return Vector3D();
}

Vector3D DummyPhysicsEngineGeometryAdapter::getAngularVelocity() const
{
    return Vector3D();
}

Vector3D DummyPhysicsEngineGeometryAdapter::getPosition() const
{
    return Vector3D();
}

Vector3D DummyPhysicsEngineGeometryAdapter::getForce() const
{
    return Vector3D();
}

Vector3D DummyPhysicsEngineGeometryAdapter::getTorque() const
{
    return Vector3D();
}

void DummyPhysicsEngineGeometryAdapter::changeState(State state)
{

}

