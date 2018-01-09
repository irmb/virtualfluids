#include "PePhysicsEngineGeometryAdapter.h"

#include <pe/basic.h>

#include "PeAdapter.h"


PePhysicsEngineGeometryAdapter::PePhysicsEngineGeometryAdapter(walberla::pe::RigidBody* peGeoObject) : peGeoObject(peGeoObject)
{
    this->id = peGeoObject->getID();
}

void PePhysicsEngineGeometryAdapter::addForce(const Vector3D& force)
{
    peGeoObject->addForce(PeConverter::convert(force));
}

void PePhysicsEngineGeometryAdapter::addTorque(const Vector3D& torque)
{
    peGeoObject->addTorque(PeConverter::convert(torque));
}

void PePhysicsEngineGeometryAdapter::setForce(const Vector3D& force)
{
    peGeoObject->setForce(PeConverter::convert(force));
}

void PePhysicsEngineGeometryAdapter::setTorque(const Vector3D& torque)
{
    peGeoObject->setTorque(PeConverter::convert(torque));
}

void PePhysicsEngineGeometryAdapter::addForceAtPosition(const Vector3D& force, const Vector3D& position)
{
    peGeoObject->addForceAtPos(PeConverter::convert(force), PeConverter::convert(position));
}

void PePhysicsEngineGeometryAdapter::setLinearVelolocity(const Vector3D& velocity)
{
    peGeoObject->setLinearVel(PeConverter::convert(velocity));
}

void PePhysicsEngineGeometryAdapter::setAngularVelocity(const Vector3D& velocity)
{
    peGeoObject->setAngularVel(PeConverter::convert(velocity));
}

void PePhysicsEngineGeometryAdapter::resetForceAndTorque()
{
    peGeoObject->resetForceAndTorque();
}

Vector3D PePhysicsEngineGeometryAdapter::getVelocityAtPosition(const Vector3D& position) const
{
    return PeConverter::convert(peGeoObject->velFromWF(PeConverter::convert(position)));
}

Vector3D PePhysicsEngineGeometryAdapter::getLinearVelocity() const
{
    return PeConverter::convert(peGeoObject->getLinearVel());
}

Vector3D PePhysicsEngineGeometryAdapter::getAngularVelocity() const
{
    return PeConverter::convert(peGeoObject->getAngularVel());
}

Vector3D PePhysicsEngineGeometryAdapter::getPosition() const
{
    return PeConverter::convert(peGeoObject->getPosition());
}

Vector3D PePhysicsEngineGeometryAdapter::getForce() const
{
    return PeConverter::convert(peGeoObject->getForce());
}

Vector3D PePhysicsEngineGeometryAdapter::getTorque() const
{
    return PeConverter::convert(peGeoObject->getTorque());
}

void PePhysicsEngineGeometryAdapter::changeState(State state)
{
    bool infinite;
    if (state == State::UNPIN)
        infinite = false;
    else
        infinite = true;
   
    peGeoObject->setMass(infinite);
}

int PePhysicsEngineGeometryAdapter::getId() const
{
    return id;
}

void PePhysicsEngineGeometryAdapter::setGeometry(walberla::pe::RigidBody* peGeoObject)
{
    this->peGeoObject = peGeoObject;
}

