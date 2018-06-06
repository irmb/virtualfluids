/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef PE_PHYSICS_ENGINE_GEOMETRY_ADAPTER_H
#define PE_PHYSICS_ENGINE_GEOMETRY_ADAPTER_H


#include "PhysicsEngineGeometryAdapter.h"

namespace walberla
{
    namespace pe
    {
        class RigidBody;
    }
}

class PePhysicsEngineGeometryAdapter : public PhysicsEngineGeometryAdapter
{
public:
    PePhysicsEngineGeometryAdapter();
    PePhysicsEngineGeometryAdapter(walberla::pe::RigidBody* peGeoObject);
    virtual ~PePhysicsEngineGeometryAdapter() {}

    void addForce(const Vector3D& force) override;
    void addTorque(const Vector3D& torque) override;

    void setForce(const Vector3D& force) override;
    void setTorque(const Vector3D& torque) override;

    void addForceAtPosition(const Vector3D& force, const Vector3D& position) override;
    void setLinearVelolocity(const Vector3D& velocity) override;
    void setAngularVelocity(const Vector3D& velocity) override;

    void resetForceAndTorque() override;

    Vector3D getVelocityAtPosition(const Vector3D& position) const override;
    Vector3D getLinearVelocity() const override;
    Vector3D getAngularVelocity() const override;
    Vector3D getPosition() const override;
    Vector3D getForce() const override;
    Vector3D getTorque() const override;

    void changeState(State state) override;

    int getId() const;
    void setId(int id);
    void setGeometry(walberla::pe::RigidBody* peGeoObject);

    void setActive();
    void setInactive();
    bool isActive();

private:
    walberla::pe::RigidBody* peGeoObject;
    //unsigned long long id;
    int id;
    bool active;
};

#endif

