/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef DUMMY_PHYSICS_ENGINE_MATERIAL_ADAPTER
#define DUMMY_PHYSICS_ENGINE_MATERIAL_ADAPTER

#include "PhysicsEngineMaterialAdapter.h"


class DummyPhysicsEngineMaterialAdapter : public PhysicsEngineMaterialAdapter
{
public:
    DummyPhysicsEngineMaterialAdapter(std::string name, double density, double restitution, double staticFriction, double dynamicFriction, double poissonRatio, double youngModul, double stiffnessInNormalDirection, double dampingoefficientNormalDirection, double dampingTangentialDirection)
        : PhysicsEngineMaterialAdapter(name, density, restitution, staticFriction, dynamicFriction, poissonRatio, youngModul, stiffnessInNormalDirection, dampingoefficientNormalDirection, dampingTangentialDirection)
    {
    }
    virtual ~DummyPhysicsEngineMaterialAdapter() {}

};

#endif

