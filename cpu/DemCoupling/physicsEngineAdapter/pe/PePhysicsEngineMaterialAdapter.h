/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef PE_PHYSICS_ENGINE_MATERIAL_ADAPTER
#define PE_PHYSICS_ENGINE_MATERIAL_ADAPTER

#include "../PhysicsEngineMaterialAdapter.h"
#include <pe/basic.h>


class PePhysicsEngineMaterialAdapter : public PhysicsEngineMaterialAdapter
{
public:
    PePhysicsEngineMaterialAdapter(std::string name, double density, double restitution, double staticFriction, double dynamicFriction, double poissonRatio, double youngModul, double stiffnessInNormalDirection, double dampingoefficientNormalDirection, double dampingTangentialDirection)
        : PhysicsEngineMaterialAdapter(name, density, restitution, staticFriction, dynamicFriction, poissonRatio, youngModul, stiffnessInNormalDirection, dampingoefficientNormalDirection, dampingTangentialDirection)
    {
    }
    virtual ~PePhysicsEngineMaterialAdapter() {}

    virtual walberla::pe::MaterialID getPeMaterial() const;

};

#endif

