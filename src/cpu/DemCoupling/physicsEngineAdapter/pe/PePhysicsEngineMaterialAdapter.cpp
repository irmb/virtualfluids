#include "PePhysicsEngineMaterialAdapter.h"


walberla::pe::MaterialID PePhysicsEngineMaterialAdapter::getPeMaterial() const
{
    if (walberla::pe::Material::find(name) != -1)
        return walberla::pe::Material::find(name);

    return walberla::pe::createMaterial(name, density, restitution, staticFriction, dynamicFriction, poissonRatio, youngModul, stiffnessInNormalDirection, dampingoefficientNormalDirection, dampingTangentialDirection);
}
