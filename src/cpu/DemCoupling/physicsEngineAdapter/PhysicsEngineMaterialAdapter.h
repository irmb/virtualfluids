/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef PHYSICS_ENGINE_MATERIAL_ADAPTER_H
#define PHYSICS_ENGINE_MATERIAL_ADAPTER_H

#include <string>

class PhysicsEngineMaterialAdapter
{
public:
    PhysicsEngineMaterialAdapter(std::string name, double density, double restitution, double staticFriction, double dynamicFriction, double poissonRatio, double youngModul, double stiffnessInNormalDirection, double dampingoefficientNormalDirection, double dampingTangentialDirection)
        : name(name), density(density), restitution(restitution), staticFriction(staticFriction), dynamicFriction(dynamicFriction), poissonRatio(poissonRatio), youngModul(youngModul), stiffnessInNormalDirection(stiffnessInNormalDirection), dampingoefficientNormalDirection(dampingoefficientNormalDirection), dampingTangentialDirection(dampingTangentialDirection)
    {}
    virtual ~PhysicsEngineMaterialAdapter() {}

protected:
    std::string name;
    double density;
    double restitution;
    double staticFriction; // Note: pe doubles the input coefficient of friction for material-material contacts.
    double dynamicFriction; //  Similar to static friction for low speed friction.
    double poissonRatio;
    double youngModul;
    double stiffnessInNormalDirection;
    double dampingoefficientNormalDirection;
    double dampingTangentialDirection;
};


#endif

