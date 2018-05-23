#ifndef BoundaryCondition_H
#define BoundaryCondition_H

#include <vector>

#include <core/PointerDefinitions.h>
#include <core/DataTypes.h>

class Side;

class BoundaryCondition
{
public:
    std::vector<uint> indices;

};

class PressureBoundaryCondition : public BoundaryCondition
{
public:
    static SPtr<PressureBoundaryCondition> make(real rho)
    {
        return SPtr<PressureBoundaryCondition>(new PressureBoundaryCondition(rho));
    }

    SPtr<Side> side;
    real rho;
private:
    PressureBoundaryCondition(real rho) : rho(rho)
    {

    }
};

class VelocityBoundaryCondition : public BoundaryCondition
{
public:
    static SPtr<VelocityBoundaryCondition> make(real vx, real vy, real vz)
    {
        return SPtr<VelocityBoundaryCondition>(new VelocityBoundaryCondition(vx, vy, vz));
    }

    SPtr<Side> side;
    real vx, vy, vz;
private:
    VelocityBoundaryCondition(real vx, real vy, real vz) : vx(vx), vy(vy), vz(vz)
    {

    }
};


#endif