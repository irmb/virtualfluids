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
    SPtr<Side> side;

    virtual void print() const = 0;

};

class PressureBoundaryCondition : public BoundaryCondition
{
public:
    static SPtr<PressureBoundaryCondition> make(real rho)
    {
        return SPtr<PressureBoundaryCondition>(new PressureBoundaryCondition(rho));
    }


    std::vector<uint> neighborIndices;

    real rho;
private:
    PressureBoundaryCondition(real rho) : rho(rho)
    {

    }

public:
    void print() const override {};
};

class VelocityBoundaryCondition : public BoundaryCondition
{
public:
    static SPtr<VelocityBoundaryCondition> make(real vx, real vy, real vz)
    {
        return SPtr<VelocityBoundaryCondition>(new VelocityBoundaryCondition(vx, vy, vz));
    }

    real vx, vy, vz;
private:
    VelocityBoundaryCondition(real vx, real vy, real vz) : vx(vx), vy(vy), vz(vz)
    {

    }

public:
    void print() const override {};
};


class GeometryBoundaryCondition : public BoundaryCondition
{
public:
    static SPtr<GeometryBoundaryCondition> make()
    {
        return SPtr<GeometryBoundaryCondition>(new GeometryBoundaryCondition());
    }

    std::vector<std::vector<real> > qs;
private:
    GeometryBoundaryCondition()
    {

    }

public:
    void print() const override {};
};


#endif