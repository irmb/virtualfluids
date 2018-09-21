#ifndef BoundaryCondition_H
#define BoundaryCondition_H

#include <vector>

#include <core/PointerDefinitions.h>
#include <core/DataTypes.h>
#include "grid/NodeValues.h"

class Side;
enum class SideType;

class BoundaryCondition
{
public:
    std::vector<uint> indices;
    SPtr<Side> side;
    std::vector<std::vector<real> > qs;

    virtual char getType() const = 0;

    bool isSide( SideType side ) const;

    real getQ( uint index, uint dir ){ return this->qs[index][dir]; }
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
    char getType() const override
    {
        return BC_PRESSURE;
    }

    real getRho()
    {
        return this->rho;
    }
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
    char getType() const override
    {
        return BC_VELOCITY;
    }

    real getVx() { return this->vx; }
    real getVy() { return this->vy; }
    real getVz() { return this->vz; }
};


class GeometryBoundaryCondition : public BoundaryCondition
{
public:
    static SPtr<GeometryBoundaryCondition> make()
    {
        return SPtr<GeometryBoundaryCondition>(new GeometryBoundaryCondition());
    }

    real vx, vy, vz;
private:
    GeometryBoundaryCondition()
    {

    }

public:
    char getType() const override
    {
        return BC_SOLID;
    }

    real getVx() { return this->vx; }
    real getVy() { return this->vy; }
    real getVz() { return this->vz; }
};


#endif