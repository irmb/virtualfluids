#ifndef BoundaryCondition_H
#define BoundaryCondition_H

#include <vector>
#include <functional>

#include "global.h"

#include "grid/NodeValues.h"

class Grid;

class Side;
enum class SideType;

class BoundaryCondition
{
public:
    virtual ~BoundaryCondition() = default;

    std::vector<uint> indices;
    SPtr<Side> side;
    std::vector<std::vector<real> > qs;

    std::vector<uint> patches;

    virtual char getType() const = 0;

    bool isSide( SideType side ) const;

    real getQ( uint index, uint dir ){ return this->qs[index][dir]; }
};

//////////////////////////////////////////////////////////////////////////

class PressureBoundaryCondition : public BoundaryCondition
{
public:
    static SPtr<PressureBoundaryCondition> make(real rho)
    {
        return SPtr<PressureBoundaryCondition>(new PressureBoundaryCondition(rho));
    }

    // matrix indices!!!
    std::vector<uint> neighborIndices;

    real rho;
protected:
    PressureBoundaryCondition(real rho) : rho(rho) { }

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

//////////////////////////////////////////////////////////////////////////

class VelocityBoundaryCondition : public BoundaryCondition
{
public:
    static SPtr<VelocityBoundaryCondition> make(real vx, real vy, real vz)
    {
        return SPtr<VelocityBoundaryCondition>(new VelocityBoundaryCondition(vx, vy, vz));
    }

    real vx, vy, vz;
    std::vector<real> vxList, vyList, vzList;
protected:
    VelocityBoundaryCondition(real vx, real vy, real vz) : vx(vx), vy(vy), vz(vz) { }

public:
    virtual char getType() const override
    {
        return BC_VELOCITY;
    }

    void fillVelocityLists()
    {
        for( uint index : this->indices ) {
            (void) index;
            this->vxList.push_back(vx);
            this->vyList.push_back(vy);
            this->vzList.push_back(vz);
        }
    }

    real getVx() { return this->vx; }
    real getVy() { return this->vy; }
    real getVz() { return this->vz; }

    real getVx(uint index) { return this->vxList[index]; }
    real getVy(uint index) { return this->vyList[index]; }
    real getVz(uint index) { return this->vzList[index]; }
    
    GRIDGENERATOR_EXPORT void setVelocityProfile( SPtr<Grid> grid, std::function<void(real,real,real,real&,real&,real&)> velocityProfile );
};

//////////////////////////////////////////////////////////////////////////


class GeometryBoundaryCondition : public VelocityBoundaryCondition
{
public:
    static SPtr<GeometryBoundaryCondition> make()
    {
        return SPtr<GeometryBoundaryCondition>(new GeometryBoundaryCondition());
    }

private:
    GeometryBoundaryCondition() : VelocityBoundaryCondition(0.0, 0.0, 0.0) { }

public:
    char getType() const override
    {
        return BC_SOLID;
    }

    void setVelocityForPatch( uint patch, real vx, real vy, real vz ){
        for( uint index = 0; index < this->indices.size(); index++ ){
            if( this->patches[index] == patch ){
                this->vxList[index] = vx;
                this->vyList[index] = vy;
                this->vzList[index] = vz;
            }
        }
    }

    GRIDGENERATOR_EXPORT void setTangentialVelocityForPatch( SPtr<Grid> grid, uint patch,
                                                  real p1x, real p1y, real p1z, 
                                                  real p2x, real p2y, real p2z, 
                                                  real v, real r );
};


#endif