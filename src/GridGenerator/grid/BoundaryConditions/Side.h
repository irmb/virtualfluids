#ifndef SIDE_H
#define SIDE_H

#include <string>
#include <vector>

#include <VirtualFluidsDefinitions.h>

#include <core/PointerDefinitions.h>
#include <core/DataTypes.h>

#define X_INDEX 0
#define Y_INDEX 1
#define Z_INDEX 2

#define POSITIVE_DIR 1
#define NEGATIVE_DIR -1

class Grid;
class BoundaryCondition;

class Side;

enum class SideType
{
    MX, PX, MY, PY, MZ, PZ, GEOMETRY
};



class Side
{
public:
    virtual void addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition) = 0;

    virtual int getCoordinate() const = 0;
    virtual int getDirection() const = 0;

protected:
    static void addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::string coord, real constant,
                           real startInner, real endInner, real startOuter, real endOuter);

    static void setPressureNeighborIndices(SPtr<BoundaryCondition> boundaryCondition, SPtr<Grid> grid, const uint index);

private:
    static uint getIndex(SPtr<Grid> grid, std::string coord, real constant, real v1, real v2);
};

class Geometry : public Side
{
public:
    void addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return X_INDEX;
    }

    int getDirection() const override
    {
        return NEGATIVE_DIR;
    }
};

class MX : public Side
{
public:
    void addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return X_INDEX;
    }

    int getDirection() const override
    {
        return NEGATIVE_DIR;
    }
};

class PX : public Side
{
public:
    void addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return X_INDEX;
    }

    int getDirection() const override
    {
        return POSITIVE_DIR;
    }
};


class MY : public Side
{
public:
    void addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return Y_INDEX;
    }

    int getDirection() const override
    {
        return NEGATIVE_DIR;
    }
};

class PY : public Side
{
public:
    void addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return Y_INDEX;
    }

    int getDirection() const override
    {
        return POSITIVE_DIR;
    }
};


class MZ : public Side
{
public:
    void addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return Z_INDEX;
    }

    int getDirection() const override
    {
        return NEGATIVE_DIR;
    }
};

class PZ : public Side
{
public:
    void addIndices(std::vector<SPtr<Grid> > grid, uint level, SPtr<BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return Z_INDEX;
    }

    int getDirection() const override
    {
        return POSITIVE_DIR;
    }
};


class SideFactory
{
public:
    static SPtr<Side> make(SideType sideType)
    {
        switch (sideType)
        {
        case SideType::MX:
            return SPtr<Side>(new MX());
        case SideType::PX:
            return SPtr<Side>(new PX());
        case SideType::MY:
            return SPtr<Side>(new MY());
        case SideType::PY:
            return SPtr<Side>(new PY());
        case SideType::MZ:
            return SPtr<Side>(new MZ());
        case SideType::PZ:
            return SPtr<Side>(new PZ());
        case SideType::GEOMETRY:
            return SPtr<Side>(new Geometry());
        }
    }
};

#endif
