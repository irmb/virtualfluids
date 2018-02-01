#ifndef MULTIPLE_GRID_BUILDER_H
#define MULTIPLE_GRID_BUILDER_H

#include "GridGenerator/global.h"

#include <vector>
#include <memory>
#include <array>
#include <exception>

class MultipleGridBuilderException : public std::exception 
{
public:
    const char* what() const noexcept override = 0;
};

class FinerGridBiggerThanCoarsestGridException : public MultipleGridBuilderException 
{
public:
    const char* what() const noexcept override
    {
        std::ostringstream getNr;
        getNr << "create two grids: second grid added is bigger than first grid but should be inside of first grid.";
        return getNr.str().c_str();
    }
};

class FirstGridMustBeCoarseException : public MultipleGridBuilderException
{
public:
    const char* what() const noexcept override
    {
        std::ostringstream getNr;
        getNr << "added a new grid without calling the method addCoarseGrid() before.";
        return getNr.str().c_str();
    }
};

class InvalidLevelException : public MultipleGridBuilderException
{
public:
    const char* what() const noexcept override
    {
        std::ostringstream getNr;
        getNr << "level is invalid.";
        return getNr.str().c_str();
    }
};

class NoIntermediateGridPossibleException : public MultipleGridBuilderException
{
public:
    const char* what() const noexcept override
    {
        std::ostringstream getNr;
        getNr << "fine grid is added. Not enough space between coarse and fine grid to create intermediate grids.";
        return getNr.str().c_str();
    }
};

template<typename Grid>
class MultipleGridBuilder
{
private:
    VF_PUBLIC MultipleGridBuilder();

public:
    VF_PUBLIC static SPtr<MultipleGridBuilder> makeShared();

    VF_PUBLIC void addCoarseGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta);
    VF_PUBLIC void addGrid(real startX, real startY, real startZ, real endX, real endY, real endZ);
    VF_PUBLIC void addFineGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, uint level);

    VF_PUBLIC uint getNumberOfLevels() const;
    VF_PUBLIC real getDelta(uint level) const;

    VF_PUBLIC real getStartX(uint level) const;
    VF_PUBLIC real getStartY(uint level) const;
    VF_PUBLIC real getStartZ(uint level) const;

    VF_PUBLIC real getEndX(uint level) const;
    VF_PUBLIC real getEndY(uint level) const;
    VF_PUBLIC real getEndZ(uint level) const;

private:
    void addGridToList(SPtr<Grid> grid);
    bool isInsideOfGrids(SPtr<Grid> grid) const;
    real calculateDelta() const;
    void checkIfCoarseGridIsMissing() const;
    void checkIfGridIsInCoarseGrid(SPtr<Grid> grid) const;
    SPtr<Grid> makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ) const;
    static std::array<real, 6> getStaggeredCoordinates(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta);
    static real getStaggeredCoordinate(real value, real offset);

    std::vector<SPtr<Grid> > grids;
};

#endif

