#ifndef MULTIPLE_GRID_BUILDER_H
#define MULTIPLE_GRID_BUILDER_H

#include "GridGenerator/global.h"

#include <vector>
#include <memory>
#include <array>
#include <exception>

template<typename Grid>
class GridFactory;

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
    VF_PUBLIC MultipleGridBuilder(SPtr<GridFactory<Grid> > gridFactory);

public:
    VF_PUBLIC static SPtr<MultipleGridBuilder> makeShared(SPtr<GridFactory<Grid> > gridFactory);

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

    VF_PUBLIC std::vector<SPtr<Grid> > getGrids() const;

private:
    void addGridToList(SPtr<Grid> grid);
    real calculateDelta(uint level) const;
    bool coarseGridExists() const;
    bool isGridInCoarseGrid(SPtr<Grid> grid) const;
    SPtr<Grid> makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, uint level) const;
    std::array<real, 6> getStaggeredCoordinates(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) const;
    SPtr<Grid> makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) const;

    static void emitNoCoarseGridExistsWarning();
    static void emitGridIsNotInCoarseGridWarning();

    std::vector<SPtr<Grid> > grids;
    SPtr<GridFactory<Grid> > gridFactory;
 
};

#endif

