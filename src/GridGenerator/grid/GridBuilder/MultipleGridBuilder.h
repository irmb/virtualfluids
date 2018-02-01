#ifndef MULTIPLE_GRID_BUILDER_H
#define MULTIPLE_GRID_BUILDER_H

#include "GridGenerator/global.h"

#include <vector>
#include <string>
#include <memory>
#include <exception>
#include "grid/GridMocks.h"

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

template<typename Grid>
class MultipleGridBuilder
{
private:
    VF_PUBLIC MultipleGridBuilder();

public:
    VF_PUBLIC static SPtr<MultipleGridBuilder> makeShared();

    VF_PUBLIC void addCoarseGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta);
    VF_PUBLIC void addGrid(real startX, real startY, real startZ, real endX, real endY, real endZ);
    VF_PUBLIC uint getNumberOfLevels() const;
    VF_PUBLIC real getDelta(int level) const;

private:
    bool isInsideOfGrids(SPtr<Grid> grid) const;
    real calculateDelta();
    void checkIfCoarseGridIsMissing() const;
    void checkIfGridIsInCoarseGrid(SPtr<Grid> grid) const;

    std::vector<SPtr<Grid> > grids;
};

#endif

