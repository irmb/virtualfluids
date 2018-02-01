#include "MultipleGridBuilder.h"
#include "utilities/math/CudaMath.cuh"
#include "../GridMocks.h"

template<typename Grid>
MultipleGridBuilder<Grid>::MultipleGridBuilder()
{

}

template<typename Grid>
SPtr<MultipleGridBuilder<Grid> > MultipleGridBuilder<Grid>::makeShared()
{
    return SPtr<MultipleGridBuilder>(new MultipleGridBuilder<Grid>());
}

template <typename Grid>
void MultipleGridBuilder<Grid>::addCoarseGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta)
{
    auto grid = Grid::makeShared(startX, startY, startZ, endX, endY, endZ, delta);
    addGridToList(grid);
}

template <typename Grid>
void MultipleGridBuilder<Grid>::addGrid(real startX, real startY, real startZ, real endX, real endY, real endZ)
{
    checkIfCoarseGridIsMissing();
    auto grid = makeGrid(startX, startY, startZ, endX, endY, endZ);
    checkIfGridIsInCoarseGrid(grid);
    addGridToList(grid);
}

template <typename Grid>
void MultipleGridBuilder<Grid>::checkIfCoarseGridIsMissing() const
{
    if (grids.empty())
        throw FirstGridMustBeCoarseException();
}

template <typename Grid>
SPtr<Grid> MultipleGridBuilder<Grid>::makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ) const
{
    const real delta = calculateDelta();

    auto staggeredCoordinates = getStaggeredCoordinates(startX, startY, startZ, endX, endY, endZ, delta);

    return Grid::makeShared(staggeredCoordinates[0], staggeredCoordinates[1], staggeredCoordinates[2], staggeredCoordinates[3], staggeredCoordinates[4], staggeredCoordinates[5], delta);
}

template <typename Grid>
real MultipleGridBuilder<Grid>::calculateDelta() const
{
    const real delta = grids[grids.size() - 1]->getDelta() / 2.0;
    return delta;
}

template <typename Grid>
std::array<real, 6> MultipleGridBuilder<Grid>::getStaggeredCoordinates(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta)
{
    const real offset = delta * 0.5;

    const real startXStaggered = getStaggeredCoordinate(startX, offset);
    const real startYStaggered = getStaggeredCoordinate(startY, offset);
    const real startZStaggered = getStaggeredCoordinate(startZ, offset);
    
    const real endXStaggered = getStaggeredCoordinate(endX, -offset);
    const real endYStaggered = getStaggeredCoordinate(endY, -offset);
    const real endZStaggered = getStaggeredCoordinate(endZ, -offset);

    return std::array<real, 6>{startXStaggered, startYStaggered, startZStaggered, endXStaggered, endYStaggered, endZStaggered};
}

template <typename Grid>
real MultipleGridBuilder<Grid>::getStaggeredCoordinate(real value, real offset)
{
    return value + offset;
}

template <typename Grid>
void MultipleGridBuilder<Grid>::checkIfGridIsInCoarseGrid(SPtr<Grid> grid) const
{
    if (!isInsideOfGrids(grid))
        throw FinerGridBiggerThanCoarsestGridException();
}

template <typename Grid>
bool MultipleGridBuilder<Grid>::isInsideOfGrids(SPtr<Grid> grid) const
{
    return
        CudaMath::greaterEqual(grid->startX, grids[0]->startX) &&
        CudaMath::greaterEqual(grid->startY, grids[0]->startY) &&
        CudaMath::greaterEqual(grid->startZ, grids[0]->startZ) &&
        CudaMath::lessEqual(grid->endX, grids[0]->endX) &&
        CudaMath::lessEqual(grid->endY, grids[0]->endY) &&
        CudaMath::lessEqual(grid->endZ, grids[0]->endZ);
}

template <typename Grid>
uint MultipleGridBuilder<Grid>::getNumberOfLevels() const
{
    return uint(grids.size());
}

template <typename Grid>
real MultipleGridBuilder<Grid>::getDelta(int level) const
{
    if (grids.empty())
        throw InvalidLevelException();
    return grids[level]->getDelta();
}

template <typename Grid>
real MultipleGridBuilder<Grid>::getStartX(uint level) const
{
    return grids[level]->startX;
}

template <typename Grid>
real MultipleGridBuilder<Grid>::getStartY(uint level) const
{
    return grids[level]->startY;
}

template <typename Grid>
real MultipleGridBuilder<Grid>::getStartZ(uint level) const
{
    return grids[level]->startZ;
}

template <typename Grid>
real MultipleGridBuilder<Grid>::getEndX(uint level) const
{
    return grids[level]->endX;
}

template <typename Grid>
real MultipleGridBuilder<Grid>::getEndY(uint level) const
{
    return grids[level]->endY;
}

template <typename Grid>
real MultipleGridBuilder<Grid>::getEndZ(uint level) const
{
    return grids[level]->endZ;
}

template <typename Grid>
void MultipleGridBuilder<Grid>::addGridToList(SPtr<Grid> grid)
{
    grids.push_back(grid);
}

template class MultipleGridBuilder<GridDummy>;
