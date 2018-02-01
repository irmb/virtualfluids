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
    auto grid = GridDummy::makeShared(startX, startY, startZ, endX, endY, endZ, delta);
    grids.push_back(grid);
}


template <typename Grid>
void MultipleGridBuilder<Grid>::addGrid(real startX, real startY, real startZ, real endX, real endY, real endZ)
{
    checkIfCoarseGridIsMissing();

    const real delta = calculateDelta();
    auto grid = GridDummy::makeShared(startX, startY, startZ, endX, endY, endZ, delta);

    checkIfGridIsInCoarseGrid(grid);

    grids.push_back(grid);
}

template <typename Grid>
void MultipleGridBuilder<Grid>::checkIfCoarseGridIsMissing() const
{
    if (grids.empty())
        throw FirstGridMustBeCoarseException();
}

template <typename Grid>
real MultipleGridBuilder<Grid>::calculateDelta()
{
    const real delta = grids[grids.size() - 1]->getDelta() / 2.0;
    return delta;
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
    return grids[level]->getDelta();
}






template class MultipleGridBuilder<GridDummy>;
