#include "MultipleGridBuilder.h"
#include "utilities/math/CudaMath.cuh"
#include "../GridMocks.h"
#include "../GridFactory.h"

#include <VirtualFluidsBasics/utilities/logger/Logger.h>

template<typename Grid>
MultipleGridBuilder<Grid>::MultipleGridBuilder(SPtr<GridFactory<Grid> > gridFactory) : gridFactory(gridFactory)
{

}

template<typename Grid>
SPtr<MultipleGridBuilder<Grid> > MultipleGridBuilder<Grid>::makeShared(SPtr<GridFactory<Grid> > gridFactory)
{
    return SPtr<MultipleGridBuilder>(new MultipleGridBuilder<Grid>(gridFactory));
}

template <typename Grid>
void MultipleGridBuilder<Grid>::addCoarseGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta)
{
    auto grid = this->makeGrid(startX, startY, startZ, endX, endY, endZ, delta);
    addGridToList(grid);
}

template <typename Grid>
void MultipleGridBuilder<Grid>::addGrid(real startX, real startY, real startZ, real endX, real endY, real endZ)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    auto grid = makeGrid(startX, startY, startZ, endX, endY, endZ, getNumberOfLevels());

    if (!isGridInCoarseGrid(grid))
        return emitGridIsNotInCoarseGridWarning();

    addGridToList(grid);
}

template <typename Grid>
void MultipleGridBuilder<Grid>::addFineGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, uint level)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    const uint nodesBetweenGrids = 8;

    const uint actualLevelSize = getNumberOfLevels();

    real offsetX = getStartX(0) - int(getStartX(0));
    real offsetY = getStartY(0) - int(getStartY(0));
    real offsetZ = getStartZ(0) - int(getStartZ(0));


    for (uint i = 1; i <= level; i++)
    {
        const real delta = calculateDelta(i);
        offsetX += delta * 0.5;
        offsetY += delta * 0.5;
        offsetZ += delta * 0.5;

        if (i >= actualLevelSize)
            addGridToList(this->makeGrid(startX + offsetX, startY + offsetY, startZ + offsetZ, endX - offsetX, endY - offsetY, endZ - offsetZ, delta));
    }

    for (auto i = level - 1; i >= actualLevelSize; i--)
    {
        real spaceBetweenInterface = nodesBetweenGrids * getDelta(i);
        real staggeredToFine = getDelta(i + 1) * 0.5;

        grids[i]->startX = grids[i + 1]->startX - staggeredToFine - spaceBetweenInterface;
        grids[i]->startY = grids[i + 1]->startY - staggeredToFine - spaceBetweenInterface;
        grids[i]->startZ = grids[i + 1]->startZ - staggeredToFine - spaceBetweenInterface;

        grids[i]->endX = grids[i + 1]->endX + staggeredToFine + spaceBetweenInterface;
        grids[i]->endY = grids[i + 1]->endY + staggeredToFine + spaceBetweenInterface;
        grids[i]->endZ = grids[i + 1]->endZ + staggeredToFine + spaceBetweenInterface;

        if (!isGridInCoarseGrid(grids[i]))
        {
            grids.erase(grids.begin() + actualLevelSize, grids.end());
            return emitGridIsNotInCoarseGridWarning();
        }
    }
}

template <typename Grid>
SPtr<Grid> MultipleGridBuilder<Grid>::makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) const
{
    return this->gridFactory->makeGrid(startX, startY, startZ, endX, endY, endZ, delta);
}

template <typename Grid>
bool MultipleGridBuilder<Grid>::coarseGridExists() const
{
    return !grids.empty();
}



template <typename Grid>
SPtr<Grid> MultipleGridBuilder<Grid>::makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, uint level) const
{
    const real delta = calculateDelta(level);

    auto staggeredCoordinates = getStaggeredCoordinates(startX, startY, startZ, endX, endY, endZ, delta);

    return this->makeGrid(staggeredCoordinates[0], staggeredCoordinates[1], staggeredCoordinates[2], staggeredCoordinates[3], staggeredCoordinates[4], staggeredCoordinates[5], delta);
}

template <typename Grid>
real MultipleGridBuilder<Grid>::calculateDelta(uint level) const
{
    real delta = this->getDelta(0);
    for (uint i = 0; i < level; i++)
        delta *= 0.5;
    return delta;
}

template <typename Grid>
std::array<real, 6> MultipleGridBuilder<Grid>::getStaggeredCoordinates(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) const
{
    const uint levelIndexCoarseGrid = getNumberOfLevels() - 1;
    const real offsetToNullStart = delta * 0.5;
    const real offsetX = getStartX(levelIndexCoarseGrid) - int(getStartX(levelIndexCoarseGrid)) + offsetToNullStart;
    const real offsetY = getStartY(levelIndexCoarseGrid) - int(getStartY(levelIndexCoarseGrid)) + offsetToNullStart;
    const real offsetZ = getStartZ(levelIndexCoarseGrid) - int(getStartZ(levelIndexCoarseGrid)) + offsetToNullStart;

    const real startXStaggered = startX + offsetX;
    const real startYStaggered = startY + offsetY;
    const real startZStaggered = startZ + offsetZ;

    const real endXStaggered = endX - offsetX;
    const real endYStaggered = endY - offsetY;
    const real endZStaggered = endZ - offsetZ;

    return std::array<real, 6>{startXStaggered, startYStaggered, startZStaggered, endXStaggered, endYStaggered, endZStaggered};
}

template <typename Grid>
bool MultipleGridBuilder<Grid>::isGridInCoarseGrid(SPtr<Grid> grid) const
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
real MultipleGridBuilder<Grid>::getDelta(uint level) const
{
    if (grids.size() <= level)
        throw std::exception("delta from invalid level was required.");
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

template <typename Grid>
std::vector<SPtr<Grid> > MultipleGridBuilder<Grid>::getGrids() const
{
    return this->grids;
}

template <typename Grid>
void MultipleGridBuilder<Grid>::emitNoCoarseGridExistsWarning()
{
    *logging::out << logging::Logger::WARNING << "No Coarse grid was added before. Actual Grid is not added, please create coarse grid before.\n";
}

template <typename Grid>
void MultipleGridBuilder<Grid>::emitGridIsNotInCoarseGridWarning()
{
    *logging::out << logging::Logger::WARNING << "Grid lies not inside of coarse grid. Actual Grid is not added.\n";
}


template class MultipleGridBuilder<GridDummy>;
