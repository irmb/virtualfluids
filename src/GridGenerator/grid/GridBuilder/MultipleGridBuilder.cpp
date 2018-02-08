#include "MultipleGridBuilder.h"
#include "utilities/math/CudaMath.cuh"
#include "../Grid.h"
#include "../GridFactory.h"

#include <VirtualFluidsBasics/utilities/logger/Logger.h>

MultipleGridBuilder::MultipleGridBuilder(SPtr<GridFactory> gridFactory, Device device, const std::string &d3qxx) :
    LevelGridBuilder(device, d3qxx), gridFactory(gridFactory)
{

}

SPtr<MultipleGridBuilder> MultipleGridBuilder::makeShared(SPtr<GridFactory> gridFactory)
{
    return SPtr<MultipleGridBuilder>(new MultipleGridBuilder(gridFactory));
}


void MultipleGridBuilder::addCoarseGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta)
{
    auto grid = this->makeGrid(startX, startY, startZ, endX, endY, endZ, delta);
    addGridToList(grid);
}


void MultipleGridBuilder::addGrid(real startX, real startY, real startZ, real endX, real endY, real endZ)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    auto grid = makeGrid(startX, startY, startZ, endX, endY, endZ, getNumberOfLevels());

    if (!isGridInCoarseGrid(grid))
        return emitGridIsNotInCoarseGridWarning();

    addGridToList(grid);
}


void MultipleGridBuilder::addFineGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, uint level)
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
        const real spaceBetweenInterface = nodesBetweenGrids * getDelta(i);
        const real staggeredToFine = getDelta(i + 1) * 0.5;

        grids[i]->setStartX(grids[i + 1]->getStartX() - staggeredToFine - spaceBetweenInterface);
        grids[i]->setStartY(grids[i + 1]->getStartY() - staggeredToFine - spaceBetweenInterface);
        grids[i]->setStartZ(grids[i + 1]->getStartZ() - staggeredToFine - spaceBetweenInterface);

        grids[i]->setEndX(grids[i + 1]->getEndX() + staggeredToFine + spaceBetweenInterface);
        grids[i]->setEndY(grids[i + 1]->getEndY() + staggeredToFine + spaceBetweenInterface);
        grids[i]->setEndZ(grids[i + 1]->getEndZ() + staggeredToFine + spaceBetweenInterface);

        if (!isGridInCoarseGrid(grids[i]))
        {
            grids.erase(grids.begin() + actualLevelSize, grids.end());
            return emitGridIsNotInCoarseGridWarning();
        }
    }
}


SPtr<Grid> MultipleGridBuilder::makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) const
{
    return gridFactory->makeGrid(startX, startY, startZ, endX, endY, endZ, delta);
}


bool MultipleGridBuilder::coarseGridExists() const
{
    return !grids.empty();
}




SPtr<Grid> MultipleGridBuilder::makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, uint level) const
{
    const real delta = calculateDelta(level);

    auto staggeredCoordinates = getStaggeredCoordinates(startX, startY, startZ, endX, endY, endZ, delta);

    return this->makeGrid(staggeredCoordinates[0], staggeredCoordinates[1], staggeredCoordinates[2], staggeredCoordinates[3], staggeredCoordinates[4], staggeredCoordinates[5], delta);
}


real MultipleGridBuilder::calculateDelta(uint level) const
{
    real delta = this->getDelta(0);
    for (uint i = 0; i < level; i++)
        delta *= 0.5;
    return delta;
}


std::array<real, 6> MultipleGridBuilder::getStaggeredCoordinates(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) const
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


bool MultipleGridBuilder::isGridInCoarseGrid(SPtr<Grid> grid) const
{
    return
        CudaMath::greaterEqual(grid->getStartX(), grids[0]->getStartX()) &&
        CudaMath::greaterEqual(grid->getStartY(), grids[0]->getStartY()) &&
        CudaMath::greaterEqual(grid->getStartZ(), grids[0]->getStartZ()) &&
        CudaMath::lessEqual(grid->getEndX(), grids[0]->getEndX()) &&
        CudaMath::lessEqual(grid->getEndY(), grids[0]->getEndY()) &&
        CudaMath::lessEqual(grid->getEndZ(), grids[0]->getEndZ());
}



uint MultipleGridBuilder::getNumberOfLevels() const
{
    return uint(grids.size());
}


real MultipleGridBuilder::getDelta(uint level) const
{
    if (grids.size() <= level)
        throw std::exception("delta from invalid level was required.");
    return grids[level]->getDelta();
}


real MultipleGridBuilder::getStartX(uint level) const
{
    return grids[level]->getStartX();
}


real MultipleGridBuilder::getStartY(uint level) const
{
    return grids[level]->getStartY();
}


real MultipleGridBuilder::getStartZ(uint level) const
{
    return grids[level]->getStartZ();
}


real MultipleGridBuilder::getEndX(uint level) const
{
    return grids[level]->getEndX();
}


real MultipleGridBuilder::getEndY(uint level) const
{
    return grids[level]->getEndY();
}


real MultipleGridBuilder::getEndZ(uint level) const
{
    return grids[level]->getEndZ();
}


void MultipleGridBuilder::addGridToList(SPtr<Grid> grid)
{
    grids.push_back(grid);
}


std::vector<SPtr<Grid> > MultipleGridBuilder::getGrids() const
{
    return this->grids;
}

void MultipleGridBuilder::createGridInterfaces()
{
    for (size_t i = 1; i < grids.size(); i++)
    {
        grids[i]->setPeriodicity(false, false, false);
    }

    for(size_t i = grids.size() - 1; i > 0; i--) {
        grids[i - 1]->removeOverlapNodes(grids[i]);
    }
}

void MultipleGridBuilder::allocateGridMemory()
{
    for (auto grid : grids)
        grid->allocateGridMemory();
}


void MultipleGridBuilder::emitNoCoarseGridExistsWarning()
{
    *logging::out << logging::Logger::WARNING << "No Coarse grid was added before. Actual Grid is not added, please create coarse grid before.\n";
}


void MultipleGridBuilder::emitGridIsNotInCoarseGridWarning()
{
    *logging::out << logging::Logger::WARNING << "Grid lies not inside of coarse grid. Actual Grid is not added.\n";
}
