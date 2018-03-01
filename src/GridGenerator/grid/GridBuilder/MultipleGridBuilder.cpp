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
    const auto grid = this->makeGrid(startX, startY, startZ, endX, endY, endZ, delta);
    addGridToList(grid);
}

void MultipleGridBuilder::addGrid(Object* gridShape)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    const auto grid = makeGrid(gridShape, getNumberOfLevels());

    addGridToListIfValid(grid);
}

void MultipleGridBuilder::addGrid(real startX, real startY, real startZ, real endX, real endY, real endZ)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    const auto grid = makeGrid(startX, startY, startZ, endX, endY, endZ, getNumberOfLevels());

    addGridToListIfValid(grid);
}

void MultipleGridBuilder::addFineGrid(real startXfine, real startYfine, real startZfine, real endXfine, real endYfine, real endZfine, uint levelFine)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    const uint nodesBetweenGrids = 8;
    const uint levelDifference = levelFine - getNumberOfLevels();
    const uint oldGridSize = this->getNumberOfLevels();

    addIntermediateGridsToList(levelDifference, levelFine, nodesBetweenGrids, startXfine, startYfine, startZfine, endXfine, endYfine, endZfine);
    addFineGridToList(levelFine, startXfine, startYfine, startZfine, endXfine, endYfine, endZfine);

    addGridsToListIfValid(oldGridSize);
}

void MultipleGridBuilder::addFineGridToList(uint level, real startXfine, real startYfine, real startZfine, real endXfine, real endYfine, real endZfine)
{
    grids.push_back(makeGrid(startXfine, startYfine, startZfine, endXfine, endYfine, endZfine, level));
}

void MultipleGridBuilder::addIntermediateGridsToList(uint levelDifference, uint levelFine, uint nodesBetweenGrids, real startXfine, real startYfine, real startZfine, real endXfine, real endYfine, real endZfine)
{
    if (levelDifference > 0) 
    {
        auto spacings = getSpacingFactors(levelDifference);

        // start = startFine - SUM(nodesBetweenGrids * 2^i * dxfine) 
        uint level = getNumberOfLevels();
        for (int i = levelDifference - 1; i >= 0; i--)
        {
            const real startX = startXfine - nodesBetweenGrids * spacings[i] * calculateDelta(levelFine);
            const real startY = startYfine - nodesBetweenGrids * spacings[i] * calculateDelta(levelFine);
            const real startZ = startZfine - nodesBetweenGrids * spacings[i] * calculateDelta(levelFine);

            const real endX = endXfine + nodesBetweenGrids * spacings[i] * calculateDelta(levelFine);
            const real endY = endYfine + nodesBetweenGrids * spacings[i] * calculateDelta(levelFine);
            const real endZ = endZfine + nodesBetweenGrids * spacings[i] * calculateDelta(levelFine);

            const auto grid = makeGrid(startX, startY, startZ, endX, endY, endZ, level++);
            grids.push_back(grid);
        }
    }
}

std::vector<uint> MultipleGridBuilder::getSpacingFactors(uint levelDifference) const
{
    std::vector<uint> spacings(levelDifference);

    spacings[0] = uint(std::pow(2, 1));
    for (uint i = 1; i < levelDifference; i++)
        spacings[i] = spacings[i - 1] + uint(std::pow(2, i + 1));

    return spacings;
}

void MultipleGridBuilder::addGridsToListIfValid(uint oldSize)
{
    if (!isGridInCoarseGrid(grids[oldSize]))
    {
        this->grids.erase(grids.begin() + oldSize, grids.end());
        emitGridIsNotInCoarseGridWarning();
    }
}

void MultipleGridBuilder::addGridToListIfValid(SPtr<Grid> grid)
{
    if (!isGridInCoarseGrid(grid))
        return emitGridIsNotInCoarseGridWarning();

    addGridToList(grid);
}


SPtr<Grid> MultipleGridBuilder::makeGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) const
{
    return gridFactory->makeGrid(new Cuboid(startX, startY, startZ, endX, endY, endZ), delta);
}

bool MultipleGridBuilder::coarseGridExists() const
{
    return !grids.empty();
}

SPtr<Grid> MultipleGridBuilder::makeGrid(Object* gridShape, uint level) const
{
    const real delta = calculateDelta(level);

    auto staggeredCoordinates = getStaggeredCoordinates(gridShape->getX1Minimum(), gridShape->getX2Minimum(), gridShape->getX3Minimum(), gridShape->getX1Maximum(), gridShape->getX2Maximum(), gridShape->getX3Maximum(), delta);
    
    gridShape->setX1Minimum(staggeredCoordinates[0]);
    gridShape->setX2Minimum(staggeredCoordinates[1]);
    gridShape->setX3Minimum(staggeredCoordinates[2]);

    gridShape->setX1Maximum(staggeredCoordinates[3]);
    gridShape->setX2Maximum(staggeredCoordinates[4]);
    gridShape->setX3Maximum(staggeredCoordinates[5]);

    return gridFactory->makeGrid(gridShape, delta);
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
    auto offset = getOffset(delta);

    const real startXStaggered = int(startX) + offset[0];
    const real startYStaggered = int(startY) + offset[1];
    const real startZStaggered = int(startZ) + offset[2];

    const real endXStaggered = int(endX) - offset[0];
    const real endYStaggered = int(endY) - offset[1];
    const real endZStaggered = int(endZ) - offset[2];

    return std::array<real, 6>{startXStaggered, startYStaggered, startZStaggered, endXStaggered, endYStaggered, endZStaggered};
}

std::array<real, 3> MultipleGridBuilder::getOffset(real delta) const
{
    const uint levelIndexCoarseGrid = getNumberOfLevels() - 1;
    const real offsetToNullStart = delta * 0.5;
    const real offsetX = getStartX(levelIndexCoarseGrid) - int(getStartX(levelIndexCoarseGrid)) + offsetToNullStart;
    const real offsetY = getStartY(levelIndexCoarseGrid) - int(getStartY(levelIndexCoarseGrid)) + offsetToNullStart;
    const real offsetZ = getStartZ(levelIndexCoarseGrid) - int(getStartZ(levelIndexCoarseGrid)) + offsetToNullStart;
    return std::array<real, 3>{offsetX, offsetY, offsetZ};
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

void MultipleGridBuilder::buildGrids()
{
    for (auto grid : grids)
        grid->allocateGridMemory();

    for (size_t i = 1; i < grids.size(); i++)
        grids[i]->setPeriodicity(false, false, false);

    for(size_t i = grids.size() - 1; i > 0; i--)
        grids[i - 1]->findGridInterface(grids[i]);

    //for (int i=0; i< grids[0]->getSize(); i++)
    //{
    //    printf("type: %d\n", grids[0]->getFieldEntry(i));
    //}
}


void MultipleGridBuilder::emitNoCoarseGridExistsWarning()
{
    *logging::out << logging::Logger::WARNING << "No Coarse grid was added before. Actual Grid is not added, please create coarse grid before.\n";
}


void MultipleGridBuilder::emitGridIsNotInCoarseGridWarning()
{
    *logging::out << logging::Logger::WARNING << "Grid lies not inside of coarse grid. Actual Grid is not added.\n";
}
