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

void MultipleGridBuilder::addFineGrid(Object* gridShape, uint levelFine)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    const uint nodesBetweenGrids = 8;
    const uint levelDifference = levelFine - getNumberOfLevels();
    const uint oldGridSize = this->getNumberOfLevels();


    addIntermediateGridsToList(levelDifference, levelFine, nodesBetweenGrids, gridShape->clone());
    addFineGridToList(levelFine, gridShape->clone());

    eraseGridsFromListIfInvalid(oldGridSize);
}

void MultipleGridBuilder::addFineGridToList(uint level, Object* gridShape)
{
    const auto grid = makeGrid(gridShape, level);
    grids.push_back(grid);
}

void MultipleGridBuilder::addIntermediateGridsToList(uint levelDifference, uint levelFine, uint nodesBetweenGrids, Object* gridShape)
{
    if (levelDifference > 0)
    {
        auto spacings = getSpacingFactors(levelDifference);

        // start = startFine - SUM(nodesBetweenGrids * 2^i * dxfine) 
        uint level = getNumberOfLevels();
        for (int i = levelDifference - 1; i >= 0; i--)
        {
            const real scalingFactor = nodesBetweenGrids * spacings[i] * calculateDelta(levelFine);
            gridShape->scale(scalingFactor);

            const auto grid = makeGrid(gridShape, level++);
            grids.push_back(grid);
        }
    }
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

    eraseGridsFromListIfInvalid(oldGridSize);
}

void MultipleGridBuilder::addFineGridToList(uint level, real startXfine, real startYfine, real startZfine, real endXfine, real endYfine, real endZfine)
{
    const auto grid = makeGrid(startXfine, startYfine, startZfine, endXfine, endYfine, endZfine, level);
    grids.push_back(grid);
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

void MultipleGridBuilder::eraseGridsFromListIfInvalid(uint oldSize)
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
    return gridFactory->makeGrid(new Cuboid(startX, startY, startZ, endX, endY, endZ), startX, startY, startZ, endX, endY, endZ, delta);
}

bool MultipleGridBuilder::coarseGridExists() const
{
    return !grids.empty();
}

SPtr<Grid> MultipleGridBuilder::makeGrid(Object* gridShape, uint level) const
{
    const real delta = calculateDelta(level);

    auto staggeredCoordinates = getStaggeredCoordinates(gridShape->getX1Minimum(), gridShape->getX2Minimum(), gridShape->getX3Minimum(), gridShape->getX1Maximum(), gridShape->getX2Maximum(), gridShape->getX3Maximum(), delta);

    return gridFactory->makeGrid(gridShape, staggeredCoordinates[0], staggeredCoordinates[1], staggeredCoordinates[2], staggeredCoordinates[3], staggeredCoordinates[4], staggeredCoordinates[5], delta);
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

    const real startXStaggered = int(startX) + offset[0] - getDelta(getNumberOfLevels() - 1);
    const real startYStaggered = int(startY) + offset[1] - getDelta(getNumberOfLevels() - 1);
    const real startZStaggered = int(startZ) + offset[2] - getDelta(getNumberOfLevels() - 1);

    const real endXStaggered = int(endX) - offset[0] + getDelta(getNumberOfLevels() - 1);
    const real endYStaggered = int(endY) - offset[1] + getDelta(getNumberOfLevels() - 1);
    const real endZStaggered = int(endZ) - offset[2] + getDelta(getNumberOfLevels() - 1);

    return std::array<real, 6>{startXStaggered, startYStaggered, startZStaggered, endXStaggered, endYStaggered, endZStaggered};
}

std::array<real, 3> MultipleGridBuilder::getOffset(real delta) const
{
    real offsetToNullStart = delta * 0.5;

    for (uint level = 1; level < getNumberOfLevels(); level++)
        offsetToNullStart += getDelta(level) * 0.5;

    const real offsetX = getStartX(0) - int(getStartX(0)) + offsetToNullStart;
    const real offsetY = getStartY(0) - int(getStartY(0)) + offsetToNullStart;
    const real offsetZ = getStartZ(0) - int(getStartZ(0)) + offsetToNullStart;
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
    for (size_t i = 1; i < grids.size(); i++)
        grids[i]->setPeriodicity(false, false, false);

    for (auto grid : grids)
        grid->inital();

    for (size_t i = 0; i < grids.size() - 1; i++)
        grids[i]->findGridInterface(grids[i + 1]);


    for (size_t i = 0; i < grids.size() - 1; i++)
        grids[i]->findSparseIndices(grids[i + 1]);

    grids[grids.size() - 1]->findSparseIndices(nullptr);
}


void MultipleGridBuilder::emitNoCoarseGridExistsWarning()
{
    *logging::out << logging::Logger::WARNING << "No Coarse grid was added before. Actual Grid is not added, please create coarse grid before.\n";
}


void MultipleGridBuilder::emitGridIsNotInCoarseGridWarning()
{
    *logging::out << logging::Logger::WARNING << "Grid lies not inside of coarse grid. Actual Grid is not added.\n";
}
