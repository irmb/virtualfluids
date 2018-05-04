#include "MultipleGridBuilder.h"
#include "utilities/math/Math.h"
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
    const auto grid = this->makeGrid(new Cuboid(startX, startY, startZ, endX, endY, endZ), startX, startY, startZ, endX, endY, endZ, delta);
    addGridToList(grid);
}

void MultipleGridBuilder::addGeometry(Object* solidObject)
{
    this->solidObject = solidObject;
}


void MultipleGridBuilder::addGrid(Object* gridShape)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    const auto grid = makeGrid(gridShape, getNumberOfLevels());

    addGridToListIfValid(grid);
}

void MultipleGridBuilder::addGrid(Object* gridShape, uint levelFine)
{
    if (!coarseGridExists())
        return emitNoCoarseGridExistsWarning();

    const uint nodesBetweenGrids = 8;
    const uint levelDifference = levelFine - getNumberOfLevels();
    const uint oldGridSize = this->getNumberOfLevels();

    addIntermediateGridsToList(levelDifference, levelFine, nodesBetweenGrids, gridShape);
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
            Object* gridShapeClone = gridShape->clone();
            gridShapeClone->scale(scalingFactor);

            const auto grid = makeGrid(gridShapeClone, level++);
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


SPtr<Grid> MultipleGridBuilder::makeGrid(Object* gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) const
{
    return gridFactory->makeGrid(gridShape, startX, startY, startZ, endX, endY, endZ, delta);
}

bool MultipleGridBuilder::coarseGridExists() const
{
    return !grids.empty();
}

SPtr<Grid> MultipleGridBuilder::makeGrid(Object* gridShape, uint level) const
{
    const real delta = calculateDelta(level);

    auto staggeredCoordinates = getStaggeredCoordinates(gridShape->getX1Minimum(), gridShape->getX2Minimum(), gridShape->getX3Minimum(), gridShape->getX1Maximum(), gridShape->getX2Maximum(), gridShape->getX3Maximum(), delta);

    return this->makeGrid(gridShape, staggeredCoordinates[0], staggeredCoordinates[1], staggeredCoordinates[2], staggeredCoordinates[3], staggeredCoordinates[4], staggeredCoordinates[5], delta);
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

    const real startXStaggered = std::floor(startX) - offset[0];
    const real startYStaggered = std::floor(startY) - offset[1];
    const real startZStaggered = std::floor(startZ) - offset[2];

    const real endXStaggered = std::ceil(endX) + offset[0];
    const real endYStaggered = std::ceil(endY) + offset[1];
    const real endZStaggered = std::ceil(endZ) + offset[2];

    auto temporaryGrid = this->makeGrid(nullptr, startXStaggered, startYStaggered, startZStaggered, endXStaggered, endYStaggered, endZStaggered, delta);
    auto startStaggered = temporaryGrid->getMinimumOnNode(Vertex(startX, startY, startZ));
    auto endStaggered = temporaryGrid->getMaximumOnNode(Vertex(endX, endY, endZ));

    return std::array<real, 6>{startStaggered.x, startStaggered.y, startStaggered.z, endStaggered.x, endStaggered.y, endStaggered.z};
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
        vf::Math::greaterEqual(grid->getStartX(), grids[0]->getStartX()) &&
        vf::Math::greaterEqual(grid->getStartY(), grids[0]->getStartY()) &&
        vf::Math::greaterEqual(grid->getStartZ(), grids[0]->getStartZ()) &&
        vf::Math::lessEqual(grid->getEndX(), grids[0]->getEndX()) &&
        vf::Math::lessEqual(grid->getEndY(), grids[0]->getEndY()) &&
        vf::Math::lessEqual(grid->getEndZ(), grids[0]->getEndZ());
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


     for (auto grid : grids)
          grid->mesh(solidObject);

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
