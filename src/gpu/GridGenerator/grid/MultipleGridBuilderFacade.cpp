#include <algorithm>
#include <climits>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>

#include "MultipleGridBuilderFacade.h"
#include "geometries/BoundingBox/BoundingBox.h"
#include "geometries/Object.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"
#include "grid/GridDimensions.h"

using namespace CommunicationDirections;

MultipleGridBuilderFacade::MultipleGridBuilderFacade(SPtr<MultipleGridBuilder> gridBuilder,
                                                     SPtr<GridDimensions> gridDimensions,
                                                     std::optional<real> overlapOfSubdomains)
    : gridBuilder(std::move(gridBuilder)), gridDimensionsDomain(std::move(gridDimensions)),
      overlapOfSubdomains(overlapOfSubdomains)
{
}

MultipleGridBuilderFacade::MultipleGridBuilderFacade(SPtr<GridDimensions> gridDimensions,
                                                     std::optional<real> overlapOfSubdomains)
    : MultipleGridBuilderFacade(std::make_shared<MultipleGridBuilder>(), std::move(gridDimensions), overlapOfSubdomains)
{
}

void MultipleGridBuilderFacade::createGrids(uint generatePart)
{
    if (createGridsHasBeenCalled)
        throw std::runtime_error("MultipleGridBuilderFacade::createGrids() has been called more than once.");
    createGridsHasBeenCalled = true;

    this->calculateNumberOfSubdomains();
    this->numberOfSubdomainsTotal =
        this->numberOfSubdomains[Axis::x] * this->numberOfSubdomains[Axis::y] * this->numberOfSubdomains[Axis::z];

    if (numberOfSubdomainsTotal > 1 && !this->overlapOfSubdomains)
        throw std::runtime_error("OverlapOfSubdomains in MultipleGridBuilderFacade is NaN.");

    if (generatePart >= numberOfSubdomainsTotal)
        throw std::runtime_error("Invalid id for subdomain: It is greater or equal to numberOfSubdomains");

    this->sortSplitLocations();

    this->checkSplitLocations(this->xSplits, this->gridDimensionsDomain->minX, this->gridDimensionsDomain->maxX);
    this->checkSplitLocations(this->ySplits, this->gridDimensionsDomain->minY, this->gridDimensionsDomain->maxY);
    this->checkSplitLocations(this->zSplits, this->gridDimensionsDomain->minZ, this->gridDimensionsDomain->maxZ);

    this->calculatedIndexOfPart(generatePart);

    this->checkForNeighbors();

    this->configureSubDomainGrids();

    this->addGeometriesToGridBuilder();

    gridBuilder->buildGrids(false); // buildGrids() has to be called before setting the BCs!!!!

    this->setUpCommunicationNeighbors();
}

void MultipleGridBuilderFacade::calculateNumberOfSubdomains()
{
    this->numberOfSubdomains[Axis::x] = (uint)(this->xSplits.size() + 1);
    this->numberOfSubdomains[Axis::y] = (uint)(this->ySplits.size() + 1);
    this->numberOfSubdomains[Axis::z] = (uint)(this->zSplits.size() + 1);
}

void MultipleGridBuilderFacade::sortSplitLocations()
{

    std::sort(this->xSplits.begin(), this->xSplits.end());
    std::sort(this->ySplits.begin(), this->ySplits.end());
    std::sort(this->zSplits.begin(), this->zSplits.end());
}

void MultipleGridBuilderFacade::calculatedIndexOfPart(uint generatePart)
{
    this->index.at(Axis::x) = this->getX3D(generatePart);
    this->index.at(Axis::y) = this->getY3D(generatePart);
    this->index.at(Axis::z) = this->getZ3D(generatePart);
}

void MultipleGridBuilderFacade::checkSplitLocations(const std::vector<real>& splits, real lowerBound, real upperBound) const
{
    if (splits.empty()) return;

    if (splits.front() < lowerBound)
        throw std::runtime_error("The domain split value " + std::to_string(splits.front()) +
                                 " is smaller than the lower bound (" + std::to_string(lowerBound) + ") of the domain.");
    if (splits.back() > upperBound)
        throw std::runtime_error("A domain split value " + std::to_string(splits.back()) +
                                 " is larger than the upper bound (" + std::to_string(upperBound) + ") of the domain.");

    auto iteratorOfDuplicate = std::adjacent_find(splits.begin(), splits.end());
    if (iteratorOfDuplicate != splits.end())
        throw std::runtime_error("The domain split value " + std::to_string(*iteratorOfDuplicate) +
                                 " was added multiple times for the same coordinate direction");
}

void MultipleGridBuilderFacade::configureSubDomainGrids()
{

    // Example: 2 subdomains in x
    // xSplits = {0}
    // xMin = {-1}
    // xMax = {1}

    //  vector xValues (bounding boxes of the subdomains)
    //
    //    /------- subdmain 0 -------/
    //    |                          |
    //  -1.0                        0.0                        1.0
    //                               |                          |
    //                               /------- subdmain 1 -------/
    //
    //  xValues = {-1.0, 0.0, 1.0}

    // create vector with the coordinates of the subdomain's bounding boxes
    std::vector<real> xValues = { this->gridDimensionsDomain->minX, this->gridDimensionsDomain->maxX };
    xValues.insert(std::prev(xValues.end()), this->xSplits.begin(), this->xSplits.end());
    std::vector<real> yValues = { this->gridDimensionsDomain->minY, this->gridDimensionsDomain->maxY };
    yValues.insert(std::prev(yValues.end()), this->ySplits.begin(), this->ySplits.end());
    std::vector<real> zValues = { this->gridDimensionsDomain->minZ, this->gridDimensionsDomain->maxZ };
    zValues.insert(std::prev(zValues.end()), this->zSplits.begin(), this->zSplits.end());

    real xMinCoarseGrid = xValues[index.at(Axis::x)];
    real yMinCoarseGrid = yValues[index.at(Axis::y)];
    real zMinCoarseGrid = zValues[index.at(Axis::z)];
    real xMaxCoarseGrid = xValues[index.at(Axis::x) + 1];
    real yMaxCoarseGrid = yValues[index.at(Axis::y) + 1];
    real zMaxCoarseGrid = zValues[index.at(Axis::z) + 1];

    // add overlap
    xMinCoarseGrid -= (hasNeighbors[CommunicationDirections::MX]) ? overlapOfSubdomains.value() : 0;
    yMinCoarseGrid -= (hasNeighbors[CommunicationDirections::MY]) ? overlapOfSubdomains.value() : 0;
    zMinCoarseGrid -= (hasNeighbors[CommunicationDirections::MZ]) ? overlapOfSubdomains.value() : 0;
    xMaxCoarseGrid += (hasNeighbors[CommunicationDirections::PX]) ? overlapOfSubdomains.value() : 0;
    yMaxCoarseGrid += (hasNeighbors[CommunicationDirections::PY]) ? overlapOfSubdomains.value() : 0;
    zMaxCoarseGrid += (hasNeighbors[CommunicationDirections::PZ]) ? overlapOfSubdomains.value() : 0;

    // add coarse grid
    gridBuilder->addCoarseGrid(xMinCoarseGrid, yMinCoarseGrid, zMinCoarseGrid, xMaxCoarseGrid, yMaxCoarseGrid,
                               zMaxCoarseGrid, gridDimensionsDomain->delta);

    // add fine grids for grid refinement
    this->addFineGridsToGridBuilder();

    // set subdomain boxes
    // subdomain boxes are only needed on multiple gpus
    if ((numberOfSubdomains[Axis::x] * numberOfSubdomains[Axis::y] * numberOfSubdomains[Axis::z]) > 1) {
        gridBuilder->setSubDomainBox(std::make_shared<BoundingBox>(xValues[index.at(Axis::x)], xValues[index.at(Axis::x) + 1],
                                                                   yValues[index.at(Axis::y)], yValues[index.at(Axis::y) + 1],
                                                                   zValues[index.at(Axis::z)], zValues[index.at(Axis::z) + 1]));
    }
}

void MultipleGridBuilderFacade::setUpCommunicationNeighbors()
{
    // Communication is only needed on multiple gpus
    if (numberOfSubdomainsTotal == 1) return;

    if (hasNeighbors.empty())
        throw std::runtime_error("checkForNeighbors() has to be called befor calling setUpCommunicationNeighbors()");

    for (const auto &[direction, hasNeighbor] : hasNeighbors) {
        if (hasNeighbor) {
            gridBuilder->findCommunicationIndices(direction);

            switch (direction) {
                case CommunicationDirections::MX:
                    gridBuilder->setCommunicationProcess(
                        direction, getIndex1D(index.at(Axis::x) - 1, index.at(Axis::y), index.at(Axis::z)));
                    break;
                case CommunicationDirections::MY:
                    gridBuilder->setCommunicationProcess(
                        direction, getIndex1D(index.at(Axis::x), index.at(Axis::y) - 1, index.at(Axis::z)));
                    break;
                case CommunicationDirections::MZ:
                    gridBuilder->setCommunicationProcess(
                        direction, getIndex1D(index.at(Axis::x), index.at(Axis::y), index.at(Axis::z) - 1));
                    break;
                case CommunicationDirections::PX:
                    gridBuilder->setCommunicationProcess(
                        direction, getIndex1D(index.at(Axis::x) + 1, index.at(Axis::y), index.at(Axis::z)));
                    break;
                case CommunicationDirections::PY:
                    gridBuilder->setCommunicationProcess(
                        direction, getIndex1D(index.at(Axis::x), index.at(Axis::y) + 1, index.at(Axis::z)));
                    break;
                case CommunicationDirections::PZ:
                    gridBuilder->setCommunicationProcess(
                        direction, getIndex1D(index.at(Axis::x), index.at(Axis::y), index.at(Axis::z) + 1));
                    break;
            }
        }
    }
}

void MultipleGridBuilderFacade::checkForNeighbors()
{
    hasNeighbors[CommunicationDirections::MX] = (index.at(Axis::x) > 0);
    hasNeighbors[CommunicationDirections::MY] = (index.at(Axis::y) > 0);
    hasNeighbors[CommunicationDirections::MZ] = (index.at(Axis::z) > 0);
    hasNeighbors[CommunicationDirections::PX] = (index.at(Axis::x) < numberOfSubdomains[Axis::x] - 1);
    hasNeighbors[CommunicationDirections::PY] = (index.at(Axis::y) < numberOfSubdomains[Axis::y] - 1);
    hasNeighbors[CommunicationDirections::PZ] = (index.at(Axis::z) < numberOfSubdomains[Axis::z] - 1);
}

void MultipleGridBuilderFacade::addFineGridsToGridBuilder()
{
    for (auto const& grid : fineGrids) {
        gridBuilder->addGrid(grid.first, grid.second);
    }
}

void MultipleGridBuilderFacade::addGeometriesToGridBuilder()
{
    for (auto const& geometry : geometries) {
        gridBuilder->addGeometry(geometry);
    }
}

void MultipleGridBuilderFacade::setOverlapOfSubdomains(real overlap)
{
    this->overlapOfSubdomains = overlap;
}

void MultipleGridBuilderFacade::addDomainSplit(real coordinate, Axis direction)
{
    if (this->createGridsHasBeenCalled)
        throw std::runtime_error("MultipleGridBuilderFacade::addSplit() should be called before createGrids().");

    switch (direction) {
        case x:
            this->xSplits.push_back(coordinate);
            break;
        case y:
            this->ySplits.push_back(coordinate);
            break;
        case z:
            this->zSplits.push_back(coordinate);
            break;
    }
}

void MultipleGridBuilderFacade::addFineGrid(std::shared_ptr<Object> gridShape, uint levelFine)
{
    if (this->createGridsHasBeenCalled)
        throw std::runtime_error("MultipleGridBuilderFacade::addFineGrid() should be called before createGrids().");

    this->fineGrids.emplace_back(gridShape, levelFine);
}

void MultipleGridBuilderFacade::addGeometry(std::shared_ptr<Object> gridShape)
{
    if (this->createGridsHasBeenCalled)
        throw std::runtime_error("MultipleGridBuilderFacade::addGeometry() should be called before createGrids().");

    this->geometries.emplace_back(gridShape);
}

void MultipleGridBuilderFacade::setNumberOfLayersForRefinement(uint numberOfLayersFine,
                                                               uint numberOfLayersBetweenLevels) const
{
    gridBuilder->setNumberOfLayers(numberOfLayersFine, numberOfLayersBetweenLevels);
}

uint MultipleGridBuilderFacade::getX3D(uint index1D) const
{
    return index1D % numberOfSubdomains[Axis::x];
}

uint MultipleGridBuilderFacade::getY3D(uint index1D) const
{
    return (index1D / numberOfSubdomains[Axis::x]) % numberOfSubdomains[Axis::y];
}

uint MultipleGridBuilderFacade::getZ3D(uint index1D) const
{
    return index1D / (numberOfSubdomains[Axis::x] * numberOfSubdomains[Axis::y]);
}

std::array<uint, 3> MultipleGridBuilderFacade::convertToIndices3D(uint index1D) const
{
    const uint xPos = getX3D(index1D);
    const uint yPos = getY3D(index1D);
    const uint zPos = getZ3D(index1D);
    return { xPos, yPos, zPos };
}

uint MultipleGridBuilderFacade::getIndex1D(uint xIndex, uint yIndex, uint zIndex) const
{
    return xIndex + yIndex * numberOfSubdomains[Axis::x] + zIndex * numberOfSubdomains[Axis::x] * numberOfSubdomains[Axis::y];
}

uint MultipleGridBuilderFacade::getIndex1D(const std::array<uint, 3>& index3D) const
{
    return getIndex1D(index3D[Axis::x], index3D[Axis::y], index3D[Axis::z]);
}

void MultipleGridBuilderFacade::setSlipBoundaryCondition(SideType sideType, real normalX, real normalY, real normalZ) const
{
    setBoundaryCondition(sideType, [&]() { gridBuilder->setSlipBoundaryCondition(sideType, normalX, normalY, normalZ); });
}

void MultipleGridBuilderFacade::setStressBoundaryCondition(SideType sideType, real normalX, real normalY, real normalZ,
                                                           uint samplingOffset, real z0, real dx) const
{
    setBoundaryCondition(sideType, [&]() {
        gridBuilder->setStressBoundaryCondition(sideType, normalX, normalY, normalZ, samplingOffset, z0, dx);
    });
}

void MultipleGridBuilderFacade::setPrecursorBoundaryCondition(SideType sideType, SPtr<FileCollection> fileCollection,
                                                              int timeStepsBetweenReads, real velocityX, real velocityY,
                                                              real velocityZ,
                                                              std::vector<uint> fileLevelToGridLevelMap) const
{
    setBoundaryCondition(sideType, [&]() {
        gridBuilder->setPrecursorBoundaryCondition(sideType, fileCollection, timeStepsBetweenReads, velocityX, velocityY,
                                                   velocityZ, fileLevelToGridLevelMap);
    });
}

void MultipleGridBuilderFacade::setVelocityBoundaryCondition(SideType sideType, real vx, real vy, real vz) const
{
    setBoundaryCondition(sideType, [&]() { gridBuilder->setVelocityBoundaryCondition(sideType, vx, vy, vz); });
}

void MultipleGridBuilderFacade::setPressureBoundaryCondition(SideType sideType, real rho) const
{
    setBoundaryCondition(sideType, [&]() { gridBuilder->setPressureBoundaryCondition(sideType, rho); });
}

void MultipleGridBuilderFacade::setNoSlipBoundaryCondition(SideType sideType) const
{
    setBoundaryCondition(sideType, [&]() { gridBuilder->setNoSlipBoundaryCondition(sideType); });
}

bool MultipleGridBuilderFacade::isFinalSubdomainInDirection(CommunicationDirection direction) const
{
    return !hasNeighbors.at(direction);
}

uint MultipleGridBuilderFacade::getIndexOfFinalSubdomainInDirection(CommunicationDirection direction) const
{
    std::array<uint, 3> resultIndex3D = index;
    const Axis axis = communicationDirectionToAxes.at(direction);

    if (isNegative(direction)) {
        resultIndex3D[axis] = 0; // first subdomain index in direction
        return getIndex1D(resultIndex3D);
    }
    if (isPositive(direction)) {
        resultIndex3D[axis] = numberOfSubdomains[axis] - 1; // last subdomain index in direction
        return getIndex1D(resultIndex3D);
    }
    return UINT_MAX;
}

void MultipleGridBuilderFacade::setPeriodicBoundaryCondition(bool periodic_X, bool periodic_Y, bool periodic_Z) const
{
    setPeriodicBoundaryCondition({ periodic_X, periodic_Y, periodic_Z });
}

void MultipleGridBuilderFacade::setPeriodicBoundaryCondition(const std::array<bool, 3>& periodicity) const
{

    std::array<bool, 3> localPeriodicity = { false, false, false };

    for (const auto coordDirection : axis::allAxes) {
        if (!periodicity[coordDirection]) {
            continue;
        }

        // only one grid in direction --> set local periodicity
        if (numberOfSubdomains[coordDirection] == 1) {
            localPeriodicity[coordDirection] = true;
            continue;
        }

        // non-local periodicity --> set communication neighbors
        const CommunicationDirection negativeDirection = getNegativeDirectionAlongAxis(coordDirection);
        const CommunicationDirection positiveDirection = getPositiveDirectionAlongAxis(coordDirection);

        if (isFinalSubdomainInDirection(negativeDirection)) {
            // set final grid in positive direction as communication neighbor
            gridBuilder->findCommunicationIndices(negativeDirection);
            gridBuilder->setCommunicationProcess(negativeDirection, getIndexOfFinalSubdomainInDirection(positiveDirection));
        } else if (isFinalSubdomainInDirection(positiveDirection)) {
            // set final grid in negative direction as communication neighbor
            gridBuilder->findCommunicationIndices(positiveDirection);
            gridBuilder->setCommunicationProcess(positiveDirection, getIndexOfFinalSubdomainInDirection(negativeDirection));
        }
    }

    gridBuilder->setPeriodicBoundaryCondition(localPeriodicity[Axis::x], localPeriodicity[Axis::y],
                                              localPeriodicity[Axis::z]);
}

SPtr<MultipleGridBuilder> MultipleGridBuilderFacade::getGridBuilder() const
{
    return gridBuilder;
}