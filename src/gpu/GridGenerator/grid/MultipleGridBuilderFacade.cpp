#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>

#include "MultipleGridBuilderFacade.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"
#include "grid/GridDimensions.h"
#include "geometries/BoundingBox/BoundingBox.h"
#include "geometries/Object.h"

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
    this->numberOfGridsTotal = this->numberGridsX * this->numberGridsY * this->numberGridsZ;

    if (numberOfGridsTotal > 1 && !this->overlapOfSubdomains)
        throw std::runtime_error("OverlapOfSubdomains in MultipleGridBuilderFacade is NaN.");

    if (generatePart >= numberOfGridsTotal)
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
    this->numberGridsX = (uint)(this->xSplits.size() + 1);
    this->numberGridsY = (uint)(this->ySplits.size() + 1);
    this->numberGridsZ = (uint)(this->zSplits.size() + 1);
}

void MultipleGridBuilderFacade::sortSplitLocations()
{

    std::sort(this->xSplits.begin(), this->xSplits.end());
    std::sort(this->ySplits.begin(), this->ySplits.end());
    std::sort(this->zSplits.begin(), this->zSplits.end());
}

void MultipleGridBuilderFacade::calculatedIndexOfPart(uint generatePart)
{
    this->xIndex = this->getX3D(generatePart);
    this->yIndex = this->getY3D(generatePart);
    this->zIndex = this->getZ3D(generatePart);
}

void MultipleGridBuilderFacade::checkSplitLocations(const std::vector<real> &splits, real lowerBound, real upperBound) const
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

    real xMinCoarseGrid = xValues[xIndex];
    real yMinCoarseGrid = yValues[yIndex];
    real zMinCoarseGrid = zValues[zIndex];
    real xMaxCoarseGrid = xValues[xIndex + 1];
    real yMaxCoarseGrid = yValues[yIndex + 1];
    real zMaxCoarseGrid = zValues[zIndex + 1];

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
    if ((numberGridsX * numberGridsY * numberGridsZ) > 1) {
        gridBuilder->setSubDomainBox(std::make_shared<BoundingBox>(xValues[xIndex], xValues[xIndex + 1],
                                                                   yValues[yIndex], yValues[yIndex + 1],
                                                                   zValues[zIndex], zValues[zIndex + 1]));
    }
}

void MultipleGridBuilderFacade::setUpCommunicationNeighbors()
{
    // Communication is only needed on multiple gpus
    if (numberOfGridsTotal == 1) return;

    if (hasNeighbors.empty())
        throw std::runtime_error("checkForNeighbors() has to be called befor calling setUpCommunicationNeighbors()");

    for (const auto &[direction, hasNeighbor] : hasNeighbors) {
        if (hasNeighbor) {
            gridBuilder->findCommunicationIndices(direction);

            switch (direction) {
                case CommunicationDirections::MX:
                    gridBuilder->setCommunicationProcess(direction, getIndex1D(xIndex - 1, yIndex, zIndex));
                    break;
                case CommunicationDirections::MY:
                    gridBuilder->setCommunicationProcess(direction, getIndex1D(xIndex, yIndex - 1, zIndex));
                    break;
                case CommunicationDirections::MZ:
                    gridBuilder->setCommunicationProcess(direction, getIndex1D(xIndex, yIndex, zIndex - 1));
                    break;
                case CommunicationDirections::PX:
                    gridBuilder->setCommunicationProcess(direction, getIndex1D(xIndex + 1, yIndex, zIndex));
                    break;
                case CommunicationDirections::PY:
                    gridBuilder->setCommunicationProcess(direction, getIndex1D(xIndex, yIndex + 1, zIndex));
                    break;
                case CommunicationDirections::PZ:
                    gridBuilder->setCommunicationProcess(direction, getIndex1D(xIndex, yIndex, zIndex + 1));
                    break;
            }
        }
    }
}

void MultipleGridBuilderFacade::checkForNeighbors()
{
    hasNeighbors[CommunicationDirections::MX] = (xIndex > 0);
    hasNeighbors[CommunicationDirections::MY] = (yIndex > 0);
    hasNeighbors[CommunicationDirections::MZ] = (zIndex > 0);
    hasNeighbors[CommunicationDirections::PX] = (xIndex < numberGridsX - 1);
    hasNeighbors[CommunicationDirections::PY] = (yIndex < numberGridsY - 1);
    hasNeighbors[CommunicationDirections::PZ] = (zIndex < numberGridsZ - 1);
}

void MultipleGridBuilderFacade::addFineGridsToGridBuilder()
{
    for (auto const &grid : fineGrids) {
        gridBuilder->addGrid(grid.first, grid.second);
    }
}

void MultipleGridBuilderFacade::addGeometriesToGridBuilder()
{
    for (auto const &geometry : geometries) {
        gridBuilder->addGeometry(geometry);
    }
}

void MultipleGridBuilderFacade::setOverlapOfSubdomains(real overlap)
{
    this->overlapOfSubdomains = overlap;
}

void MultipleGridBuilderFacade::addDomainSplit(real coordinate, MultipleGridBuilderFacade::CoordDirection direction)
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
    return index1D % numberGridsX;
}

uint MultipleGridBuilderFacade::getY3D(uint index1D) const
{
    return (index1D / numberGridsX) % numberGridsY;
}

uint MultipleGridBuilderFacade::getZ3D(uint index1D) const
{
    return index1D / (numberGridsX * numberGridsY);
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
    return xIndex + yIndex * numberGridsX + zIndex * numberGridsX * numberGridsY;
}

void MultipleGridBuilderFacade::setSlipBoundaryCondition(SideType sideType, real normalX, real normalY, real normalZ)
{
    setBoundaryCondition(sideType, [&]() { gridBuilder->setSlipBoundaryCondition(sideType, normalX, normalY, normalZ); });
}

void MultipleGridBuilderFacade::setStressBoundaryCondition(SideType sideType, real normalX, real normalY, real normalZ,
                                                           uint samplingOffset, real z0, real dx)
{
    setBoundaryCondition(sideType, [&]() {
        gridBuilder->setStressBoundaryCondition(sideType, normalX, normalY, normalZ, samplingOffset, z0, dx);
    });
}

void MultipleGridBuilderFacade::setPrecursorBoundaryCondition(SideType sideType, SPtr<FileCollection> fileCollection,
                                                              int timeStepsBetweenReads, real velocityX, real velocityY,
                                                              real velocityZ, std::vector<uint> fileLevelToGridLevelMap)
{
    setBoundaryCondition(sideType, [&]() {
        gridBuilder->setPrecursorBoundaryCondition(sideType, fileCollection, timeStepsBetweenReads, velocityX, velocityY,
                                                   velocityZ, fileLevelToGridLevelMap);
    });
}

void MultipleGridBuilderFacade::setVelocityBoundaryCondition(SideType sideType, real vx, real vy, real vz)
{
    setBoundaryCondition(sideType, [&]() { gridBuilder->setVelocityBoundaryCondition(sideType, vx, vy, vz); });
}

void MultipleGridBuilderFacade::setPressureBoundaryCondition(SideType sideType, real rho)
{
    setBoundaryCondition(sideType, [&]() { gridBuilder->setPressureBoundaryCondition(sideType, rho); });
}

void MultipleGridBuilderFacade::setNoSlipBoundaryCondition(SideType sideType)
{
    setBoundaryCondition(sideType, [&]() { gridBuilder->setNoSlipBoundaryCondition(sideType); });
}

void MultipleGridBuilderFacade::setPeriodicBoundaryCondition(bool periodic_X, bool periodic_Y, bool periodic_Z)
{
    bool localPeriodicityX = false;
    bool localPeriodicityY = false;
    bool localPeriodicityZ = false;

    if (periodic_X) {
        if (numberGridsX == 1) {
            localPeriodicityX = true;
        }
        if (numberGridsX > 1 && !hasNeighbors[CommunicationDirections::MX]) {
            // set last grid in x-direction as communication neighbor
            gridBuilder->findCommunicationIndices(CommunicationDirections::MX);
            gridBuilder->setCommunicationProcess(CommunicationDirections::MX, getIndex1D(numberGridsX - 1, yIndex, zIndex));
        } else if (numberGridsX > 1 && !hasNeighbors[CommunicationDirections::PX]) {
            // set first grid in x-direction as communication neighbor
            gridBuilder->findCommunicationIndices(CommunicationDirections::PX);
            gridBuilder->setCommunicationProcess(CommunicationDirections::PX, getIndex1D(0, yIndex, zIndex));
        }
    }

    if (periodic_Y) {
        if (numberGridsY == 1) {
            localPeriodicityY = true;
        }
        if (numberGridsY > 1 && !hasNeighbors[CommunicationDirections::MY]) {
            // set last grid in x-direction as communication neighbor
            gridBuilder->findCommunicationIndices(CommunicationDirections::MY);
            gridBuilder->setCommunicationProcess(CommunicationDirections::MY, getIndex1D(xIndex, numberGridsY - 1, zIndex));
        } else if (numberGridsY > 1 && !hasNeighbors[CommunicationDirections::PY]) {
            // set first grid in x-direction as communication neighbor
            gridBuilder->findCommunicationIndices(CommunicationDirections::PY);
            gridBuilder->setCommunicationProcess(CommunicationDirections::PY, getIndex1D(xIndex, 0, zIndex));
        }
    }

    if (periodic_Z) {
        if (numberGridsZ == 1) {
            localPeriodicityZ = true;
        }
        if (numberGridsZ > 1 && !hasNeighbors[CommunicationDirections::MZ]) {
            // set last grid in x-direction as communication neighbor
            gridBuilder->findCommunicationIndices(CommunicationDirections::MZ);
            gridBuilder->setCommunicationProcess(CommunicationDirections::MZ, getIndex1D(xIndex, yIndex, numberGridsZ - 1));
        } else if (numberGridsZ > 1 && !hasNeighbors[CommunicationDirections::PZ]) {
            // set first grid in x-direction as communication neighbor
            gridBuilder->findCommunicationIndices(CommunicationDirections::PZ);
            gridBuilder->setCommunicationProcess(CommunicationDirections::PZ, getIndex1D(xIndex, yIndex, 0));
        }
    }

    gridBuilder->setPeriodicBoundaryCondition(localPeriodicityX, localPeriodicityY, localPeriodicityZ);
}

SPtr<MultipleGridBuilder> MultipleGridBuilderFacade::getGridBuilder()
{
    return gridBuilder;
}