#include "GridInterface.cuh"

#include "Grid.cuh"
#include "NodeValues.h"
#include "utilities/math/CudaMath.cuh"


GridInterface::GridInterface(const Grid* fineGrid)
{
    const uint sizeCF = fineGrid->nx * fineGrid->ny + fineGrid->ny * fineGrid->nz + fineGrid->nx * fineGrid->nz;
    initalCoarseToFine(sizeCF, fineGrid);
    initalFineToCoarse(sizeCF, fineGrid);
}

GridInterface::~GridInterface()
{
    delete[] cf.coarse;
    delete[] cf.fine;

    delete[] fc.coarse;
    delete[] fc.fine;
}

void GridInterface::initalCoarseToFine(const uint sizeCF, const Grid* fineGrid)
{
    cf.coarse = new uint[sizeCF];
    cf.fine = new uint[sizeCF];

    cf.startCoarseX = fineGrid->startX - fineGrid->delta * 0.5;
    cf.startCoarseY = fineGrid->startY - fineGrid->delta * 0.5;
    cf.startCoarseZ = fineGrid->startZ - fineGrid->delta * 0.5;

    cf.endCoarseX = fineGrid->endX - fineGrid->delta * 1.5;
    cf.endCoarseY = fineGrid->endY - fineGrid->delta * 1.5;
    cf.endCoarseZ = fineGrid->endZ - fineGrid->delta * 1.5;

    cf.coarseEntry = CFC;
    cf.fineEntry = CFF;
}

void GridInterface::initalFineToCoarse(const uint sizeCF, const Grid* fineGrid)
{
    fc.coarse = new uint[sizeCF];
    fc.fine = new uint[sizeCF];

    fc.startCoarseX = cf.startCoarseX + 4 * fineGrid->delta;
    fc.startCoarseY = cf.startCoarseY + 4 * fineGrid->delta;
    fc.startCoarseZ = cf.startCoarseZ + 4 * fineGrid->delta;

    fc.endCoarseX = cf.endCoarseX - 2 * fineGrid->delta;
    fc.endCoarseY = cf.endCoarseY - 2 * fineGrid->delta;
    fc.endCoarseZ = cf.endCoarseZ - 2 * fineGrid->delta;

    fc.coarseEntry = FCC;
    fc.fineEntry = FCF;
}

void GridInterface::findCF(const uint& index, const Grid* coarseGrid, const Grid* fineGrid)
{
    findInterface(cf, 1, index, coarseGrid, fineGrid);
}

void GridInterface::findFC(const uint& index, const Grid* coarseGrid, const Grid* fineGrid)
{
    findInterface(fc, -1, index, coarseGrid, fineGrid);
}

void GridInterface::findInterface(Interface& interface, const int& factor, const uint& index, const Grid* coarseGrid, const Grid* fineGrid)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(index, x, y, z);

    if (isOnInterface(interface, x, y, z))
    {
        const uint indexFinerGrid = getIndexOnFinerGrid(factor, fineGrid, x, y, z);

        interface.coarse[interface.numberOfEntries] = index;
        interface.fine[interface.numberOfEntries] = indexFinerGrid;

        interface.numberOfEntries++;

        coarseGrid->field[index] = interface.coarseEntry;
        fineGrid->field[indexFinerGrid] = interface.fineEntry;
    }
}

HOSTDEVICE bool isOn(const real coord, const real plane1, const real plane2)
{
    return  CudaMath::equal(coord, plane1) || CudaMath::equal(coord, plane2);
}

HOSTDEVICE bool isBetween(const real coord, const real start, const real end)
{
    return  CudaMath::greaterEqual(coord, start) && CudaMath::lessEqual(coord, end);
}

bool GridInterface::isOnInterface(Interface& interface, const real& x, const real& y, const real& z)
{
    const bool isOnXYPlanes = isOn(z, interface.startCoarseZ, interface.endCoarseZ) && isBetween(y, interface.startCoarseY, interface.endCoarseY) && isBetween(x, interface.startCoarseX, interface.endCoarseX);
    const bool isOnXZPlanes = isOn(y, interface.startCoarseY, interface.endCoarseY) && isBetween(x, interface.startCoarseX, interface.endCoarseX) && isBetween(z, interface.startCoarseZ, interface.endCoarseZ);
    const bool isOnYZPlanes = isOn(x, interface.startCoarseX, interface.endCoarseX) && isBetween(y, interface.startCoarseY, interface.endCoarseY) && isBetween(z, interface.startCoarseZ, interface.endCoarseZ);

    return isOnXYPlanes || isOnXZPlanes || isOnYZPlanes;
}

uint GridInterface::getIndexOnFinerGrid(const int& factor, const Grid* fineGrid, const real& x, const real& y, const real& z)
{
    const real xFine = x + factor * (fineGrid->delta * 0.5);
    const real yFine = y + factor * (fineGrid->delta * 0.5);
    const real zFine = z + factor * (fineGrid->delta * 0.5);

    return fineGrid->transCoordToIndex(xFine, yFine, zFine);
}
