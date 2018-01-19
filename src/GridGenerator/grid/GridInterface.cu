#include "GridInterface.cuh"

#include "Grid.cuh"
#include "NodeValues.h"
#include "utilities/math/CudaMath.cuh"

GridInterface::GridInterface(const Grid* fineGrid)
{
    uint sizeCF = fineGrid->nx * fineGrid->ny + fineGrid->ny * fineGrid->nz + fineGrid->nx * fineGrid->nz;
    cfc = new uint[sizeCF];
    cff = new uint[sizeCF];

    fcf = new uint[sizeCF];
    fcc = new uint[sizeCF];

    startCFCx = fineGrid->startX - fineGrid->delta * 0.5;
    startCFCy = fineGrid->startY - fineGrid->delta * 0.5;
    startCFCz = fineGrid->startZ - fineGrid->delta * 0.5;

    endCFCx = fineGrid->endX - fineGrid->delta * 1.5;
    endCFCy = fineGrid->endY - fineGrid->delta * 1.5;
    endCFCz = fineGrid->endZ - fineGrid->delta * 1.5;
}

GridInterface::GridInterface()
{
    startCFCx = 0.0;
    startCFCy = 0.0;
    startCFCz = 0.0;

    endCFCx = 0.0;
    endCFCy = 0.0;
    endCFCz = 0.0;
}

void GridInterface::findCF(const uint index, const Grid* coarseGrid, const Grid* fineGrid)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(index, x, y, z);

    const bool isOnXYPlanes = 
        (CudaMath::equal(z, startCFCz) || CudaMath::equal(z, endCFCz)) &&
        CudaMath::lessEqual(x, endCFCx) &&
        CudaMath::greaterEqual(x, startCFCx) &&
        CudaMath::lessEqual(y, endCFCy) &&
        CudaMath::greaterEqual(y, startCFCy);

    const bool isOnXZPlanes =
        (CudaMath::equal(y, startCFCy) || CudaMath::equal(y, endCFCy)) &&
        CudaMath::lessEqual(x, endCFCx) &&
        CudaMath::greaterEqual(x, startCFCx) &&
        CudaMath::lessEqual(z, endCFCz) &&
        CudaMath::greaterEqual(z, startCFCz);

    const bool isOnYZPlanes =
        (CudaMath::equal(x, startCFCx) || CudaMath::equal(x, endCFCx)) &&
        CudaMath::lessEqual(z, endCFCz) &&
        CudaMath::greaterEqual(z, startCFCz) &&
        CudaMath::lessEqual(y, endCFCy) &&
        CudaMath::greaterEqual(y, startCFCy);

    const bool isCF = isOnXYPlanes || isOnXZPlanes || isOnYZPlanes;
    if (isCF)
    {
        const real xFine = x + fineGrid->delta * 0.5;
        const real yFine = y + fineGrid->delta * 0.5;
        const real zFine = z + fineGrid->delta * 0.5;
        const uint indexFinerGrid = fineGrid->transCoordToIndex(xFine, yFine, zFine);
        cfc[numberOfEntriesInCF++] = index;
        cff[numberOfEntriesInCF++] = indexFinerGrid;

        coarseGrid->field[index] = CFC;
        fineGrid->field[indexFinerGrid] = CFF;
    }
}

void GridInterface::findFC(uint index, const Grid* coarseGrid, const Grid* fineGrid)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(index, x, y, z);

    const real startFCCx = this->startCFCx + 2 * coarseGrid->delta;
    const real startFCCy = this->startCFCy + 2 * coarseGrid->delta;
    const real startFCCz = this->startCFCz + 2 * coarseGrid->delta;
    
    const real endFCCx = this->endCFCx - 2 * coarseGrid->delta;
    const real endFCCy = this->endCFCy - 2 * coarseGrid->delta;
    const real endFCCz = this->endCFCz - 2 * coarseGrid->delta;

    const bool isOnXYPlanes =
        (CudaMath::equal(z, startFCCz) || CudaMath::equal(z, endFCCz)) &&
        CudaMath::lessEqual(x, endFCCx) &&
        CudaMath::greaterEqual(x, startFCCx) &&
        CudaMath::lessEqual(y, endFCCy) &&
        CudaMath::greaterEqual(y, startFCCy);

    const bool isOnXZPlanes =
        (CudaMath::equal(y, startFCCy) || CudaMath::equal(y, endFCCy)) &&
        CudaMath::lessEqual(x, endFCCx) &&
        CudaMath::greaterEqual(x, startFCCx) &&
        CudaMath::lessEqual(z, endFCCz) &&
        CudaMath::greaterEqual(z, startFCCz);

    const bool isOnYZPlanes =
        (CudaMath::equal(x, startFCCx) || CudaMath::equal(x, endFCCx)) &&
        CudaMath::lessEqual(z, endFCCz) &&
        CudaMath::greaterEqual(z, startFCCz) &&
        CudaMath::lessEqual(y, endFCCy) &&
        CudaMath::greaterEqual(y, startFCCy);

    const bool isCF = isOnXYPlanes || isOnXZPlanes || isOnYZPlanes;
    if (isCF)
    {
        const real xFine = x - fineGrid->delta * 0.5;
        const real yFine = y - fineGrid->delta * 0.5;
        const real zFine = z - fineGrid->delta * 0.5;
        const uint indexFinerGrid = fineGrid->transCoordToIndex(xFine, yFine, zFine);
        fcc[numberOfEntriesInFC++] = index;
        fcf[numberOfEntriesInFC++] = indexFinerGrid;

        coarseGrid->field[index] = FCC;
        fineGrid->field[indexFinerGrid] = FCF;
    }
}
