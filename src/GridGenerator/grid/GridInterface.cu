#include "GridInterface.cuh"

#include "Grid.cuh"
#include "utilities/math/CudaMath.cuh"

GridInterface::GridInterface(const Grid* fineGrid)
{
    uint sizeCF = fineGrid->nx * fineGrid->ny + fineGrid->ny * fineGrid->nz + fineGrid->nx * fineGrid->nz;
    cfc = new uint[sizeCF];
    cff = new uint[sizeCF];

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

void GridInterface::findCF(uint index, const Grid* coarseGrid, const Grid* fineGrid)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(index, x, y, z);

    bool isOnXYPlanes = 
        (CudaMath::equal(z, startCFCz) || CudaMath::equal(z, endCFCz)) &&
        CudaMath::lessEqual(x, endCFCx) &&
        CudaMath::greaterEqual(x, startCFCx) &&
        CudaMath::lessEqual(y, endCFCy) &&
        CudaMath::greaterEqual(y, startCFCy);

    bool isOnXZPlanes =
        (CudaMath::equal(y, startCFCy) || CudaMath::equal(y, endCFCy)) &&
        CudaMath::lessEqual(x, endCFCx) &&
        CudaMath::greaterEqual(x, startCFCx) &&
        CudaMath::lessEqual(z, endCFCz) &&
        CudaMath::greaterEqual(z, startCFCz);

    bool isOnYZPlanes =
        (CudaMath::equal(x, startCFCx) || CudaMath::equal(x, endCFCx)) &&
        CudaMath::lessEqual(z, endCFCz) &&
        CudaMath::greaterEqual(z, startCFCz) &&
        CudaMath::lessEqual(y, endCFCy) &&
        CudaMath::greaterEqual(y, startCFCy);

    const bool isCF = isOnXYPlanes || isOnXZPlanes || isOnYZPlanes;
    if (isCF)
    {
        coarseGrid->field[index] = 77;

        const real xFine = x + fineGrid->delta * 0.5;
        const real yFine = y + fineGrid->delta * 0.5;
        const real zFine = z + fineGrid->delta * 0.5;
        const uint indexFinerGrid = fineGrid->transCoordToIndex(xFine, yFine, zFine);
        cfc[numberOfEntriesInCF++] = index;
        cff[numberOfEntriesInCF++] = indexFinerGrid;
    }
}
