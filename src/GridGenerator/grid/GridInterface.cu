#include "GridInterface.cuh"

#include "GridImp.cuh"
#include "NodeValues.h"
#include "utilities/math/CudaMath.cuh"
#include <iostream>

GridInterface::GridInterface()
{
}

GridInterface::~GridInterface()
{

}

void GridInterface::initalGridInterface(const GridImp* fineGrid)
{
    initalCoarseToFine(fineGrid);
    initalFineToCoarse(fineGrid);
}

void GridInterface::initalCoarseToFine(const GridImp* fineGrid)
{
    cf.coarseEntry = CFC;
    cf.fineEntry = CFF;

    cf.startOffset = -0.5 * fineGrid->delta;
    cf.endOffset = -1.5 * fineGrid->delta;
}

void GridInterface::initalFineToCoarse(const GridImp* fineGrid)
{
    fc.coarseEntry = FCC;
    fc.fineEntry = FCF;

    fc.startOffset = cf.startOffset + 4 * fineGrid->delta;
    fc.endOffset   = cf.endOffset - 2 * fineGrid->delta;
}




void GridInterface::findCF(const uint& index, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    //findInterface(cf, 1, index, coarseGrid, fineGrid);
}

void GridInterface::findFC(const uint& index, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    //findInterface(fc, -1, index, coarseGrid, fineGrid);
}

int GridInterface::getIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);
    const real xFine = x + (fineGrid->delta * 0.5);
    const real yFine = y + (fineGrid->delta * 0.5);
    const real zFine = z + (fineGrid->delta * 0.5);

    return fineGrid->transCoordToIndex(xFine, yFine, zFine);
}

void GridInterface::findInterface(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsFluid = coarseGrid->isFluid(indexOnCoarseGrid);
    if (!nodeOnCoarseGridIsFluid)
        return;

    const int indexOnFineGridCF = getIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridCF == -1)
        return;
    //const int matrixIndexOnFineGrid = fineGrid->matrixIndex[indexOnFineGrid];
    //if (matrixIndexOnFineGrid == -1)
    //    return;

    const bool fineGridNodeIsFluid = fineGrid->isFluid(indexOnFineGridCF);
    if (!fineGridNodeIsFluid)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    bool fineGridNodeIsFluidX = false;
    bool fineGridNodeIsFluidX2 = false;
    bool fineGridNodeIsFluidY = false;
    bool fineGridNodeIsFluidY2 = false;
    bool fineGridNodeIsFluidZ = false;
    bool fineGridNodeIsFluidZ2 = false;

    const int xNeighbor = coarseGrid->transCoordToIndex(x + coarseGrid->delta, y, z);
    int indexOnFineGrid = getIndexOnFineGrid(xNeighbor, coarseGrid, fineGrid);
    if (indexOnFineGrid != -1) {
            fineGridNodeIsFluidX = fineGrid->isFluid(indexOnFineGrid);
    }

    const int xNeighbor2 = coarseGrid->transCoordToIndex(x - coarseGrid->delta, y, z);
    indexOnFineGrid = getIndexOnFineGrid(xNeighbor2, coarseGrid, fineGrid);
    if (indexOnFineGrid != -1) {
            fineGridNodeIsFluidX2 = fineGrid->isFluid(indexOnFineGrid);
    }

    const int yNeighbor = coarseGrid->transCoordToIndex(x, y + coarseGrid->delta, z);
    indexOnFineGrid = getIndexOnFineGrid(yNeighbor, coarseGrid, fineGrid);
    if (indexOnFineGrid != -1) {
            fineGridNodeIsFluidY = fineGrid->isFluid(indexOnFineGrid);
    }

    const int yNeighbor2 = coarseGrid->transCoordToIndex(x, y - coarseGrid->delta, z);
    indexOnFineGrid = getIndexOnFineGrid(yNeighbor2, coarseGrid, fineGrid);
    if (indexOnFineGrid != -1) {
            fineGridNodeIsFluidY2 = fineGrid->isFluid(indexOnFineGrid);
    }

    const int zNeighbor = coarseGrid->transCoordToIndex(x, y, z + coarseGrid->delta);
    indexOnFineGrid = getIndexOnFineGrid(zNeighbor, coarseGrid, fineGrid);
    if (indexOnFineGrid != -1) {
            fineGridNodeIsFluidZ = fineGrid->isFluid(indexOnFineGrid);
    }

    const int zNeighbor2 = coarseGrid->transCoordToIndex(x, y, z - coarseGrid->delta);
    indexOnFineGrid = getIndexOnFineGrid(zNeighbor2, coarseGrid, fineGrid);
    if (indexOnFineGrid != -1) {
            fineGridNodeIsFluidZ2 = fineGrid->isFluid(indexOnFineGrid);
    }

   if (!fineGridNodeIsFluidX || !fineGridNodeIsFluidX2 || !fineGridNodeIsFluidY || !fineGridNodeIsFluidY2 || !fineGridNodeIsFluidZ || !fineGridNodeIsFluidZ2) {
        cf.coarse[cf.numberOfEntries] = indexOnCoarseGrid;
        cf.fine[cf.numberOfEntries] = fineGrid->matrixIndex[indexOnFineGridCF];

        cf.numberOfEntries++;

        //coarseGrid->field[indexOnCoarseGrid] = cf.coarseEntry;
        //fineGrid->field[indexOnFineGridCF] = cf.fineEntry;
    }
}


//void GridInterface::findInterface(Interface& gridInterface, const int& factor, const uint& index, const GridImp* coarseGrid, const GridImp* fineGrid)
//{
//    real x, y, z;
//    coarseGrid->transIndexToCoords(index, x, y, z);
//
//    if (fineGrid->isOnInterface(x, y, z, gridInterface.startOffset, gridInterface.endOffset))
//    {
//        const uint indexFinerGrid = getIndexOnFinerGrid(factor, fineGrid, x, y, z);
//
//        gridInterface.coarse[gridInterface.numberOfEntries] = index;
//        gridInterface.fine[gridInterface.numberOfEntries] = indexFinerGrid;
//
//        gridInterface.numberOfEntries++;
//
//        coarseGrid->field[index] = gridInterface.coarseEntry;
//        fineGrid->field[indexFinerGrid] = gridInterface.fineEntry;
//    }
//}

//HOSTDEVICE uint GridInterface::getIndexOnFinerGrid(const real& factor, const GridImp* fineGrid, const real& x, const real& y, const real& z)
//{
//    const real xFine = x + factor * (fineGrid->delta * 0.5);
//    const real yFine = y + factor * (fineGrid->delta * 0.5);
//    const real zFine = z + factor * (fineGrid->delta * 0.5);
//
//    return fineGrid->matrixIndex[fineGrid->transCoordToIndex(xFine, yFine, zFine)];
//}

void GridInterface::print() const
{
    printf("offset cf: (%2.2f, %2.2f); offset fc: (%2.2f, %2.2f); ", cf.startOffset, cf.endOffset, cf.startOffset, cf.endOffset);
    printf("Grid Interface - CF nodes: %d, FC nodes: %d\n", cf.numberOfEntries, fc.numberOfEntries);
}
