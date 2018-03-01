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
    cf.coarseEntry = FLUID_CFC;
    cf.fineEntry = FLUID_CFF;

    cf.startOffset = -0.5 * fineGrid->delta;
    cf.endOffset = -1.5 * fineGrid->delta;
}

void GridInterface::initalFineToCoarse(const GridImp* fineGrid)
{
    fc.coarseEntry = FLUID_FCC;
    fc.fineEntry = FLUID_FCF;

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

int GridInterface::getCoarseToFineIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);
    const real xFine = x + (fineGrid->delta * 0.5);
    const real yFine = y + (fineGrid->delta * 0.5);
    const real zFine = z + (fineGrid->delta * 0.5);

    return fineGrid->transCoordToIndex(xFine, yFine, zFine);
}



bool GridInterface::isNeighborFineFluid(real x, real y, real z, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    const int neighbor = coarseGrid->transCoordToIndex(x, y, z);
    const int indexOnFineGrid = getCoarseToFineIndexOnFineGrid(neighbor, coarseGrid, fineGrid);
    if (indexOnFineGrid == -1)
        return false;
    return fineGrid->isFluid(indexOnFineGrid);
}

void GridInterface::markCellTo(const GridImp* coarseGrid, uint index, char type)
{
    coarseGrid->field[index] = type;

    real x, y, z;
    coarseGrid->transIndexToCoords(index, x, y, z);

    const int neighborX = coarseGrid->transCoordToIndex(x + coarseGrid->getDelta(), y, z);
    const int neighborY = coarseGrid->transCoordToIndex(x, y + coarseGrid->getDelta(), z);
    const int neighborZ = coarseGrid->transCoordToIndex(x , y, z + coarseGrid->getDelta());

    coarseGrid->field[neighborX] = type;
    coarseGrid->field[neighborY] = type;
    coarseGrid->field[neighborZ] = type;

    const int neighborYZ = coarseGrid->transCoordToIndex(x, y + coarseGrid->getDelta(), z + coarseGrid->getDelta());
    const int neighborXZ = coarseGrid->transCoordToIndex(x + coarseGrid->getDelta(), y , z + coarseGrid->getDelta());
    const int neighborXY = coarseGrid->transCoordToIndex(x + coarseGrid->getDelta(), y + coarseGrid->getDelta(), z);

    coarseGrid->field[neighborYZ] = type;
    coarseGrid->field[neighborXZ] = type;
    coarseGrid->field[neighborXY] = type;

    const int neighborXYZ = coarseGrid->transCoordToIndex(x + coarseGrid->getDelta(), y + coarseGrid->getDelta(), z + coarseGrid->getDelta());

    coarseGrid->field[neighborXYZ] = type;
}

void GridInterface::findInterfaceCF(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsFluid = coarseGrid->isFluid(indexOnCoarseGrid);
    if (!nodeOnCoarseGridIsFluid)
        return;

    const int indexOnFineGridCF = getCoarseToFineIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridCF == -1)
        return;

    const bool fineGridNodeIsFluid = fineGrid->isFluid(indexOnFineGridCF);
    if (!fineGridNodeIsFluid)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    const bool fineGridNodeIsFluidX = isNeighborFineFluid(x + coarseGrid->delta, y, z, coarseGrid, fineGrid);
    const bool fineGridNodeIsFluidX2 = isNeighborFineFluid(x - coarseGrid->delta, y, z, coarseGrid, fineGrid);
    const bool fineGridNodeIsFluidY = isNeighborFineFluid(x, y + coarseGrid->delta, z, coarseGrid, fineGrid);
    const bool fineGridNodeIsFluidY2 = isNeighborFineFluid(x, y - coarseGrid->delta, z, coarseGrid, fineGrid);
    const bool fineGridNodeIsFluidZ = isNeighborFineFluid(x, y, z + coarseGrid->delta, coarseGrid, fineGrid);
    const bool fineGridNodeIsFluidZ2 = isNeighborFineFluid(x, y, z - coarseGrid->delta, coarseGrid, fineGrid);

   if (!fineGridNodeIsFluidX || !fineGridNodeIsFluidX2 || !fineGridNodeIsFluidY || !fineGridNodeIsFluidY2 || !fineGridNodeIsFluidZ || !fineGridNodeIsFluidZ2) {
        cf.coarse[cf.numberOfEntries] = indexOnCoarseGrid;
        cf.fine[cf.numberOfEntries] = indexOnFineGridCF;

        cf.numberOfEntries++;

        this->markCellTo(coarseGrid, indexOnCoarseGrid, cf.coarseEntry);

        //fineGrid->field[indexOnFineGridCF] = cf.fineEntry;
    }
}

HOSTDEVICE void GridInterface::findInterfaceFC(const uint& indexOnCoarseGrid, GridImp* coarseGrid, const GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsFluid = coarseGrid->isFluid(indexOnCoarseGrid);
    const bool nodeOnCoarseGridIsCoarseToFine = coarseGrid->isCoarseToFineNode(indexOnCoarseGrid);
    if (!nodeOnCoarseGridIsFluid || nodeOnCoarseGridIsCoarseToFine)
        return;

    const int indexOnFineGridFC = getFineToCoarseIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridFC == -1)
        return;

    const bool fineGridNodeIsFluid = fineGrid->isFluid(indexOnFineGridFC);
    if (!fineGridNodeIsFluid)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    const bool neighborBelongsToCoarseToFineInterpolationCellX  = belongsNeighborToCoarseToFineInterpolationCell(x + coarseGrid->delta, y, z, coarseGrid, fineGrid);
    const bool neighborBelongsToCoarseToFineInterpolationCellX2 = belongsNeighborToCoarseToFineInterpolationCell(x - coarseGrid->delta, y, z, coarseGrid, fineGrid);
    const bool neighborBelongsToCoarseToFineInterpolationCellY  = belongsNeighborToCoarseToFineInterpolationCell(x, y + coarseGrid->delta, z, coarseGrid, fineGrid);
    const bool neighborBelongsToCoarseToFineInterpolationCellY2 = belongsNeighborToCoarseToFineInterpolationCell(x, y - coarseGrid->delta, z, coarseGrid, fineGrid);
    const bool neighborBelongsToCoarseToFineInterpolationCellZ  = belongsNeighborToCoarseToFineInterpolationCell(x, y, z + coarseGrid->delta, coarseGrid, fineGrid);
    const bool neighborBelongsToCoarseToFineInterpolationCellZ2 = belongsNeighborToCoarseToFineInterpolationCell(x, y, z - coarseGrid->delta, coarseGrid, fineGrid);

    const bool secondNeighborBelongsToCoarseToFineInterpolationCellX2 = belongsNeighborToCoarseToFineInterpolationCell(x - 2 * coarseGrid->delta, y, z, coarseGrid, fineGrid);
    const bool secondNeighborBelongsToCoarseToFineInterpolationCellY2 = belongsNeighborToCoarseToFineInterpolationCell(x, y - 2 * coarseGrid->delta, z, coarseGrid, fineGrid);
    const bool secondNeighborBelongsToCoarseToFineInterpolationCellZ2 = belongsNeighborToCoarseToFineInterpolationCell(x, y, z - 2 * coarseGrid->delta, coarseGrid, fineGrid);

    if (neighborBelongsToCoarseToFineInterpolationCellX ||
        neighborBelongsToCoarseToFineInterpolationCellX2 ||
        neighborBelongsToCoarseToFineInterpolationCellY ||
        neighborBelongsToCoarseToFineInterpolationCellY2 ||
        neighborBelongsToCoarseToFineInterpolationCellZ ||
        neighborBelongsToCoarseToFineInterpolationCellZ2) 
    {
        fc.coarse[fc.numberOfEntries] = indexOnCoarseGrid;
        fc.fine[fc.numberOfEntries] = indexOnFineGridFC;

        fc.numberOfEntries++;

        //this->markCellTo(coarseGrid, indexOnCoarseGrid, cf.coarseEntry);

        //fineGrid->field[indexOnFineGridCF] = cf.fineEntry;
    } else if( // isStopper
        secondNeighborBelongsToCoarseToFineInterpolationCellX2 ||
        secondNeighborBelongsToCoarseToFineInterpolationCellY2 ||
        secondNeighborBelongsToCoarseToFineInterpolationCellZ2)
    {
        coarseGrid->setFieldEntryToStopperOverlapGrid(indexOnCoarseGrid);
    } 
    else //should be inside of fine grid and can be deleted
    {
        coarseGrid->setFieldEntryToInvalid(indexOnCoarseGrid);
    }
}

int GridInterface::getFineToCoarseIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);
    const real xFine = x - (fineGrid->delta * 0.5);
    const real yFine = y - (fineGrid->delta * 0.5);
    const real zFine = z - (fineGrid->delta * 0.5);

    return fineGrid->transCoordToIndex(xFine, yFine, zFine);
}

bool GridInterface::belongsNeighborToCoarseToFineInterpolationCell(real x, real y, real z, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    const int neighborIndex = coarseGrid->transCoordToIndex(x, y, z);
    return coarseGrid->isCoarseToFineNode(neighborIndex);
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
