#include "GridInterface.cuh"

#include "GridImp.cuh"
#include "NodeValues.h"
#include "utilities/math/CudaMath.cuh"


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
    findInterface(cf, 1, index, coarseGrid, fineGrid);
}

void GridInterface::findFC(const uint& index, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    findInterface(fc, -1, index, coarseGrid, fineGrid);
}

void GridInterface::findInterface(Interface& interface, const int& factor, const uint& index, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(index, x, y, z);

    if (fineGrid->isOnInterface(x, y, z, interface.startOffset, interface.endOffset))
    {
        const uint indexFinerGrid = getIndexOnFinerGrid(factor, fineGrid, x, y, z);

        interface.coarse[interface.numberOfEntries] = index;
        interface.fine[interface.numberOfEntries] = indexFinerGrid;

        interface.numberOfEntries++;

        coarseGrid->field[index] = interface.coarseEntry;
        fineGrid->field[indexFinerGrid] = interface.fineEntry;
    }
}

HOSTDEVICE uint GridInterface::getIndexOnFinerGrid(const real& factor, const GridImp* fineGrid, const real& x, const real& y, const real& z)
{
    const real xFine = x + factor * (fineGrid->delta * 0.5);
    const real yFine = y + factor * (fineGrid->delta * 0.5);
    const real zFine = z + factor * (fineGrid->delta * 0.5);

    return fineGrid->matrixIndex[fineGrid->transCoordToIndex(xFine, yFine, zFine)];
}

void GridInterface::print() const
{
    printf("offset cf: (%2.2f, %2.2f); offset fc: (%2.2f, %2.2f); ", cf.startOffset, cf.endOffset, cf.startOffset, cf.endOffset);
    printf("Grid Interface - CF nodes: %d, FC nodes: %d\n", cf.numberOfEntries, fc.numberOfEntries);
}
