#include "GridInterface.h"

#include "GridImp.h"
#include "Field.h"
#include "NodeValues.h"
#include "utilities/math/CudaMath.cuh"

GridInterface::GridInterface()
{

}

GridInterface::~GridInterface()
{

}

void GridInterface::findInterfaceCF(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsFluid = coarseGrid->getField().isFluid(indexOnCoarseGrid);
    if (!nodeOnCoarseGridIsFluid)
        return;

    const int indexOnFineGridCF = getCoarseToFineIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridCF == -1)
        return;

    const bool fineGridNodeIsFluid = fineGrid->getField().isFluid(indexOnFineGridCF);
    if (!fineGridNodeIsFluid)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    int dirX, dirY, dirZ;

    for(const auto dir : coarseGrid->distribution)
    {
        const bool isFineGridNeighborFluid = isNeighborFineFluid(x + dir[0] * coarseGrid->delta, y + dir[1] * coarseGrid->delta, z + dir[2] * coarseGrid->delta, coarseGrid, fineGrid);
        if(!isFineGridNeighborFluid)
        {
            cf.coarse[cf.numberOfEntries] = indexOnCoarseGrid;
            cf.fine[cf.numberOfEntries] = fineGrid->getSparseIndex(indexOnFineGridCF);

            cf.numberOfEntries++;

            coarseGrid->setCellTo(indexOnCoarseGrid, FLUID_CFC);
            fineGrid->setCellTo(indexOnFineGridCF, FLUID_CFF);
            break;
        }
    }
}



HOSTDEVICE void GridInterface::findInterfaceFC(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsFluid = coarseGrid->getField().isFluid(indexOnCoarseGrid);
    const bool nodeOnCoarseGridIsCoarseToFine = coarseGrid->getField().isCoarseToFineNode(indexOnCoarseGrid);
    if (!nodeOnCoarseGridIsFluid || nodeOnCoarseGridIsCoarseToFine)
        return;

    const int indexOnFineGridFC = getFineToCoarseIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridFC == -1)
        return;

    const bool fineGridNodeIsFluid = fineGrid->getField().isFluid(indexOnFineGridFC);
    if (!fineGridNodeIsFluid)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    for (const auto dir : coarseGrid->distribution)
    {
        const int neighborIndex = coarseGrid->transCoordToIndex(x + dir[0] * coarseGrid->delta, y + dir[1] * coarseGrid->delta, z + dir[2] * coarseGrid->delta);
        const bool neighborBelongsToCoarseToFineInterpolationCell = coarseGrid->getField().isCoarseToFineNode(neighborIndex);
        if (neighborBelongsToCoarseToFineInterpolationCell)
        {
            fc.coarse[fc.numberOfEntries] = indexOnCoarseGrid;
            fc.fine[fc.numberOfEntries] = fineGrid->getSparseIndex(indexOnFineGridFC);

            fc.numberOfEntries++;

            fineGrid->setCellTo(indexOnFineGridFC, FLUID_FCF);
            coarseGrid->getField().setFieldEntry(indexOnCoarseGrid, FLUID_FCC);
            break;
        }
    }
}

void GridInterface::findOverlapStopper(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid)
{
    const bool nodeOnCoarseGridIsFluid = coarseGrid->getField().isFluid(indexOnCoarseGrid);
    const bool nodeOnCoarseGridIsCoarseToFine = coarseGrid->getField().isCoarseToFineNode(indexOnCoarseGrid);
    const bool nodeOnCoarseGridIsFineToCoarse = coarseGrid->getField().isFineToCoarseNode(indexOnCoarseGrid);
    if (!nodeOnCoarseGridIsFluid || nodeOnCoarseGridIsCoarseToFine || nodeOnCoarseGridIsFineToCoarse)
        return;

    const int indexOnFineGridFC = getFineToCoarseIndexOnFineGrid(indexOnCoarseGrid, coarseGrid, fineGrid);
    if (indexOnFineGridFC == -1)
        return;

    const bool fineGridNodeIsFluid = fineGrid->getField().isFluid(indexOnFineGridFC);
    if (!fineGridNodeIsFluid)
        return;

    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    bool neighborBelongsToFineToCoarseInterpolationCell = false;
    for (const auto dir : coarseGrid->distribution)
    {
        //if (dir[0] > 0 || dir[1] > 0 || dir[2] > 0)
        //    continue;

        const int neighborIndex = coarseGrid->transCoordToIndex(x + dir[0] * coarseGrid->delta, y + dir[1] * coarseGrid->delta, z + dir[2] * coarseGrid->delta);
        neighborBelongsToFineToCoarseInterpolationCell = coarseGrid->getField().isFineToCoarseNode(neighborIndex);
        if (neighborBelongsToFineToCoarseInterpolationCell)
        {
            coarseGrid->getField().setFieldEntryToStopperOverlapGrid(indexOnCoarseGrid);
            break;
        }

    }
    if(!neighborBelongsToFineToCoarseInterpolationCell) //should be inside of fine grid and can be deleted
        coarseGrid->getField().setFieldEntryToInvalid(indexOnCoarseGrid);
}

bool GridInterface::isNeighborFineFluid(real x, real y, real z, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    const int neighbor = coarseGrid->transCoordToIndex(x, y, z);
    const int indexOnFineGrid = getCoarseToFineIndexOnFineGrid(neighbor, coarseGrid, fineGrid);
    if (indexOnFineGrid == -1)
        return false;
    return fineGrid->getField().isFluid(indexOnFineGrid);
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

int GridInterface::getFineToCoarseIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);
    const real xFine = x - (fineGrid->delta * 0.5);
    const real yFine = y - (fineGrid->delta * 0.5);
    const real zFine = z - (fineGrid->delta * 0.5);

    return fineGrid->transCoordToIndex(xFine, yFine, zFine);
}

void GridInterface::print() const
{
    printf("Grid Interface - CF nodes: %d, FC nodes: %d\n", cf.numberOfEntries, fc.numberOfEntries);
}
