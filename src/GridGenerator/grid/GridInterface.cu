#include "GridInterface.h"

#include "GridImp.h"
#include "Field.h"
#include "NodeValues.h"

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

    for(const auto dir : coarseGrid->distribution)
    {
        const bool isFineGridNeighborInvalid = isNeighborFineInvalid(x + dir[0] * coarseGrid->getDelta(), y + dir[1] * coarseGrid->getDelta(), z + dir[2] * coarseGrid->getDelta(), coarseGrid, fineGrid);
        if(isFineGridNeighborInvalid)
        {
            cf.coarse[cf.numberOfEntries] = indexOnCoarseGrid;
            cf.fine[cf.numberOfEntries] = indexOnFineGridCF;

            cf.numberOfEntries++;

            coarseGrid->setCellTo(indexOnCoarseGrid, FLUID_CFC);
            fineGrid->setCellTo(indexOnFineGridCF, FLUID_CFF);
            break;
        }
    }
}

void GridInterface::findInterfaceCF_GKS(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid)
{
	const bool nodeOnCoarseGridIsFluid = coarseGrid->getField().isFluid(indexOnCoarseGrid);
	if (!nodeOnCoarseGridIsFluid)
		return;

	real x, y, z;
	coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

	for (const auto dir : coarseGrid->distribution)
	{
		const uint indexOnFineGrid = fineGrid->transCoordToIndex(x + 0.25 * dir[0] * coarseGrid->getDelta(),
																 y + 0.25 * dir[1] * coarseGrid->getDelta(),
																 z + 0.25 * dir[2] * coarseGrid->getDelta());

		if (indexOnFineGrid != -1 && fineGrid->getField().is(indexOnFineGrid, STOPPER_OUT_OF_GRID)) 
		{
			coarseGrid->getField().setFieldEntry(indexOnCoarseGrid, FLUID_CFC);
			break;
		}
	}
}

void GridInterface::findInterfaceFC(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid)
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
        const int neighborIndex = coarseGrid->transCoordToIndex(x + dir[0] * coarseGrid->getDelta(), y + dir[1] * coarseGrid->getDelta(), z + dir[2] * coarseGrid->getDelta());
        const bool neighborBelongsToCoarseToFineInterpolationCell = coarseGrid->getField().isCoarseToFineNode(neighborIndex);
        if (neighborBelongsToCoarseToFineInterpolationCell)
        {
            fc.coarse[fc.numberOfEntries] = indexOnCoarseGrid;
            fc.fine[fc.numberOfEntries] = indexOnFineGridFC;

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

	real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);

    bool neighborBelongsToFineToCoarseInterpolationCell = false;
    for (const auto dir : coarseGrid->distribution)
    {
        //if (dir[0] > 0 || dir[1] > 0 || dir[2] > 0)  //only Esoteric Twist stopper, not perfectly implemented
        //    continue;								   //should not be here, should be made conditional

        const int neighborIndex = coarseGrid->transCoordToIndex(x + dir[0] * coarseGrid->getDelta(), y + dir[1] * coarseGrid->getDelta(), z + dir[2] * coarseGrid->getDelta());
        neighborBelongsToFineToCoarseInterpolationCell = coarseGrid->getField().isFineToCoarseNode(neighborIndex);
        if (neighborBelongsToFineToCoarseInterpolationCell)
        {
            coarseGrid->getField().setFieldEntryToStopperCoarseUnderFine(indexOnCoarseGrid);
            break;
        }
    }

	//should be inside of fine grid and can be deleted
    if(!neighborBelongsToFineToCoarseInterpolationCell && (fineGrid->getField().isInvalidSolid(indexOnFineGridFC) || 
	                                                       fineGrid->getField().isFluid(indexOnFineGridFC) ||
	                                                       fineGrid->getField().is(indexOnFineGridFC, STOPPER_SOLID) || 
	                                                       fineGrid->getField().is(indexOnFineGridFC, BC_GEOMETRY)))
        coarseGrid->getField().setFieldEntryToInvalidCoarseUnderFine(indexOnCoarseGrid);
}

bool GridInterface::isNeighborFineInvalid(real x, real y, real z, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    const int neighbor = coarseGrid->transCoordToIndex(x, y, z);
    const int indexOnFineGrid = getCoarseToFineIndexOnFineGrid(neighbor, coarseGrid, fineGrid);
    if (indexOnFineGrid == -1)
        return true;
    return fineGrid->getField().isInvalidOutOfGrid(indexOnFineGrid) || fineGrid->getField().isStopperOutOfGrid(indexOnFineGrid);
}

int GridInterface::getCoarseToFineIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);
    const real xFine = x + (fineGrid->getDelta() * 0.5);
    const real yFine = y + (fineGrid->getDelta() * 0.5);
    const real zFine = z + (fineGrid->getDelta() * 0.5);

    return fineGrid->transCoordToIndex(xFine, yFine, zFine);
}

int GridInterface::getFineToCoarseIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid)
{
    real x, y, z;
    coarseGrid->transIndexToCoords(indexOnCoarseGrid, x, y, z);
    const real xFine = x - (fineGrid->getDelta() * 0.5);
    const real yFine = y - (fineGrid->getDelta() * 0.5);
    const real zFine = z - (fineGrid->getDelta() * 0.5);

    return fineGrid->transCoordToIndex(xFine, yFine, zFine);
}

void GridInterface::findForGridInterfaceSparseIndexCF(GridImp* coarseGrid, GridImp* fineGrid, uint index)
{
    findSparseIndex(cf.coarse, coarseGrid, index);
    findSparseIndex(cf.fine, fineGrid, index);
}

void GridInterface::findForGridInterfaceSparseIndexFC(GridImp* coarseGrid, GridImp* fineGrid, uint index)
{
    findSparseIndex(fc.coarse, coarseGrid, index);
    findSparseIndex(fc.fine, fineGrid, index);
}

void GridInterface::findSparseIndex(uint* indices, GridImp* grid, uint index)
{
    const uint matrixIndex = indices[index];
    const uint sparseIndex = grid->getSparseIndex(matrixIndex);
    indices[index] = sparseIndex;
}

void GridInterface::print() const
{
    printf("Grid Interface - CF nodes: %d, FC nodes: %d\n", cf.numberOfEntries, fc.numberOfEntries);
}
