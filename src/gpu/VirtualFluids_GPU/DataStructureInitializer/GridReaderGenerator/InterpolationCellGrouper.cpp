#include "InterpolationCellGrouper.h"

#include "Parameter/Parameter.h"
#include <GridGenerator/grid/Grid.h>
#include <GridGenerator/grid/GridBuilder/GridBuilder.h>

InterpolationCellGrouper::InterpolationCellGrouper(const LBMSimulationParameters &parHs,
                                                   const LBMSimulationParameters &parDs, SPtr<GridBuilder> builder)
    : parHs(parHs), parDs(parDs), builder(builder)
{
}

void InterpolationCellGrouper::splitFineToCoarseIntoBorderAndBulk(uint level) const
{
    this->reorderFineToCoarseIntoBorderAndBulk(level);

    parDs[level]->fineToCoarseBorder.numberOfCells = parHs[level]->fineToCoarseBorder.numberOfCells;
    parDs[level]->fineToCoarseBulk.numberOfCells = parHs[level]->fineToCoarseBulk.numberOfCells;
    parDs[level]->fineToCoarseBorder.coarseCellIndices = parDs[level]->fineToCoarse.coarseCellIndices;
    parDs[level]->fineToCoarseBulk.coarseCellIndices = parDs[level]->fineToCoarseBorder.coarseCellIndices + parDs[level]->fineToCoarseBorder.numberOfCells;
    parDs[level]->fineToCoarseBorder.fineCellIndices = parDs[level]->fineToCoarse.fineCellIndices;
    parDs[level]->fineToCoarseBulk.fineCellIndices = parDs[level]->fineToCoarseBorder.fineCellIndices + parDs[level]->fineToCoarseBorder.numberOfCells;
    parDs[level]->neighborFineToCoarseBulk.x = parDs[level]->neighborFineToCoarse.x + parDs[level]->fineToCoarseBorder.numberOfCells;
    parDs[level]->neighborFineToCoarseBulk.y = parDs[level]->neighborFineToCoarse.y + parDs[level]->fineToCoarseBorder.numberOfCells;
    parDs[level]->neighborFineToCoarseBulk.z = parDs[level]->neighborFineToCoarse.z + parDs[level]->fineToCoarseBorder.numberOfCells;
}

void InterpolationCellGrouper::reorderFineToCoarseIntoBorderAndBulk(uint level) const
{
    // create some local variables for better readability
    uint *fineToCoarseCoarseAll = parHs[level]->fineToCoarse.coarseCellIndices;
    uint *fineToCoarseFineAll = parHs[level]->fineToCoarse.fineCellIndices;
    auto grid = this->builder->getGrid(level);

    std::vector<uint> fineToCoarseCoarseBorderVector;
    std::vector<uint> fineToCoarseCoarseBulkVector;
    std::vector<uint> fineToCoarseFineBorderVector;
    std::vector<uint> fineToCoarseFineBulkVector;
    std::vector<real> neighborXBorder;
    std::vector<real> neighborYBorder;
    std::vector<real> neighborZBorder;
    std::vector<real> neighborXBulk;
    std::vector<real> neighborYBulk;
    std::vector<real> neighborZBulk;

    // fill border and bulk vectors with interpolation cells fine to coarse
    for (uint i = 0; i < parHs[level]->fineToCoarse.numberOfCells; i++)
        if (grid->isSparseIndexInFluidNodeIndicesBorder(fineToCoarseCoarseAll[i])) {
            fineToCoarseCoarseBorderVector.push_back(fineToCoarseCoarseAll[i]);
            fineToCoarseFineBorderVector.push_back(fineToCoarseFineAll[i]);
            neighborXBorder.push_back(parHs[level]->neighborFineToCoarse.x[i]);
            neighborYBorder.push_back(parHs[level]->neighborFineToCoarse.y[i]);
            neighborZBorder.push_back(parHs[level]->neighborFineToCoarse.z[i]);
        } else {
            fineToCoarseCoarseBulkVector.push_back(fineToCoarseCoarseAll[i]);
            fineToCoarseFineBulkVector.push_back(fineToCoarseFineAll[i]);
            neighborXBulk.push_back(parHs[level]->neighborFineToCoarse.x[i]);
            neighborYBulk.push_back(parHs[level]->neighborFineToCoarse.y[i]);
            neighborZBulk.push_back(parHs[level]->neighborFineToCoarse.z[i]);
        }

    // set new sizes and pointers
    parHs[level]->fineToCoarseBorder.coarseCellIndices = fineToCoarseCoarseAll;
    parHs[level]->fineToCoarseBorder.fineCellIndices = fineToCoarseFineAll;
    parHs[level]->fineToCoarseBorder.numberOfCells = (uint)fineToCoarseCoarseBorderVector.size();
    parHs[level]->fineToCoarseBulk.numberOfCells = (uint)fineToCoarseCoarseBulkVector.size();
    parHs[level]->fineToCoarseBulk.coarseCellIndices = fineToCoarseCoarseAll + parHs[level]->fineToCoarseBorder.numberOfCells;
    parHs[level]->fineToCoarseBulk.fineCellIndices = fineToCoarseFineAll + parHs[level]->fineToCoarseBorder.numberOfCells;
    parHs[level]->neighborFineToCoarseBulk.x = parHs[level]->neighborFineToCoarse.x + parHs[level]->fineToCoarseBorder.numberOfCells;
    parHs[level]->neighborFineToCoarseBulk.y = parHs[level]->neighborFineToCoarse.y + parHs[level]->fineToCoarseBorder.numberOfCells;
    parHs[level]->neighborFineToCoarseBulk.z = parHs[level]->neighborFineToCoarse.z + parHs[level]->fineToCoarseBorder.numberOfCells;

    // copy the created vectors to the memory addresses of the old arrays
    // this is inefficient :(
    for (uint i = 0; i < (uint)fineToCoarseCoarseBorderVector.size(); i++) {
        fineToCoarseCoarseAll[i] = fineToCoarseCoarseBorderVector[i];
        fineToCoarseFineAll[i] = fineToCoarseFineBorderVector[i];
        parHs[level]->neighborFineToCoarse.x[i] = neighborXBorder[i];
        parHs[level]->neighborFineToCoarse.y[i] = neighborYBorder[i];
        parHs[level]->neighborFineToCoarse.z[i] = neighborZBorder[i];
    }
    for (uint i = 0; i < (uint)fineToCoarseCoarseBulkVector.size(); i++) {
        parHs[level]->fineToCoarseBulk.coarseCellIndices[i] = fineToCoarseCoarseBulkVector[i];
        parHs[level]->fineToCoarseBulk.fineCellIndices[i] = fineToCoarseFineBulkVector[i];
        parHs[level]->neighborFineToCoarseBulk.x[i] = neighborXBulk[i];
        parHs[level]->neighborFineToCoarseBulk.y[i] = neighborYBulk[i];
        parHs[level]->neighborFineToCoarseBulk.z[i] = neighborZBulk[i];
    }
}

void InterpolationCellGrouper::splitCoarseToFineIntoBorderAndBulk(uint level) const
{
    this->reorderCoarseToFineIntoBorderAndBulk(level);

    parDs[level]->coarseToFineBorder.numberOfCells = parHs[level]->coarseToFineBorder.numberOfCells;
    parDs[level]->coarseToFineBulk.numberOfCells = parHs[level]->coarseToFineBulk.numberOfCells;
    parDs[level]->coarseToFineBorder.coarseCellIndices = parDs[level]->coarseToFine.coarseCellIndices;
    parDs[level]->coarseToFineBulk.coarseCellIndices = parDs[level]->coarseToFineBorder.coarseCellIndices + parDs[level]->coarseToFineBorder.numberOfCells;
    parDs[level]->coarseToFineBorder.fineCellIndices = parDs[level]->coarseToFine.fineCellIndices;
    parDs[level]->coarseToFineBulk.fineCellIndices = parDs[level]->coarseToFineBorder.fineCellIndices + parDs[level]->coarseToFineBorder.numberOfCells;
    parDs[level]->neighborCoarseToFineBulk.x = parDs[level]->neighborCoarseToFine.x + parDs[level]->coarseToFineBorder.numberOfCells;
    parDs[level]->neighborCoarseToFineBulk.y = parDs[level]->neighborCoarseToFine.y + parDs[level]->coarseToFineBorder.numberOfCells;
    parDs[level]->neighborCoarseToFineBulk.z = parDs[level]->neighborCoarseToFine.z + parDs[level]->coarseToFineBorder.numberOfCells;
}

void InterpolationCellGrouper::reorderCoarseToFineIntoBorderAndBulk(uint level) const
{
    // create some local variables for better readability
    uint *coarseToFineCoarseAll = parHs[level]->coarseToFine.coarseCellIndices;
    uint *coarseToFineFineAll = parHs[level]->coarseToFine.fineCellIndices;
    uint *neighborX = this->parHs[level]->neighborX;
    uint *neighborY = this->parHs[level]->neighborY;
    uint *neighborZ = this->parHs[level]->neighborZ;
    auto grid = this->builder->getGrid(level);

    std::vector<uint> coarseToFineCoarseBorderVector;
    std::vector<uint> coarseToFineCoarseBulkVector;
    std::vector<uint> coarseToFineFineBorderVector;
    std::vector<uint> coarseToFineFineBulkVector;
    std::vector<real> neighborXBorder;
    std::vector<real> neighborYBorder;
    std::vector<real> neighborZBorder;
    std::vector<real> neighborXBulk;
    std::vector<real> neighborYBulk;
    std::vector<real> neighborZBulk;
    uint sparseIndexOfICellBSW;

    // fill border and bulk vectors with interpolation cells coarse to fine
    for (uint i = 0; i < parHs[level]->coarseToFine.numberOfCells; i++) {
        sparseIndexOfICellBSW = coarseToFineCoarseAll[i];

        if (grid->isSparseIndexInFluidNodeIndicesBorder(sparseIndexOfICellBSW) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborX[sparseIndexOfICellBSW]) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborY[sparseIndexOfICellBSW]) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborZ[sparseIndexOfICellBSW]) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborY[neighborX[sparseIndexOfICellBSW]]) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborZ[neighborX[sparseIndexOfICellBSW]]) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborZ[neighborY[sparseIndexOfICellBSW]]) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborZ[neighborY[neighborX[sparseIndexOfICellBSW]]])) {

            coarseToFineCoarseBorderVector.push_back(coarseToFineCoarseAll[i]);
            coarseToFineFineBorderVector.push_back(coarseToFineFineAll[i]);
            neighborXBorder.push_back(parHs[level]->neighborCoarseToFine.x[i]);
            neighborYBorder.push_back(parHs[level]->neighborCoarseToFine.y[i]);
            neighborZBorder.push_back(parHs[level]->neighborCoarseToFine.z[i]);
        } else {
            coarseToFineCoarseBulkVector.push_back(coarseToFineCoarseAll[i]);
            coarseToFineFineBulkVector.push_back(coarseToFineFineAll[i]);
            neighborXBulk.push_back(parHs[level]->neighborCoarseToFine.x[i]);
            neighborYBulk.push_back(parHs[level]->neighborCoarseToFine.y[i]);
            neighborZBulk.push_back(parHs[level]->neighborCoarseToFine.z[i]);
        }
    }

    // set new sizes and pointers
    parHs[level]->coarseToFineBorder.coarseCellIndices = parHs[level]->coarseToFine.coarseCellIndices;
    parHs[level]->coarseToFineBorder.fineCellIndices = parHs[level]->coarseToFine.fineCellIndices;
    parHs[level]->coarseToFineBorder.numberOfCells = (uint)coarseToFineCoarseBorderVector.size();
    parHs[level]->coarseToFineBulk.numberOfCells = (uint)coarseToFineCoarseBulkVector.size();
    parHs[level]->coarseToFineBulk.coarseCellIndices = parHs[level]->coarseToFine.coarseCellIndices + parHs[level]->coarseToFineBorder.numberOfCells;
    parHs[level]->coarseToFineBulk.fineCellIndices = parHs[level]->coarseToFine.fineCellIndices + parHs[level]->coarseToFineBorder.numberOfCells;
    parHs[level]->neighborCoarseToFineBulk.x = parHs[level]->neighborCoarseToFine.x + parHs[level]->coarseToFineBorder.numberOfCells;
    parHs[level]->neighborCoarseToFineBulk.y = parHs[level]->neighborCoarseToFine.y + parHs[level]->coarseToFineBorder.numberOfCells;
    parHs[level]->neighborCoarseToFineBulk.z = parHs[level]->neighborCoarseToFine.z + parHs[level]->coarseToFineBorder.numberOfCells;

    // copy the created vectors to the memory addresses of the old arrays
    // this is inefficient :(
    for (uint i = 0; i < (uint)coarseToFineCoarseBorderVector.size(); i++) {
        parHs[level]->coarseToFineBorder.coarseCellIndices[i] = coarseToFineCoarseBorderVector[i];
        parHs[level]->coarseToFineBorder.fineCellIndices[i] = coarseToFineFineBorderVector[i];
        parHs[level]->neighborCoarseToFine.x[i] = neighborXBorder[i];
        parHs[level]->neighborCoarseToFine.y[i] = neighborYBorder[i];
        parHs[level]->neighborCoarseToFine.z[i] = neighborZBorder[i];
    }
    for (uint i = 0; i < (uint)coarseToFineCoarseBulkVector.size(); i++) {
        parHs[level]->coarseToFineBulk.coarseCellIndices[i] = coarseToFineCoarseBulkVector[i];
        parHs[level]->coarseToFineBulk.fineCellIndices[i] = coarseToFineFineBulkVector[i];
        parHs[level]->neighborCoarseToFineBulk.x[i] = neighborXBulk[i];
        parHs[level]->neighborCoarseToFineBulk.y[i] = neighborYBulk[i];
        parHs[level]->neighborCoarseToFineBulk.z[i] = neighborZBulk[i];
    }
}
