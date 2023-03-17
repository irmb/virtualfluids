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
    uint *iCellFccAll = parHs[level]->fineToCoarse.coarseCellIndices;
    uint *iCellFcfAll = parHs[level]->fineToCoarse.fineCellIndices;
    auto grid = this->builder->getGrid(level);

    std::vector<uint> iCellFccBorderVector;
    std::vector<uint> iCellFccBulkVector;
    std::vector<uint> iCellFcfBorderVector;
    std::vector<uint> iCellFcfBulkVector;
    std::vector<real> xBorderVector;
    std::vector<real> yBorderVector;
    std::vector<real> zBorderVector;
    std::vector<real> xBulkVector;
    std::vector<real> yBulkVector;
    std::vector<real> zBulkVector;

    // fill border and bulk vectors with iCellFCs
    for (uint i = 0; i < parHs[level]->fineToCoarse.numberOfCells; i++)
        if (grid->isSparseIndexInFluidNodeIndicesBorder(iCellFccAll[i])) {
            iCellFccBorderVector.push_back(iCellFccAll[i]);
            iCellFcfBorderVector.push_back(iCellFcfAll[i]);
            xBorderVector.push_back(parHs[level]->neighborFineToCoarse.x[i]);
            yBorderVector.push_back(parHs[level]->neighborFineToCoarse.y[i]);
            zBorderVector.push_back(parHs[level]->neighborFineToCoarse.z[i]);
        } else {
            iCellFccBulkVector.push_back(iCellFccAll[i]);
            iCellFcfBulkVector.push_back(iCellFcfAll[i]);
            xBulkVector.push_back(parHs[level]->neighborFineToCoarse.x[i]);
            yBulkVector.push_back(parHs[level]->neighborFineToCoarse.y[i]);
            zBulkVector.push_back(parHs[level]->neighborFineToCoarse.z[i]);
        }

    // set new sizes and pointers
    parHs[level]->fineToCoarseBorder.coarseCellIndices = iCellFccAll;
    parHs[level]->fineToCoarseBorder.fineCellIndices = iCellFcfAll;
    parHs[level]->fineToCoarseBorder.numberOfCells = (uint)iCellFccBorderVector.size();
    parHs[level]->fineToCoarseBulk.numberOfCells = (uint)iCellFccBulkVector.size();
    parHs[level]->fineToCoarseBulk.coarseCellIndices = iCellFccAll + parHs[level]->fineToCoarseBorder.numberOfCells;
    parHs[level]->fineToCoarseBulk.fineCellIndices = iCellFcfAll + parHs[level]->fineToCoarseBorder.numberOfCells;
    parHs[level]->neighborFineToCoarseBulk.x = parHs[level]->neighborFineToCoarse.x + parHs[level]->fineToCoarseBorder.numberOfCells;
    parHs[level]->neighborFineToCoarseBulk.y = parHs[level]->neighborFineToCoarse.y + parHs[level]->fineToCoarseBorder.numberOfCells;
    parHs[level]->neighborFineToCoarseBulk.z = parHs[level]->neighborFineToCoarse.z + parHs[level]->fineToCoarseBorder.numberOfCells;

    // copy the created vectors to the memory addresses of the old arrays
    // this is inefficient :(
    for (uint i = 0; i < (uint)iCellFccBorderVector.size(); i++) {
        iCellFccAll[i] = iCellFccBorderVector[i];
        iCellFcfAll[i] = iCellFcfBorderVector[i];
        parHs[level]->neighborFineToCoarse.x[i] = xBorderVector[i];
        parHs[level]->neighborFineToCoarse.y[i] = yBorderVector[i];
        parHs[level]->neighborFineToCoarse.z[i] = zBorderVector[i];
    }
    for (uint i = 0; i < (uint)iCellFccBulkVector.size(); i++) {
        parHs[level]->fineToCoarseBulk.coarseCellIndices[i] = iCellFccBulkVector[i];
        parHs[level]->fineToCoarseBulk.fineCellIndices[i] = iCellFcfBulkVector[i];
        parHs[level]->neighborFineToCoarseBulk.x[i] = xBulkVector[i];
        parHs[level]->neighborFineToCoarseBulk.y[i] = yBulkVector[i];
        parHs[level]->neighborFineToCoarseBulk.z[i] = zBulkVector[i];
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
    uint *iCellCfcAll = parHs[level]->coarseToFine.coarseCellIndices;
    uint *iCellCffAll = parHs[level]->coarseToFine.fineCellIndices;
    uint *neighborX = this->parHs[level]->neighborX;
    uint *neighborY = this->parHs[level]->neighborY;
    uint *neighborZ = this->parHs[level]->neighborZ;
    auto grid = this->builder->getGrid(level);

    std::vector<uint> iCellCfcBorderVector;
    std::vector<uint> iCellCfcBulkVector;
    std::vector<uint> iCellCffBorderVector;
    std::vector<uint> iCellCffBulkVector;
    std::vector<real> xBorderVector;
    std::vector<real> yBorderVector;
    std::vector<real> zBorderVector;
    std::vector<real> xBulkVector;
    std::vector<real> yBulkVector;
    std::vector<real> zBulkVector;
    uint sparseIndexOfICellBSW;

    // fill border and bulk vectors with iCellCFs
    for (uint i = 0; i < parHs[level]->coarseToFine.numberOfCells; i++) {
        sparseIndexOfICellBSW = iCellCfcAll[i];

        if (grid->isSparseIndexInFluidNodeIndicesBorder(sparseIndexOfICellBSW) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborX[sparseIndexOfICellBSW]) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborY[sparseIndexOfICellBSW]) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborZ[sparseIndexOfICellBSW]) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborY[neighborX[sparseIndexOfICellBSW]]) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborZ[neighborX[sparseIndexOfICellBSW]]) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborZ[neighborY[sparseIndexOfICellBSW]]) ||
            grid->isSparseIndexInFluidNodeIndicesBorder(neighborZ[neighborY[neighborX[sparseIndexOfICellBSW]]])) {

            iCellCfcBorderVector.push_back(iCellCfcAll[i]);
            iCellCffBorderVector.push_back(iCellCffAll[i]);
            xBorderVector.push_back(parHs[level]->neighborCoarseToFine.x[i]);
            yBorderVector.push_back(parHs[level]->neighborCoarseToFine.y[i]);
            zBorderVector.push_back(parHs[level]->neighborCoarseToFine.z[i]);
        } else {
            iCellCfcBulkVector.push_back(iCellCfcAll[i]);
            iCellCffBulkVector.push_back(iCellCffAll[i]);
            xBulkVector.push_back(parHs[level]->neighborCoarseToFine.x[i]);
            yBulkVector.push_back(parHs[level]->neighborCoarseToFine.y[i]);
            zBulkVector.push_back(parHs[level]->neighborCoarseToFine.z[i]);
        }
    }

    // set new sizes and pointers
    parHs[level]->coarseToFineBorder.coarseCellIndices = parHs[level]->coarseToFine.coarseCellIndices;
    parHs[level]->coarseToFineBorder.fineCellIndices = parHs[level]->coarseToFine.fineCellIndices;
    parHs[level]->coarseToFineBorder.numberOfCells = (uint)iCellCfcBorderVector.size();
    parHs[level]->coarseToFineBulk.numberOfCells = (uint)iCellCfcBulkVector.size();
    parHs[level]->coarseToFineBulk.coarseCellIndices = parHs[level]->coarseToFine.coarseCellIndices + parHs[level]->coarseToFineBorder.numberOfCells;
    parHs[level]->coarseToFineBulk.fineCellIndices = parHs[level]->coarseToFine.fineCellIndices + parHs[level]->coarseToFineBorder.numberOfCells;
    parHs[level]->neighborCoarseToFineBulk.x = parHs[level]->neighborCoarseToFine.x + parHs[level]->coarseToFineBorder.numberOfCells;
    parHs[level]->neighborCoarseToFineBulk.y = parHs[level]->neighborCoarseToFine.y + parHs[level]->coarseToFineBorder.numberOfCells;
    parHs[level]->neighborCoarseToFineBulk.z = parHs[level]->neighborCoarseToFine.z + parHs[level]->coarseToFineBorder.numberOfCells;

    // copy the created vectors to the memory addresses of the old arrays
    // this is inefficient :(
    for (uint i = 0; i < (uint)iCellCfcBorderVector.size(); i++) {
        parHs[level]->coarseToFineBorder.coarseCellIndices[i] = iCellCfcBorderVector[i];
        parHs[level]->coarseToFineBorder.fineCellIndices[i] = iCellCffBorderVector[i];
        parHs[level]->neighborCoarseToFine.x[i] = xBorderVector[i];
        parHs[level]->neighborCoarseToFine.y[i] = yBorderVector[i];
        parHs[level]->neighborCoarseToFine.z[i] = zBorderVector[i];
    }
    for (uint i = 0; i < (uint)iCellCfcBulkVector.size(); i++) {
        parHs[level]->coarseToFineBulk.coarseCellIndices[i] = iCellCfcBulkVector[i];
        parHs[level]->coarseToFineBulk.fineCellIndices[i] = iCellCffBulkVector[i];
        parHs[level]->neighborCoarseToFineBulk.x[i] = xBulkVector[i];
        parHs[level]->neighborCoarseToFineBulk.y[i] = yBulkVector[i];
        parHs[level]->neighborCoarseToFineBulk.z[i] = zBulkVector[i];
    }
}
