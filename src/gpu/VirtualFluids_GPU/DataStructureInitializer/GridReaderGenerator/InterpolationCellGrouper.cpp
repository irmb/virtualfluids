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
    parDs[level]->neighborFCBulk.x = parDs[level]->neighborFC.x + parDs[level]->fineToCoarseBorder.numberOfCells;
    parDs[level]->neighborFCBulk.y = parDs[level]->neighborFC.y + parDs[level]->fineToCoarseBorder.numberOfCells;
    parDs[level]->neighborFCBulk.z = parDs[level]->neighborFC.z + parDs[level]->fineToCoarseBorder.numberOfCells;
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
            xBorderVector.push_back(parHs[level]->neighborFC.x[i]);
            yBorderVector.push_back(parHs[level]->neighborFC.y[i]);
            zBorderVector.push_back(parHs[level]->neighborFC.z[i]);
        } else {
            iCellFccBulkVector.push_back(iCellFccAll[i]);
            iCellFcfBulkVector.push_back(iCellFcfAll[i]);
            xBulkVector.push_back(parHs[level]->neighborFC.x[i]);
            yBulkVector.push_back(parHs[level]->neighborFC.y[i]);
            zBulkVector.push_back(parHs[level]->neighborFC.z[i]);
        }

    // set new sizes and pointers
    parHs[level]->fineToCoarseBorder.coarseCellIndices = iCellFccAll;
    parHs[level]->fineToCoarseBorder.fineCellIndices = iCellFcfAll;
    parHs[level]->fineToCoarseBorder.numberOfCells = (uint)iCellFccBorderVector.size();
    parHs[level]->fineToCoarseBulk.numberOfCells = (uint)iCellFccBulkVector.size();
    parHs[level]->fineToCoarseBulk.coarseCellIndices = iCellFccAll + parHs[level]->fineToCoarseBorder.numberOfCells;
    parHs[level]->fineToCoarseBulk.fineCellIndices = iCellFcfAll + parHs[level]->fineToCoarseBorder.numberOfCells;
    parHs[level]->neighborFCBulk.x = parHs[level]->neighborFC.x + parHs[level]->fineToCoarseBorder.numberOfCells;
    parHs[level]->neighborFCBulk.y = parHs[level]->neighborFC.y + parHs[level]->fineToCoarseBorder.numberOfCells;
    parHs[level]->neighborFCBulk.z = parHs[level]->neighborFC.z + parHs[level]->fineToCoarseBorder.numberOfCells;

    // copy the created vectors to the memory addresses of the old arrays
    // this is inefficient :(
    for (uint i = 0; i < (uint)iCellFccBorderVector.size(); i++) {
        iCellFccAll[i] = iCellFccBorderVector[i];
        iCellFcfAll[i] = iCellFcfBorderVector[i];
        parHs[level]->neighborFC.x[i] = xBorderVector[i];
        parHs[level]->neighborFC.y[i] = yBorderVector[i];
        parHs[level]->neighborFC.z[i] = zBorderVector[i];
    }
    for (uint i = 0; i < (uint)iCellFccBulkVector.size(); i++) {
        parHs[level]->fineToCoarseBulk.coarseCellIndices[i] = iCellFccBulkVector[i];
        parHs[level]->fineToCoarseBulk.fineCellIndices[i] = iCellFcfBulkVector[i];
        parHs[level]->neighborFCBulk.x[i] = xBulkVector[i];
        parHs[level]->neighborFCBulk.y[i] = yBulkVector[i];
        parHs[level]->neighborFCBulk.z[i] = zBulkVector[i];
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
    parDs[level]->neighborCFBulk.x = parDs[level]->neighborCF.x + parDs[level]->coarseToFineBorder.numberOfCells;
    parDs[level]->neighborCFBulk.y = parDs[level]->neighborCF.y + parDs[level]->coarseToFineBorder.numberOfCells;
    parDs[level]->neighborCFBulk.z = parDs[level]->neighborCF.z + parDs[level]->coarseToFineBorder.numberOfCells;
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
            xBorderVector.push_back(parHs[level]->neighborCF.x[i]);
            yBorderVector.push_back(parHs[level]->neighborCF.y[i]);
            zBorderVector.push_back(parHs[level]->neighborCF.z[i]);
        } else {
            iCellCfcBulkVector.push_back(iCellCfcAll[i]);
            iCellCffBulkVector.push_back(iCellCffAll[i]);
            xBulkVector.push_back(parHs[level]->neighborCF.x[i]);
            yBulkVector.push_back(parHs[level]->neighborCF.y[i]);
            zBulkVector.push_back(parHs[level]->neighborCF.z[i]);
        }
    }

    // set new sizes and pointers
    parHs[level]->coarseToFineBorder.coarseCellIndices = parHs[level]->coarseToFine.coarseCellIndices;
    parHs[level]->coarseToFineBorder.fineCellIndices = parHs[level]->coarseToFine.fineCellIndices;
    parHs[level]->coarseToFineBorder.numberOfCells = (uint)iCellCfcBorderVector.size();
    parHs[level]->coarseToFineBulk.numberOfCells = (uint)iCellCfcBulkVector.size();
    parHs[level]->coarseToFineBulk.coarseCellIndices = parHs[level]->coarseToFine.coarseCellIndices + parHs[level]->coarseToFineBorder.numberOfCells;
    parHs[level]->coarseToFineBulk.fineCellIndices = parHs[level]->coarseToFine.fineCellIndices + parHs[level]->coarseToFineBorder.numberOfCells;
    parHs[level]->neighborCFBulk.x = parHs[level]->neighborCF.x + parHs[level]->coarseToFineBorder.numberOfCells;
    parHs[level]->neighborCFBulk.y = parHs[level]->neighborCF.y + parHs[level]->coarseToFineBorder.numberOfCells;
    parHs[level]->neighborCFBulk.z = parHs[level]->neighborCF.z + parHs[level]->coarseToFineBorder.numberOfCells;

    // copy the created vectors to the memory addresses of the old arrays
    // this is inefficient :(
    for (uint i = 0; i < (uint)iCellCfcBorderVector.size(); i++) {
        parHs[level]->coarseToFineBorder.coarseCellIndices[i] = iCellCfcBorderVector[i];
        parHs[level]->coarseToFineBorder.fineCellIndices[i] = iCellCffBorderVector[i];
        parHs[level]->neighborCF.x[i] = xBorderVector[i];
        parHs[level]->neighborCF.y[i] = yBorderVector[i];
        parHs[level]->neighborCF.z[i] = zBorderVector[i];
    }
    for (uint i = 0; i < (uint)iCellCfcBulkVector.size(); i++) {
        parHs[level]->coarseToFineBulk.coarseCellIndices[i] = iCellCfcBulkVector[i];
        parHs[level]->coarseToFineBulk.fineCellIndices[i] = iCellCffBulkVector[i];
        parHs[level]->neighborCFBulk.x[i] = xBulkVector[i];
        parHs[level]->neighborCFBulk.y[i] = yBulkVector[i];
        parHs[level]->neighborCFBulk.z[i] = zBulkVector[i];
    }
}
