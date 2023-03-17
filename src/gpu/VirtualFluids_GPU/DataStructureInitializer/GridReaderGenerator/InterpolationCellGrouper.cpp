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

    parDs[level]->intFCBorder.numberOfCells = parHs[level]->intFCBorder.numberOfCells;
    parDs[level]->intFCBulk.numberOfCells = parHs[level]->intFCBulk.numberOfCells;
    parDs[level]->intFCBorder.coarseCellIndices = parDs[level]->fineToCoarse.coarseCellIndices;
    parDs[level]->intFCBulk.coarseCellIndices = parDs[level]->intFCBorder.coarseCellIndices + parDs[level]->intFCBorder.numberOfCells;
    parDs[level]->intFCBorder.fineCellIndices = parDs[level]->fineToCoarse.fineCellIndices;
    parDs[level]->intFCBulk.fineCellIndices = parDs[level]->intFCBorder.fineCellIndices + parDs[level]->intFCBorder.numberOfCells;
    parDs[level]->neighborFCBulk.x = parDs[level]->neighborFC.x + parDs[level]->intFCBorder.numberOfCells;
    parDs[level]->neighborFCBulk.y = parDs[level]->neighborFC.y + parDs[level]->intFCBorder.numberOfCells;
    parDs[level]->neighborFCBulk.z = parDs[level]->neighborFC.z + parDs[level]->intFCBorder.numberOfCells;
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
    parHs[level]->intFCBorder.coarseCellIndices = iCellFccAll;
    parHs[level]->intFCBorder.fineCellIndices = iCellFcfAll;
    parHs[level]->intFCBorder.numberOfCells = (uint)iCellFccBorderVector.size();
    parHs[level]->intFCBulk.numberOfCells = (uint)iCellFccBulkVector.size();
    parHs[level]->intFCBulk.coarseCellIndices = iCellFccAll + parHs[level]->intFCBorder.numberOfCells;
    parHs[level]->intFCBulk.fineCellIndices = iCellFcfAll + parHs[level]->intFCBorder.numberOfCells;
    parHs[level]->neighborFCBulk.x = parHs[level]->neighborFC.x + parHs[level]->intFCBorder.numberOfCells;
    parHs[level]->neighborFCBulk.y = parHs[level]->neighborFC.y + parHs[level]->intFCBorder.numberOfCells;
    parHs[level]->neighborFCBulk.z = parHs[level]->neighborFC.z + parHs[level]->intFCBorder.numberOfCells;

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
        parHs[level]->intFCBulk.coarseCellIndices[i] = iCellFccBulkVector[i];
        parHs[level]->intFCBulk.fineCellIndices[i] = iCellFcfBulkVector[i];
        parHs[level]->neighborFCBulk.x[i] = xBulkVector[i];
        parHs[level]->neighborFCBulk.y[i] = yBulkVector[i];
        parHs[level]->neighborFCBulk.z[i] = zBulkVector[i];
    }
}

void InterpolationCellGrouper::splitCoarseToFineIntoBorderAndBulk(uint level) const
{
    this->reorderCoarseToFineIntoBorderAndBulk(level);

    parDs[level]->intCFBorder.numberOfCells = parHs[level]->intCFBorder.numberOfCells;
    parDs[level]->intCFBulk.numberOfCells = parHs[level]->intCFBulk.numberOfCells;
    parDs[level]->intCFBorder.coarseCellIndices = parDs[level]->coarseToFine.coarseCellIndices;
    parDs[level]->intCFBulk.coarseCellIndices = parDs[level]->intCFBorder.coarseCellIndices + parDs[level]->intCFBorder.numberOfCells;
    parDs[level]->intCFBorder.fineCellIndices = parDs[level]->coarseToFine.fineCellIndices;
    parDs[level]->intCFBulk.fineCellIndices = parDs[level]->intCFBorder.fineCellIndices + parDs[level]->intCFBorder.numberOfCells;
    parDs[level]->neighborCFBulk.x = parDs[level]->neighborCF.x + parDs[level]->intCFBorder.numberOfCells;
    parDs[level]->neighborCFBulk.y = parDs[level]->neighborCF.y + parDs[level]->intCFBorder.numberOfCells;
    parDs[level]->neighborCFBulk.z = parDs[level]->neighborCF.z + parDs[level]->intCFBorder.numberOfCells;
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
    parHs[level]->intCFBorder.coarseCellIndices = parHs[level]->coarseToFine.coarseCellIndices;
    parHs[level]->intCFBorder.fineCellIndices = parHs[level]->coarseToFine.fineCellIndices;
    parHs[level]->intCFBorder.numberOfCells = (uint)iCellCfcBorderVector.size();
    parHs[level]->intCFBulk.numberOfCells = (uint)iCellCfcBulkVector.size();
    parHs[level]->intCFBulk.coarseCellIndices = parHs[level]->coarseToFine.coarseCellIndices + parHs[level]->intCFBorder.numberOfCells;
    parHs[level]->intCFBulk.fineCellIndices = parHs[level]->coarseToFine.fineCellIndices + parHs[level]->intCFBorder.numberOfCells;
    parHs[level]->neighborCFBulk.x = parHs[level]->neighborCF.x + parHs[level]->intCFBorder.numberOfCells;
    parHs[level]->neighborCFBulk.y = parHs[level]->neighborCF.y + parHs[level]->intCFBorder.numberOfCells;
    parHs[level]->neighborCFBulk.z = parHs[level]->neighborCF.z + parHs[level]->intCFBorder.numberOfCells;

    // copy the created vectors to the memory addresses of the old arrays
    // this is inefficient :(
    for (uint i = 0; i < (uint)iCellCfcBorderVector.size(); i++) {
        parHs[level]->intCFBorder.coarseCellIndices[i] = iCellCfcBorderVector[i];
        parHs[level]->intCFBorder.fineCellIndices[i] = iCellCffBorderVector[i];
        parHs[level]->neighborCF.x[i] = xBorderVector[i];
        parHs[level]->neighborCF.y[i] = yBorderVector[i];
        parHs[level]->neighborCF.z[i] = zBorderVector[i];
    }
    for (uint i = 0; i < (uint)iCellCfcBulkVector.size(); i++) {
        parHs[level]->intCFBulk.coarseCellIndices[i] = iCellCfcBulkVector[i];
        parHs[level]->intCFBulk.fineCellIndices[i] = iCellCffBulkVector[i];
        parHs[level]->neighborCFBulk.x[i] = xBulkVector[i];
        parHs[level]->neighborCFBulk.y[i] = yBulkVector[i];
        parHs[level]->neighborCFBulk.z[i] = zBulkVector[i];
    }
}
