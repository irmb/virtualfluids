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

    parDs[level]->intFCBorder.kFC = parHs[level]->intFCBorder.kFC;
    parDs[level]->intFCBulk.kFC = parHs[level]->intFCBulk.kFC;
    parDs[level]->intFCBorder.ICellFCC = parDs[level]->intFC.ICellFCC;
    parDs[level]->intFCBulk.ICellFCC = parDs[level]->intFCBorder.ICellFCC + parDs[level]->intFCBorder.kFC;
    parDs[level]->intFCBorder.ICellFCF = parDs[level]->intFC.ICellFCF;
    parDs[level]->intFCBulk.ICellFCF = parDs[level]->intFCBorder.ICellFCF + parDs[level]->intFCBorder.kFC;
    parDs[level]->neighborFCBulk.x = parDs[level]->neighborFC.x + parDs[level]->intFCBorder.kFC;
    parDs[level]->neighborFCBulk.y = parDs[level]->neighborFC.y + parDs[level]->intFCBorder.kFC;
    parDs[level]->neighborFCBulk.z = parDs[level]->neighborFC.z + parDs[level]->intFCBorder.kFC;
}

void InterpolationCellGrouper::reorderFineToCoarseIntoBorderAndBulk(uint level) const
{
    // create some local variables for better readability
    uint *iCellFccAll = parHs[level]->intFC.ICellFCC;
    uint *iCellFcfAll = parHs[level]->intFC.ICellFCF;
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
    for (uint i = 0; i < parHs[level]->intFC.kFC; i++)
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
    parHs[level]->intFCBorder.ICellFCC = iCellFccAll;
    parHs[level]->intFCBorder.ICellFCF = iCellFcfAll;
    parHs[level]->intFCBorder.kFC = (uint)iCellFccBorderVector.size();
    parHs[level]->intFCBulk.kFC = (uint)iCellFccBulkVector.size();
    parHs[level]->intFCBulk.ICellFCC = iCellFccAll + parHs[level]->intFCBorder.kFC;
    parHs[level]->intFCBulk.ICellFCF = iCellFcfAll + parHs[level]->intFCBorder.kFC;
    parHs[level]->neighborFCBulk.x = parHs[level]->neighborFC.x + parHs[level]->intFCBorder.kFC;
    parHs[level]->neighborFCBulk.y = parHs[level]->neighborFC.y + parHs[level]->intFCBorder.kFC;
    parHs[level]->neighborFCBulk.z = parHs[level]->neighborFC.z + parHs[level]->intFCBorder.kFC;

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
        parHs[level]->intFCBulk.ICellFCC[i] = iCellFccBulkVector[i];
        parHs[level]->intFCBulk.ICellFCF[i] = iCellFcfBulkVector[i];
        parHs[level]->neighborFCBulk.x[i] = xBulkVector[i];
        parHs[level]->neighborFCBulk.y[i] = yBulkVector[i];
        parHs[level]->neighborFCBulk.z[i] = zBulkVector[i];
    }
}

void InterpolationCellGrouper::splitCoarseToFineIntoBorderAndBulk(uint level) const
{
    this->reorderCoarseToFineIntoBorderAndBulk(level);

    parDs[level]->intCFBorder.kCF = parHs[level]->intCFBorder.kCF;
    parDs[level]->intCFBulk.kCF = parHs[level]->intCFBulk.kCF;
    parDs[level]->intCFBorder.ICellCFC = parDs[level]->intCF.ICellCFC;
    parDs[level]->intCFBulk.ICellCFC = parDs[level]->intCFBorder.ICellCFC + parDs[level]->intCFBorder.kCF;
    parDs[level]->intCFBorder.ICellCFF = parDs[level]->intCF.ICellCFF;
    parDs[level]->intCFBulk.ICellCFF = parDs[level]->intCFBorder.ICellCFF + parDs[level]->intCFBorder.kCF;
    parDs[level]->neighborCFBulk.x = parDs[level]->neighborCF.x + parDs[level]->intCFBorder.kCF;
    parDs[level]->neighborCFBulk.y = parDs[level]->neighborCF.y + parDs[level]->intCFBorder.kCF;
    parDs[level]->neighborCFBulk.z = parDs[level]->neighborCF.z + parDs[level]->intCFBorder.kCF;
}

void InterpolationCellGrouper::reorderCoarseToFineIntoBorderAndBulk(uint level) const
{
    // create some local variables for better readability
    uint *iCellCfcAll = parHs[level]->intCF.ICellCFC;
    uint *iCellCffAll = parHs[level]->intCF.ICellCFF;
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
    for (uint i = 0; i < parHs[level]->intCF.kCF; i++) {
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
    parHs[level]->intCFBorder.ICellCFC = parHs[level]->intCF.ICellCFC;
    parHs[level]->intCFBorder.ICellCFF = parHs[level]->intCF.ICellCFF;
    parHs[level]->intCFBorder.kCF = (uint)iCellCfcBorderVector.size();
    parHs[level]->intCFBulk.kCF = (uint)iCellCfcBulkVector.size();
    parHs[level]->intCFBulk.ICellCFC = parHs[level]->intCF.ICellCFC + parHs[level]->intCFBorder.kCF;
    parHs[level]->intCFBulk.ICellCFF = parHs[level]->intCF.ICellCFF + parHs[level]->intCFBorder.kCF;
    parHs[level]->neighborCFBulk.x = parHs[level]->neighborCF.x + parHs[level]->intCFBorder.kCF;
    parHs[level]->neighborCFBulk.y = parHs[level]->neighborCF.y + parHs[level]->intCFBorder.kCF;
    parHs[level]->neighborCFBulk.z = parHs[level]->neighborCF.z + parHs[level]->intCFBorder.kCF;

    // copy the created vectors to the memory addresses of the old arrays
    // this is inefficient :(
    for (uint i = 0; i < (uint)iCellCfcBorderVector.size(); i++) {
        parHs[level]->intCFBorder.ICellCFC[i] = iCellCfcBorderVector[i];
        parHs[level]->intCFBorder.ICellCFF[i] = iCellCffBorderVector[i];
        parHs[level]->neighborCF.x[i] = xBorderVector[i];
        parHs[level]->neighborCF.y[i] = yBorderVector[i];
        parHs[level]->neighborCF.z[i] = zBorderVector[i];
    }
    for (uint i = 0; i < (uint)iCellCfcBulkVector.size(); i++) {
        parHs[level]->intCFBulk.ICellCFC[i] = iCellCfcBulkVector[i];
        parHs[level]->intCFBulk.ICellCFF[i] = iCellCffBulkVector[i];
        parHs[level]->neighborCFBulk.x[i] = xBulkVector[i];
        parHs[level]->neighborCFBulk.y[i] = yBulkVector[i];
        parHs[level]->neighborCFBulk.z[i] = zBulkVector[i];
    }
}
