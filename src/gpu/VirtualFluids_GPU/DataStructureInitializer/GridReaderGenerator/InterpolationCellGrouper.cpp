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
    parDs[level]->offFCBulk.xOffFC = parDs[level]->offFC.xOffFC + parDs[level]->intFCBorder.kFC;
    parDs[level]->offFCBulk.yOffFC = parDs[level]->offFC.yOffFC + parDs[level]->intFCBorder.kFC;
    parDs[level]->offFCBulk.zOffFC = parDs[level]->offFC.zOffFC + parDs[level]->intFCBorder.kFC;
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
    std::vector<real> xOffFCBorderVector;
    std::vector<real> yOffFCBorderVector;
    std::vector<real> zOffFCBorderVector;
    std::vector<real> xOffFCBulkVector;
    std::vector<real> yOffFCBulkVector;
    std::vector<real> zOffFCBulkVector;

    // fill border and bulk vectors with iCellFCs
    for (uint i = 0; i < parHs[level]->intFC.kFC; i++)
        if (grid->isSparseIndexInFluidNodeIndicesBorder(iCellFccAll[i])) {
            iCellFccBorderVector.push_back(iCellFccAll[i]);
            iCellFcfBorderVector.push_back(iCellFcfAll[i]);
            xOffFCBorderVector.push_back(parHs[level]->offFC.xOffFC[i]);
            yOffFCBorderVector.push_back(parHs[level]->offFC.yOffFC[i]);
            zOffFCBorderVector.push_back(parHs[level]->offFC.zOffFC[i]);
        } else {
            iCellFccBulkVector.push_back(iCellFccAll[i]);
            iCellFcfBulkVector.push_back(iCellFcfAll[i]);
            xOffFCBulkVector.push_back(parHs[level]->offFC.xOffFC[i]);
            yOffFCBulkVector.push_back(parHs[level]->offFC.yOffFC[i]);
            zOffFCBulkVector.push_back(parHs[level]->offFC.zOffFC[i]);
        }

    // set new sizes and pointers
    parHs[level]->intFCBorder.ICellFCC = iCellFccAll;
    parHs[level]->intFCBorder.ICellFCF = iCellFcfAll;
    parHs[level]->intFCBorder.kFC = (uint)iCellFccBorderVector.size();
    parHs[level]->intFCBulk.kFC = (uint)iCellFccBulkVector.size();
    parHs[level]->intFCBulk.ICellFCC = iCellFccAll + parHs[level]->intFCBorder.kFC;
    parHs[level]->intFCBulk.ICellFCF = iCellFcfAll + parHs[level]->intFCBorder.kFC;
    parHs[level]->offFCBulk.xOffFC = parHs[level]->offFC.xOffFC + parHs[level]->intFCBorder.kFC;
    parHs[level]->offFCBulk.yOffFC = parHs[level]->offFC.yOffFC + parHs[level]->intFCBorder.kFC;
    parHs[level]->offFCBulk.zOffFC = parHs[level]->offFC.zOffFC + parHs[level]->intFCBorder.kFC;

    // copy the created vectors to the memory addresses of the old arrays
    // this is inefficient :(
    for (uint i = 0; i < (uint)iCellFccBorderVector.size(); i++) {
        iCellFccAll[i] = iCellFccBorderVector[i];
        iCellFcfAll[i] = iCellFcfBorderVector[i];
        parHs[level]->offFC.xOffFC[i] = xOffFCBorderVector[i];
        parHs[level]->offFC.yOffFC[i] = yOffFCBorderVector[i];
        parHs[level]->offFC.zOffFC[i] = zOffFCBorderVector[i];
    }
    for (uint i = 0; i < (uint)iCellFccBulkVector.size(); i++) {
        parHs[level]->intFCBulk.ICellFCC[i] = iCellFccBulkVector[i];
        parHs[level]->intFCBulk.ICellFCF[i] = iCellFcfBulkVector[i];
        parHs[level]->offFCBulk.xOffFC[i] = xOffFCBulkVector[i];
        parHs[level]->offFCBulk.yOffFC[i] = yOffFCBulkVector[i];
        parHs[level]->offFCBulk.zOffFC[i] = zOffFCBulkVector[i];
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
    parDs[level]->offCFBulk.xOffCF = parDs[level]->offCF.xOffCF + parDs[level]->intCFBorder.kCF;
    parDs[level]->offCFBulk.yOffCF = parDs[level]->offCF.yOffCF + parDs[level]->intCFBorder.kCF;
    parDs[level]->offCFBulk.zOffCF = parDs[level]->offCF.zOffCF + parDs[level]->intCFBorder.kCF;
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
    std::vector<real> xOffCFBorderVector;
    std::vector<real> yOffCFBorderVector;
    std::vector<real> zOffCFBorderVector;
    std::vector<real> xOffCFBulkVector;
    std::vector<real> yOffCFBulkVector;
    std::vector<real> zOffCFBulkVector;
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
            xOffCFBorderVector.push_back(parHs[level]->offCF.xOffCF[i]);
            yOffCFBorderVector.push_back(parHs[level]->offCF.yOffCF[i]);
            zOffCFBorderVector.push_back(parHs[level]->offCF.zOffCF[i]);
        } else {
            iCellCfcBulkVector.push_back(iCellCfcAll[i]);
            iCellCffBulkVector.push_back(iCellCffAll[i]);
            xOffCFBulkVector.push_back(parHs[level]->offCF.xOffCF[i]);
            yOffCFBulkVector.push_back(parHs[level]->offCF.yOffCF[i]);
            zOffCFBulkVector.push_back(parHs[level]->offCF.zOffCF[i]);
        }
    }

    // set new sizes and pointers
    parHs[level]->intCFBorder.ICellCFC = parHs[level]->intCF.ICellCFC;
    parHs[level]->intCFBorder.ICellCFF = parHs[level]->intCF.ICellCFF;
    parHs[level]->intCFBorder.kCF = (uint)iCellCfcBorderVector.size();
    parHs[level]->intCFBulk.kCF = (uint)iCellCfcBulkVector.size();
    parHs[level]->intCFBulk.ICellCFC = parHs[level]->intCF.ICellCFC + parHs[level]->intCFBorder.kCF;
    parHs[level]->intCFBulk.ICellCFF = parHs[level]->intCF.ICellCFF + parHs[level]->intCFBorder.kCF;
    parHs[level]->offCFBulk.xOffCF = parHs[level]->offCF.xOffCF + parHs[level]->intCFBorder.kCF;
    parHs[level]->offCFBulk.yOffCF = parHs[level]->offCF.yOffCF + parHs[level]->intCFBorder.kCF;
    parHs[level]->offCFBulk.zOffCF = parHs[level]->offCF.zOffCF + parHs[level]->intCFBorder.kCF;

    // copy the created vectors to the memory addresses of the old arrays
    // this is inefficient :(
    for (uint i = 0; i < (uint)iCellCfcBorderVector.size(); i++) {
        parHs[level]->intCFBorder.ICellCFC[i] = iCellCfcBorderVector[i];
        parHs[level]->intCFBorder.ICellCFF[i] = iCellCffBorderVector[i];
        parHs[level]->offCF.xOffCF[i] = xOffCFBorderVector[i];
        parHs[level]->offCF.yOffCF[i] = yOffCFBorderVector[i];
        parHs[level]->offCF.zOffCF[i] = zOffCFBorderVector[i];
    }
    for (uint i = 0; i < (uint)iCellCfcBulkVector.size(); i++) {
        parHs[level]->intCFBulk.ICellCFC[i] = iCellCfcBulkVector[i];
        parHs[level]->intCFBulk.ICellCFF[i] = iCellCffBulkVector[i];
        parHs[level]->offCFBulk.xOffCF[i] = xOffCFBulkVector[i];
        parHs[level]->offCFBulk.yOffCF[i] = yOffCFBulkVector[i];
        parHs[level]->offCFBulk.zOffCF[i] = zOffCFBulkVector[i];
    }
}
