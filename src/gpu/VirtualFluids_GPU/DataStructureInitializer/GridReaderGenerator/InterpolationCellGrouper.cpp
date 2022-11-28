#include "InterpolationCellGrouper.h"

#include "Parameter/Parameter.h"
#include <GridGenerator/grid/Grid.h>
#include <GridGenerator/grid/GridBuilder/GridBuilder.h>

InterpolationCellGrouper::InterpolationCellGrouper(std::shared_ptr<Parameter> para,
                                                           std::shared_ptr<GridBuilder> builder)
    : para(para), builder(builder)
{
}

void InterpolationCellGrouper::splitFineToCoarseIntoBorderAndBulk(uint level)
{
    this->reorderFineToCoarseIntoBorderAndBulk(level);

    para->getParD(level)->intFCBorder.kFC      = para->getParH(level)->intFCBorder.kFC;
    para->getParD(level)->intFCBulk.kFC        = para->getParH(level)->intFCBulk.kFC;
    para->getParD(level)->intFCBorder.ICellFCC = para->getParD(level)->intFC.ICellFCC;
    para->getParD(level)->intFCBulk.ICellFCC =
        para->getParD(level)->intFCBorder.ICellFCC + para->getParD(level)->intFCBorder.kFC;
    para->getParD(level)->intFCBorder.ICellFCF = para->getParD(level)->intFC.ICellFCF;
    para->getParD(level)->intFCBulk.ICellFCF =
        para->getParD(level)->intFCBorder.ICellFCF + para->getParD(level)->intFCBorder.kFC;
    para->getParD(level)->offFCBulk.xOffFC = para->getParD(level)->offFC.xOffFC + para->getParD(level)->intFCBorder.kFC;
    para->getParD(level)->offFCBulk.yOffFC = para->getParD(level)->offFC.yOffFC + para->getParD(level)->intFCBorder.kFC;
    para->getParD(level)->offFCBulk.zOffFC = para->getParD(level)->offFC.zOffFC + para->getParD(level)->intFCBorder.kFC;
}

void InterpolationCellGrouper::reorderFineToCoarseIntoBorderAndBulk(int level)
{
    // create some local variables for better readability
    uint *iCellFccAll = para->getParH(level)->intFC.ICellFCC;
    uint *iCellFcfAll = para->getParH(level)->intFC.ICellFCF;
    auto grid         = this->builder->getGrid((uint)level);

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
    for (uint i = 0; i < para->getParH(level)->intFC.kFC; i++)
        if (grid->isSparseIndexInFluidNodeIndicesBorder(iCellFccAll[i])) {
            iCellFccBorderVector.push_back(iCellFccAll[i]);
            iCellFcfBorderVector.push_back(iCellFcfAll[i]);
            xOffFCBorderVector.push_back(para->getParH(level)->offFC.xOffFC[i]);
            yOffFCBorderVector.push_back(para->getParH(level)->offFC.yOffFC[i]);
            zOffFCBorderVector.push_back(para->getParH(level)->offFC.zOffFC[i]);
        } else {
            iCellFccBulkVector.push_back(iCellFccAll[i]);
            iCellFcfBulkVector.push_back(iCellFcfAll[i]);
            xOffFCBulkVector.push_back(para->getParH(level)->offFC.xOffFC[i]);
            yOffFCBulkVector.push_back(para->getParH(level)->offFC.yOffFC[i]);
            zOffFCBulkVector.push_back(para->getParH(level)->offFC.zOffFC[i]);
        }

    // set new sizes and pointers
    para->getParH(level)->intFCBorder.ICellFCC = iCellFccAll;
    para->getParH(level)->intFCBorder.ICellFCF = iCellFcfAll;
    para->getParH(level)->intFCBorder.kFC      = (uint)iCellFccBorderVector.size();
    para->getParH(level)->intFCBulk.kFC        = (uint)iCellFccBulkVector.size();
    para->getParH(level)->intFCBulk.ICellFCC   = iCellFccAll + para->getParH(level)->intFCBorder.kFC;
    para->getParH(level)->intFCBulk.ICellFCF   = iCellFcfAll + para->getParH(level)->intFCBorder.kFC;
    para->getParH(level)->offFCBulk.xOffFC = para->getParH(level)->offFC.xOffFC + para->getParH(level)->intFCBorder.kFC;
    para->getParH(level)->offFCBulk.yOffFC = para->getParH(level)->offFC.yOffFC + para->getParH(level)->intFCBorder.kFC;
    para->getParH(level)->offFCBulk.zOffFC = para->getParH(level)->offFC.zOffFC + para->getParH(level)->intFCBorder.kFC;


    // copy the created vectors to the memory addresses of the old arrays
    // this is inefficient :(
    for (uint i = 0; i < (uint)iCellFccBorderVector.size(); i++) {
        iCellFccAll[i] = iCellFccBorderVector[i];
        iCellFcfAll[i] = iCellFcfBorderVector[i];
        para->getParH(level)->offFC.xOffFC[i] = xOffFCBorderVector[i];
        para->getParH(level)->offFC.yOffFC[i] = yOffFCBorderVector[i];
        para->getParH(level)->offFC.zOffFC[i] = zOffFCBorderVector[i];
    }
    for (uint i = 0; i < (uint)iCellFccBulkVector.size(); i++) {
        para->getParH(level)->intFCBulk.ICellFCC[i] = iCellFccBulkVector[i];
        para->getParH(level)->intFCBulk.ICellFCF[i] = iCellFcfBulkVector[i];
        para->getParH(level)->offFCBulk.xOffFC[i]   = xOffFCBulkVector[i];
        para->getParH(level)->offFCBulk.yOffFC[i]   = yOffFCBulkVector[i];
        para->getParH(level)->offFCBulk.zOffFC[i]   = zOffFCBulkVector[i];
    }
}

void InterpolationCellGrouper::splitCoarseToFineIntoBorderAndBulk(uint level)
{
    this->reorderCoarseToFineIntoBorderAndBulk(level);

    para->getParD(level)->intCFBorder.kCF      = para->getParH(level)->intCFBorder.kCF;
    para->getParD(level)->intCFBulk.kCF        = para->getParH(level)->intCFBulk.kCF;
    para->getParD(level)->intCFBorder.ICellCFC = para->getParD(level)->intCF.ICellCFC;
    para->getParD(level)->intCFBulk.ICellCFC =
        para->getParD(level)->intCFBorder.ICellCFC + para->getParD(level)->intCFBorder.kCF;
    para->getParD(level)->intCFBorder.ICellCFF = para->getParD(level)->intCF.ICellCFF;
    para->getParD(level)->intCFBulk.ICellCFF =
        para->getParD(level)->intCFBorder.ICellCFF + para->getParD(level)->intCFBorder.kCF;
    para->getParD(level)->offCFBulk.xOffCF = para->getParD(level)->offCF.xOffCF + para->getParD(level)->intCFBorder.kCF;
    para->getParD(level)->offCFBulk.yOffCF = para->getParD(level)->offCF.yOffCF + para->getParD(level)->intCFBorder.kCF;
    para->getParD(level)->offCFBulk.zOffCF = para->getParD(level)->offCF.zOffCF + para->getParD(level)->intCFBorder.kCF;
}

void InterpolationCellGrouper::reorderCoarseToFineIntoBorderAndBulk(int level)
{
    // create some local variables for better readability
    uint *iCellCfcAll  = para->getParH(level)->intCF.ICellCFC;
    uint *iCellCffAll  = para->getParH(level)->intCF.ICellCFF;
    uint *neighborX = this->para->getParH(level)->neighborX;
    uint *neighborY = this->para->getParH(level)->neighborY;
    uint *neighborZ = this->para->getParH(level)->neighborZ;
    auto grid          = this->builder->getGrid((uint)level);

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
    for (uint i = 0; i < para->getParH(level)->intCF.kCF; i++) {
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
            xOffCFBorderVector.push_back(para->getParH(level)->offCF.xOffCF[i]);
            yOffCFBorderVector.push_back(para->getParH(level)->offCF.yOffCF[i]);
            zOffCFBorderVector.push_back(para->getParH(level)->offCF.zOffCF[i]);
        } else {
            iCellCfcBulkVector.push_back(iCellCfcAll[i]);
            iCellCffBulkVector.push_back(iCellCffAll[i]);
            xOffCFBulkVector.push_back(para->getParH(level)->offCF.xOffCF[i]);
            yOffCFBulkVector.push_back(para->getParH(level)->offCF.yOffCF[i]);
            zOffCFBulkVector.push_back(para->getParH(level)->offCF.zOffCF[i]);
        }
    }

    // set new sizes and pointers
    para->getParH(level)->intCFBorder.ICellCFC = para->getParH(level)->intCF.ICellCFC;
    para->getParH(level)->intCFBorder.ICellCFF = para->getParH(level)->intCF.ICellCFF;
    para->getParH(level)->intCFBorder.kCF      = (uint)iCellCfcBorderVector.size();
    para->getParH(level)->intCFBulk.kCF        = (uint)iCellCfcBulkVector.size();
    para->getParH(level)->intCFBulk.ICellCFC =
        para->getParH(level)->intCF.ICellCFC + para->getParH(level)->intCFBorder.kCF;
    para->getParH(level)->intCFBulk.ICellCFF =
        para->getParH(level)->intCF.ICellCFF + para->getParH(level)->intCFBorder.kCF;
    para->getParH(level)->offCFBulk.xOffCF = para->getParH(level)->offCF.xOffCF + para->getParH(level)->intCFBorder.kCF;
    para->getParH(level)->offCFBulk.yOffCF = para->getParH(level)->offCF.yOffCF + para->getParH(level)->intCFBorder.kCF;
    para->getParH(level)->offCFBulk.zOffCF = para->getParH(level)->offCF.zOffCF + para->getParH(level)->intCFBorder.kCF;

    // copy the created vectors to the memory addresses of the old arrays
    // this is inefficient :(
    for (uint i = 0; i < (uint)iCellCfcBorderVector.size(); i++) {
        para->getParH(level)->intCFBorder.ICellCFC[i] = iCellCfcBorderVector[i];
        para->getParH(level)->intCFBorder.ICellCFF[i] = iCellCffBorderVector[i];
        para->getParH(level)->offCF.xOffCF[i]         = xOffCFBorderVector[i];
        para->getParH(level)->offCF.yOffCF[i]         = yOffCFBorderVector[i];
        para->getParH(level)->offCF.zOffCF[i]         = zOffCFBorderVector[i];
    }
    for (uint i = 0; i < (uint)iCellCfcBulkVector.size(); i++) {
        para->getParH(level)->intCFBulk.ICellCFC[i] = iCellCfcBulkVector[i];
        para->getParH(level)->intCFBulk.ICellCFF[i] = iCellCffBulkVector[i];
        para->getParH(level)->offCFBulk.xOffCF[i]   = xOffCFBulkVector[i];
        para->getParH(level)->offCFBulk.yOffCF[i]   = yOffCFBulkVector[i];
        para->getParH(level)->offCFBulk.zOffCF[i]   = zOffCFBulkVector[i];
    }
}
