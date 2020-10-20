#include "FineToCoarseNodeSetBlock3DConnector.h"
#include "BCProcessor.h"
#include "DataSet3D.h"

//////////////////////////////////////////////////////////////////////////
FineToCoarseNodeSetBlock3DConnector::FineToCoarseNodeSetBlock3DConnector(SPtr<Block3D> block,
                                                                         VectorTransmitterPtr sender,
                                                                         VectorTransmitterPtr receiver, int sendDir,
                                                                         InterpolationProcessorPtr iprocessor,
                                                                         CFconnectorType connType)
    : FineToCoarseBlock3DConnector(block, sender, receiver, sendDir, iprocessor, connType)
{
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseNodeSetBlock3DConnector::init()
{
    bMaxX1 = (int)FineToCoarseBlock3DConnector::block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1();
    bMaxX2 = (int)FineToCoarseBlock3DConnector::block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2();
    bMaxX3 = (int)FineToCoarseBlock3DConnector::block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3();

    minX1 = 0;
    minX2 = 0;
    minX3 = 0;
    maxX1 = bMaxX1 - 1;
    maxX2 = bMaxX2 - 1;
    maxX3 = bMaxX3 - 1;

    minOffX1 = 0;
    minOffX2 = 0;
    minOffX3 = 0;

    maxOffX1 = 0;
    maxOffX2 = 0;
    maxOffX3 = 0;

    if (Utilities::isEven(bMaxX1)) {
        minOffX1 = 0;
        maxOffX1 = 0;
    } else if (Utilities::isOdd(bMaxX1)) {
        minOffX1 = 1;
        maxOffX1 = -1;
    }

    if (Utilities::isEven(bMaxX2)) {
        minOffX2 = 0;
        maxOffX2 = 0;
    } else if (Utilities::isOdd(bMaxX2)) {
        minOffX2 = 1;
        maxOffX2 = -1;
    }

    if (Utilities::isEven(bMaxX3)) {
        minOffX3 = 0;
        maxOffX3 = 0;
    } else if (Utilities::isOdd(bMaxX3)) {
        minOffX3 = 1;
        maxOffX3 = -1;
    }

    // int       sendSize = 0;
    LBMReal initValue = -999.0;

    int sendDataPerNode = 27 /*f*/;
    int iCellSize       = 1; // size of interpolation cell

    findFCCells();
    findCFCells();

    //////////////////////////////////////////////////////
    // Debug
    //////////////////////////////////////////////////////
    //   if (FineToCoarseBlock3DConnector::block.lock()->getGlobalID() == 2183)
    //   {
    //      int test = 0;
    //   }

    if (FineToCoarseBlock3DConnector::sender)
        FineToCoarseBlock3DConnector::sender->getData().resize(iNodeSetSender.size() * iCellSize * sendDataPerNode,
                                                               initValue);
    else
        FineToCoarseBlock3DConnector::sender = VectorTransmitterPtr(new TbLocalTransmitter<CbVector<LBMReal>>());

    if (!FineToCoarseBlock3DConnector::receiver)
        FineToCoarseBlock3DConnector::receiver = VectorTransmitterPtr(new TbLocalTransmitter<CbVector<LBMReal>>());
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseNodeSetBlock3DConnector::findFCCells(int lMinX1, int lMinX2, int lMinX3, int lMaxX1, int lMaxX2,
                                                      int lMaxX3, INodeSet &inodes)
{
    //////////////////////////////////////////////////////
    // Debug
    //////////////////////////////////////////////////////
    //   if (FineToCoarseBlock3DConnector::block.lock()->getGlobalID() == 2183)
    //   {
    //      int test = 0;
    //   }

    int ix1, ix2, ix3;
    LBMReal x1off, x2off, x3off;

    SPtr<DistributionArray3D> fFrom =
        FineToCoarseBlock3DConnector::block.lock()->getKernel()->getDataSet()->getFdistributions();
    SPtr<BCArray3D> bcArray = FineToCoarseBlock3DConnector::block.lock()->getKernel()->getBCProcessor()->getBCArray();

    for (ix3 = lMinX3; ix3 <= lMaxX3; ix3 += 2) {
        for (ix2 = lMinX2; ix2 <= lMaxX2; ix2 += 2) {
            for (ix1 = lMinX1; ix1 <= lMaxX1; ix1 += 2) {
                D3Q27ICell icellC;

                int howManySolids =
                    FineToCoarseBlock3DConnector::iprocessor->iCellHowManySolids(bcArray, ix1, ix2, ix3);

                if (howManySolids == 0 || howManySolids == 8) {
                    x1off = 0.0;
                    x2off = 0.0;
                    x3off = 0.0;
                } else {
                    if (!iprocessor->findNeighborICell(bcArray, fFrom, icellC, bMaxX1, bMaxX2, bMaxX3, ix1, ix2, ix3,
                                                       x1off, x2off, x3off)) {
                        std::string err = "For " + FineToCoarseBlock3DConnector::block.lock()->toString() +
                                          " x1=" + UbSystem::toString(ix1) + ", x2=" + UbSystem::toString(ix2) +
                                          ", x3=" + UbSystem::toString(ix3) +
                                          " interpolation is not implemented for other direction" +
                                          " by using in: " + (std::string) typeid(*this).name() +
                                          " or maybe you have a solid on the block boundary";
                        UB_THROW(UbException(UB_EXARGS, err));
                    }
                }

                INodeVector inv;
                inv.push_back(ix1 + (int)x1off);
                inv.push_back(ix2 + (int)x2off);
                inv.push_back(ix3 + (int)x3off);
                inv.push_back((int)x1off);
                inv.push_back((int)x2off);
                inv.push_back((int)x3off);
                // inodes.insert(inv);
                inodes.push_back(inv);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////

void FineToCoarseNodeSetBlock3DConnector::findFCCells()
{
    using namespace D3Q27System;

    int lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3;
    int lMin2X1, lMin2X2, lMin2X3, lMax2X1, lMax2X2, lMax2X3;
    int lMin3X1, lMin3X2, lMin3X3, lMax3X1, lMax3X2, lMax3X3;

    // lMin1X1 = minX1+1; lMin1X2 = minX2+1; lMin1X3 = minX3+1;
    // lMax1X1 = maxX1-1; lMax1X2 = maxX2-1; lMax1X3 = maxX3-1;
    // getLocalMinMax(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3);

    // lMin2X1 = minX1+1; lMin2X2 = minX2+1; lMin2X3 = minX3+1;
    // lMax2X1 = maxX1-1; lMax2X2 = maxX2-1; lMax2X3 = maxX3-1;
    // getLocalMinMax(lMin2X1, lMin2X2, lMin2X3, lMax2X1, lMax2X2, lMax2X3);

    // lMin3X1 = minX1+1; lMin3X2 = minX2+1; lMin3X3 = minX3+1;
    // lMax3X1 = maxX1-1; lMax3X2 = maxX2-1; lMax3X3 = maxX3-1;
    // getLocalMinMax(lMin3X1, lMin3X2, lMin3X3, lMax3X1, lMax3X2, lMax3X3);

    switch (sendDir) {
            // faces
        case E:
        case W:
            if (sendDir == E) {
                lMin1X1 = maxX1 - 6;
                lMax1X1 = lMin1X1;
            } else if (sendDir == W) {
                lMin1X1 = 5;
                lMax1X1 = lMin1X1;
            }

            if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type00) {
                lMin1X2 = minX2;
                lMax1X2 = maxX2 + maxOffX2 - 1;
                lMin1X3 = minX3;
                lMax1X3 = maxX3 + maxOffX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type10) {
                lMin1X2 = minX2 + minOffX2;
                lMax1X2 = maxX2 - 1;
                lMin1X3 = minX3;
                lMax1X3 = maxX3 + maxOffX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type01) {
                lMin1X2 = minX2;
                lMax1X2 = maxX2 + maxOffX2 - 1;
                lMin1X3 = minX3 + minOffX3;
                lMax1X3 = maxX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type11) {
                lMin1X2 = minX2 + minOffX2;
                lMax1X2 = maxX2 - 1;
                lMin1X3 = minX3 + minOffX3;
                lMax1X3 = maxX3 - 1;
            }
            findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetSender);
            break;
        case N:
        case S:
            if (sendDir == N) {
                lMin1X2 = maxX2 - 6;
                lMax1X2 = lMin1X2;
            } else if (sendDir == S) {
                lMin1X2 = 5;
                lMax1X2 = lMin1X2;
            }

            if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type00) {
                lMin1X1 = minX1;
                lMax1X1 = maxX1 + maxOffX1 - 1;
                lMin1X3 = minX3;
                lMax1X3 = maxX3 + maxOffX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type10) {
                lMin1X1 = minX1 + minOffX1;
                lMax1X1 = maxX1 - 1;
                lMin1X3 = minX3;
                lMax1X3 = maxX3 + maxOffX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type01) {
                lMin1X1 = minX1;
                lMax1X1 = maxX1 + maxOffX1 - 1;
                lMin1X3 = minX3 + minOffX3;
                lMax1X3 = maxX3;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type11) {
                lMin1X1 = minX1 + minOffX1;
                lMax1X1 = maxX1 - 1;
                lMin1X3 = minX3 + minOffX3;
                lMax1X3 = maxX3 - 1;
            }
            findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetSender);
            break;
        case T:
        case B:
            if (sendDir == T) {
                lMin1X3 = maxX3 - 6;
                lMax1X3 = lMin1X3;
            } else if (sendDir == B) {
                lMin1X3 = 5;
                lMax1X3 = lMin1X3;
            }

            if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type00) {
                lMin1X1 = minX1;
                lMax1X1 = maxX1 + maxOffX1 - 1;
                lMin1X2 = minX2;
                lMax1X2 = maxX2 + maxOffX2 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type10) {
                lMin1X1 = minX1 + minOffX1;
                lMax1X1 = maxX1 - 1;
                lMin1X2 = minX2;
                lMax1X2 = maxX2 + maxOffX2 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type01) {
                lMin1X1 = minX1;
                lMax1X1 = maxX1 + maxOffX1 - 1;
                lMin1X2 = minX2 + minOffX2;
                lMax1X2 = maxX2 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type11) {
                lMin1X1 = minX1 + minOffX1;
                lMax1X1 = maxX1 - 1;
                lMin1X2 = minX2 + minOffX2;
                lMax1X2 = maxX2 - 1;
            }
            findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetSender);
            break;
            // edges
            // N-S-E-W
        case NE:
        case SW:
        case SE:
        case NW:
            if (sendDir == NE) {
                lMin1X1 = maxX1 - 6;
                lMax1X1 = lMin1X1 + 4;
                lMin1X2 = maxX2 - 6;
                lMax1X2 = lMin1X2;

                lMin2X1 = maxX1 - 6;
                lMax2X1 = lMin2X1;
                lMin2X2 = maxX2 - 6;
                lMax2X2 = lMin2X2 + 4;
            } else if (sendDir == SW) {
                lMin1X1 = 1;
                lMax1X1 = lMin1X1 + 4;
                lMin1X2 = 5;
                lMax1X2 = lMin1X2;

                lMin2X1 = 5;
                lMax2X1 = lMin2X1;
                lMin2X2 = 1;
                lMax2X2 = lMin2X2 + 4;
            } else if (sendDir == SE) {
                lMin1X1 = maxX1 - 6;
                lMax1X1 = lMin1X1 + 4;
                lMin1X2 = 5;
                lMax1X2 = lMin1X2;

                lMin2X1 = maxX1 - 6;
                lMax2X1 = lMin2X1;
                lMin2X2 = 1;
                lMax2X2 = lMin2X2 + 4;
            } else if (sendDir == NW) {
                lMin1X1 = 1;
                lMax1X1 = lMin1X1 + 4;
                lMin1X2 = maxX2 - 6;
                lMax1X2 = lMin1X2;

                lMin2X1 = 5;
                lMax2X1 = lMin2X1;
                lMin2X2 = maxX2 - 6;
                lMax2X2 = lMin2X2 + 4;
            }

            if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type00) {
                lMin1X3 = minX3;
                lMax1X3 = maxX3 + maxOffX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type10) {
                lMin1X3 = minX3 + minOffX3;
                lMax1X3 = maxX3 - 1;
            }

            findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetSender);
            findFCCells(lMin2X1, lMin2X2, lMin1X3, lMax2X1, lMax2X2, lMax1X3, iNodeSetSender);

            break;
            // T-B-E-W
        case TE:
        case BW:
        case BE:
        case TW:
            if (sendDir == TE) {
                lMin1X1 = maxX1 - 6;
                lMax1X1 = lMin1X1 + 4;
                lMin1X3 = maxX3 - 6;
                lMax1X3 = lMin1X3;

                lMin2X1 = maxX1 - 6;
                lMax2X1 = lMin2X1;
                lMin2X3 = maxX3 - 6;
                lMax2X3 = lMin2X3 + 4;
            } else if (sendDir == BW) {
                lMin1X1 = 1;
                lMax1X1 = lMin1X1 + 4;
                lMin1X3 = 5;
                lMax1X3 = lMin1X3;

                lMin2X1 = 5;
                lMax2X1 = lMin2X1;
                lMin2X3 = 1;
                lMax2X3 = lMin2X3 + 4;
            } else if (sendDir == BE) {
                lMin1X1 = maxX1 - 6;
                lMax1X1 = lMin1X1 + 4;
                lMin1X3 = 5;
                lMax1X3 = lMin1X3;

                lMin2X1 = maxX1 - 6;
                lMax2X1 = lMin2X1;
                lMin2X3 = 1;
                lMax2X3 = lMin2X3 + 4;
            } else if (sendDir == TW) {
                lMin1X1 = 1;
                lMax1X1 = lMin1X1 + 5;
                lMin1X3 = maxX3 - 6;
                lMax1X3 = lMin1X3;

                lMin2X1 = 5;
                lMax2X1 = lMin2X1;
                lMin2X3 = maxX3 - 6;
                lMax2X3 = lMin2X3 + 4;
            }

            if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type00) {
                lMin1X2 = minX2;
                lMax1X2 = maxX2 + maxOffX2 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type10) {
                lMin1X2 = minX2 + minOffX2;
                lMax1X2 = maxX2 - 1;
            }

            findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetSender);
            findFCCells(lMin2X1, lMin1X2, lMin2X3, lMax2X1, lMax1X2, lMax2X3, iNodeSetSender);
            break;
            // T-B-N-S
        case TN:
        case BS:
        case BN:
        case TS:
            if (sendDir == TN) {
                lMin1X2 = maxX2 - 6;
                lMax1X2 = lMin1X2 + 4;
                lMin1X3 = maxX3 - 6;
                lMax1X3 = lMin1X3;

                lMin2X2 = maxX2 - 6;
                lMax2X2 = lMin2X2;
                lMin2X3 = maxX3 - 6;
                lMax2X3 = lMin2X3 + 4;
            } else if (sendDir == BS) {
                lMin1X2 = 1;
                lMax1X2 = lMin1X2 + 4;
                lMin1X3 = 5;
                lMax1X3 = lMin1X3;

                lMin2X2 = 5;
                lMax2X2 = lMin2X2;
                lMin2X3 = 1;
                lMax2X3 = lMin2X3 + 4;
            } else if (sendDir == BN) {
                lMin1X2 = maxX2 - 6;
                lMax1X2 = lMin1X2 + 4;
                lMin1X3 = 5;
                lMax1X3 = lMin1X3;

                lMin2X2 = maxX2 - 6;
                lMax2X2 = lMin2X2;
                lMin2X3 = 1;
                lMax2X3 = lMin2X3 + 4;
            } else if (sendDir == TS) {
                lMin1X2 = 1;
                lMax1X2 = lMin1X2 + 4;
                lMin1X3 = maxX3 - 6;
                lMax1X3 = lMin1X3;

                lMin2X2 = 5;
                lMax2X2 = lMin2X2;
                lMin2X3 = maxX3 - 6;
                lMax2X3 = lMin2X3 + 4;
            }

            if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type00) {
                lMin1X1 = minX1;
                lMax1X1 = maxX1 + maxOffX1 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type10) {
                lMin1X1 = minX1 + minOffX1;
                lMax1X1 = maxX1 - 1;
            }

            findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetSender);
            findFCCells(lMin1X1, lMin2X2, lMin2X3, lMax1X1, lMax2X2, lMax2X3, iNodeSetSender);
            break;
            // corners
        case TNE:
        case TNW:
        case TSE:
        case TSW:
        case BNE:
        case BNW:
        case BSE:
        case BSW:
            if (sendDir == TNE) {
                lMin1X1 = maxX1 - 6;
                lMax1X1 = maxX1 - 6;
                lMin1X2 = maxX2 - 6;
                lMax1X2 = maxX2 - 2;
                lMin1X3 = maxX3 - 6;
                lMax1X3 = maxX3 - 2;

                lMin2X1 = maxX1 - 6;
                lMax2X1 = maxX1 - 2;
                lMin2X2 = maxX2 - 6;
                lMax2X2 = maxX2 - 6;
                lMin2X3 = maxX3 - 6;
                lMax2X3 = maxX3 - 1;

                lMin3X1 = maxX1 - 6;
                lMax3X1 = maxX1 - 2;
                lMin3X2 = maxX2 - 6;
                lMax3X2 = maxX2 - 2;
                lMin3X3 = maxX3 - 6;
                lMax3X3 = maxX3 - 5;
            } else if (sendDir == TNW) {
                lMin1X1 = 5;
                lMax1X1 = 5;
                lMin1X2 = maxX2 - 6;
                lMax1X2 = maxX2 - 2;
                lMin1X3 = maxX3 - 6;
                lMax1X3 = maxX3 - 2;

                lMin2X1 = 1;
                lMax2X1 = 5;
                lMin2X2 = maxX2 - 6;
                lMax2X2 = maxX2 - 6;
                lMin2X3 = maxX3 - 6;
                lMax2X3 = maxX3 - 2;

                lMin3X1 = 1;
                lMax3X1 = 5;
                lMin3X2 = maxX2 - 6;
                lMax3X2 = maxX2 - 2;
                lMin3X3 = maxX3 - 6;
                lMax3X3 = maxX3 - 6;
            } else if (sendDir == TSE) {
                lMin1X1 = maxX1 - 6;
                lMax1X1 = maxX1 - 6;
                lMin1X2 = 1;
                lMax1X2 = 5;
                lMin1X3 = maxX3 - 6;
                lMax1X3 = maxX3 - 2;

                lMin2X1 = maxX1 - 6;
                lMax2X1 = maxX1 - 2;
                lMin2X2 = 5;
                lMax2X2 = 5;
                lMin2X3 = maxX3 - 6;
                lMax2X3 = maxX3 - 2;

                lMin3X1 = maxX1 - 6;
                lMax3X1 = maxX1 - 2;
                lMin3X2 = 1;
                lMax3X2 = 5;
                lMin3X3 = maxX3 - 6;
                lMax3X3 = maxX3 - 6;
            } else if (sendDir == TSW) {
                lMin1X1 = 5;
                lMax1X1 = 5;
                lMin1X2 = 1;
                lMax1X2 = 5;
                lMin1X3 = maxX3 - 6;
                lMax1X3 = maxX3 - 2;

                lMin2X1 = 1;
                lMax2X1 = 5;
                lMin2X2 = 5;
                lMax2X2 = 5;
                lMin2X3 = maxX3 - 6;
                lMax2X3 = maxX3 - 2;

                lMin3X1 = 1;
                lMax3X1 = 5;
                lMin3X2 = 1;
                lMax3X2 = 5;
                lMin3X3 = maxX3 - 6;
                lMax3X3 = maxX3 - 6;
            } else if (sendDir == BNE) {
                lMin1X1 = maxX1 - 6;
                lMax1X1 = maxX1 - 6;
                lMin1X2 = maxX2 - 6;
                lMax1X2 = maxX2 - 2;
                lMin1X3 = 1;
                lMax1X3 = 5;

                lMin2X1 = maxX1 - 6;
                lMax2X1 = maxX1 - 2;
                lMin2X2 = maxX2 - 6;
                lMax2X2 = maxX2 - 6;
                lMin2X3 = 1;
                lMax2X3 = 5;

                lMin3X1 = maxX1 - 6;
                lMax3X1 = maxX1 - 2;
                lMin3X2 = maxX2 - 6;
                lMax3X2 = maxX2 - 2;
                lMin3X3 = 5;
                lMax3X3 = 5;
            } else if (sendDir == BNW) {
                lMin1X1 = 5;
                lMax1X1 = 5;
                lMin1X2 = maxX2 - 6;
                lMax1X2 = maxX2 - 2;
                lMin1X3 = 1;
                lMax1X3 = 5;

                lMin2X1 = 1;
                lMax2X1 = 5;
                lMin2X2 = maxX2 - 6;
                lMax2X2 = maxX2 - 6;
                lMin2X3 = 1;
                lMax2X3 = 5;

                lMin3X1 = 1;
                lMax3X1 = 5;
                lMin3X2 = maxX2 - 6;
                lMax3X2 = maxX2 - 2;
                lMin3X3 = 5;
                lMax3X3 = 5;
            } else if (sendDir == BSE) {
                lMin1X1 = maxX1 - 6;
                lMax1X1 = maxX1 - 6;
                lMin1X2 = 1;
                lMax1X2 = 5;
                lMin1X3 = 1;
                lMax1X3 = 5;

                lMin2X1 = maxX1 - 6;
                lMax2X1 = maxX1 - 2;
                lMin2X2 = 5;
                lMax2X2 = 5;
                lMin2X3 = 1;
                lMax2X3 = 5;

                lMin3X1 = maxX1 - 5;
                lMax3X1 = maxX1 - 2;
                lMin3X2 = 1;
                lMax3X2 = 5;
                lMin3X3 = 5;
                lMax3X3 = 5;
            } else if (sendDir == BSW) {
                lMin1X1 = 5;
                lMax1X1 = 5;
                lMin1X2 = 1;
                lMax1X2 = 5;
                lMin1X3 = 1;
                lMax1X3 = 5;

                lMin2X1 = 1;
                lMax2X1 = 5;
                lMin2X2 = 5;
                lMax2X2 = 5;
                lMin2X3 = 1;
                lMax2X3 = 5;

                lMin3X1 = 1;
                lMax3X1 = 5;
                lMin3X2 = 1;
                lMax3X2 = 5;
                lMin3X3 = 5;
                lMax3X3 = 5;
            }
            findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetSender);
            findFCCells(lMin2X1, lMin2X2, lMin2X3, lMax2X1, lMax2X2, lMax2X3, iNodeSetSender);
            findFCCells(lMin3X1, lMin3X2, lMin3X3, lMax3X1, lMax3X2, lMax3X3, iNodeSetSender);
            break;
    }
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseNodeSetBlock3DConnector::fillSendVectors()
{
    using namespace D3Q27System;

    SPtr<DistributionArray3D> fFrom = block.lock()->getKernel()->getDataSet()->getFdistributions();

    int index = 0;

    vector_type &data = this->sender->getData();

    for (INodeVector inode : iNodeSetSender) {
        LBMReal icellC[27];
        D3Q27ICell icellF;
        iprocessor->readICell(fFrom, icellF, inode[0], inode[1], inode[2]);
        iprocessor->interpolateFineToCoarse(icellF, icellC, inode[3], inode[4], inode[5]);
        writeICellCtoData(data, index, icellC);
    }
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseNodeSetBlock3DConnector::readICellFfromData(vector_type &data, int &index, D3Q27ICell &icellF)
{
    readNodeFromVector(data, index, icellF.BSW);
    readNodeFromVector(data, index, icellF.BSE);
    readNodeFromVector(data, index, icellF.BNW);
    readNodeFromVector(data, index, icellF.BNE);
    readNodeFromVector(data, index, icellF.TSW);
    readNodeFromVector(data, index, icellF.TSE);
    readNodeFromVector(data, index, icellF.TNW);
    readNodeFromVector(data, index, icellF.TNE);
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseNodeSetBlock3DConnector::readNodeFromVector(vector_type &data, int &index, LBMReal *inode)
{
    for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF + 1; i++) {
        inode[i] = data[index++];
    }
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseNodeSetBlock3DConnector::writeICellCtoData(vector_type &data, int &index, LBMReal *icellC)
{
    for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF + 1; i++) {
        data[index++] = icellC[i];
    }
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseNodeSetBlock3DConnector::findCFCells(int lMinX1, int lMinX2, int lMinX3, int lMaxX1, int lMaxX2,
                                                      int lMaxX3, INodeSet &inodes)
{
    int ix1, ix2, ix3;

    for (ix3 = lMinX3; ix3 <= lMaxX3; ix3 += 2) {
        for (ix2 = lMinX2; ix2 <= lMaxX2; ix2 += 2) {
            for (ix1 = lMinX1; ix1 <= lMaxX1; ix1 += 2) {
                INodeVector inv;
                inv.push_back(ix1);
                inv.push_back(ix2);
                inv.push_back(ix3);
                // inodes.insert(inv);
                inodes.push_back(inv);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseNodeSetBlock3DConnector::findCFCells()
{
    using namespace D3Q27System;

    int lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3;

    //////////////////////////////////////////////////////
    // Debug
    //////////////////////////////////////////////////////
    //   if (block.lock()->getGlobalID() == 2183)
    //   {
    //      int test = 0;
    //   }

    switch (sendDir) {
            // faces
        case E:
        case W:
            if (sendDir == E) {
                lMinX1 = maxX1 - 3;
                lMaxX1 = lMinX1;
            } else if (sendDir == W) {
                lMinX1 = 2;
                lMaxX1 = lMinX1;
            }

            if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type00) {
                lMinX2 = minX2;
                lMaxX2 = maxX2 + maxOffX2 - 1;
                lMinX3 = minX3;
                lMaxX3 = maxX3 + maxOffX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type10) {
                lMinX2 = minX2 + minOffX2;
                lMaxX2 = maxX2 - 1;
                lMinX3 = minX3;
                lMaxX3 = maxX3 + maxOffX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type01) {
                lMinX2 = minX2;
                lMaxX2 = maxX2 + maxOffX2 - 1;
                lMinX3 = minX3 + minOffX3;
                lMaxX3 = maxX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type11) {
                lMinX2 = minX2 + minOffX2;
                lMaxX2 = maxX2 - 1;
                lMinX3 = minX3 + minOffX3;
                lMaxX3 = maxX3 - 1;
            }
            findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetReceiver);
            break;
        case N:
        case S:
            if (sendDir == N) {
                lMinX2 = maxX2 - 3;
                lMaxX2 = lMinX2;
            } else if (sendDir == S) {
                lMinX2 = 2;
                lMaxX2 = lMinX2;
            }

            if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type00) {
                lMinX1 = minX1;
                lMaxX1 = maxX1 + maxOffX1 - 1;
                lMinX3 = minX3;
                lMaxX3 = maxX3 + maxOffX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type10) {
                lMinX1 = minX1 + minOffX1;
                lMaxX1 = maxX1;
                lMinX3 = minX3;
                lMaxX3 = maxX3 + maxOffX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type01) {
                lMinX1 = minX1;
                lMaxX1 = maxX1 + maxOffX1 - 1;
                lMinX3 = minX3 + minOffX3;
                lMaxX3 = maxX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type11) {
                lMinX1 = minX1 + minOffX1;
                lMaxX1 = maxX1 - 1;
                lMinX3 = minX3 + minOffX3;
                lMaxX3 = maxX3 - 1;
            }
            findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetReceiver);
            break;
        case T:
        case B:
            if (sendDir == T) {
                lMinX3 = maxX3 - 3;
                lMaxX3 = lMinX3;
            } else if (sendDir == B) {
                lMinX3 = 2;
                lMaxX3 = lMinX3;
            }

            if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type00) {
                lMinX1 = minX1;
                lMaxX1 = maxX1 + maxOffX1 - 1;
                lMinX2 = minX2;
                lMaxX2 = maxX2 + maxOffX2 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type10) {
                lMinX1 = minX1 + minOffX1;
                lMaxX1 = maxX1 - 1;
                lMinX2 = minX2;
                lMaxX2 = maxX2 + maxOffX2 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type01) {
                lMinX1 = minX1;
                lMaxX1 = maxX1 + maxOffX1 - 1;
                lMinX2 = minX2 + minOffX2;
                lMaxX2 = maxX2 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type11) {
                lMinX1 = minX1 + minOffX1;
                lMaxX1 = maxX1 - 1;
                lMinX2 = minX2 + minOffX2;
                lMaxX2 = maxX2 - 1;
            }
            findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetReceiver);
            break;

            // edges
            // N-S-E-W
        case NE:
        case SW:
        case SE:
        case NW:
            if (sendDir == NE) {
                lMinX1 = maxX1 - 3;
                lMaxX1 = lMinX1 + 2;
                lMinX2 = maxX2 - 3;
                lMaxX2 = lMinX2 + 2;
            } else if (sendDir == SW) {
                lMinX1 = 0;
                lMaxX1 = lMinX1 + 3;
                lMinX2 = 0;
                lMaxX2 = lMinX2 + 3;
            } else if (sendDir == SE) {
                lMinX1 = maxX1 - 3;
                lMaxX1 = lMinX1 + 2;
                lMinX2 = 0;
                lMaxX2 = lMinX2 + 2;
            } else if (sendDir == NW) {
                lMinX1 = 0;
                lMaxX1 = lMinX1 + 2;
                lMinX2 = maxX2 - 3;
                lMaxX2 = lMinX2 + 2;
            }

            if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type00) {
                lMinX3 = minX3;
                lMaxX3 = maxX3 + maxOffX3 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type10) {
                lMinX3 = minX3 + minOffX3;
                lMaxX3 = maxX3 - 1;
            }

            findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetReceiver);
            break;
            // T-B-E-W
        case TE:
        case BW:
        case BE:
        case TW:
            if (sendDir == TE) {
                lMinX1 = maxX1 - 3;
                lMaxX1 = lMinX1 + 2;
                lMinX3 = maxX3 - 3;
                lMaxX3 = lMinX3 + 2;
            } else if (sendDir == BW) {
                lMinX1 = 0;
                lMaxX1 = lMinX1 + 2;
                lMinX3 = 0;
                lMaxX3 = lMinX3 + 2;
            } else if (sendDir == BE) {
                lMinX1 = maxX1 - 3;
                lMaxX1 = lMinX1 + 2;
                lMinX3 = 0;
                lMaxX3 = lMinX3 + 2;
            } else if (sendDir == TW) {
                lMinX1 = 0;
                lMaxX1 = lMinX1 + 2;
                lMinX3 = maxX3 - 3;
                lMaxX3 = lMinX3 + 2;
            }

            if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type00) {
                lMinX2 = minX2;
                lMaxX2 = maxX2 + maxOffX2 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type10) {
                lMinX2 = minX2 + minOffX2;
                lMaxX2 = maxX2 - 1;
            }

            findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetReceiver);
            break;
            // T-B-N-S
        case TN:
        case BS:
        case BN:
        case TS:
            if (sendDir == TN) {
                lMinX2 = maxX2 - 3;
                lMaxX2 = lMinX2 + 2;
                lMinX3 = maxX3 - 3;
                lMaxX3 = lMinX3 + 2;
            } else if (sendDir == BS) {
                lMinX2 = 0;
                lMaxX2 = lMinX2 + 2;
                lMinX3 = 0;
                lMaxX3 = lMinX3 + 2;
            } else if (sendDir == BN) {
                lMinX2 = maxX2 - 3;
                lMaxX2 = lMinX2 + 2;
                lMinX3 = 0;
                lMaxX3 = lMinX3 + 2;
            } else if (sendDir == TS) {
                lMinX2 = 0;
                lMaxX2 = lMinX2 + 2;
                lMinX3 = maxX3 - 3;
                lMaxX3 = lMinX3 + 2;
            }

            if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type00) {
                lMinX1 = minX1;
                lMaxX1 = maxX1 + maxOffX1 - 1;
            } else if (FineToCoarseBlock3DConnector::connType == FineToCoarseBlock3DConnector::Type10) {
                lMinX1 = minX1 + minOffX1;
                lMaxX1 = maxX1 - 1;
            }

            findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetReceiver);
            break;
            // corners
        case TNE:
        case TNW:
        case TSE:
        case TSW:
        case BNE:
        case BNW:
        case BSE:
        case BSW:
            if (sendDir == TNE) {
                lMinX1 = maxX1 - 3;
                lMaxX1 = maxX1 - 1;
                lMinX2 = maxX2 - 3;
                lMaxX2 = maxX2 - 1;
                lMinX3 = maxX3 - 3;
                lMaxX3 = maxX3 - 1;
            } else if (sendDir == TNW) {
                lMinX1 = 0;
                lMaxX1 = 2;
                lMinX2 = maxX2 - 3;
                lMaxX2 = maxX2 - 1;
                lMinX3 = maxX3 - 3;
                lMaxX3 = maxX3 - 1;
            } else if (sendDir == TSE) {
                lMinX1 = maxX1 - 3;
                lMaxX1 = maxX1 - 1;
                lMinX2 = 0;
                lMaxX2 = 2;
                lMinX3 = maxX3 - 3;
                lMaxX3 = maxX3 - 1;
            } else if (sendDir == TSW) {
                lMinX1 = 0;
                lMaxX1 = 2;
                lMinX2 = 0;
                lMaxX2 = 2;
                lMinX3 = maxX3 - 3;
                lMaxX3 = maxX3 - 1;
            } else if (sendDir == BNE) {
                lMinX1 = maxX1 - 3;
                lMaxX1 = maxX1 - 1;
                lMinX2 = maxX2 - 3;
                lMaxX2 = maxX2 - 1;
                lMinX3 = 0;
                lMaxX3 = 2;
            } else if (sendDir == BNW) {
                lMinX1 = 0;
                lMaxX1 = 2;
                lMinX2 = maxX2 - 3;
                lMaxX2 = maxX2 - 1;
                lMinX3 = 0;
                lMaxX3 = 2;
            } else if (sendDir == BSE) {
                lMinX1 = maxX1 - 3;
                lMaxX1 = maxX1 - 1;
                lMinX2 = 0;
                lMaxX2 = 2;
                lMinX3 = 0;
                lMaxX3 = 2;
            } else if (sendDir == BSW) {
                lMinX1 = 0;
                lMaxX1 = 2;
                lMinX2 = 0;
                lMaxX2 = 2;
                lMinX3 = 0;
                lMaxX3 = 2;
            }
            findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetReceiver);
            break;
    }
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseNodeSetBlock3DConnector::distributeReceiveVectors()
{
    using namespace D3Q27System;

    SPtr<DistributionArray3D> fTo =
        FineToCoarseBlock3DConnector::block.lock()->getKernel()->getDataSet()->getFdistributions();

    int index = 0;

    vector_type &data = this->receiver->getData();

    for (INodeVector inode : iNodeSetReceiver) {
        D3Q27ICell icellF;
        this->readICellFfromData(data, index, icellF);
        iprocessor->writeICellInv(fTo, icellF, inode[0], inode[1], inode[2]);
    }
}
//////////////////////////////////////////////////////////////////////////
//
// void FineToCoarseNodeSetBlock3DConnector::getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2,
// int& maxX3)
//{
//   using namespace D3Q27System;
//   int TminX1 = minX1; int TminX2 = minX2; int TminX3 = minX3; int TmaxX1 = maxX1; int TmaxX2 = maxX2; int TmaxX3 =
//   maxX3;
//
//   if (block.lock()->hasInterpolationFlagFC(E))
//   {
//      if (maxX1==TmaxX1) maxX1 -= 2;
//   }
//   if (block.lock()->hasInterpolationFlagFC(W))
//   {
//      if (minX1==TminX1) minX1 += 4;
//   }
//   if (block.lock()->hasInterpolationFlagFC(N))
//   {
//      if (maxX2==TmaxX2) maxX2 -= 2;
//   }
//   if (block.lock()->hasInterpolationFlagFC(S))
//   {
//      if (minX2==TminX2) minX2 += 4;
//   }
//   if (block.lock()->hasInterpolationFlagFC(T))
//   {
//      if (maxX3==TmaxX3) maxX3 -= 2;
//   }
//   if (block.lock()->hasInterpolationFlagFC(B))
//   {
//      if (minX3==TminX3) minX3 += 4;
//   }
//
//   ////////////
//   /////E-W-N-S
//   if (block.lock()->hasInterpolationFlagFC(NE)&& !block.lock()->hasInterpolationFlagFC(N) &&
//   !block.lock()->hasInterpolationFlagFC(E))
//   {
//      if (maxX1==TmaxX1) maxX1 -= 2;
//      if (maxX2==TmaxX2) maxX2 -= 2;
//   }
//   if (block.lock()->hasInterpolationFlagFC(SW)&& !block.lock()->hasInterpolationFlagFC(W) &&
//   !block.lock()->hasInterpolationFlagFC(S))
//   {
//      if (minX1==TminX1) minX1 += 4;
//      if (minX2==TminX2) minX2 += 4;
//   }
//   if (block.lock()->hasInterpolationFlagFC(SE)&& !block.lock()->hasInterpolationFlagFC(E) &&
//   !block.lock()->hasInterpolationFlagFC(S))
//   {
//      if (maxX1==TmaxX1) maxX1 -= 2;
//      if (minX2==TminX2) minX2 += 4;
//   }
//   if (block.lock()->hasInterpolationFlagFC(NW)&& !block.lock()->hasInterpolationFlagFC(N) &&
//   !block.lock()->hasInterpolationFlagFC(W))
//   {
//      if (minX1==TminX1) minX1 += 4;
//      if (maxX2==TmaxX2) maxX2 -= 2;
//   }
//
//   //////T-B-E-W
//   if (block.lock()->hasInterpolationFlagFC(TE) && !block.lock()->hasInterpolationFlagFC(E) &&
//   !block.lock()->hasInterpolationFlagFC(T))
//   {
//      if (maxX1==TmaxX1) maxX1 -= 2;
//      if (maxX3==TmaxX3) maxX3 -= 2;
//   }
//   if (block.lock()->hasInterpolationFlagFC(BW)&& !block.lock()->hasInterpolationFlagFC(W) &&
//   !block.lock()->hasInterpolationFlagFC(B))
//   {
//      if (minX1==TminX1) minX1 += 4;
//      if (minX3==TminX3) minX3 += 4;
//   }
//   if (block.lock()->hasInterpolationFlagFC(BE)&& !block.lock()->hasInterpolationFlagFC(E) &&
//   !block.lock()->hasInterpolationFlagFC(B))
//   {
//      if (maxX1==TmaxX1) maxX1 -= 2;
//      if (minX3==TminX3) minX3 += 4;
//   }
//   if (block.lock()->hasInterpolationFlagFC(TW)&& !block.lock()->hasInterpolationFlagFC(W) &&
//   !block.lock()->hasInterpolationFlagFC(T))
//   {
//      if (minX1==TminX1) minX1 += 4;
//      if (maxX3==TmaxX3) maxX3 -= 2;
//   }
//
//
//   ////T-B-N-S
//   if (block.lock()->hasInterpolationFlagFC(TN)&& !block.lock()->hasInterpolationFlagFC(N) &&
//   !block.lock()->hasInterpolationFlagFC(T))
//   {
//      if (maxX2==TmaxX2) maxX2 -= 2;
//      if (maxX3==TmaxX3) maxX3 -= 2;
//   }
//   if (block.lock()->hasInterpolationFlagFC(BS)&& !block.lock()->hasInterpolationFlagFC(S) &&
//   !block.lock()->hasInterpolationFlagFC(B))
//   {
//      if (minX2==TminX2) minX2 += 4;
//      if (minX3==TminX3) minX3 += 4;
//   }
//   if (block.lock()->hasInterpolationFlagFC(BN)&& !block.lock()->hasInterpolationFlagFC(N) &&
//   !block.lock()->hasInterpolationFlagFC(B))
//   {
//      if (maxX2==TmaxX2) maxX2 -= 2;
//      if (minX3==TminX3) minX3 += 4;
//   }
//   if (block.lock()->hasInterpolationFlagFC(TS) && !block.lock()->hasInterpolationFlagFC(S) &&
//   !block.lock()->hasInterpolationFlagFC(T))
//   {
//      if (minX2==TminX2) minX2 += 4;
//      if (maxX3==TmaxX3) maxX3 -= 2;
//   }
//}
//////////////////////////////////////////////////////////////////////////
