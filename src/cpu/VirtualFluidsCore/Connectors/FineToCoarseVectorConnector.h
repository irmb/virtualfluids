//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file FineToCoarseVectorConnector.h
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef FineToCoarseVectorConnector_H
#define FineToCoarseVectorConnector_H

#include <vector>

#include "Block3D.h"
#include "Block3DConnector.h"
#include "D3Q27System.h"
#include "Grid3D.h"
#include "Interpolator.h"
#include "LBMKernel.h"
#include "MathUtil.hpp"
#include "basics/transmitter/TbTransmitter.h"
#include <PointerDefinitions.h>
#include "basics/constants/NumericConstants.h"

#include "BCSet.h"
#include "DataSet3D.h"

class Block3D;

enum CFconnectorType { EvenOddNW, EvenEvenSW, OddEvenSE, OddOddNE };

//! \brief Interpolation from fine level to coarse.
//! \details Data are copied into a vector (this is located in the transmitter) the vector is transmitted via transmitter.
//!  Transmitter can be a local, MPI, RCG, CTL or whatever be a transmitter that is derived from the transmitter ;-)
template <typename VectorTransmitter>
class FineToCoarseVectorConnector : public Block3DConnector
{
public:
protected:
    using vector_type          = typename VectorTransmitter::value_type;
    using VectorTransmitterPtr = SPtr<VectorTransmitter>;

public:
    FineToCoarseVectorConnector(SPtr<Block3D> block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver,
                                int sendDir, InterpolationProcessorPtr iprocessor, CFconnectorType connType);

    bool isLocalConnector() override;
    bool isRemoteConnector() override;
    void init() override;

    void sendTransmitterDataSize() override;
    void receiveTransmitterDataSize() override;

    void prepareForSend() override;
    void sendVectors() override;

    void prepareForReceive() override;
    void receiveVectors() override;

    void fillSendVectors() override;
    void distributeReceiveVectors() override;

    bool isInterpolationConnectorCF() override { return false; }
    bool isInterpolationConnectorFC() override { return true; }

    real getSendRecieveTime();

    void prepareForSendX1() override {}
    void prepareForSendX2() override {}
    void prepareForSendX3() override {}

    void sendVectorsX1() override {}
    void sendVectorsX2() override {}
    void sendVectorsX3() override {}

    void prepareForReceiveX1() override {}
    void prepareForReceiveX2() override {}
    void prepareForReceiveX3() override {}

    void receiveVectorsX1() override {}
    void receiveVectorsX2() override {}
    void receiveVectorsX3() override {}
    void getLocalMinMax(int &minX1, int &minX2, int &minX3, int &maxX1, int &maxX2, int &maxX3);

protected:
    WPtr<Block3D> block; // dieser nvd sendet daten und die empfangenen werden diesem nvd zugeordnet
    // gegenstelle muss "inversen" connector besitzen
    VectorTransmitterPtr sender, receiver;

    InterpolationProcessorPtr iprocessor;

    CFconnectorType connType;

    void writeICellCtoData(vector_type &data, int &index, real *icellC);
    void writeNodeToVector(vector_type &data, int &index, real *inode);
    //void getLocalMinMax(int &minX1, int &minX2, int &minX3, int &maxX1, int &maxX2, int &maxX3);
    void getLocalMinMax(int &minX1, int &minX2, int &minX3, int &maxX1, int &maxX2, int &maxX3,
                        CFconnectorType connType);
    void getLocalMinMaxCF(int gMax, int &lMin, int &lMax);
    void fillSendVector(SPtr<DistributionArray3D> fFrom, const int &lMinX1, const int &lMinX2, const int &lMinX3,
                        const int &lMaxX1, const int &lMaxX2, const int &lMaxX3, vector_type &data, int &index);

    void distributeReceiveVector(SPtr<DistributionArray3D> fTo, const int &lMinX1, const int &lMinX2, const int &lMinX3,
                                 const int &lMaxX1, const int &lMaxX2, const int &lMaxX3, vector_type &data,
                                 int &index);
    void readICellFfromData(vector_type &data, int &index, D3Q27ICell &icellF);
    void readNodeFromVector(vector_type &data, int &index, real *inode);
    void getLocalOffsets(const int &gMax, int &oMin);
    void getLocalMins(int &minX1, int &minX2, int &minX3, const int &oMinX1, const int &oMinX2, const int &oMinX3);

    int bMaxX1, bMaxX2, bMaxX3;
};
////////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
FineToCoarseVectorConnector<VectorTransmitter>::FineToCoarseVectorConnector(SPtr<Block3D> block,
                                                                            VectorTransmitterPtr sender,
                                                                            VectorTransmitterPtr receiver, int sendDir,
                                                                            InterpolationProcessorPtr iprocessor,
                                                                            CFconnectorType connType)
    : Block3DConnector(sendDir), block(block), sender(sender), receiver(receiver), iprocessor(iprocessor),
      connType(connType)
{
    using namespace vf::lbm::dir;

    if (!(sendDir == DIR_P00 || sendDir == DIR_M00 || sendDir == DIR_0P0 ||
          sendDir == DIR_0M0 || sendDir == DIR_00P || sendDir == DIR_00M ||
          sendDir == DIR_PP0 || sendDir == DIR_MM0 || sendDir == DIR_PM0 ||
          sendDir == DIR_MP0 || sendDir == DIR_P0P || sendDir == DIR_M0M ||
          sendDir == DIR_P0M || sendDir == DIR_M0P || sendDir == DIR_0PP ||
          sendDir == DIR_0MM || sendDir == DIR_0PM || sendDir == DIR_0MP

          || sendDir == DIR_PPP || sendDir == DIR_MPP || sendDir == DIR_PMP ||
          sendDir == DIR_MMP || sendDir == DIR_PPM || sendDir == DIR_MPM ||
          sendDir == DIR_PMM || sendDir == DIR_MMM

          )) {
        throw UbException(UB_EXARGS, "invalid constructor for this direction");
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
bool FineToCoarseVectorConnector<VectorTransmitter>::isLocalConnector()
{
    return !this->isRemoteConnector();
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
bool FineToCoarseVectorConnector<VectorTransmitter>::isRemoteConnector()
{
    return sender->isRemoteTransmitter() || receiver->isRemoteTransmitter();
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::sendTransmitterDataSize()
{
    if (sender) {
        UBLOG(logDEBUG5, "FineToCoarseVectorConnector<VectorTransmitter>::sendTransmitterDataSize()"
                             << block.lock()->toString() + "sendDir=" << sendDir);
        sender->sendDataSize();
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()
{
    if (receiver) {
        UBLOG(logDEBUG5, "FineToCoarseVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()"
                             << block.lock()->toString() << "sendDir=" << sendDir);
        receiver->receiveDataSize();
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::prepareForSend()
{
    if (sender)
        sender->prepareForSend();
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::sendVectors()
{
    if (sender)
        sender->sendData();
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::prepareForReceive()
{
    if (receiver)
        receiver->prepareForReceive();
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::receiveVectors()
{
    if (receiver)
        receiver->receiveData();
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::init()
{
    using namespace D3Q27System;
    using namespace vf::lbm::dir;

    bMaxX1 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1();
    bMaxX2 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2();
    bMaxX3 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3();

    int sendSize      = 0;
    real initValue = -999.0;

    int sendDataPerNode = 27 /*f*/;
    int iCellSize       = 1; // size of interpolation cell

    switch (this->sendDir) {
        case DIR_P00:
        case DIR_M00:
            sendSize = (bMaxX2 - 1) / 2 * (bMaxX3 - 1) / 2 * sendDataPerNode * iCellSize;
            break;
        case DIR_0P0:
        case DIR_0M0:
            sendSize = (bMaxX1 - 1) / 2 * (bMaxX3 - 1) / 2 * sendDataPerNode * iCellSize;
            break;
        case DIR_00P:
        case DIR_00M:
            sendSize = (bMaxX1 - 1) / 2 * (bMaxX2 - 1) / 2 * sendDataPerNode * iCellSize;
            break;
        case DIR_PP0:
        case DIR_MM0:
        case DIR_PM0:
        case DIR_MP0:
            sendSize = (3 * bMaxX3 - 3) * sendDataPerNode * iCellSize;
            break; // buffer overhead, should be (3*bMaxX3-6) for even bMax3
        case DIR_P0P:
        case DIR_M0M:
        case DIR_P0M:
        case DIR_M0P:
            sendSize = (3 * bMaxX2 - 3) * sendDataPerNode * iCellSize;
            break;
        case DIR_0PP:
        case DIR_0MM:
        case DIR_0PM:
        case DIR_0MP:
            sendSize = (3 * bMaxX1 - 3) * sendDataPerNode * iCellSize;
            break;
        case DIR_PPP:
        case DIR_MPP:
        case DIR_PMP:
        case DIR_MMP:
        case DIR_PPM:
        case DIR_MPM:
        case DIR_PMM:
        case DIR_MMM:
            sendSize = 3 * (3 * bMaxX1 - 3) * sendDataPerNode * iCellSize;
            break;
        default:
            throw UbException(UB_EXARGS, "direction not allowed in this constructor");
    }
    sender->getData().resize(sendSize, initValue);
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::fillSendVectors()
{
    using namespace D3Q27System;
    using namespace vf::lbm::dir;

    SPtr<DistributionArray3D> fFrom = block.lock()->getKernel()->getDataSet()->getFdistributions();
    int maxX1                       = (int)fFrom->getNX1();
    int maxX2                       = (int)fFrom->getNX2();
    int maxX3                       = (int)fFrom->getNX3();
    int minX1                       = 0;
    int minX2                       = 0;
    int minX3                       = 0;

    int oMinX1, oMinX2, oMinX3;
    getLocalOffsets(maxX1, oMinX1);
    getLocalOffsets(maxX2, oMinX2);
    getLocalOffsets(maxX3, oMinX3);

    int lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3;
    int index         = 0;
    vector_type &data = sender->getData();

    lMinX1 = minX1 + 1;
    lMinX2 = minX2 + 1;
    lMinX3 = minX3 + 1;
    lMaxX1 = maxX1 - 2;
    lMaxX2 = maxX2 - 2;
    lMaxX3 = maxX3 - 2;

    ///////////////////////////////////////
    /// DEBUG
#ifdef _DEBUG
    // if (block.lock()->getGlobalID() == 2558)
    //{
    //   int test = 0;
    //}
#endif
    //////////////

    switch (sendDir) {
        case DIR_P00:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX1 = maxX1 - 7;
            lMaxX1 = lMinX1 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
        case DIR_M00:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX1 = 5;
            lMaxX1 = lMinX1 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
        case DIR_0P0:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX2 = maxX2 - 7;
            lMaxX2 = lMinX2 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
        case DIR_0M0:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX2 = 5;
            lMaxX2 = lMinX2 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
        case DIR_00P:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX3 = maxX3 - 7;
            lMaxX3 = lMinX3 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
        case DIR_00M:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX3 = 5;
            lMaxX3 = lMinX3 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        //	////N-S-E-W
        case DIR_PP0:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX1 = maxX1 - 7;
            lMaxX1 = lMinX1 + 5;
            lMinX2 = maxX2 - 7;
            lMaxX2 = lMinX2 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = maxX1 - 7;
            lMaxX1 = lMinX1 + 1;
            lMinX2 = maxX2 - 7;
            lMaxX2 = lMinX2 + 5;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
        case DIR_MM0:

            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX1 = 1;
            lMaxX1 = lMinX1 + 5;
            lMinX2 = 5;
            lMaxX2 = lMinX2 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = 5;
            lMaxX1 = lMinX1 + 1;
            lMinX2 = 1;
            lMaxX2 = lMinX2 + 5;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_PM0:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX1 = maxX1 - 7;
            lMaxX1 = lMinX1 + 5;
            lMinX2 = 5;
            lMaxX2 = lMinX2 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = maxX1 - 7;
            lMaxX1 = lMinX1 + 1;
            lMinX2 = 1;
            lMaxX2 = lMinX2 + 5;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            break;

        case DIR_MP0:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX1 = 1;
            lMaxX1 = lMinX1 + 5;
            lMinX2 = maxX2 - 7;
            lMaxX2 = lMinX2 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = 5;
            lMaxX1 = lMinX1 + 1;
            lMinX2 = maxX2 - 7;
            lMaxX2 = lMinX2 + 5;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
            //////T-B-E-W
        case DIR_P0P:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX1 = maxX1 - 7;
            lMaxX1 = lMinX1 + 5;
            lMinX3 = maxX3 - 7;
            lMaxX3 = lMinX3 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = maxX1 - 7;
            lMaxX1 = lMinX1 + 1;
            lMinX3 = maxX3 - 7;
            lMaxX3 = lMinX3 + 5;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_M0M:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX1 = 1;
            lMaxX1 = lMinX1 + 5;
            lMinX3 = 5;
            lMaxX3 = lMinX3 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = 5;
            lMaxX1 = lMinX1 + 1;
            lMinX3 = 1;
            lMaxX3 = lMinX3 + 5;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_P0M:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX1 = maxX1 - 7;
            lMaxX1 = lMinX1 + 5;
            lMinX3 = 5;
            lMaxX3 = lMinX3 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = maxX1 - 7;
            lMaxX1 = lMinX1 + 1;
            lMinX3 = 1;
            lMaxX3 = lMinX3 + 5;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_M0P:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX1 = 1;
            lMaxX1 = lMinX1 + 5;
            lMinX3 = maxX3 - 7;
            lMaxX3 = lMinX3 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = 5;
            lMaxX1 = lMinX1 + 1;
            lMinX3 = maxX3 - 7;
            lMaxX3 = lMinX3 + 5;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
            ///////////////T-B-N-S
            //
        case DIR_0PP:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX2 = maxX2 - 7;
            lMaxX2 = lMinX2 + 5;
            lMinX3 = maxX3 - 7;
            lMaxX3 = lMinX3 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX2 = maxX2 - 7;
            lMaxX2 = lMinX2 + 1;
            lMinX3 = maxX3 - 7;
            lMaxX3 = lMinX3 + 5;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_0MM:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX2 = 1;
            lMaxX2 = lMinX2 + 5;
            lMinX3 = 5;
            lMaxX3 = lMinX3 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX2 = 5;
            lMaxX2 = lMinX2 + 1;
            lMinX3 = 1;
            lMaxX3 = lMinX3 + 5;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_0PM:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX2 = maxX2 - 7;
            lMaxX2 = lMinX2 + 5;
            lMinX3 = 5;
            lMaxX3 = lMinX3 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX2 = maxX2 - 7;
            lMaxX2 = lMinX2 + 1;
            lMinX3 = 1;
            lMaxX3 = lMinX3 + 5;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_0MP:
            getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
            getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
            lMinX2 = 1;
            lMaxX2 = lMinX2 + 5;
            lMinX3 = maxX3 - 7;
            lMaxX3 = lMinX3 + 1;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX2 = 5;
            lMaxX2 = lMinX2 + 1;
            lMinX3 = maxX3 - 7;
            lMaxX3 = lMinX3 + 5;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        // TNE
        case DIR_PPP:
            lMinX1 = maxX1 - 7;
            lMaxX1 = maxX1 - 6;
            lMinX2 = maxX2 - 7;
            lMaxX2 = maxX2 - 2;
            lMinX3 = maxX3 - 7;
            lMaxX3 = maxX3 - 2;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = maxX1 - 7;
            lMaxX1 = maxX1 - 2;
            lMinX2 = maxX2 - 7;
            lMaxX2 = maxX2 - 6;
            lMinX3 = maxX3 - 7;
            lMaxX3 = maxX3 - 2;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = maxX1 - 7;
            lMaxX1 = maxX1 - 2;
            lMinX2 = maxX2 - 7;
            lMaxX2 = maxX2 - 2;
            lMinX3 = maxX3 - 7;
            lMaxX3 = maxX3 - 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        // TNW
        case DIR_MPP:
            lMinX1 = 5;
            lMaxX1 = 6;
            lMinX2 = maxX2 - 7;
            lMaxX2 = maxX2 - 2;
            lMinX3 = maxX3 - 7;
            lMaxX3 = maxX3 - 2;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = 1;
            lMaxX1 = 6;
            lMinX2 = maxX2 - 7;
            lMaxX2 = maxX2 - 6;
            lMinX3 = maxX3 - 7;
            lMaxX3 = maxX3 - 2;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = 1;
            lMaxX1 = 6;
            lMinX2 = maxX2 - 7;
            lMaxX2 = maxX2 - 2;
            lMinX3 = maxX3 - 7;
            lMaxX3 = maxX3 - 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

            break;

        //      TSE
        case DIR_PMP:
            lMinX1 = maxX1 - 7;
            lMaxX1 = maxX1 - 6;
            lMinX2 = 1;
            lMaxX2 = 6;
            lMinX3 = maxX3 - 7;
            lMaxX3 = maxX3 - 2;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = maxX1 - 7;
            lMaxX1 = maxX1 - 2;
            lMinX2 = 5;
            lMaxX2 = 6;
            lMinX3 = maxX3 - 7;
            lMaxX3 = maxX3 - 2;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = maxX1 - 7;
            lMaxX1 = maxX1 - 2;
            lMinX2 = 1;
            lMaxX2 = 6;
            lMinX3 = maxX3 - 7;
            lMaxX3 = maxX3 - 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            break;
        //      TSW
        case DIR_MMP:
            lMinX1 = 5;
            lMaxX1 = 6;
            lMinX2 = 1;
            lMaxX2 = 6;
            lMinX3 = maxX3 - 7;
            lMaxX3 = maxX3 - 2;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = 1;
            lMaxX1 = 6;
            lMinX2 = 5;
            lMaxX2 = 6;
            lMinX3 = maxX3 - 7;
            lMaxX3 = maxX3 - 2;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = 1;
            lMaxX1 = 6;
            lMinX2 = 1;
            lMaxX2 = 6;
            lMinX3 = maxX3 - 7;
            lMaxX3 = maxX3 - 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            break;
        //      BNE
        case DIR_PPM:
            lMinX1 = maxX1 - 7;
            lMaxX1 = maxX1 - 6;
            lMinX2 = maxX2 - 7;
            lMaxX2 = maxX2 - 2;
            lMinX3 = 1;
            lMaxX3 = 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = maxX1 - 7;
            lMaxX1 = maxX1 - 2;
            lMinX2 = maxX2 - 7;
            lMaxX2 = maxX2 - 6;
            lMinX3 = 1;
            lMaxX3 = 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = maxX1 - 7;
            lMaxX1 = maxX1 - 2;
            lMinX2 = maxX2 - 7;
            lMaxX2 = maxX2 - 2;
            lMinX3 = 5;
            lMaxX3 = 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            break;
        //      BNW
        case DIR_MPM:
            lMinX1 = 5;
            lMaxX1 = 6;
            lMinX2 = maxX2 - 7;
            lMaxX2 = maxX2 - 2;
            lMinX3 = 1;
            lMaxX3 = 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = 1;
            lMaxX1 = 6;
            lMinX2 = maxX2 - 7;
            lMaxX2 = maxX2 - 6;
            lMinX3 = 1;
            lMaxX3 = 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = 1;
            lMaxX1 = 6;
            lMinX2 = maxX2 - 7;
            lMaxX2 = maxX2 - 2;
            lMinX3 = 5;
            lMaxX3 = 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            break;

        //      BSE
        case DIR_PMM:
            lMinX1 = maxX1 - 7;
            lMaxX1 = maxX1 - 6;
            lMinX2 = 1;
            lMaxX2 = 6;
            lMinX3 = 1;
            lMaxX3 = 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = maxX1 - 7;
            lMaxX1 = maxX1 - 2;
            lMinX2 = 5;
            lMaxX2 = 6;
            lMinX3 = 1;
            lMaxX3 = 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = maxX1 - 7;
            lMaxX1 = maxX1 - 2;
            lMinX2 = 1;
            lMaxX2 = 6;
            lMinX3 = 5;
            lMaxX3 = 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            break;

        // BSW
        case DIR_MMM:
            lMinX1 = 5;
            lMaxX1 = 6;
            lMinX2 = 1;
            lMaxX2 = 6;
            lMinX3 = 1;
            lMaxX3 = 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = 1;
            lMaxX1 = 6;
            lMinX2 = 5;
            lMaxX2 = 6;
            lMinX3 = 1;
            lMaxX3 = 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            lMinX1 = 1;
            lMaxX1 = 6;
            lMinX2 = 1;
            lMaxX2 = 6;
            lMinX3 = 5;
            lMaxX3 = 6;
            fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

            break;
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::fillSendVector(SPtr<DistributionArray3D> fFrom, const int &lMinX1,
                                                                    const int &lMinX2, const int &lMinX3,
                                                                    const int &lMaxX1, const int &lMaxX2,
                                                                    const int &lMaxX3, vector_type &data, int &index)
{
    using namespace vf::basics::constant;

    int ix1, ix2, ix3;
    real xoff, yoff, zoff;
    SPtr<BCArray3D> bcArray = block.lock()->getKernel()->getBCSet()->getBCArray();

    for (ix3 = lMinX3; ix3 < lMaxX3; ix3 += 2) {
        for (ix2 = lMinX2; ix2 < lMaxX2; ix2 += 2) {
            for (ix1 = lMinX1; ix1 < lMaxX1; ix1 += 2) {
                real icellC[27];
                D3Q27ICell icellF;

                int howManySolids = iprocessor->iCellHowManySolids(bcArray, ix1, ix2, ix3);

                if (howManySolids == 0 || howManySolids == 8) {
                    iprocessor->readICell(fFrom, icellF, ix1, ix2, ix3);
                    xoff = c0o1;
                    yoff = c0o1;
                    zoff = c0o1;
                } else {
                    if (!iprocessor->findNeighborICell(bcArray, fFrom, icellF, bMaxX1, bMaxX2, bMaxX3, ix1, ix2, ix3,
                                                       xoff, yoff, zoff)) {
                        std::string err = "For " + block.lock()->toString() + " x1=" + UbSystem::toString(ix1) +
                                          ", x2=" + UbSystem::toString(ix2) + ", x3=" + UbSystem::toString(ix3) +
                                          " interpolation is not implemented for other direction" +
                                          " by using in: " + (std::string) typeid(*this).name() +
                                          " or maybe you have a solid on the block boundary";
                        // UBLOG(logINFO, err);
                        UB_THROW(UbException(UB_EXARGS, err));
                    }
                }

                iprocessor->interpolateFineToCoarse(icellF, icellC, xoff, yoff, zoff);
                this->writeICellCtoData(data, index, icellC);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::writeICellCtoData(vector_type &data, int &index, real *icellC)
{
    for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF + 1; i++) {
        data[index++] = icellC[i];
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::getLocalMinMaxCF(int gMax, int &lMin, int &lMax)
{
    if (Utilities::isOdd(gMax)) {
        if (connType == OddEvenSE || connType == OddOddNE) {
            lMin = 1;
            lMax = gMax;
        }
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::distributeReceiveVectors()
{
    using namespace D3Q27System;
    using namespace vf::lbm::dir;

    SPtr<DistributionArray3D> fTo = block.lock()->getKernel()->getDataSet()->getFdistributions();
    int maxX1                     = (int)fTo->getNX1();
    int maxX2                     = (int)fTo->getNX2();
    int maxX3                     = (int)fTo->getNX3();
    int minX1                     = 0;
    int minX2                     = 0;
    int minX3                     = 0;

    int lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3;
    int index         = 0;
    vector_type &data = receiver->getData();

    lMinX1 = minX1;
    lMinX2 = minX2;
    lMinX3 = minX3;
    lMaxX1 = maxX1 - 1;
    lMaxX2 = maxX2 - 1;
    lMaxX3 = maxX3 - 1;

    switch (sendDir) {
        case DIR_P00:
            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 1;
            getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
            getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
        case DIR_M00:
            lMinX1 = 2;
            lMaxX1 = lMinX1 + 1;
            getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
            getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
        case DIR_0P0:
            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 1;
            getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
            getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
        case DIR_0M0:
            lMinX2 = 2;
            lMaxX2 = lMinX2 + 1;
            getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
            getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
        case DIR_00P:
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 1;
            getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
            getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
        case DIR_00M:
            lMinX3 = 2;
            lMaxX3 = lMinX3 + 1;
            getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
            getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

            /////E-W-N-S
        case DIR_PP0:
            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 3;
            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 3;
            getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_MM0:
            lMinX1 = 0;
            lMaxX1 = lMinX1 + 3;
            lMinX2 = 0;
            lMaxX2 = lMinX2 + 3;
            getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_PM0:
            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 3;
            lMinX2 = 0;
            lMaxX2 = lMinX2 + 3;
            getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_MP0:
            lMinX1 = 0;
            lMaxX1 = lMinX1 + 3;
            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 3;
            getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
        //
        //	/////T-B-E-W
        case DIR_P0P:
            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 3;
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_M0M:
            lMinX1 = 0;
            lMaxX1 = lMinX1 + 3;
            lMinX3 = 0;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_P0M:
            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 3;
            lMinX3 = 0;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_M0P:
            lMinX1 = 0;
            lMaxX1 = lMinX1 + 3;
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        //	////////////////T-B-N-S
        //
        case DIR_0PP:
            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 3;
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_0MM:
            lMinX2 = 0;
            lMaxX2 = lMinX2 + 3;
            lMinX3 = 0;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_0PM:
            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 3;
            lMinX3 = 0;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        case DIR_0MP:
            lMinX2 = 0;
            lMaxX2 = lMinX2 + 3;
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;

        //   //TNE
        case DIR_PPP:
            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 1;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 1;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 1;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
            //   TNW
        case DIR_MPP:
            lMinX1 = 0;
            lMaxX1 = 3;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 1;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 1;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
            //   TSE
        case DIR_PMP:
            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 1;
            lMinX2 = 0;
            lMaxX2 = 3;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 1;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
            //   TSW
        case DIR_MMP:
            lMinX1 = 0;
            lMaxX1 = 3;
            lMinX2 = 0;
            lMaxX2 = 3;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 1;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
            //   BNE
        case DIR_PPM:
            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 1;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 1;
            lMinX3 = 0;
            lMaxX3 = 3;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
            //   BNW
        case DIR_MPM:
            lMinX1 = 0;
            lMaxX1 = 3;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 1;
            lMinX3 = 0;
            lMaxX3 = 3;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
            //   BSE
        case DIR_PMM:
            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 1;
            lMinX2 = 0;
            lMaxX2 = 3;
            lMinX3 = 0;
            lMaxX3 = 3;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
            // BSW
        case DIR_MMM:
            lMinX1 = 0;
            lMaxX1 = 3;
            lMinX2 = 0;
            lMaxX2 = 3;
            lMinX3 = 0;
            lMaxX3 = 3;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
            break;
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::distributeReceiveVector(SPtr<DistributionArray3D> fTo,
                                                                             const int &lMinX1, const int &lMinX2,
                                                                             const int &lMinX3, const int &lMaxX1,
                                                                             const int &lMaxX2, const int &lMaxX3,
                                                                             vector_type &data, int &index)
{
    if (data.size() == 0)
        return;

    int ix1, ix2, ix3;
    for (ix3 = lMinX3; ix3 < lMaxX3; ix3 += 2) {
        for (ix2 = lMinX2; ix2 < lMaxX2; ix2 += 2) {
            for (ix1 = lMinX1; ix1 < lMaxX1; ix1 += 2) {
                D3Q27ICell icellF;
                this->readICellFfromData(data, index, icellF);
                iprocessor->writeICellInv(fTo, icellF, ix1, ix2, ix3);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::readICellFfromData(vector_type &data, int &index,
                                                                        D3Q27ICell &icellF)
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
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::readNodeFromVector(vector_type &data, int &index, real *inode)
{
    for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF + 1; i++) {
        inode[i] = data[index++];
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::getLocalMinMax(int &minX1, int &minX2, int &minX3, int &maxX1,
                                                                    int &maxX2, int &maxX3)
{
    using namespace D3Q27System;
    using namespace vf::lbm::dir;
    
    int TminX1 = minX1;
    int TminX2 = minX2;
    int TminX3 = minX3;
    int TmaxX1 = maxX1;
    int TmaxX2 = maxX2;
    int TmaxX3 = maxX3;

    if (block.lock()->hasInterpolationFlagFC(DIR_P00)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 3;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_M00)) {
        if (minX1 == TminX1)
            minX1 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_0P0)) {
        if (maxX2 == TmaxX2)
            maxX2 -= 3;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_0M0)) {
        if (minX2 == TminX2)
            minX2 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_00P)) {
        if (maxX3 == TmaxX3)
            maxX3 -= 3;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_00M)) {
        if (minX3 == TminX3)
            minX3 += 4;
    }

    ////////////
    /////E-W-N-S
    if (block.lock()->hasInterpolationFlagFC(DIR_PP0) && !block.lock()->hasInterpolationFlagFC(DIR_0P0) &&
        !block.lock()->hasInterpolationFlagFC(DIR_P00)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 3;
        if (maxX2 == TmaxX2)
            maxX2 -= 3;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_MM0) && !block.lock()->hasInterpolationFlagFC(DIR_M00) &&
        !block.lock()->hasInterpolationFlagFC(DIR_0M0)) {
        if (minX1 == TminX1)
            minX1 += 4;
        if (minX2 == TminX2)
            minX2 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_PM0) && !block.lock()->hasInterpolationFlagFC(DIR_P00) &&
        !block.lock()->hasInterpolationFlagFC(DIR_0M0)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 3;
        if (minX2 == TminX2)
            minX2 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_MP0) && !block.lock()->hasInterpolationFlagFC(DIR_0P0) &&
        !block.lock()->hasInterpolationFlagFC(DIR_M00)) {
        if (minX1 == TminX1)
            minX1 += 4;
        if (maxX2 == TmaxX2)
            maxX2 -= 3;
    }

    //////T-B-E-W
    if (block.lock()->hasInterpolationFlagFC(DIR_P0P) && !block.lock()->hasInterpolationFlagFC(DIR_P00) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00P)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 3;
        if (maxX3 == TmaxX3)
            maxX3 -= 3;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_M0M) && !block.lock()->hasInterpolationFlagFC(DIR_M00) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00M)) {
        if (minX1 == TminX1)
            minX1 += 4;
        if (minX3 == TminX3)
            minX3 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_P0M) && !block.lock()->hasInterpolationFlagFC(DIR_P00) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00M)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 3;
        if (minX3 == TminX3)
            minX3 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_M0P) && !block.lock()->hasInterpolationFlagFC(DIR_M00) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00P)) {
        if (minX1 == TminX1)
            minX1 += 4;
        if (maxX3 == TmaxX3)
            maxX3 -= 3;
    }

    ////T-B-N-S
    if (block.lock()->hasInterpolationFlagFC(DIR_0PP) && !block.lock()->hasInterpolationFlagFC(DIR_0P0) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00P)) {
        if (maxX2 == TmaxX2)
            maxX2 -= 3;
        if (maxX3 == TmaxX3)
            maxX3 -= 3;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_0MM) && !block.lock()->hasInterpolationFlagFC(DIR_0M0) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00M)) {
        if (minX2 == TminX2)
            minX2 += 4;
        if (minX3 == TminX3)
            minX3 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_0PM) && !block.lock()->hasInterpolationFlagFC(DIR_0P0) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00M)) {
        if (maxX2 == TmaxX2)
            maxX2 -= 3;
        if (minX3 == TminX3)
            minX3 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_0MP) && !block.lock()->hasInterpolationFlagFC(DIR_0M0) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00P)) {
        if (minX2 == TminX2)
            minX2 += 4;
        if (maxX3 == TmaxX3)
            maxX3 -= 3;
    }

    // if
    // (block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_PPP)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_P0P)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_0PP)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_PP0)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_00P)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_0P0)
    // && !block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_P00)) if
    // (!block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_P0P)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_00P) &&
    // !block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_P00))
    //{
    //   if (maxX1==TmaxX1) maxX1 -= 3;
    //   if (maxX2==TmaxX2) maxX2 -= 3;
    //   if (maxX3==TmaxX3) maxX3 -= 3;
    //}
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::getLocalMinMax(int &minX1, int &minX2, int &minX3, int &maxX1,
                                                                    int &maxX2, int &maxX3,
                                                                    CFconnectorType /*connType*/)
{
    using namespace D3Q27System;
    using namespace vf::lbm::dir;

    int TminX1 = minX1;
    int TminX2 = minX2;
    int TminX3 = minX3;
    int TmaxX1 = maxX1;
    int TmaxX2 = maxX2;
    int TmaxX3 = maxX3;

    if (block.lock()->hasInterpolationFlagFC(DIR_P00)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 3;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_M00)) {
        if (minX1 == TminX1)
            minX1 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_0P0)) {
        if (maxX2 == TmaxX2)
            maxX2 -= 3;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_0M0)) {
        if (minX2 == TminX2)
            minX2 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_00P)) {
        if (maxX3 == TmaxX3)
            maxX3 -= 3;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_00M)) {
        if (minX3 == TminX3)
            minX3 += 4;
    }

    ////////////
    /////E-W-N-S
    if (block.lock()->hasInterpolationFlagFC(DIR_PP0) && !block.lock()->hasInterpolationFlagFC(DIR_0P0) &&
        !block.lock()->hasInterpolationFlagFC(DIR_P00)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 3;
        if (maxX2 == TmaxX2)
            maxX2 -= 3;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_MM0) && !block.lock()->hasInterpolationFlagFC(DIR_M00) &&
        !block.lock()->hasInterpolationFlagFC(DIR_0M0)) {
        if (minX1 == TminX1)
            minX1 += 4;
        if (minX2 == TminX2)
            minX2 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_PM0) && !block.lock()->hasInterpolationFlagFC(DIR_P00) &&
        !block.lock()->hasInterpolationFlagFC(DIR_0M0)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 3;
        if (minX2 == TminX2)
            minX2 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_MP0) && !block.lock()->hasInterpolationFlagFC(DIR_0P0) &&
        !block.lock()->hasInterpolationFlagFC(DIR_M00)) {
        if (minX1 == TminX1)
            minX1 += 4;
        if (maxX2 == TmaxX2)
            maxX2 -= 3;
    }

    //////T-B-E-W
    if (block.lock()->hasInterpolationFlagFC(DIR_P0P) && !block.lock()->hasInterpolationFlagFC(DIR_P00) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00P)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 3;
        if (maxX3 == TmaxX3)
            maxX3 -= 3;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_M0M) && !block.lock()->hasInterpolationFlagFC(DIR_M00) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00M)) {
        if (minX1 == TminX1)
            minX1 += 4;
        if (minX3 == TminX3)
            minX3 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_P0M) && !block.lock()->hasInterpolationFlagFC(DIR_P00) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00M)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 3;
        if (minX3 == TminX3)
            minX3 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_M0P) && !block.lock()->hasInterpolationFlagFC(DIR_M00) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00P)) {
        if (minX1 == TminX1)
            minX1 += 4;
        if (maxX3 == TmaxX3)
            maxX3 -= 3;
    }

    ////T-B-N-S
    if (block.lock()->hasInterpolationFlagFC(DIR_0PP) && !block.lock()->hasInterpolationFlagFC(DIR_0P0) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00P)) {
        if (maxX2 == TmaxX2)
            maxX2 -= 3;
        if (maxX3 == TmaxX3)
            maxX3 -= 3;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_0MM) && !block.lock()->hasInterpolationFlagFC(DIR_0M0) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00M)) {
        if (minX2 == TminX2)
            minX2 += 4;
        if (minX3 == TminX3)
            minX3 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_0PM) && !block.lock()->hasInterpolationFlagFC(DIR_0P0) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00M)) {
        if (maxX2 == TmaxX2)
            maxX2 -= 3;
        if (minX3 == TminX3)
            minX3 += 4;
    }
    if (block.lock()->hasInterpolationFlagFC(DIR_0MP) && !block.lock()->hasInterpolationFlagFC(DIR_0M0) &&
        !block.lock()->hasInterpolationFlagFC(DIR_00P)) {
        if (minX2 == TminX2)
            minX2 += 4;
        if (maxX3 == TmaxX3)
            maxX3 -= 3;
    }

    // if
    // (block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_PPP)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_P0P)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_0PP)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_PP0)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_00P)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_0P0)
    // && !block.lock()->hasInterpolationFlagFC(D3Q27System::DIR_P00))
    //{
    //   if (maxX1==TmaxX1) maxX1 -= 3;
    //   if (maxX2==TmaxX2) maxX2 -= 3;
    //   if (maxX3==TmaxX3) maxX3 -= 3;
    //}
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::getLocalOffsets(const int &gMax, int &oMin)
{
    if (Utilities::isEven(gMax)) {
        oMin = 0;
    }
    if (Utilities::isOdd(gMax)) {
        oMin = -1;
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void FineToCoarseVectorConnector<VectorTransmitter>::getLocalMins(int &minX1, int &minX2, int &minX3, const int &oMinX1,
                                                                  const int &oMinX2, const int &oMinX3)
{
    using namespace D3Q27System;
    using namespace vf::lbm::dir;

    switch (sendDir) {
        case DIR_P00:
        case DIR_M00:
            if (connType == OddEvenSE)
                minX2 += oMinX2;
            if (connType == OddOddNE) {
                minX2 += oMinX2;
                minX3 += oMinX3;
            }
            if (connType == EvenOddNW)
                minX3 += oMinX3;
            break;
        case DIR_0P0:
        case DIR_0M0:
            if (connType == OddEvenSE)
                minX1 += oMinX1;
            if (connType == OddOddNE) {
                minX1 += oMinX1;
                minX3 += oMinX3;
            }
            if (connType == EvenOddNW)
                minX3 += oMinX3;
            break;
        case DIR_00P:
        case DIR_00M:
            if (connType == OddEvenSE)
                minX1 += oMinX1;
            if (connType == OddOddNE) {
                minX1 += oMinX1;
                minX2 += oMinX2;
            }
            if (connType == EvenOddNW)
                minX2 += oMinX2;
            break;

            /////
        case DIR_PP0:
        case DIR_MM0:
        case DIR_PM0:
        case DIR_MP0:
            // case SW:
            if (connType == OddEvenSE)
                // minX2 += oMinX2;
                if (connType == OddOddNE) {
                    // minX2 += oMinX2;
                    minX3 += oMinX3;
                }
            if (connType == EvenOddNW)
                minX3 += oMinX3;
            break;

            //////
        case DIR_P0P:
        case DIR_M0M:
        case DIR_P0M:
        case DIR_M0P:
            if (connType == OddEvenSE)
                //		minX1 += oMinX1;
                if (connType == OddOddNE) {
                    //		minX1 += oMinX1;
                    minX2 += oMinX2;
                }
            if (connType == EvenOddNW)
                minX2 += oMinX2;
            break;

        //	//////
        case DIR_0PP:
        case DIR_0MM:
        case DIR_0PM:
        case DIR_0MP:
            if (connType == OddEvenSE)
                minX1 += oMinX1;
            if (connType == OddOddNE) {
                minX1 += oMinX1;
                // minX3 += oMinX3;
            }
            if (connType == EvenOddNW)
                // minX3 += oMinX3;
                break;

            //	/////
            //	case TNE: case TNW: case TSE: case TSW: case BNE: case BNW: case BSE: case BSW:
            //	if(connType == OddEvenSE)
            //	//	minX1 += oMinX1;
            //	if(connType == OddOddNE)
            //	{
            //		//minX1 += oMinX1;
            //		//minX3 += oMinX3;
            //	}
            //	if(connType == EvenOddNW)
            //		//minX3 += oMinX3;
            //	break;
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
real FineToCoarseVectorConnector<VectorTransmitter>::getSendRecieveTime()
{
    return 0;
}

#endif
