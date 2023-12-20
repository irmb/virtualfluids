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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_Connectors Connectors
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef CoarseToFineVectorConnector_H
#define CoarseToFineVectorConnector_H

#include <vector>

#include "Block3D.h"
#include "Block3DConnector.h"
#include "D3Q27System.h"
#include "Grid3D.h"
#include "Interpolator.h"
#include "LBMKernel.h"
#include "MathUtil.hpp"
#include "basics/container/CbVector.h"
#include "parallel/transmitter/TbTransmitter.h"
#include "parallel/transmitter/TbTransmitterLocal.h"
#include <PointerDefinitions.h>
#include "basics/constants/NumericConstants.h"

#include "BCSet.h"
#include "FineToCoarseVectorConnector.h"

class Block3D;

//! \brief Interpolation from coarse level to fine.
//! \details Data are copied into a vector (this is located in the transmitter) the vector is transmitted via transmitter.
//!  Transmitter can be a local, MPI, RCG, CTL or whatever be a transmitter that is derived from the transmitter ;-)
//! send direction:    E<->W     N<->S    T<->B
//!  ---------       x3       x3        x2
//! | NW | NE |      ^        ^         ^
//! |----+----|      +-> x2   +->x1     +->x1
//! | SW | SE |
//!  ---------
//! NW==even-odd, SW==even-even, SE==odd-even, NE==odd-odd
template <typename VectorTransmitter>
class CoarseToFineVectorConnector : public Block3DConnector
{
public:
    using vector_type          = typename VectorTransmitter::value_type;
    using VectorTransmitterPtr = SPtr<VectorTransmitter>;

public:
    CoarseToFineVectorConnector(SPtr<Block3D> block, VectorTransmitterPtr senderEvenEvenSW,
                                VectorTransmitterPtr receiverEvenEvenSW, VectorTransmitterPtr senderEvenOddNW,
                                VectorTransmitterPtr receiverEvenOddNW, VectorTransmitterPtr senderOddEvenSE,
                                VectorTransmitterPtr receiverOddEvenSE, VectorTransmitterPtr senderOddOddNE,
                                VectorTransmitterPtr receiverOddOddNE, int sendDir,
                                InterpolationProcessorPtr iprocessor);

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

    bool isInterpolationConnectorCF() override { return true; }
    bool isInterpolationConnectorFC() override { return false; }

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

protected:
    WPtr<Block3D> block; // dieser nvd sendet daten und die empfangenen werden diesem nvd zugeordnet
    VectorTransmitterPtr senderEvenEvenSW, receiverEvenEvenSW, senderEvenOddNW, receiverEvenOddNW, senderOddEvenSE,
        receiverOddEvenSE, senderOddOddNE, receiverOddOddNE;

    InterpolationProcessorPtr iprocessor;

    void writeICellFtoData(vector_type &data, int &index, D3Q27ICell &icellF);
    void writeNodeToVector(vector_type &data, int &index, real *inode);
    void getLocalMinMax(const int &gMin, const int &gMax, const bool &even, int &lMin, int &lMax,
                        const bool &dataDistribution);
    void getLocalMinMax(int &minX1, int &minX2, int &minX3, int &maxX1, int &maxX2, int &maxX3);
    void getLocalMinMax(int &minX1, int &minX2, int &minX3, int &maxX1, int &maxX2, int &maxX3,
                        CFconnectorType connType);
    void fillSendVectorExt(SPtr<DistributionArray3D> fFrom, const int &lMinX1, const int &lMinX2, const int &lMinX3,
                           const int &lMaxX1, const int &lMaxX2, const int &lMaxX3, vector_type &data, int &index);

    void distributeReceiveVector(SPtr<DistributionArray3D> fTo, const int &lMinX1, const int &lMinX2, const int &lMinX3,
                                 const int &lMaxX1, const int &lMaxX2, const int &lMaxX3, vector_type &data,
                                 int &index);
    void readICellCfromData(vector_type &data, int &index, real *icellC);

    void findCFnodes();
    void findCFnodes(SPtr<DistributionArray3D> fFrom, const int &lMinX1, const int &lMinX2, const int &lMinX3,
                     const int &lMaxX1, const int &lMaxX2, const int &lMaxX3, vector_type &data, int &index);

    int bMaxX1, bMaxX2, bMaxX3;
};

////////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
CoarseToFineVectorConnector<VectorTransmitter>::CoarseToFineVectorConnector(
    SPtr<Block3D> block, VectorTransmitterPtr senderEvenEvenSW, VectorTransmitterPtr receiverEvenEvenSW,
    VectorTransmitterPtr senderEvenOddNW, VectorTransmitterPtr receiverEvenOddNW, VectorTransmitterPtr senderOddEvenSE,
    VectorTransmitterPtr receiverOddEvenSE, VectorTransmitterPtr senderOddOddNE, VectorTransmitterPtr receiverOddOddNE,
    int sendDir, InterpolationProcessorPtr iprocessor)
    : Block3DConnector(sendDir), block(block), senderEvenEvenSW(senderEvenEvenSW), senderEvenOddNW(senderEvenOddNW),
      senderOddEvenSE(senderOddEvenSE), senderOddOddNE(senderOddOddNE), receiverEvenEvenSW(receiverEvenEvenSW),
      receiverEvenOddNW(receiverEvenOddNW), receiverOddEvenSE(receiverOddEvenSE), receiverOddOddNE(receiverOddOddNE),
      iprocessor(iprocessor)
{
    using namespace vf::lbm::dir;

    if (!(sendDir == dP00 || sendDir == dM00 || sendDir == d0P0 ||
          sendDir == d0M0 || sendDir == d00P || sendDir == d00M ||
          sendDir == dPP0 || sendDir == dMM0 || sendDir == dPM0 ||
          sendDir == dMP0 || sendDir == dP0P || sendDir == dM0M ||
          sendDir == dP0M || sendDir == dM0P || sendDir == d0PP ||
          sendDir == d0MM || sendDir == d0PM || sendDir == d0MP ||
          sendDir == dPPP || sendDir == dMPP || sendDir == dPMP ||
          sendDir == dMMP || sendDir == dPPM || sendDir == dMPM ||
          sendDir == dPMM || sendDir == dMMM)) {
        throw UbException(UB_EXARGS, "invalid constructor for this direction");
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
bool CoarseToFineVectorConnector<VectorTransmitter>::isLocalConnector()
{
    return !this->isRemoteConnector();
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
bool CoarseToFineVectorConnector<VectorTransmitter>::isRemoteConnector()
{
    return ((senderOddOddNE && senderOddOddNE->isRemoteTransmitter()) ||
            (receiverOddOddNE && receiverOddOddNE->isRemoteTransmitter()) ||
            (senderEvenEvenSW && senderEvenEvenSW->isRemoteTransmitter()) ||
            (receiverEvenEvenSW && receiverEvenEvenSW->isRemoteTransmitter()) ||
            (senderEvenOddNW && senderEvenOddNW->isRemoteTransmitter()) ||
            (receiverEvenOddNW && receiverEvenOddNW->isRemoteTransmitter()) ||
            (senderOddEvenSE && senderOddEvenSE->isRemoteTransmitter()) ||
            (receiverOddEvenSE && receiverOddEvenSE->isRemoteTransmitter()));
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::sendTransmitterDataSize()
{
    if (senderEvenEvenSW) {
        UBLOG(logDEBUG5, "CoarseToFineVectorConnector<VectorTransmitter>::sendTransmitterDataSize()-senderEvenEvenSW "
                             << block.lock()->toString() << " sendDir=" << sendDir);
        senderEvenEvenSW->sendDataSize();
    }
    if (senderEvenOddNW) {
        UBLOG(logDEBUG5, "CoarseToFineVectorConnector<VectorTransmitter>::sendTransmitterDataSize()-senderEvenOddNW "
                             << block.lock()->toString() << "sendDir=" << sendDir);
        senderEvenOddNW->sendDataSize();
    }
    if (senderOddEvenSE) {
        UBLOG(logDEBUG5, "CoarseToFineVectorConnector<VectorTransmitter>::sendTransmitterDataSize()-senderOddEvenSE "
                             << block.lock()->toString() + "sendDir=" << sendDir);
        senderOddEvenSE->sendDataSize();
    }
    if (senderOddOddNE) {
        UBLOG(logDEBUG5, "CoarseToFineVectorConnector<VectorTransmitter>::sendTransmitterDataSize()-senderOddOddNE "
                             << block.lock()->toString() << "sendDir=" << sendDir);
        senderOddOddNE->sendDataSize();
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()
{
    if (receiverEvenEvenSW) {
        UBLOG(logDEBUG5,
              "CoarseToFineVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()-receiverEvenEvenSW "
                  << block.lock()->toString() << "sendDir=" << sendDir);
        receiverEvenEvenSW->receiveDataSize();
    }
    if (receiverEvenOddNW) {
        UBLOG(logDEBUG5,
              "CoarseToFineVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()-receiverEvenOddNW "
                  << block.lock()->toString() << "sendDir=" << sendDir);
        receiverEvenOddNW->receiveDataSize();
    }
    if (receiverOddEvenSE) {
        UBLOG(logDEBUG5,
              "CoarseToFineVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()-receiverOddEvenSE "
                  << block.lock()->toString() << "sendDir=" << sendDir);
        receiverOddEvenSE->receiveDataSize();
    }
    if (receiverOddOddNE) {
        UBLOG(logDEBUG5,
              "CoarseToFineVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()-receiverOddOddNE "
                  << block.lock()->toString() << "sendDir=" << sendDir);
        receiverOddOddNE->receiveDataSize();
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::prepareForSend()
{
    if (senderEvenEvenSW)
        senderEvenEvenSW->prepareForSend();
    if (senderEvenOddNW)
        senderEvenOddNW->prepareForSend();
    if (senderOddEvenSE)
        senderOddEvenSE->prepareForSend();
    if (senderOddOddNE)
        senderOddOddNE->prepareForSend();
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::sendVectors()
{
    if (senderEvenEvenSW)
        senderEvenEvenSW->sendData();
    if (senderEvenOddNW)
        senderEvenOddNW->sendData();
    if (senderOddEvenSE)
        senderOddEvenSE->sendData();
    if (senderOddOddNE)
        senderOddOddNE->sendData();
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::prepareForReceive()
{
    if (receiverEvenEvenSW)
        receiverEvenEvenSW->prepareForReceive();
    if (receiverEvenOddNW)
        receiverEvenOddNW->prepareForReceive();
    if (receiverOddEvenSE)
        receiverOddEvenSE->prepareForReceive();
    if (receiverOddOddNE)
        receiverOddOddNE->prepareForReceive();
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::receiveVectors()
{
    if (receiverEvenEvenSW)
        receiverEvenEvenSW->receiveData();
    if (receiverEvenOddNW)
        receiverEvenOddNW->receiveData();
    if (receiverOddEvenSE)
        receiverOddEvenSE->receiveData();
    if (receiverOddOddNE)
        receiverOddOddNE->receiveData();
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::init()
{
    using namespace D3Q27System;
    using namespace vf::lbm::dir;

    bMaxX1 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1();
    bMaxX2 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2();
    bMaxX3 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3();

    int sendSize      = 0;
    real initValue = -999.0;

    int sendDataPerNode = 27 /*f*/;
    int iCellSize       = 8; // size of interpolation cell

    switch (this->sendDir) {
        case dP00:
        case dM00:
            sendSize = bMaxX2 * bMaxX3 * sendDataPerNode * iCellSize;
            break;
        case d0P0:
        case d0M0:
            sendSize = bMaxX1 * bMaxX3 * sendDataPerNode * iCellSize;
            break;
        case d00P:
        case d00M:
            sendSize = bMaxX1 * bMaxX2 * sendDataPerNode * iCellSize;
            break;
        case dPP0:
        case dMM0:
        case dPM0:
        case dMP0:
            sendSize = 2 * bMaxX3 * sendDataPerNode * iCellSize;
            break;
        case dP0P:
        case dM0M:
        case dP0M:
        case dM0P:
            sendSize = 2 * bMaxX2 * sendDataPerNode * iCellSize;
            break;
        case d0PP:
        case d0MM:
        case d0PM:
        case d0MP:
            sendSize = 2 * bMaxX1 * sendDataPerNode * iCellSize;
            break;
        case dPPP:
        case dMPP:
        case dPMP:
        case dMMP:
        case dPPM:
        case dMPM:
        case dPMM:
        case dMMM:
            sendSize = 6 * bMaxX1 * sendDataPerNode * iCellSize;
            break;
        default:
            throw UbException(UB_EXARGS, "direction not allowed in this constructor");
    }
    if (senderEvenEvenSW)
        senderEvenEvenSW->getData().resize(sendSize, initValue);
    else
        senderEvenEvenSW = VectorTransmitterPtr(new TbLocalTransmitter<CbVector<real>>());
    if (senderEvenOddNW)
        senderEvenOddNW->getData().resize(sendSize, initValue);
    else
        senderEvenOddNW = VectorTransmitterPtr(new TbLocalTransmitter<CbVector<real>>());
    if (senderOddEvenSE)
        senderOddEvenSE->getData().resize(sendSize, initValue);
    else
        senderOddEvenSE = VectorTransmitterPtr(new TbLocalTransmitter<CbVector<real>>());
    if (senderOddOddNE)
        senderOddOddNE->getData().resize(sendSize, initValue);
    else
        senderOddOddNE = VectorTransmitterPtr(new TbLocalTransmitter<CbVector<real>>());

    if (!receiverEvenEvenSW)
        receiverEvenEvenSW = VectorTransmitterPtr(new TbLocalTransmitter<CbVector<real>>());
    if (!receiverEvenOddNW)
        receiverEvenOddNW = VectorTransmitterPtr(new TbLocalTransmitter<CbVector<real>>());
    if (!receiverOddEvenSE)
        receiverOddEvenSE = VectorTransmitterPtr(new TbLocalTransmitter<CbVector<real>>());
    if (!receiverOddOddNE)
        receiverOddOddNE = VectorTransmitterPtr(new TbLocalTransmitter<CbVector<real>>());

    // findCFnodes();
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::fillSendVectors()
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

    int indexEvEv = 0;
    int indexEvOd = 0;
    int indexOdEv = 0;
    int indexOdOd = 0;

    vector_type &dataEvEv = this->senderEvenEvenSW->getData();
    vector_type &dataEvOd = this->senderEvenOddNW->getData();
    vector_type &dataOdEv = this->senderOddEvenSE->getData();
    vector_type &dataOdOd = this->senderOddOddNE->getData();

    int lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3;

    switch (sendDir) {
        case dP00:
            lMinX1 = maxX1 - 3;
            lMaxX1 = lMinX1 + 1;

            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
            break;
        case dM00:
            lMinX1 = 1;
            lMaxX1 = lMinX1 + 1;

            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
            break;
        case d0P0:
            lMinX2 = maxX2 - 3;
            lMaxX2 = lMinX2 + 1;

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
            break;
        case d0M0:
            lMinX2 = 1;
            lMaxX2 = lMinX2 + 1;

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
            break;
        case d00P:
            lMinX3 = maxX3 - 3;
            lMaxX3 = lMinX3 + 1;

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
            break;
        case d00M:
            lMinX3 = 1;
            lMaxX3 = lMinX3 + 1;

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
            break;
            /// N-S-E-W
        case dPP0:
            lMinX1 = maxX1 - 3;
            lMaxX1 = lMinX1 + 2;
            lMinX2 = maxX2 - 3;
            lMaxX2 = lMinX2 + 2;

            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case dMM0:
            lMinX1 = 0;
            lMaxX1 = lMinX1 + 2;
            lMinX2 = 0;
            lMaxX2 = lMinX2 + 2;

            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case dPM0:
            lMinX1 = maxX1 - 3;
            lMaxX1 = lMinX1 + 2;
            lMinX2 = 0;
            lMaxX2 = lMinX2 + 2;

            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case dMP0:
            lMinX1 = 0;
            lMaxX1 = lMinX1 + 2;
            lMinX2 = maxX2 - 3;
            lMaxX2 = lMinX2 + 2;

            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;
            /////T-B-E-W
        case dP0P:
            lMinX1 = maxX1 - 3;
            lMaxX1 = lMinX1 + 2;
            lMinX3 = maxX3 - 3;
            lMaxX3 = lMinX3 + 2;

            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case dM0M:
            lMinX1 = 0;
            lMaxX1 = lMinX1 + 2;
            lMinX3 = 0;
            lMaxX3 = lMinX3 + 2;

            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case dP0M:
            lMinX1 = maxX1 - 3;
            lMaxX1 = lMinX1 + 2;
            lMinX3 = 0;
            lMaxX3 = lMinX3 + 2;

            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case dM0P:
            lMinX1 = 0;
            lMaxX1 = lMinX1 + 2;
            lMinX3 = maxX3 - 3;
            lMaxX3 = lMinX3 + 2;

            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;
            ////
            /////T-B-N-S
        case d0PP:
            lMinX2 = maxX2 - 3;
            lMaxX2 = lMinX2 + 2;
            lMinX3 = maxX3 - 3;
            lMaxX3 = lMinX3 + 2;

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case d0MM:
            lMinX2 = 0;
            lMaxX2 = lMinX2 + 2;
            lMinX3 = 0;
            lMaxX3 = lMinX3 + 2;

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case d0PM:
            lMinX2 = maxX2 - 3;
            lMaxX2 = lMinX2 + 2;
            lMinX3 = 0;
            lMaxX3 = lMinX3 + 2;

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case d0MP:
            lMinX2 = 0;
            lMaxX2 = lMinX2 + 2;
            lMinX3 = maxX3 - 3;
            lMaxX3 = lMinX3 + 2;

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

            // TNE
        case dPPP:
            lMinX1 = maxX1 - 3;
            lMaxX1 = maxX1 - 1;
            lMinX2 = maxX2 - 3;
            lMaxX2 = maxX2 - 1;
            lMinX3 = maxX3 - 3;
            lMaxX3 = maxX3 - 1;

            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
            break;
            //   TNW
        case dMPP:
            lMinX1 = 0;
            lMaxX1 = 2;
            lMinX2 = maxX2 - 3;
            lMaxX2 = maxX2 - 1;
            lMinX3 = maxX3 - 3;
            lMaxX3 = maxX3 - 1;

            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
            break;
            //   TSE
        case dPMP:
            lMinX1 = maxX1 - 3;
            lMaxX1 = maxX1 - 1;
            lMinX2 = 0;
            lMaxX2 = 2;
            lMinX3 = maxX3 - 3;
            lMaxX3 = maxX3 - 1;

            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
            break;
            //   TSW
        case dMMP:
            lMinX1 = 0;
            lMaxX1 = 2;
            lMinX2 = 0;
            lMaxX2 = 2;
            lMinX3 = maxX3 - 3;
            lMaxX3 = maxX3 - 1;

            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
            break;
            //   BNE
        case dPPM:
            lMinX1 = maxX1 - 3;
            lMaxX1 = maxX1 - 1;
            lMinX2 = maxX2 - 3;
            lMaxX2 = maxX2 - 1;
            lMinX3 = 0;
            lMaxX3 = 2;

            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
            break;
            //   BNW
        case dMPM:
            lMinX1 = 0;
            lMaxX1 = 2;
            lMinX2 = maxX2 - 3;
            lMaxX2 = maxX2 - 1;
            lMinX3 = 0;
            lMaxX3 = 2;

            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
            break;
            //   BSE
        case dPMM:
            lMinX1 = maxX1 - 3;
            lMaxX1 = maxX1 - 1;
            lMinX2 = 0;
            lMaxX2 = 2;
            lMinX3 = 0;
            lMaxX3 = 2;

            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
            break;
            //   BSW
        case dMMM:
            lMinX1 = 0;
            lMaxX1 = 2;
            lMinX2 = 0;
            lMaxX2 = 2;
            lMinX3 = 0;
            lMaxX3 = 2;

            fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
            break;
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::getLocalMinMax(const int &gMin, const int &gMax, const bool &even,
                                                                    int &lMin, int &lMax, const bool &dataDistribution)
{
    int halfEven = 0;
    int halfOdd  = 0;
    int dCoef    = 0;

    if (dataDistribution)
        dCoef = 1;

    if (Utilities::isOdd(gMax)) {
        halfEven = gMax / 2;
        halfOdd  = gMax / 2;
    }
    if (Utilities::isEven(gMax)) {
        halfEven = gMax / 2;
        halfOdd  = gMax / 2 - 1 + dCoef;
    }

    if (even) {
        lMin = gMin + dCoef;
        lMax = lMin + halfEven - dCoef;
    } else {
        lMin = gMin + halfOdd;
        lMax = gMax - 1;
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::fillSendVectorExt(SPtr<DistributionArray3D> fFrom,
                                                                       const int &lMinX1, const int &lMinX2,
                                                                       const int &lMinX3, const int &lMaxX1,
                                                                       const int &lMaxX2, const int &lMaxX3,
                                                                       vector_type &data, int &index)
{
    using namespace vf::basics::constant;

    if (data.size() == 0)
        return;
    int ix1, ix2, ix3;
    real xoff, yoff, zoff;
    SPtr<BCArray3D> bcArray = block.lock()->getKernel()->getBCSet()->getBCArray();

    for (ix3 = lMinX3; ix3 < lMaxX3; ix3++) {
        for (ix2 = lMinX2; ix2 < lMaxX2; ix2++) {
            for (ix1 = lMinX1; ix1 < lMaxX1; ix1++) {
                D3Q27ICell icellC;
                D3Q27ICell icellF;

                int howManySolids = iprocessor->iCellHowManySolids(bcArray, ix1, ix2, ix3);

                if (howManySolids == 0 || howManySolids == 8) {
                    iprocessor->readICell(fFrom, icellC, ix1, ix2, ix3);
                    xoff = c0o1;
                    yoff = c0o1;
                    zoff = c0o1;
                } else {
                    if (!iprocessor->findNeighborICell(bcArray, fFrom, icellC, bMaxX1, bMaxX2, bMaxX3, ix1, ix2, ix3,
                                                       xoff, yoff, zoff)) {
                        std::string err = "For " + block.lock()->toString() + " x1=" + UbSystem::toString(ix1) +
                                          ", x2=" + UbSystem::toString(ix2) + ", x3=" + UbSystem::toString(ix3) +
                                          " interpolation is not implemented for other direction" +
                                          " by using in: " + (std::string) typeid(*this).name() +
                                          " or maybe you have a solid on the block boundary";
                        UB_THROW(UbException(UB_EXARGS, err));
                    }
                }

                iprocessor->interpolateCoarseToFine(icellC, icellF, xoff, yoff, zoff);
                this->writeICellFtoData(data, index, icellF);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::writeICellFtoData(vector_type &data, int &index,
                                                                       D3Q27ICell &icellF)
{
    writeNodeToVector(data, index, icellF.BSW);
    writeNodeToVector(data, index, icellF.BSE);
    writeNodeToVector(data, index, icellF.BNW);
    writeNodeToVector(data, index, icellF.BNE);
    writeNodeToVector(data, index, icellF.TSW);
    writeNodeToVector(data, index, icellF.TSE);
    writeNodeToVector(data, index, icellF.TNW);
    writeNodeToVector(data, index, icellF.TNE);
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::writeNodeToVector(vector_type &data, int &index, real *inode)
{
    for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF + 1; i++) {
        data[index++] = inode[i];
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::distributeReceiveVectors()
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

    int indexEvEv = 0;
    int indexEvOd = 0;
    int indexOdEv = 0;
    int indexOdOd = 0;

    vector_type &dataEvEv = this->receiverEvenEvenSW->getData();
    vector_type &dataEvOd = this->receiverEvenOddNW->getData();
    vector_type &dataOdEv = this->receiverOddEvenSE->getData();
    vector_type &dataOdOd = this->receiverOddOddNE->getData();

    int lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3;
    int dummy;

    switch (sendDir) {
        case dP00:
            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 1;
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, lMinX2, lMinX3, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, lMinX3, dummy, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
            break;
        case dM00:
            lMinX1 = 3;
            lMaxX1 = lMinX1 + 1;
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, lMinX2, lMinX3, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, lMinX3, dummy, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
            break;
        case d0P0:
            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 1;
            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(lMinX1, dummy, lMinX3, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, lMinX3, lMaxX1, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
            break;
        case d0M0:
            lMinX2 = 3;
            lMaxX2 = lMinX2 + 1;
            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(lMinX1, dummy, lMinX3, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, lMinX3, lMaxX1, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
            break;
        case d00P:
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 1;
            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(lMinX1, lMinX2, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(lMinX1, dummy, dummy, dummy, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, lMinX2, dummy, lMaxX1, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, dummy, dummy, lMaxX1, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
            break;
        case d00M:
            lMinX3 = 3;
            lMaxX3 = lMinX3 + 1;
            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(lMinX1, lMinX2, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(lMinX1, dummy, dummy, dummy, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, lMinX2, dummy, lMaxX1, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, dummy, dummy, lMaxX1, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
            break;

            //    /////E-W-N-S
        case dPP0:
            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 3;
            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 1;
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 1;
            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 3;
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case dMM0:
            lMinX1 = 1;
            lMaxX1 = lMinX1 + 3;
            lMinX2 = 3;
            lMaxX2 = lMinX2 + 1;
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            lMinX1 = 3;
            lMaxX1 = lMinX1 + 1;
            lMinX2 = 1;
            lMaxX2 = lMinX2 + 3;
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case dPM0:
            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 3;
            lMinX2 = 3;
            lMaxX2 = lMinX2 + 1;
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 1;
            lMinX2 = 1;
            lMaxX2 = lMinX2 + 3;
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case dMP0:
            lMinX1 = 1;
            lMaxX1 = lMinX1 + 3;
            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 1;
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            lMinX1 = 3;
            lMaxX1 = lMinX1 + 1;
            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 3;
            getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;
            //        /////T-B-E-W
        case dP0P:
            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 3;
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 1;
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 1;
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case dM0M:
            lMinX1 = 1;
            lMaxX1 = lMinX1 + 3;
            lMinX3 = 3;
            lMaxX3 = lMinX3 + 1;
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            lMinX1 = 3;
            lMaxX1 = lMinX1 + 1;
            lMinX3 = 1;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case dP0M:
            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 3;
            lMinX3 = 3;
            lMaxX3 = lMinX3 + 1;
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            lMinX1 = maxX1 - 4;
            lMaxX1 = lMinX1 + 1;
            lMinX3 = 1;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case dM0P:
            lMinX1 = 1;
            lMaxX1 = lMinX1 + 3;
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 1;
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            lMinX1 = 3;
            lMaxX1 = lMinX1 + 1;
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
            getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

            /////////////////////////T-N-B-S
        case d0PP:
            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 3;
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 1;
            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 1;
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case d0MM:
            lMinX2 = 1;
            lMaxX2 = lMinX2 + 3;
            lMinX3 = 3;
            lMaxX3 = lMinX3 + 1;
            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            lMinX2 = 3;
            lMaxX2 = lMinX2 + 1;
            lMinX3 = 1;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case d0PM:
            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 3;
            lMinX3 = 3;
            lMaxX3 = lMinX3 + 1;
            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            lMinX2 = maxX2 - 4;
            lMaxX2 = lMinX2 + 1;
            lMinX3 = 1;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

        case d0MP:
            lMinX2 = 1;
            lMaxX2 = lMinX2 + 3;
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 1;
            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            lMinX2 = 3;
            lMaxX2 = lMinX2 + 1;
            lMinX3 = maxX3 - 4;
            lMaxX3 = lMinX3 + 3;
            getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
            getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
            getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

            break;

            // TNE
        case dPPP:
            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 3;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 1;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 1;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 1;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 3;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 1;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 1;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 1;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 3;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            break;
            //   TNW
        case dMPP:
            lMinX1 = 3;
            lMaxX1 = 4;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 1;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 1;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = 1;
            lMaxX1 = 4;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 3;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 1;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = 1;
            lMaxX1 = 4;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 1;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 3;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            break;
            //   TSE
        case dPMP:
            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 3;
            lMinX2 = 1;
            lMaxX2 = 4;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 1;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 1;
            lMinX2 = 3;
            lMaxX2 = 4;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 1;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 1;
            lMinX2 = 1;
            lMaxX2 = 4;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 3;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
            break;
            //   TSW
        case dMMP:
            lMinX1 = 3;
            lMaxX1 = 4;
            lMinX2 = 1;
            lMaxX2 = 4;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 1;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = 1;
            lMaxX1 = 4;
            lMinX2 = 3;
            lMaxX2 = 4;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 1;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = 1;
            lMaxX1 = 4;
            lMinX2 = 1;
            lMaxX2 = 4;
            lMinX3 = maxX3 - 4;
            lMaxX3 = maxX3 - 3;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
            break;
            //   BNE
        case dPPM:
            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 3;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 1;
            lMinX3 = 1;
            lMaxX3 = 4;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 1;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 3;
            lMinX3 = 1;
            lMaxX3 = 4;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 1;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 1;
            lMinX3 = 3;
            lMaxX3 = 4;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            break;
            //   BNW
        case dMPM:
            lMinX1 = 3;
            lMaxX1 = 4;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 1;
            lMinX3 = 1;
            lMaxX3 = 4;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = 1;
            lMaxX1 = 4;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 3;
            lMinX3 = 1;
            lMaxX3 = 4;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = 1;
            lMaxX1 = 4;
            lMinX2 = maxX2 - 4;
            lMaxX2 = maxX2 - 1;
            lMinX3 = 3;
            lMaxX3 = 4;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
            break;
            //   BSE
        case dPMM:
            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 3;
            lMinX2 = 1;
            lMaxX2 = 4;
            lMinX3 = 1;
            lMaxX3 = 4;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 1;
            lMinX2 = 3;
            lMaxX2 = 4;
            lMinX3 = 1;
            lMaxX3 = 4;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = maxX1 - 4;
            lMaxX1 = maxX1 - 1;
            lMinX2 = 1;
            lMaxX2 = 4;
            lMinX3 = 3;
            lMaxX3 = 4;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
            break;
            //   BSW
        case dMMM:
            lMinX1 = 3;
            lMaxX1 = 4;
            lMinX2 = 1;
            lMaxX2 = 4;
            lMinX3 = 1;
            lMaxX3 = 4;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = 1;
            lMaxX1 = 4;
            lMinX2 = 3;
            lMaxX2 = 4;
            lMinX3 = 1;
            lMaxX3 = 4;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            lMinX1 = 1;
            lMaxX1 = 4;
            lMinX2 = 1;
            lMaxX2 = 4;
            lMinX3 = 3;
            lMaxX3 = 4;
            distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

            break;
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::distributeReceiveVector(SPtr<DistributionArray3D> fTo,
                                                                             const int &lMinX1, const int &lMinX2,
                                                                             const int &lMinX3, const int &lMaxX1,
                                                                             const int &lMaxX2, const int &lMaxX3,
                                                                             vector_type &data, int &index)
{
    if (data.size() == 0)
        return;

    int ix1, ix2, ix3;
    for (ix3 = lMinX3; ix3 < lMaxX3; ix3++) {
        for (ix2 = lMinX2; ix2 < lMaxX2; ix2++) {
            for (ix1 = lMinX1; ix1 < lMaxX1; ix1++) {
                real icellC[27];
                this->readICellCfromData(data, index, icellC);
                iprocessor->writeINodeInv(fTo, icellC, ix1, ix2, ix3);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::readICellCfromData(vector_type &data, int &index, real *icellC)
{
    for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF + 1; i++) {
        icellC[i] = data[index++];
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::getLocalMinMax(int &minX1, int &minX2, int &minX3, int &maxX1,
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

    if (block.lock()->hasInterpolationFlagCF(dP00)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dM00)) {
        if (minX1 == TminX1)
            minX1 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d0P0)) {
        if (maxX2 == TmaxX2)
            maxX2 -= 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d0M0)) {
        if (minX2 == TminX2)
            minX2 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d00P)) {
        if (maxX3 == TmaxX3)
            maxX3 -= 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d00M)) {
        if (minX3 == TminX3)
            minX3 += 2;
    }

    // E-W-N-S
    if (block.lock()->hasInterpolationFlagCF(dPP0) && !block.lock()->hasInterpolationFlagCF(d0P0) &&
        !block.lock()->hasInterpolationFlagCF(dP00)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 2;
        if (maxX2 == TmaxX2)
            maxX2 -= 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dMM0) && !block.lock()->hasInterpolationFlagCF(dM00) &&
        !block.lock()->hasInterpolationFlagCF(d0M0)) {
        if (minX1 == TminX1)
            minX1 += 2;
        if (minX2 == TminX2)
            minX2 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dPM0) && !block.lock()->hasInterpolationFlagCF(dP00) &&
        !block.lock()->hasInterpolationFlagCF(d0M0)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 2;
        if (minX2 == TminX2)
            minX2 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dMP0) && !block.lock()->hasInterpolationFlagCF(d0P0) &&
        !block.lock()->hasInterpolationFlagCF(dM00)) {
        if (minX1 == TminX1)
            minX1 += 2;
        if (maxX2 == TmaxX2)
            maxX2 -= 2;
    }

    //    ////T-B-E-W
    if (block.lock()->hasInterpolationFlagCF(dP0P) && !block.lock()->hasInterpolationFlagCF(dP00) &&
        !block.lock()->hasInterpolationFlagCF(d00P)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 2;
        if (maxX3 == TmaxX3)
            maxX3 -= 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dM0M) && !block.lock()->hasInterpolationFlagCF(dM00) &&
        !block.lock()->hasInterpolationFlagCF(d00M)) {
        if (minX1 == TminX1)
            minX1 += 2;
        if (minX3 == TminX3)
            minX3 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dP0M) && !block.lock()->hasInterpolationFlagCF(dP00) &&
        !block.lock()->hasInterpolationFlagCF(d00M)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 2;
        if (minX3 == TminX3)
            minX3 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dM0P) && !block.lock()->hasInterpolationFlagCF(dM00) &&
        !block.lock()->hasInterpolationFlagCF(d00P)) {
        if (minX1 == TminX1)
            minX1 += 2;
        if (maxX3 == TmaxX3)
            maxX3 -= 2;
    }

    ////T-B-N-S
    if (block.lock()->hasInterpolationFlagCF(d0PP) && !block.lock()->hasInterpolationFlagCF(d0P0) &&
        !block.lock()->hasInterpolationFlagCF(d00P)) {
        if (maxX2 == TmaxX2)
            maxX2 -= 2;
        if (maxX3 == TmaxX3)
            maxX3 -= 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d0MM) && !block.lock()->hasInterpolationFlagCF(d0M0) &&
        !block.lock()->hasInterpolationFlagCF(d00M)) {
        if (minX2 == TminX2)
            minX2 += 2;
        if (minX3 == TminX3)
            minX3 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d0PM) && !block.lock()->hasInterpolationFlagCF(d0P0) &&
        !block.lock()->hasInterpolationFlagCF(d00M)) {
        if (maxX2 == TmaxX2)
            maxX2 -= 2;
        if (minX3 == TminX3)
            minX3 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d0MP) && !block.lock()->hasInterpolationFlagCF(d0M0) &&
        !block.lock()->hasInterpolationFlagCF(d00P)) {
        if (minX2 == TminX2)
            minX2 += 2;
        if (maxX3 == TmaxX3)
            maxX3 -= 2;
    }

    // if
    // (block.lock()->hasInterpolationFlagCF(D3Q27System::dPPP)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::dP0P)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::d0PP)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::dPP0)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::d00P)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::d0P0)
    // && !block.lock()->hasInterpolationFlagCF(D3Q27System::dP00)) if
    // (!block.lock()->hasInterpolationFlagCF(D3Q27System::dP0P)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::d00P)&&
    // !block.lock()->hasInterpolationFlagCF(D3Q27System::dP00))
    //{
    //   if (maxX1==TmaxX1) maxX1 -= 2;
    //   if (maxX2==TmaxX2) maxX2 -= 2;
    //   if (maxX3==TmaxX3) maxX3 -= 2;
    //}
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::getLocalMinMax(int &minX1, int &minX2, int &minX3, int &maxX1,
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

    if (block.lock()->hasInterpolationFlagCF(dP00)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dM00)) {
        if (minX1 == TminX1)
            minX1 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d0P0)) {
        if (maxX2 == TmaxX2)
            maxX2 -= 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d0M0)) {
        if (minX2 == TminX2)
            minX2 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d00P)) {
        if (maxX3 == TmaxX3)
            maxX3 -= 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d00M)) {
        if (minX3 == TminX3)
            minX3 += 2;
    }

    // E-W-N-S
    if (block.lock()->hasInterpolationFlagCF(dPP0) && !block.lock()->hasInterpolationFlagCF(d0P0) &&
        !block.lock()->hasInterpolationFlagCF(dP00)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 2;
        if (maxX2 == TmaxX2)
            maxX2 -= 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dMM0) && !block.lock()->hasInterpolationFlagCF(dM00) &&
        !block.lock()->hasInterpolationFlagCF(d0M0)) {
        if (minX1 == TminX1)
            minX1 += 2;
        if (minX2 == TminX2)
            minX2 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dPM0) && !block.lock()->hasInterpolationFlagCF(dP00) &&
        !block.lock()->hasInterpolationFlagCF(d0M0)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 2;
        if (minX2 == TminX2)
            minX2 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dMP0) && !block.lock()->hasInterpolationFlagCF(d0P0) &&
        !block.lock()->hasInterpolationFlagCF(dM00)) {
        if (minX1 == TminX1)
            minX1 += 2;
        if (maxX2 == TmaxX2)
            maxX2 -= 2;
    }

    //    ////T-B-E-W
    if (block.lock()->hasInterpolationFlagCF(dP0P) && !block.lock()->hasInterpolationFlagCF(dP00) &&
        !block.lock()->hasInterpolationFlagCF(d00P)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 2;
        if (maxX3 == TmaxX3)
            maxX3 -= 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dM0M) && !block.lock()->hasInterpolationFlagCF(dM00) &&
        !block.lock()->hasInterpolationFlagCF(d00M)) {
        if (minX1 == TminX1)
            minX1 += 2;
        if (minX3 == TminX3)
            minX3 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dP0M) && !block.lock()->hasInterpolationFlagCF(dP00) &&
        !block.lock()->hasInterpolationFlagCF(d00M)) {
        if (maxX1 == TmaxX1)
            maxX1 -= 2;
        if (minX3 == TminX3)
            minX3 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(dM0P) && !block.lock()->hasInterpolationFlagCF(dM00) &&
        !block.lock()->hasInterpolationFlagCF(d00P)) {
        if (minX1 == TminX1)
            minX1 += 2;
        if (maxX3 == TmaxX3)
            maxX3 -= 2;
    }

    ////T-B-N-S
    if (block.lock()->hasInterpolationFlagCF(d0PP) && !block.lock()->hasInterpolationFlagCF(d0P0) &&
        !block.lock()->hasInterpolationFlagCF(d00P)) {
        if (maxX2 == TmaxX2)
            maxX2 -= 2;
        if (maxX3 == TmaxX3)
            maxX3 -= 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d0MM) && !block.lock()->hasInterpolationFlagCF(d0M0) &&
        !block.lock()->hasInterpolationFlagCF(d00M)) {
        if (minX2 == TminX2)
            minX2 += 2;
        if (minX3 == TminX3)
            minX3 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d0PM) && !block.lock()->hasInterpolationFlagCF(d0P0) &&
        !block.lock()->hasInterpolationFlagCF(d00M)) {
        if (maxX2 == TmaxX2)
            maxX2 -= 2;
        if (minX3 == TminX3)
            minX3 += 2;
    }
    if (block.lock()->hasInterpolationFlagCF(d0MP) && !block.lock()->hasInterpolationFlagCF(d0M0) &&
        !block.lock()->hasInterpolationFlagCF(d00P)) {
        if (minX2 == TminX2)
            minX2 += 2;
        if (maxX3 == TmaxX3)
            maxX3 -= 2;
    }

    // if
    // (block.lock()->hasInterpolationFlagCF(D3Q27System::dPPP)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::dP0P)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::d0PP)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::dPP0)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::d00P)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::d0P0)
    // && !block.lock()->hasInterpolationFlagCF(D3Q27System::dP00))
    //{
    //   if (maxX1==TmaxX1) maxX1 -= 2;
    //   if (maxX2==TmaxX2) maxX2 -= 2;
    //   if (maxX3==TmaxX3) maxX3 -= 2;
    //}
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::findCFnodes()
{
    SPtr<DistributionArray3D> fFrom = block.lock()->getKernel()->getDataSet()->getFdistributions();
    int maxX1                       = (int)fFrom->getNX1();
    int maxX2                       = (int)fFrom->getNX2();
    int maxX3                       = (int)fFrom->getNX3();
    int minX1                       = 0;
    int minX2                       = 0;
    int minX3                       = 0;

    int indexEvEv = 0;
    int indexEvOd = 0;
    int indexOdEv = 0;
    int indexOdOd = 0;

    vector_type &dataEvEv = this->senderEvenEvenSW->getData();
    vector_type &dataEvOd = this->senderEvenOddNW->getData();
    vector_type &dataOdEv = this->senderOddEvenSE->getData();
    vector_type &dataOdOd = this->senderOddOddNE->getData();

    int lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3;

    using namespace D3Q27System;
    using namespace vf::lbm::dir;

    if (block.lock()->hasInterpolationFlagCF(dM00)) {
        lMinX1 = 1;
        lMaxX1 = lMinX1 + 1;

        getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
        getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
        findCFnodes(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

        getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
        getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
        findCFnodes(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

        getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
        getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
        findCFnodes(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

        getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
        getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
        findCFnodes(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
    }
    if (block.lock()->hasInterpolationFlagCF(d0PP) && !block.lock()->hasInterpolationFlagCF(d0P0) &&
        !block.lock()->hasInterpolationFlagCF(d00P)) {
        lMinX2 = maxX2 - 3;
        lMaxX2 = lMinX2 + 1;
        lMinX3 = maxX3 - 3;
        lMaxX3 = lMinX3 + 1;

        getLocalMinMax(minX1 + 1, maxX1, true, lMinX1, lMaxX1, false);
        findCFnodes(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

        getLocalMinMax(minX1 + 1, maxX1, false, lMinX1, lMaxX1, false);
        findCFnodes(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
void CoarseToFineVectorConnector<VectorTransmitter>::findCFnodes(SPtr<DistributionArray3D> fFrom, const int &lMinX1,
                                                                 const int &lMinX2, const int &lMinX3,
                                                                 const int &lMaxX1, const int &lMaxX2,
                                                                 const int &lMaxX3, vector_type &data, int &index)
{
    using namespace vf::basics::constant;
    if (data.size() == 0)
        return;
    int ix1, ix2, ix3;
    real xoff, yoff, zoff;
    SPtr<BCArray3D> bcArray = block.lock()->getKernel()->getBCSet()->getBCArray();

    for (ix3 = lMinX3; ix3 < lMaxX3; ix3++) {
        for (ix2 = lMinX2; ix2 < lMaxX2; ix2++) {
            for (ix1 = lMinX1; ix1 < lMaxX1; ix1++) {
                D3Q27ICell icellC;
                D3Q27ICell icellF;

                int howManySolids = iprocessor->iCellHowManySolids(bcArray, ix1, ix2, ix3);

                if (howManySolids == 0 || howManySolids == 8) {
                    iprocessor->readICell(fFrom, icellC, ix1, ix2, ix3);
                    xoff = c0o1;
                    yoff = c0o1;
                    zoff = c0o1;
                } else {
                    if (!iprocessor->findNeighborICell(bcArray, fFrom, icellC, bMaxX1, bMaxX2, bMaxX3, ix1, ix2, ix3,
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

                iprocessor->interpolateCoarseToFine(icellC, icellF, xoff, yoff, zoff);
                this->writeICellFtoData(data, index, icellF);
                // for (int iix3 = ix3; iix3<=ix3+1; iix3++)
                //{
                //   for (int iix2 = ix2; iix2<=ix2+1; iix2++)
                //   {
                //      for (int iix1 = ix1; iix1<=ix1+1; iix1++)
                //      {
                //         bcArray->setInterfaceCF(iix1, iix2, iix3);
                //      }
                //   }
                //}
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
template <typename VectorTransmitter>
real CoarseToFineVectorConnector<VectorTransmitter>::getSendRecieveTime()
{
    return 0;
}

#endif

//! \}
