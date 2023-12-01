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
//! \file CreateTransmittersHelper.cpp
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================

#include "CreateTransmittersHelper.h"
#include <parallel/Communicator.h>
#include <D3Q27System.h>
#include <string>

#ifdef VF_FETOL
#include <FETOLTransmitterBondPool.h>
#endif
#include <MathUtil.hpp>

unsigned CreateTransmittersHelper::vKey = 0;

using namespace std;

//////////////////////////////////////////////////////////////////////////
CreateTransmittersHelper::CreateTransmittersHelper() = default;
//////////////////////////////////////////////////////////////////////////
void CreateTransmittersHelper::createTransmitters(SPtr<Block3D> sblock, SPtr<Block3D> tblock, int dir, IBlock ib,
                                                  TransmitterPtr &sender, TransmitterPtr &receiver,
                                                  std::shared_ptr<vf::parallel::Communicator> comm, TransmitterType tType)
{
    // SourceBlock
    int srcLevel = sblock->getLevel();
    //   int srcID    = sblock->getGlobalID();

    // TargetBlock
    int tgtLevel = tblock->getLevel();
    //   int tgtID    = tblock->getGlobalID();

    int invDir = D3Q27System::INVDIR[dir];

    if (srcLevel != tgtLevel)
        invDir = dir;

    int srcRank = 0;
    int tgtRank = 0;

    if (tType == MPI) {
        srcRank = sblock->getRank();
        tgtRank = tblock->getRank();
    }
#ifdef VF_FETOL
    else if (tType == MPI2BOND) {
        srcRank = sblock->getLocalRank();
        tgtRank = tblock->getLocalRank();
    }
#endif

    if (tType == MPI
#ifdef VF_FETOL
        || tType == MPI2BOND
#endif
    ) {
        string sendPoolKey    = generatePoolKey(srcRank, srcLevel, tgtRank, tgtLevel);
        string receivePoolKey = generatePoolKey(tgtRank, tgtLevel, srcRank, srcLevel);

        TbCbVectorMpiPool<real>::MpiPoolPtr sendPool = TbCbVectorMpiPool<real>::getTbCbVectorMpiPool(sendPoolKey);
        TbCbVectorMpiPool<real>::MpiPoolPtr recvPool =
            TbCbVectorMpiPool<real>::getTbCbVectorMpiPool(receivePoolKey);

        MPI_Comm mpi_comm = *((MPI_Comm *)comm->getNativeCommunicator());

        if (!sendPool)
            sendPool = TbCbVectorMpiPool<real>::createTbCbVectorMpiPool(
                sendPoolKey, tgtRank, generateMPITag(srcLevel, tgtLevel), mpi_comm);
        if (!recvPool)
            recvPool = TbCbVectorMpiPool<real>::createTbCbVectorMpiPool(
                receivePoolKey, tgtRank, generateMPITag(tgtLevel, srcLevel), mpi_comm);

        TbCbVectorMpiPool<real>::CbVectorKey keyOfSendCbVectorKey =
            generateVectorKey(sblock->getX1(), sblock->getX2(), sblock->getX3() /*tgtID*/, dir, ib);
        TbCbVectorMpiPool<real>::CbVectorKey keyOfRecvCbVectorKey =
            generateVectorKey(tblock->getX1(), tblock->getX2(), tblock->getX3() /*srcID*/, invDir, ib);

        ////////////////////////////////////////////////////////
        // DEBUG
        // int myid = comm->getProcessID();
        // FILE * file;
        ////char * name = "d:/temp/sendPoolKey.csv";
        // std::string name = "d:/temp/VectorKey" + UbSystem::toString(myid) + ".csv";
        // file = fopen(name.c_str(), "a");
        // fprintf(file, "%d;%d%;%d;%d;%d;%u;%d;%d%;%d;%d;%d;%u\n", sblock->getX1(), sblock->getX2(),
        // sblock->getX3()/*tgtID*/, dir, ib, keyOfSendCbVectorKey, tblock->getX1(), tblock->getX2(),
        // tblock->getX3()/*srcID*/, invDir, ib, keyOfRecvCbVectorKey); fclose(file);
        ////////////////////////////////////////////////////////

        // create sender-/receiver
        sender   = TransmitterPtr(new TbCbVectorSenderMpiPool<real>(keyOfSendCbVectorKey, sendPool.get()));
        receiver = TransmitterPtr(new TbCbVectorReceiverMpiPool<real>(keyOfRecvCbVectorKey, recvPool.get()));
    }
#ifdef VF_FETOL
    if (tType == BOND) {
        int srcBondRank = sblock->getRank();
        int tgtBondRank = tblock->getRank();

        int sendBondPoolKey    = generatePoolKey(srcBondRank, srcLevel, tgtBondRank, tgtLevel);
        int receiveBondPoolKey = generatePoolKey(tgtBondRank, tgtLevel, srcBondRank, srcLevel);

        TbCbVectorBondPool<real>::BondPoolPtr sendPool =
            TbCbVectorBondPool<real>::getTbCbVectorBondPool(sendBondPoolKey);
        TbCbVectorBondPool<real>::BondPoolPtr recvPool =
            TbCbVectorBondPool<real>::getTbCbVectorBondPool(receiveBondPoolKey);

        if (!sendPool)
            sendPool = TbCbVectorBondPool<real>::createTbCbVectorBondPool(sendBondPoolKey, tgtBondRank,
                                                                             generateMPITag(srcLevel, tgtLevel));
        if (!recvPool)
            recvPool = TbCbVectorBondPool<real>::createTbCbVectorBondPool(receiveBondPoolKey, tgtBondRank,
                                                                             generateMPITag(tgtLevel, srcLevel));

        TbCbVectorBondPool<real>::CbVectorKey keyOfSendCbVectorKey = generateVectorKey(tgtID, dir, ib);
        TbCbVectorBondPool<real>::CbVectorKey keyOfRecvCbVectorKey = generateVectorKey(srcID, invDir, ib);

        // create sender-/receiver
        sender   = TransmitterPtr(new TbCbVectorSenderBondPool<real>(keyOfSendCbVectorKey, sendPool.get()));
        receiver = TransmitterPtr(new TbCbVectorReceiverBondPool<real>(keyOfRecvCbVectorKey, recvPool.get()));
    }
#endif
}
//////////////////////////////////////////////////////////////////////////
int CreateTransmittersHelper::generateMPITag(int srcLevel, int tgtLevel)
{
    // The MPI standard guarantees that integers 0-32767 can be used as tags
    if (srcLevel == tgtLevel) {
        return srcLevel;
    } else {
        srcLevel++;
        tgtLevel++;
        std::string str = UbSystem::toString<int>(srcLevel) + UbSystem::toString<int>(tgtLevel);
        int r           = UbSystem::stringTo<int>(str);
        return r;
    }
}
//////////////////////////////////////////////////////////////////////////
string CreateTransmittersHelper::generatePoolKey(int srcRank, int srcLevel, int tgtRank, int tgtLevel)
{
    std::string str;
    str = UbSystem::toString<int>(srcLevel);
    str += "#";
    str += UbSystem::toString<int>(tgtLevel);
    str += "#";
    str += UbSystem::toString<int>(srcRank);
    str += "#";
    str += UbSystem::toString<int>(tgtRank);

    return str;
}
//////////////////////////////////////////////////////////////////////////
string CreateTransmittersHelper::generateVectorKey(int x1, int x2, int x3, int dir, IBlock ib)
{
    std::string str;

    str += UbSystem::toString<int>(x1);
    str += "#";
    str += UbSystem::toString<int>(x2);
    str += "#";
    str += UbSystem::toString<int>(x3);
    str += "#";
    str += UbSystem::toString<int>(dir);
    str += "#";
    str += UbSystem::toString<int>(ib);

    return str;
}
