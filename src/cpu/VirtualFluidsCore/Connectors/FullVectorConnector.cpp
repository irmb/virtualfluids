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
//! \file FullVectorConnector.cpp
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#include "FullVectorConnector.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "DataSet3D.h"
#include "LBMKernel.h"
//////////////////////////////////////////////////////////////////////////
FullVectorConnector::FullVectorConnector(SPtr<Block3D> block, VectorTransmitterPtr sender,
                                                       VectorTransmitterPtr receiver, int sendDir)
    : RemoteBlock3DConnector(block, sender, receiver, sendDir)
{
    if (!block || !sender || !receiver)
        UB_THROW(UbException(UB_EXARGS, "sender or receiver == NULL!!"));
}
//////////////////////////////////////////////////////////////////////////
void FullVectorConnector::init()
{
    maxX1 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1() - 1;
    maxX2 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2() - 1;
    maxX3 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3() - 1;
}
//////////////////////////////////////////////////////////////////////////
void FullVectorConnector::fillSendVectors() 
{ 
    updatePointers();
    fillData();
}
////////////////////////////////////////////////////////////////////////
void FullVectorConnector::fillData()
{
    using namespace vf::lbm::dir;

    vector_type &sdata = sender->getData();

    int index = 0;
    // EAST
    if (sendDir == DIR_P00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                fillData(sdata, index, maxX1 - 1, x2, x3);
            }
        }
    }
    // WEST
    else if (sendDir == DIR_M00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                fillData(sdata, index, 1, x2, x3);
            }
        }
    }
    // NORTH
    else if (sendDir == DIR_0P0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                fillData(sdata, index, x1, maxX2 - 1, x3);
            }
        }
    }
    // SOUTH
    else if (sendDir == DIR_0M0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                fillData(sdata, index, x1, 1, x3);
            }
        }
    }
    // TOP
    else if (sendDir == DIR_00P) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                fillData(sdata, index, x1, x2, maxX3 - 1);
            }
        }
    }
    // BOTTOM
    else if (sendDir == DIR_00M) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                fillData(sdata, index, x1, x2, 1);
            }
        }
    }
    // NE NW SW SE
    else if (sendDir == DIR_PP0 || sendDir == DIR_MP0 || sendDir == DIR_MM0 ||
             sendDir == DIR_PM0) {
        int x1 = 0;
        int x2 = 0;
        switch (sendDir) {
            case DIR_PP0:
                x1 = maxX1 - 1;
                x2 = maxX2 - 1;
                break;
            case DIR_MP0:
                x1 = 1;
                x2 = maxX2 - 1;
                break;
            case DIR_MM0:
                x1 = 1;
                x2 = 1;
                break;
            case DIR_PM0:
                x1 = maxX1 - 1;
                x2 = 1;
                break;
        }
        for (int x3 = 1; x3 < maxX3; x3++) {
            fillData(sdata, index, x1, x2, x3);
        }
    }
    // TE TW BW BE
    else if (sendDir == DIR_P0P || sendDir == DIR_M0P || sendDir == DIR_M0M ||
             sendDir == DIR_P0M) {
        int x1 = 0;
        int x3 = 0;
        switch (sendDir) {
            case DIR_P0P:
                x1 = maxX1 - 1;
                x3 = maxX3 - 1;
                break;
            case DIR_M0P:
                x1 = 1;
                x3 = maxX3 - 1;
                break;
            case DIR_M0M:
                x1 = 1;
                x3 = 1;
                break;
            case DIR_P0M:
                x1 = maxX1 - 1;
                x3 = 1;
                break;
        }
        for (int x2 = 1; x2 < maxX2; x2++) {
            fillData(sdata, index, x1, x2, x3);
        }
    }
    // TN BN BS TS
    else if (sendDir == DIR_0PP || sendDir == DIR_0PM || sendDir == DIR_0MM ||
             sendDir == DIR_0MP) {
        int x2 = 0;
        int x3 = 0;
        switch (sendDir) {
            case DIR_0PP:
                x3 = maxX3 - 1;
                x2 = maxX2 - 1;
                break;
            case DIR_0PM:
                x3 = 1;
                x2 = maxX2 - 1;
                break;
            case DIR_0MM:
                x3 = 1;
                x2 = 1;
                break;
            case DIR_0MP:
                x3 = maxX3 - 1;
                x2 = 1;
                break;
        }
        for (int x1 = 1; x1 < maxX1; x1++) {
            fillData(sdata, index, x1, x2, x3);
        }
    }
    // TNE TNW TSW TSE BNE BNW BSW BSE
    else if (sendDir == DIR_PPP || sendDir == DIR_MPP || sendDir == DIR_MMP ||
             sendDir == DIR_PMP || sendDir == DIR_PPM || sendDir == DIR_MPM ||
             sendDir == DIR_MMM || sendDir == DIR_PMM) {
        int x1 = 0;
        int x2 = 0;
        int x3 = 0;
        switch (sendDir) {
            case DIR_PPP:
                x1 = maxX1 - 1;
                x2 = maxX2 - 1;
                x3 = maxX3 - 1;
                break;
            case DIR_MPP:
                x1 = 1;
                x2 = maxX2 - 1;
                x3 = maxX3 - 1;
                break;
            case DIR_MMP:
                x1 = 1;
                x2 = 1;
                x3 = maxX3 - 1;
                break;
            case DIR_PMP:
                x1 = maxX1 - 1;
                x2 = 1;
                x3 = maxX3 - 1;
                break;
            case DIR_PPM:
                x1 = maxX1 - 1;
                x2 = maxX2 - 1;
                x3 = 1;
                break;
            case DIR_MPM:
                x1 = 1;
                x2 = maxX2 - 1;
                x3 = 1;
                break;
            case DIR_MMM:
                x1 = 1;
                x2 = 1;
                x3 = 1;
                break;
            case DIR_PMM:
                x1 = maxX1 - 1;
                x2 = 1;
                x3 = 1;
                break;
        }
        fillData(sdata, index, x1, x2, x3);
    } else
        UB_THROW(UbException(UB_EXARGS, "unknown dir"));
}
////////////////////////////////////////////////////////////////////////
void FullVectorConnector::distributeReceiveVectors() 
{
    updatePointers();
    distributeData();
}
////////////////////////////////////////////////////////////////////////
void FullVectorConnector::distributeData()
{
    using namespace vf::lbm::dir;

    vector_type &rdata = receiver->getData();

    int index = 0;

    if (sendDir == DIR_M00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                distributeData(rdata, index, 0, x2, x3);
            }
        }
    } else if (sendDir == DIR_P00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                distributeData(rdata, index, maxX1, x2, x3);
            }
        }
    } else if (sendDir == DIR_0M0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                distributeData(rdata, index, x1, 0, x3);
            }
        }
    } else if (sendDir == DIR_0P0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                distributeData(rdata, index, x1, maxX2, x3);
            }
        }
    } else if (sendDir == DIR_00M) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                distributeData(rdata, index, x1, x2, 0);
            }
        }
    } else if (sendDir == DIR_00P) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                distributeData(rdata, index, x1, x2, maxX3);
            }
        }
    }
    // NE NW SW SE
    else if (sendDir == DIR_PP0 || sendDir == DIR_MP0 || sendDir == DIR_MM0 ||
             sendDir == DIR_PM0) {
        int x1 = 0;
        int x2 = 0;
        switch (sendDir) // wenn sendir NE dann kommen werte von SW
        {
            case DIR_PP0:
                x1 = maxX1;
                x2 = maxX2;
                break;
            case DIR_MP0:
                x1 = 0;
                x2 = maxX2;
                break;
            case DIR_MM0:
                x1 = 0;
                x2 = 0;
                break;
            case DIR_PM0:
                x1 = maxX1;
                x2 = 0;
                break;
        }
        for (int x3 = 1; x3 < maxX3; x3++) {
            distributeData(rdata, index, x1, x2, x3);
        }

    }
    // TE TW BW BE
    else if (sendDir == DIR_P0P || sendDir == DIR_M0P || sendDir == DIR_M0M ||
             sendDir == DIR_P0M)

    {
        int x1 = 0;
        int x3 = 0;
        switch (sendDir) // wenn sendir NE dann kommen werte von SW
        {
            case DIR_P0P:
                x1 = maxX1;
                x3 = maxX3;
                break;
            case DIR_M0P:
                x1 = 0;
                x3 = maxX3;
                break;
            case DIR_M0M:
                x1 = 0;
                x3 = 0;
                break;
            case DIR_P0M:
                x1 = maxX1;
                x3 = 0;
                break;
        }
        for (int x2 = 1; x2 < maxX2; x2++) {
            distributeData(rdata, index, x1, x2, x3);
        }
    }
    // TN BN BS TS
    else if (sendDir == DIR_0PP || sendDir == DIR_0PM || sendDir == DIR_0MM ||
             sendDir == DIR_0MP) {
        int x2 = 0;
        int x3 = 0;
        switch (sendDir) {
            case DIR_0PP:
                x3 = maxX3;
                x2 = maxX2;
                break;
            case DIR_0PM:
                x3 = 0;
                x2 = maxX2;
                break;
            case DIR_0MM:
                x3 = 0;
                x2 = 0;
                break;
            case DIR_0MP:
                x3 = maxX3;
                x2 = 0;
                break;
        }
        for (int x1 = 1; x1 < maxX1; x1++) {
            distributeData(rdata, index, x1, x2, x3);
        }
    }
    // TNE TNW TSW TSE BNE BNW BSW BSE
    else if (sendDir == DIR_PPP || sendDir == DIR_MPP || sendDir == DIR_MMP ||
             sendDir == DIR_PMP || sendDir == DIR_PPM || sendDir == DIR_MPM ||
             sendDir == DIR_MMM || sendDir == DIR_PMM) {
        int x1 = 0;
        int x2 = 0;
        int x3 = 0;

        switch (sendDir) {
            case DIR_PPP:
                x1 = maxX1;
                x2 = maxX2;
                x3 = maxX3;
                break;
            case DIR_MPP:
                x1 = 0;
                x2 = maxX2;
                x3 = maxX3;
                break;
            case DIR_MMP:
                x1 = 0;
                x2 = 0;
                x3 = maxX3;
                break;
            case DIR_PMP:
                x1 = maxX1;
                x2 = 0;
                x3 = maxX3;
                break;
            case DIR_PPM:
                x1 = maxX1;
                x2 = maxX2;
                x3 = 0;
                break;
            case DIR_MPM:
                x1 = 0;
                x2 = maxX2;
                x3 = 0;
                break;
            case DIR_MMM:
                x1 = 0;
                x2 = 0;
                x3 = 0;
                break;
            case DIR_PMM:
                x1 = maxX1;
                x2 = 0;
                x3 = 0;
                break;
        }
        distributeData(rdata, index, x1, x2, x3);
    } else
        UB_THROW(UbException(UB_EXARGS, "unknown dir"));
}
//////////////////////////////////////////////////////////////////////////
