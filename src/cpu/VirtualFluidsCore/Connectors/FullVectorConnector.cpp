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
    vector_type &sdata = sender->getData();

    int index = 0;
    // EAST
    if (sendDir == D3Q27System::DIR_P00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                fillData(sdata, index, maxX1 - 1, x2, x3);
            }
        }
    }
    // WEST
    else if (sendDir == D3Q27System::DIR_M00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                fillData(sdata, index, 1, x2, x3);
            }
        }
    }
    // NORTH
    else if (sendDir == D3Q27System::DIR_0P0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                fillData(sdata, index, x1, maxX2 - 1, x3);
            }
        }
    }
    // SOUTH
    else if (sendDir == D3Q27System::DIR_0M0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                fillData(sdata, index, x1, 1, x3);
            }
        }
    }
    // TOP
    else if (sendDir == D3Q27System::DIR_00P) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                fillData(sdata, index, x1, x2, maxX3 - 1);
            }
        }
    }
    // BOTTOM
    else if (sendDir == D3Q27System::DIR_00M) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                fillData(sdata, index, x1, x2, 1);
            }
        }
    }
    // NE NW SW SE
    else if (sendDir == D3Q27System::DIR_PP0 || sendDir == D3Q27System::DIR_MP0 || sendDir == D3Q27System::DIR_MM0 ||
             sendDir == D3Q27System::DIR_PM0) {
        int x1 = 0;
        int x2 = 0;
        switch (sendDir) {
            case D3Q27System::DIR_PP0:
                x1 = maxX1 - 1;
                x2 = maxX2 - 1;
                break;
            case D3Q27System::DIR_MP0:
                x1 = 1;
                x2 = maxX2 - 1;
                break;
            case D3Q27System::DIR_MM0:
                x1 = 1;
                x2 = 1;
                break;
            case D3Q27System::DIR_PM0:
                x1 = maxX1 - 1;
                x2 = 1;
                break;
        }
        for (int x3 = 1; x3 < maxX3; x3++) {
            fillData(sdata, index, x1, x2, x3);
        }
    }
    // TE TW BW BE
    else if (sendDir == D3Q27System::DIR_P0P || sendDir == D3Q27System::DIR_M0P || sendDir == D3Q27System::DIR_M0M ||
             sendDir == D3Q27System::DIR_P0M) {
        int x1 = 0;
        int x3 = 0;
        switch (sendDir) {
            case D3Q27System::DIR_P0P:
                x1 = maxX1 - 1;
                x3 = maxX3 - 1;
                break;
            case D3Q27System::DIR_M0P:
                x1 = 1;
                x3 = maxX3 - 1;
                break;
            case D3Q27System::DIR_M0M:
                x1 = 1;
                x3 = 1;
                break;
            case D3Q27System::DIR_P0M:
                x1 = maxX1 - 1;
                x3 = 1;
                break;
        }
        for (int x2 = 1; x2 < maxX2; x2++) {
            fillData(sdata, index, x1, x2, x3);
        }
    }
    // TN BN BS TS
    else if (sendDir == D3Q27System::DIR_0PP || sendDir == D3Q27System::DIR_0PM || sendDir == D3Q27System::DIR_0MM ||
             sendDir == D3Q27System::DIR_0MP) {
        int x2 = 0;
        int x3 = 0;
        switch (sendDir) {
            case D3Q27System::DIR_0PP:
                x3 = maxX3 - 1;
                x2 = maxX2 - 1;
                break;
            case D3Q27System::DIR_0PM:
                x3 = 1;
                x2 = maxX2 - 1;
                break;
            case D3Q27System::DIR_0MM:
                x3 = 1;
                x2 = 1;
                break;
            case D3Q27System::DIR_0MP:
                x3 = maxX3 - 1;
                x2 = 1;
                break;
        }
        for (int x1 = 1; x1 < maxX1; x1++) {
            fillData(sdata, index, x1, x2, x3);
        }
    }
    // TNE TNW TSW TSE BNE BNW BSW BSE
    else if (sendDir == D3Q27System::DIR_PPP || sendDir == D3Q27System::DIR_MPP || sendDir == D3Q27System::DIR_MMP ||
             sendDir == D3Q27System::DIR_PMP || sendDir == D3Q27System::DIR_PPM || sendDir == D3Q27System::DIR_MPM ||
             sendDir == D3Q27System::DIR_MMM || sendDir == D3Q27System::DIR_PMM) {
        int x1 = 0;
        int x2 = 0;
        int x3 = 0;
        switch (sendDir) {
            case D3Q27System::DIR_PPP:
                x1 = maxX1 - 1;
                x2 = maxX2 - 1;
                x3 = maxX3 - 1;
                break;
            case D3Q27System::DIR_MPP:
                x1 = 1;
                x2 = maxX2 - 1;
                x3 = maxX3 - 1;
                break;
            case D3Q27System::DIR_MMP:
                x1 = 1;
                x2 = 1;
                x3 = maxX3 - 1;
                break;
            case D3Q27System::DIR_PMP:
                x1 = maxX1 - 1;
                x2 = 1;
                x3 = maxX3 - 1;
                break;
            case D3Q27System::DIR_PPM:
                x1 = maxX1 - 1;
                x2 = maxX2 - 1;
                x3 = 1;
                break;
            case D3Q27System::DIR_MPM:
                x1 = 1;
                x2 = maxX2 - 1;
                x3 = 1;
                break;
            case D3Q27System::DIR_MMM:
                x1 = 1;
                x2 = 1;
                x3 = 1;
                break;
            case D3Q27System::DIR_PMM:
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
    vector_type &rdata = receiver->getData();

    int index = 0;

    if (sendDir == D3Q27System::DIR_M00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                distributeData(rdata, index, 0, x2, x3);
            }
        }
    } else if (sendDir == D3Q27System::DIR_P00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                distributeData(rdata, index, maxX1, x2, x3);
            }
        }
    } else if (sendDir == D3Q27System::DIR_0M0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                distributeData(rdata, index, x1, 0, x3);
            }
        }
    } else if (sendDir == D3Q27System::DIR_0P0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                distributeData(rdata, index, x1, maxX2, x3);
            }
        }
    } else if (sendDir == D3Q27System::DIR_00M) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                distributeData(rdata, index, x1, x2, 0);
            }
        }
    } else if (sendDir == D3Q27System::DIR_00P) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                distributeData(rdata, index, x1, x2, maxX3);
            }
        }
    }
    // NE NW SW SE
    else if (sendDir == D3Q27System::DIR_PP0 || sendDir == D3Q27System::DIR_MP0 || sendDir == D3Q27System::DIR_MM0 ||
             sendDir == D3Q27System::DIR_PM0) {
        int x1 = 0;
        int x2 = 0;
        switch (sendDir) // wenn sendir NE dann kommen werte von SW
        {
            case D3Q27System::DIR_PP0:
                x1 = maxX1;
                x2 = maxX2;
                break;
            case D3Q27System::DIR_MP0:
                x1 = 0;
                x2 = maxX2;
                break;
            case D3Q27System::DIR_MM0:
                x1 = 0;
                x2 = 0;
                break;
            case D3Q27System::DIR_PM0:
                x1 = maxX1;
                x2 = 0;
                break;
        }
        for (int x3 = 1; x3 < maxX3; x3++) {
            distributeData(rdata, index, x1, x2, x3);
        }

    }
    // TE TW BW BE
    else if (sendDir == D3Q27System::DIR_P0P || sendDir == D3Q27System::DIR_M0P || sendDir == D3Q27System::DIR_M0M ||
             sendDir == D3Q27System::DIR_P0M)

    {
        int x1 = 0;
        int x3 = 0;
        switch (sendDir) // wenn sendir NE dann kommen werte von SW
        {
            case D3Q27System::DIR_P0P:
                x1 = maxX1;
                x3 = maxX3;
                break;
            case D3Q27System::DIR_M0P:
                x1 = 0;
                x3 = maxX3;
                break;
            case D3Q27System::DIR_M0M:
                x1 = 0;
                x3 = 0;
                break;
            case D3Q27System::DIR_P0M:
                x1 = maxX1;
                x3 = 0;
                break;
        }
        for (int x2 = 1; x2 < maxX2; x2++) {
            distributeData(rdata, index, x1, x2, x3);
        }
    }
    // TN BN BS TS
    else if (sendDir == D3Q27System::DIR_0PP || sendDir == D3Q27System::DIR_0PM || sendDir == D3Q27System::DIR_0MM ||
             sendDir == D3Q27System::DIR_0MP) {
        int x2 = 0;
        int x3 = 0;
        switch (sendDir) {
            case D3Q27System::DIR_0PP:
                x3 = maxX3;
                x2 = maxX2;
                break;
            case D3Q27System::DIR_0PM:
                x3 = 0;
                x2 = maxX2;
                break;
            case D3Q27System::DIR_0MM:
                x3 = 0;
                x2 = 0;
                break;
            case D3Q27System::DIR_0MP:
                x3 = maxX3;
                x2 = 0;
                break;
        }
        for (int x1 = 1; x1 < maxX1; x1++) {
            distributeData(rdata, index, x1, x2, x3);
        }
    }
    // TNE TNW TSW TSE BNE BNW BSW BSE
    else if (sendDir == D3Q27System::DIR_PPP || sendDir == D3Q27System::DIR_MPP || sendDir == D3Q27System::DIR_MMP ||
             sendDir == D3Q27System::DIR_PMP || sendDir == D3Q27System::DIR_PPM || sendDir == D3Q27System::DIR_MPM ||
             sendDir == D3Q27System::DIR_MMM || sendDir == D3Q27System::DIR_PMM) {
        int x1 = 0;
        int x2 = 0;
        int x3 = 0;

        switch (sendDir) {
            case D3Q27System::DIR_PPP:
                x1 = maxX1;
                x2 = maxX2;
                x3 = maxX3;
                break;
            case D3Q27System::DIR_MPP:
                x1 = 0;
                x2 = maxX2;
                x3 = maxX3;
                break;
            case D3Q27System::DIR_MMP:
                x1 = 0;
                x2 = 0;
                x3 = maxX3;
                break;
            case D3Q27System::DIR_PMP:
                x1 = maxX1;
                x2 = 0;
                x3 = maxX3;
                break;
            case D3Q27System::DIR_PPM:
                x1 = maxX1;
                x2 = maxX2;
                x3 = 0;
                break;
            case D3Q27System::DIR_MPM:
                x1 = 0;
                x2 = maxX2;
                x3 = 0;
                break;
            case D3Q27System::DIR_MMM:
                x1 = 0;
                x2 = 0;
                x3 = 0;
                break;
            case D3Q27System::DIR_PMM:
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
