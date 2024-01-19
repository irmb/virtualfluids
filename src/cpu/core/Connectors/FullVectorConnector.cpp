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

#include "FullVectorConnector.h"
#include "EsoSplit.h"
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
    if (sendDir == dP00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                fillData(sdata, index, maxX1 - 1, x2, x3);
            }
        }
    }
    // WEST
    else if (sendDir == dM00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                fillData(sdata, index, 1, x2, x3);
            }
        }
    }
    // NORTH
    else if (sendDir == d0P0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                fillData(sdata, index, x1, maxX2 - 1, x3);
            }
        }
    }
    // SOUTH
    else if (sendDir == d0M0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                fillData(sdata, index, x1, 1, x3);
            }
        }
    }
    // TOP
    else if (sendDir == d00P) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                fillData(sdata, index, x1, x2, maxX3 - 1);
            }
        }
    }
    // BOTTOM
    else if (sendDir == d00M) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                fillData(sdata, index, x1, x2, 1);
            }
        }
    }
    // NE NW SW SE
    else if (sendDir == dPP0 || sendDir == dMP0 || sendDir == dMM0 ||
             sendDir == dPM0) {
        int x1 = 0;
        int x2 = 0;
        switch (sendDir) {
            case dPP0:
                x1 = maxX1 - 1;
                x2 = maxX2 - 1;
                break;
            case dMP0:
                x1 = 1;
                x2 = maxX2 - 1;
                break;
            case dMM0:
                x1 = 1;
                x2 = 1;
                break;
            case dPM0:
                x1 = maxX1 - 1;
                x2 = 1;
                break;
        }
        for (int x3 = 1; x3 < maxX3; x3++) {
            fillData(sdata, index, x1, x2, x3);
        }
    }
    // TE TW BW BE
    else if (sendDir == dP0P || sendDir == dM0P || sendDir == dM0M ||
             sendDir == dP0M) {
        int x1 = 0;
        int x3 = 0;
        switch (sendDir) {
            case dP0P:
                x1 = maxX1 - 1;
                x3 = maxX3 - 1;
                break;
            case dM0P:
                x1 = 1;
                x3 = maxX3 - 1;
                break;
            case dM0M:
                x1 = 1;
                x3 = 1;
                break;
            case dP0M:
                x1 = maxX1 - 1;
                x3 = 1;
                break;
        }
        for (int x2 = 1; x2 < maxX2; x2++) {
            fillData(sdata, index, x1, x2, x3);
        }
    }
    // TN BN BS TS
    else if (sendDir == d0PP || sendDir == d0PM || sendDir == d0MM ||
             sendDir == d0MP) {
        int x2 = 0;
        int x3 = 0;
        switch (sendDir) {
            case d0PP:
                x3 = maxX3 - 1;
                x2 = maxX2 - 1;
                break;
            case d0PM:
                x3 = 1;
                x2 = maxX2 - 1;
                break;
            case d0MM:
                x3 = 1;
                x2 = 1;
                break;
            case d0MP:
                x3 = maxX3 - 1;
                x2 = 1;
                break;
        }
        for (int x1 = 1; x1 < maxX1; x1++) {
            fillData(sdata, index, x1, x2, x3);
        }
    }
    // TNE TNW TSW TSE BNE BNW BSW BSE
    else if (sendDir == dPPP || sendDir == dMPP || sendDir == dMMP ||
             sendDir == dPMP || sendDir == dPPM || sendDir == dMPM ||
             sendDir == dMMM || sendDir == dPMM) {
        int x1 = 0;
        int x2 = 0;
        int x3 = 0;
        switch (sendDir) {
            case dPPP:
                x1 = maxX1 - 1;
                x2 = maxX2 - 1;
                x3 = maxX3 - 1;
                break;
            case dMPP:
                x1 = 1;
                x2 = maxX2 - 1;
                x3 = maxX3 - 1;
                break;
            case dMMP:
                x1 = 1;
                x2 = 1;
                x3 = maxX3 - 1;
                break;
            case dPMP:
                x1 = maxX1 - 1;
                x2 = 1;
                x3 = maxX3 - 1;
                break;
            case dPPM:
                x1 = maxX1 - 1;
                x2 = maxX2 - 1;
                x3 = 1;
                break;
            case dMPM:
                x1 = 1;
                x2 = maxX2 - 1;
                x3 = 1;
                break;
            case dMMM:
                x1 = 1;
                x2 = 1;
                x3 = 1;
                break;
            case dPMM:
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

    if (sendDir == dM00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                distributeData(rdata, index, 0, x2, x3);
            }
        }
    } else if (sendDir == dP00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                distributeData(rdata, index, maxX1, x2, x3);
            }
        }
    } else if (sendDir == d0M0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                distributeData(rdata, index, x1, 0, x3);
            }
        }
    } else if (sendDir == d0P0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                distributeData(rdata, index, x1, maxX2, x3);
            }
        }
    } else if (sendDir == d00M) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                distributeData(rdata, index, x1, x2, 0);
            }
        }
    } else if (sendDir == d00P) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                distributeData(rdata, index, x1, x2, maxX3);
            }
        }
    }
    // NE NW SW SE
    else if (sendDir == dPP0 || sendDir == dMP0 || sendDir == dMM0 ||
             sendDir == dPM0) {
        int x1 = 0;
        int x2 = 0;
        switch (sendDir) // wenn sendir NE dann kommen werte von SW
        {
            case dPP0:
                x1 = maxX1;
                x2 = maxX2;
                break;
            case dMP0:
                x1 = 0;
                x2 = maxX2;
                break;
            case dMM0:
                x1 = 0;
                x2 = 0;
                break;
            case dPM0:
                x1 = maxX1;
                x2 = 0;
                break;
        }
        for (int x3 = 1; x3 < maxX3; x3++) {
            distributeData(rdata, index, x1, x2, x3);
        }

    }
    // TE TW BW BE
    else if (sendDir == dP0P || sendDir == dM0P || sendDir == dM0M ||
             sendDir == dP0M)

    {
        int x1 = 0;
        int x3 = 0;
        switch (sendDir) // wenn sendir NE dann kommen werte von SW
        {
            case dP0P:
                x1 = maxX1;
                x3 = maxX3;
                break;
            case dM0P:
                x1 = 0;
                x3 = maxX3;
                break;
            case dM0M:
                x1 = 0;
                x3 = 0;
                break;
            case dP0M:
                x1 = maxX1;
                x3 = 0;
                break;
        }
        for (int x2 = 1; x2 < maxX2; x2++) {
            distributeData(rdata, index, x1, x2, x3);
        }
    }
    // TN BN BS TS
    else if (sendDir == d0PP || sendDir == d0PM || sendDir == d0MM ||
             sendDir == d0MP) {
        int x2 = 0;
        int x3 = 0;
        switch (sendDir) {
            case d0PP:
                x3 = maxX3;
                x2 = maxX2;
                break;
            case d0PM:
                x3 = 0;
                x2 = maxX2;
                break;
            case d0MM:
                x3 = 0;
                x2 = 0;
                break;
            case d0MP:
                x3 = maxX3;
                x2 = 0;
                break;
        }
        for (int x1 = 1; x1 < maxX1; x1++) {
            distributeData(rdata, index, x1, x2, x3);
        }
    }
    // TNE TNW TSW TSE BNE BNW BSW BSE
    else if (sendDir == dPPP || sendDir == dMPP || sendDir == dMMP ||
             sendDir == dPMP || sendDir == dPPM || sendDir == dMPM ||
             sendDir == dMMM || sendDir == dPMM) {
        int x1 = 0;
        int x2 = 0;
        int x3 = 0;

        switch (sendDir) {
            case dPPP:
                x1 = maxX1;
                x2 = maxX2;
                x3 = maxX3;
                break;
            case dMPP:
                x1 = 0;
                x2 = maxX2;
                x3 = maxX3;
                break;
            case dMMP:
                x1 = 0;
                x2 = 0;
                x3 = maxX3;
                break;
            case dPMP:
                x1 = maxX1;
                x2 = 0;
                x3 = maxX3;
                break;
            case dPPM:
                x1 = maxX1;
                x2 = maxX2;
                x3 = 0;
                break;
            case dMPM:
                x1 = 0;
                x2 = maxX2;
                x3 = 0;
                break;
            case dMMM:
                x1 = 0;
                x2 = 0;
                x3 = 0;
                break;
            case dPMM:
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

//! \}
