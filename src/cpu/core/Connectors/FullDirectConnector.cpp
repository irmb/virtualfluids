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

#include "FullDirectConnector.h"
#include "EsoSplit.h"
#include "DataSet3D.h"
#include "LBMKernel.h"

using namespace std;

FullDirectConnector::FullDirectConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir)
    : LocalBlock3DConnector(from, to, sendDir)

{
}
//////////////////////////////////////////////////////////////////////////
void FullDirectConnector::init()
{ 
    maxX1 = (int)this->from.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1() - 1;
    maxX2 = (int)this->from.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2() - 1;
    maxX3 = (int)this->from.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3() - 1;
}
//////////////////////////////////////////////////////////////////////////
void FullDirectConnector::sendVectors()
{
    updatePointers();
    exchangeData();
}
//////////////////////////////////////////////////////////////////////////
void FullDirectConnector::exchangeData()
{
    using namespace vf::lbm::dir;

    // EAST
    if (sendDir == dP00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                exchangeData(maxX1 - 1, x2, x3, 0, x2, x3);
            }
        }
    }
    // WEST
    else if (sendDir == dM00) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                exchangeData(1, x2, x3, maxX1, x2, x3);
            }
        }
    }
    // NORTH
    else if (sendDir == d0P0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                exchangeData(x1, maxX2 - 1, x3, x1, 0, x3);
            }
        }
    }
    // SOUTH
    else if (sendDir == d0M0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                exchangeData(x1, 1, x3, x1, maxX2, x3);
            }
        }
    }

    // TOP
    else if (sendDir == d00P) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                exchangeData(x1, x2, maxX3 - 1, x1, x2, 0);
            }
        }
    }
    // BOTTOM
    else if (sendDir == d00M) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                exchangeData(x1, x2, 1, x1, x2, maxX3);
            }
        }
    }
    // NORTHEAST
    else if (sendDir == dPP0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            exchangeData(maxX1 - 1, maxX2 - 1, x3, 0, 0, x3);
        }
    }
    // NORTHWEST
    else if (sendDir == dMP0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            exchangeData(1, maxX2 - 1, x3, maxX1, 0, x3);
        }
    }
    // SOUTHWEST
    else if (sendDir == dMM0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            exchangeData(1, 1, x3, maxX1, maxX2, x3);
        }
    }
    // SOUTHEAST
    else if (sendDir == dPM0) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            exchangeData(maxX1 - 1, 1, x3, 0, maxX2, x3);
        }
    } else if (sendDir == dP0P)
        for (int x2 = 1; x2 < maxX2; x2++) {
            exchangeData(maxX1 - 1, x2, maxX3 - 1, 0, x2, 0);
        }
    else if (sendDir == dM0M)
        for (int x2 = 1; x2 < maxX2; x2++) {
            exchangeData(1, x2, 1, maxX1, x2, maxX3);
        }
    else if (sendDir == dP0M)
        for (int x2 = 1; x2 < maxX2; x2++) {
            exchangeData(maxX1 - 1, x2, 1, 0, x2, maxX3);
        }
    else if (sendDir == dM0P)
        for (int x2 = 1; x2 < maxX2; x2++) {
            exchangeData(1, x2, maxX3 - 1, maxX1, x2, 0);
        }
    else if (sendDir == d0PP)
        for (int x1 = 1; x1 < maxX1; x1++) {
            exchangeData(x1, maxX2 - 1, maxX3 - 1, x1, 0, 0);
        }
    else if (sendDir == d0MM)
        for (int x1 = 1; x1 < maxX1; x1++) {
            exchangeData(x1, 1, 1, x1, maxX2, maxX3);
        }
    else if (sendDir == d0PM)
        for (int x1 = 1; x1 < maxX1; x1++) {
            exchangeData(x1, maxX2 - 1, 1, x1, 0, maxX3);
        }

    else if (sendDir == d0MP)
        for (int x1 = 1; x1 < maxX1; x1++) {
            exchangeData(x1, 1, maxX3 - 1, x1, maxX2, 0);
        }

    else if (sendDir == dMMP) {
        exchangeData(1, 1, maxX3 - 1, maxX1, maxX2, 0);
    } else if (sendDir == dPMP) {
        exchangeData(maxX1 - 1, 1, maxX3 - 1, 0, maxX2, 0);
    } else if (sendDir == dMPP) {
        exchangeData(1, maxX2 - 1, maxX3 - 1, maxX1, 0, 0);
    } else if (sendDir == dPPP) {
        exchangeData(maxX1 - 1, maxX2 - 1, maxX3 - 1, 0, 0, 0);
    } else if (sendDir == dMMM) {
        exchangeData(1, 1, 1, maxX1, maxX2, maxX3);
    } else if (sendDir == dPMM) {
        exchangeData(maxX1 - 1, 1, 1, 0, maxX2, maxX3);
    } else if (sendDir == dMPM) {
        exchangeData(1, maxX2 - 1, 1, maxX1, 0, maxX3);
    } else if (sendDir == dPPM) {
        exchangeData(maxX1 - 1, maxX2 - 1, 1, 0, 0, maxX3);
    } else
        UB_THROW(UbException(UB_EXARGS, "unknown dir"));
}

//! \}
