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

#include "TwoDistributionsDoubleGhostLayerFullVectorConnector.h"
#include "Block3D.h"
#include "LBMKernel.h"
#include "EsoTwist3D.h"
#include "DataSet3D.h"

//////////////////////////////////////////////////////////////////////////
TwoDistributionsDoubleGhostLayerFullVectorConnector::TwoDistributionsDoubleGhostLayerFullVectorConnector(SPtr<Block3D> block,
                                                                         VectorTransmitterPtr sender,
                                                                         VectorTransmitterPtr receiver, int sendDir)
    : FullVectorConnector(block, sender, receiver, sendDir)
{
   if (!block || !sender || !receiver)
      UB_THROW(UbException(UB_EXARGS, "sender or receiver == NULL!!"));

}
//////////////////////////////////////////////////////////////////////////
void TwoDistributionsDoubleGhostLayerFullVectorConnector::init()
{
   using namespace vf::lbm::dir;
   using namespace vf::basics::constant;

   FullVectorConnector::init();

   fDis = dynamicPointerCast<EsoTwist3D>(block.lock()->getKernel()->getDataSet()->getFdistributions());
   hDis = dynamicPointerCast<EsoTwist3D>(block.lock()->getKernel()->getDataSet()->getHdistributions());
   pressure   = block.lock()->getKernel()->getDataSet()->getPressureField();

   int anz = 2*27+1;
   switch (sendDir)
   {
   case d000: UB_THROW(UbException(UB_EXARGS, "ZERO not allowed")); break;
   case dP00:
   case dM00: sender->getData().resize(maxX2*maxX3*anz*2, c0o1);   break;
   case d0P0:
   case d0M0: sender->getData().resize(maxX1*maxX3*anz*2, c0o1);   break;
   case d00P:
   case d00M: sender->getData().resize(maxX1*maxX2*anz*2, c0o1);   break;

   case dPP0:
   case dMM0:
   case dPM0:
   case dMP0:  sender->getData().resize(maxX3*anz*4, c0o1);   break;

   case dP0P:
   case dM0M:
   case dP0M:
   case dM0P:  sender->getData().resize(maxX2*anz*4, c0o1);   break;

   case d0PP:
   case d0MM:
   case d0PM:
   case d0MP:  sender->getData().resize(maxX1*anz*4, c0o1);   break;

   case dPPP:
   case dMMM:
   case dPPM:
   case dMMP:
   case dPMP:
   case dMPM:
   case dPMM:
   case dMPP:  sender->getData().resize(anz*8, c0o1);   break;

   default: UB_THROW(UbException(UB_EXARGS, "unknown sendDir"));
   }
}
//////////////////////////////////////////////////////////////////////////
void TwoDistributionsDoubleGhostLayerFullVectorConnector::fillSendVectors() 
{ 
    updatePointers();
    fillData();
}
////////////////////////////////////////////////////////////////////////
void TwoDistributionsDoubleGhostLayerFullVectorConnector::fillData()
{
    using namespace vf::lbm::dir;

    ////////////////////////////////////////////////////////////
    // relation between ghost layer and regular nodes
    // maxX1m3 maxX1m2 ... minX1p2 minX1p3 - regular nodes
    // minX1   minX1p1 ... maxX1m1 maxX1   - ghost layer
    ////////////////////////////////////////////////////////////

    int minX1   = 0;
    //int minX1p1 = minX1 + 1;
    int minX1p2 = minX1 + 2;
    int minX1p3 = minX1 + 3;
    //int maxX1m1 = maxX1 - 1;
    int maxX1m2 = maxX1 - 2;
    int maxX1m3 = maxX1 - 3;

    int minX2   = 0;
    //int minX2p1 = minX2 + 1;
    int minX2p2 = minX2 + 2;
    int minX2p3 = minX2 + 3;
    //int maxX2m1 = maxX2 - 1;
    int maxX2m2 = maxX2 - 2;
    int maxX2m3 = maxX2 - 3;

    int minX3   = 0;
    //int minX3p1 = minX3 + 1;
    int minX3p2 = minX3 + 2;
    int minX3p3 = minX3 + 3;
    //int maxX3m1 = maxX3 - 1;
    int maxX3m2 = maxX3 - 2;
    int maxX3m3 = maxX3 - 3;

    vector_type &sdata = sender->getData();

    int index = 0;
    // EAST
    if (sendDir == dP00) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
                fillData(sdata, index, maxX1m3, x2, x3);
                fillData(sdata, index, maxX1m2, x2, x3);
            }
        }
    }
    // WEST
    else if (sendDir == dM00) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
                fillData(sdata, index, minX1p3, x2, x3);
                fillData(sdata, index, minX1p2, x2, x3);
            }
        }
    }
    // NORTH
    else if (sendDir == d0P0) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
                fillData(sdata, index, x1, maxX2m3, x3);
                fillData(sdata, index, x1, maxX2m2, x3);
            }
        }
    }
    // SOUTH
    else if (sendDir == d0M0) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
                fillData(sdata, index, x1, minX2p3, x3);
                fillData(sdata, index, x1, minX2p2, x3);
            }
        }
    }

    // TOP
    else if (sendDir == d00P) {
        for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
            for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
                fillData(sdata, index, x1, x2, maxX3m3);
                fillData(sdata, index, x1, x2, maxX3m2);
            }
        }
    }
    // BOTTOM
    else if (sendDir == d00M) {
        for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
            for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
                fillData(sdata, index, x1, x2, minX3p3);
                fillData(sdata, index, x1, x2, minX3p2);
            }
        }
    }
    // NORTHEAST
    else if (sendDir == dPP0) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            fillData(sdata, index, maxX1m3, maxX2m3, x3);
            fillData(sdata, index, maxX1m2, maxX2m2, x3);
            fillData(sdata, index, maxX1m3, maxX2m2, x3);
            fillData(sdata, index, maxX1m2, maxX2m3, x3);
        }
    }
    // NORTHWEST
    else if (sendDir == dMP0) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            fillData(sdata, index, minX1p3, maxX2m3, x3);
            fillData(sdata, index, minX1p2, maxX2m2, x3);
            fillData(sdata, index, minX1p3, maxX2m2, x3);
            fillData(sdata, index, minX1p2, maxX2m3, x3);
        }
    }
    // SOUTHWEST
    else if (sendDir == dMM0) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            fillData(sdata, index, minX1p3, minX2p3, x3);
            fillData(sdata, index, minX1p2, minX2p2, x3);
            fillData(sdata, index, minX1p3, minX2p2, x3);
            fillData(sdata, index, minX1p2, minX2p3, x3);
        }
    }
    // SOUTHEAST
    else if (sendDir == dPM0) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            fillData(sdata, index, maxX1m3, minX2p3, x3);
            fillData(sdata, index, maxX1m2, minX2p2, x3);
            fillData(sdata, index, maxX1m3, minX2p2, x3);
            fillData(sdata, index, maxX1m2, minX2p3, x3);
        }
    } else if (sendDir == dP0P)
        for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
            fillData(sdata, index, maxX1m3, x2, maxX3m3);
            fillData(sdata, index, maxX1m2, x2, maxX3m2);
            fillData(sdata, index, maxX1m3, x2, maxX3m2);
            fillData(sdata, index, maxX1m2, x2, maxX3m3);
        }
    else if (sendDir == dM0M)
        for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
            fillData(sdata, index, minX1p3, x2, minX3p3);
            fillData(sdata, index, minX1p2, x2, minX3p2);
            fillData(sdata, index, minX1p3, x2, minX3p2);
            fillData(sdata, index, minX1p2, x2, minX3p3);
        }
    else if (sendDir == dP0M)
        for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
            fillData(sdata, index, maxX1m3, x2, minX3p3);
            fillData(sdata, index, maxX1m2, x2, minX3p2);
            fillData(sdata, index, maxX1m3, x2, minX3p2);
            fillData(sdata, index, maxX1m2, x2, minX3p3);
        }
    else if (sendDir == dM0P)
        for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
            fillData(sdata, index, minX1p3, x2, maxX3m3);
            fillData(sdata, index, minX1p2, x2, maxX3m2);
            fillData(sdata, index, minX1p3, x2, maxX3m2);
            fillData(sdata, index, minX1p2, x2, maxX3m3);
        }
    else if (sendDir == d0PP)
        for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
            fillData(sdata, index, x1, maxX2m3, maxX3m3);
            fillData(sdata, index, x1, maxX2m2, maxX3m2);
            fillData(sdata, index, x1, maxX2m3, maxX3m2);
            fillData(sdata, index, x1, maxX2m2, maxX3m3);
        }
    else if (sendDir == d0MM)
        for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
            fillData(sdata, index, x1, minX2p3, minX3p3);
            fillData(sdata, index, x1, minX2p2, minX3p2);
            fillData(sdata, index, x1, minX2p3, minX3p2);
            fillData(sdata, index, x1, minX2p2, minX3p3);
        }
    else if (sendDir == d0PM)
        for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
            fillData(sdata, index, x1, maxX2m3, minX3p3);
            fillData(sdata, index, x1, maxX2m2, minX3p2);
            fillData(sdata, index, x1, maxX2m3, minX3p2);
            fillData(sdata, index, x1, maxX2m2, minX3p3);
        }
    else if (sendDir == d0MP)
        for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
            fillData(sdata, index, x1, minX2p3, maxX3m3);
            fillData(sdata, index, x1, minX2p2, maxX3m2);
            fillData(sdata, index, x1, minX2p3, maxX3m2);
            fillData(sdata, index, x1, minX2p2, maxX3m3);
        }
    else if (sendDir == dMMP) {
        fillData(sdata, index, minX1p3, minX2p3, maxX3m3);
        fillData(sdata, index, minX1p2, minX2p2, maxX3m2);
        fillData(sdata, index, minX1p3, minX2p2, maxX3m2);
        fillData(sdata, index, minX1p2, minX2p3, maxX3m2);
        fillData(sdata, index, minX1p2, minX2p2, maxX3m3);
        fillData(sdata, index, minX1p3, minX2p3, maxX3m2);
        fillData(sdata, index, minX1p3, minX2p2, maxX3m3);
        fillData(sdata, index, minX1p2, minX2p3, maxX3m3);
    } else if (sendDir == dPMP) {
        fillData(sdata, index, maxX1m3, minX1p3, maxX3m3);
        fillData(sdata, index, maxX1m2, minX1p2, maxX3m2);
        fillData(sdata, index, maxX1m3, minX1p2, maxX3m2);
        fillData(sdata, index, maxX1m2, minX1p3, maxX3m2);
        fillData(sdata, index, maxX1m2, minX1p2, maxX3m3);
        fillData(sdata, index, maxX1m3, minX1p3, maxX3m2);
        fillData(sdata, index, maxX1m3, minX1p2, maxX3m3);
        fillData(sdata, index, maxX1m2, minX1p3, maxX3m3);
    } else if (sendDir == dMPP) {
        fillData(sdata, index, minX1p3, maxX2m3, maxX3m3);
        fillData(sdata, index, minX1p2, maxX2m2, maxX3m2);
        fillData(sdata, index, minX1p3, maxX2m2, maxX3m2);
        fillData(sdata, index, minX1p2, maxX2m3, maxX3m2);
        fillData(sdata, index, minX1p2, maxX2m2, maxX3m3);
        fillData(sdata, index, minX1p3, maxX2m3, maxX3m2);
        fillData(sdata, index, minX1p3, maxX2m2, maxX3m3);
        fillData(sdata, index, minX1p2, maxX2m3, maxX3m3);
    } else if (sendDir == dPPP) {
        fillData(sdata, index, maxX1m3, maxX2m3, maxX3m3);
        fillData(sdata, index, maxX1m2, maxX2m2, maxX3m2);
        fillData(sdata, index, maxX1m3, maxX2m2, maxX3m2);
        fillData(sdata, index, maxX1m2, maxX2m3, maxX3m2);
        fillData(sdata, index, maxX1m2, maxX2m2, maxX3m3);
        fillData(sdata, index, maxX1m3, maxX2m3, maxX3m2);
        fillData(sdata, index, maxX1m3, maxX2m2, maxX3m3);
        fillData(sdata, index, maxX1m2, maxX2m3, maxX3m3);
    } else if (sendDir == dMMM) {
        fillData(sdata, index, minX1p3, minX2p3, minX3p3);
        fillData(sdata, index, minX1p2, minX2p2, minX3p2);
        fillData(sdata, index, minX1p3, minX2p2, minX3p2);
        fillData(sdata, index, minX1p2, minX2p3, minX3p2);
        fillData(sdata, index, minX1p2, minX2p2, minX3p3);
        fillData(sdata, index, minX1p3, minX2p3, minX3p2);
        fillData(sdata, index, minX1p3, minX2p2, minX3p3);
        fillData(sdata, index, minX1p2, minX2p3, minX3p3);
    } else if (sendDir == dPMM) {
        fillData(sdata, index, maxX1m3, minX2p3, minX3p3);
        fillData(sdata, index, maxX1m2, minX2p2, minX3p2);
        fillData(sdata, index, maxX1m3, minX2p2, minX3p2);
        fillData(sdata, index, maxX1m2, minX2p3, minX3p2);
        fillData(sdata, index, maxX1m2, minX2p2, minX3p3);
        fillData(sdata, index, maxX1m3, minX2p3, minX3p2);
        fillData(sdata, index, maxX1m3, minX2p2, minX3p3);
        fillData(sdata, index, maxX1m2, minX2p3, minX3p3);
    } else if (sendDir == dMPM) {
        fillData(sdata, index, minX1p3, maxX2m3, minX3p3);
        fillData(sdata, index, minX1p2, maxX2m2, minX3p2);
        fillData(sdata, index, minX1p3, maxX2m2, minX3p2);
        fillData(sdata, index, minX1p2, maxX2m3, minX3p2);
        fillData(sdata, index, minX1p2, maxX2m2, minX3p3);
        fillData(sdata, index, minX1p3, maxX2m3, minX3p2);
        fillData(sdata, index, minX1p3, maxX2m2, minX3p3);
        fillData(sdata, index, minX1p2, maxX2m3, minX3p3);
    } else if (sendDir == dPPM) {
        fillData(sdata, index, maxX1m3, maxX2m3, minX3p3);
        fillData(sdata, index, maxX1m2, maxX2m2, minX3p2);
        fillData(sdata, index, maxX1m3, maxX2m2, minX3p2);
        fillData(sdata, index, maxX1m2, maxX2m3, minX3p2);
        fillData(sdata, index, maxX1m2, maxX2m2, minX3p3);
        fillData(sdata, index, maxX1m3, maxX2m3, minX3p2);
        fillData(sdata, index, maxX1m3, maxX2m2, minX3p3);
        fillData(sdata, index, maxX1m2, maxX2m3, minX3p3);
    } else
        UB_THROW(UbException(UB_EXARGS, "unknown dir"));
}
////////////////////////////////////////////////////////////////////////
void TwoDistributionsDoubleGhostLayerFullVectorConnector::distributeReceiveVectors() 
{
    updatePointers();
    distributeData();
}
////////////////////////////////////////////////////////////////////////
void TwoDistributionsDoubleGhostLayerFullVectorConnector::distributeData()
{
    using namespace vf::lbm::dir;

    vector_type &rdata = receiver->getData();

    int index = 0;
    ////////////////////////////////////////////////////////////
    // relation between ghost layer and regular nodes
    // maxX1m3 maxX1m2 ... minX1p2 minX1p3 - regular nodes
    // minX1   minX1p1 ... maxX1m1 maxX1   - ghost layer
    ////////////////////////////////////////////////////////////

    int minX1   = 0;
    int minX1p1 = minX1 + 1;
    int minX1p2 = minX1 + 2;
    //int minX1p3 = minX1 + 3;
    int maxX1m1 = maxX1 - 1;
    int maxX1m2 = maxX1 - 2;
    //int maxX1m3 = maxX1 - 3;

    int minX2   = 0;
    int minX2p1 = minX2 + 1;
    int minX2p2 = minX2 + 2;
    //int minX2p3 = minX2 + 3;
    int maxX2m1 = maxX2 - 1;
    int maxX2m2 = maxX2 - 2;
    //int maxX2m3 = maxX2 - 3;

    int minX3   = 0;
    int minX3p1 = minX3 + 1;
    int minX3p2 = minX3 + 2;
    //int minX3p3 = minX3 + 3;
    int maxX3m1 = maxX3 - 1;
    int maxX3m2 = maxX3 - 2;
    //int maxX3m3 = maxX3 - 3;

    if (sendDir == dM00) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
                distributeData(rdata, index, minX1, x2, x3);
                distributeData(rdata, index, minX1p1, x2, x3);
            }
        }
    }
    else if (sendDir == dP00) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
                distributeData(rdata, index, maxX1, x2, x3);
                distributeData(rdata, index, maxX1m1, x2, x3);
            }
        }
    }
    else if (sendDir == d0M0) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
                distributeData(rdata, index, x1, minX2, x3);
                distributeData(rdata, index, x1, minX2p1, x3);
            }
        }
    }
    else if (sendDir == d0P0) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
                distributeData(rdata, index, x1, maxX2, x3);
                distributeData(rdata, index, x1, maxX2m1, x3);
            }
        }
    }
    else if (sendDir == d00M) {
        for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
            for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
                distributeData(rdata, index, x1, x2, minX3);
                distributeData(rdata, index, x1, x2, minX3p1);
            }
        }
    }
    else if (sendDir == d00P) {
        for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
            for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
                distributeData(rdata, index, x1, x2, maxX3);
                distributeData(rdata, index, x1, x2, maxX3m1);
            }
        }
    }
    else if (sendDir == dMM0) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            distributeData(rdata, index, minX1, minX2, x3);
            distributeData(rdata, index, minX1p1, minX2p1, x3);
            distributeData(rdata, index, minX1, minX2p1, x3);
            distributeData(rdata, index, minX1p1, minX2, x3);
        }
    }
    else if (sendDir == dPM0) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            distributeData(rdata, index, maxX1, minX2, x3);
            distributeData(rdata, index, maxX1m1, minX2p1, x3);
            distributeData(rdata, index, maxX1, minX2p1, x3);
            distributeData(rdata, index, maxX1m1, minX2, x3);
        }
    }
    else if (sendDir == dPP0) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            distributeData(rdata, index, maxX1, maxX2, x3);
            distributeData(rdata, index, maxX1m1, maxX2m1, x3);
            distributeData(rdata, index, maxX1, maxX2m1, x3);
            distributeData(rdata, index, maxX1m1, maxX2, x3);
        }
    }
    else if (sendDir == dMP0) {
        for (int x3 = minX3p2; x3 <= maxX3m2; x3++) {
            distributeData(rdata, index, minX1, maxX2, x3);
            distributeData(rdata, index, minX1p1, maxX2m1, x3);
            distributeData(rdata, index, minX1, maxX2m1, x3);
            distributeData(rdata, index, minX1p1, maxX2, x3);
        }
    } else if (sendDir == dM0M)
        for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
            distributeData(rdata, index, minX1, x2, minX3);
            distributeData(rdata, index, minX1p1, x2, minX3p1);
            distributeData(rdata, index, minX1, x2, minX3p1);
            distributeData(rdata, index, minX1p1, x2, minX3);
        }
    else if (sendDir == dP0P)
        for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
            distributeData(rdata, index, maxX1, x2, maxX3);
            distributeData(rdata, index, maxX1m1, x2, maxX3m1);
            distributeData(rdata, index, maxX1, x2, maxX3m1);
            distributeData(rdata, index, maxX1m1, x2, maxX3);
        }
    else if (sendDir == dM0P)
        for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
            distributeData(rdata, index, minX1, x2, maxX3);
            distributeData(rdata, index, minX1p1, x2, maxX3m1);
            distributeData(rdata, index, minX1, x2, maxX3m1);
            distributeData(rdata, index, minX1p1, x2, maxX3);
        }
    else if (sendDir == dP0M)
        for (int x2 = minX2p2; x2 <= maxX2m2; x2++) {
            distributeData(rdata, index, maxX1, x2, minX3);
            distributeData(rdata, index, maxX1m1, x2, minX3p1);
            distributeData(rdata, index, maxX1, x2, minX3p1);
            distributeData(rdata, index, maxX1m1, x2, minX3);
        }
    else if (sendDir == d0MM)
        for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
            distributeData(rdata, index, x1, minX2, minX3);
            distributeData(rdata, index, x1, minX2p1, minX3p1);
            distributeData(rdata, index, x1, minX2, minX3p1);
            distributeData(rdata, index, x1, minX2p1, minX3);
        }
    else if (sendDir == d0PP)
        for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
            distributeData(rdata, index, x1, maxX2, maxX3);
            distributeData(rdata, index, x1, maxX2m1, maxX3m1);
            distributeData(rdata, index, x1, maxX2, maxX3m1);
            distributeData(rdata, index, x1, maxX2m1, maxX3);
        }
    else if (sendDir == d0MP)
        for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
            distributeData(rdata, index, x1, minX2, maxX3);
            distributeData(rdata, index, x1, minX2p1, maxX3m1);
            distributeData(rdata, index, x1, minX2, maxX3m1);
            distributeData(rdata, index, x1, minX2p1, maxX3);
        }
    else if (sendDir == d0PM)
        for (int x1 = minX1p2; x1 <= maxX1m2; x1++) {
            distributeData(rdata, index, x1, maxX2, minX3);
            distributeData(rdata, index, x1, maxX2m1, minX3p1);
            distributeData(rdata, index, x1, maxX2, minX3p1);
            distributeData(rdata, index, x1, maxX2m1, minX3);
        }
    else if (sendDir == dPPM) {
        distributeData(rdata, index, maxX1, maxX2, minX3);
        distributeData(rdata, index, maxX1m1, maxX2m1, minX3p1);
        distributeData(rdata, index, maxX1, maxX2m1, minX3p1);
        distributeData(rdata, index, maxX1m1, maxX2, minX3p1);
        distributeData(rdata, index, maxX1m1, maxX2m1, minX3);
        distributeData(rdata, index, maxX1, maxX2, minX3p1);
        distributeData(rdata, index, maxX1, maxX2m1, minX3);
        distributeData(rdata, index, maxX1m1, maxX2, minX3);
    } else if (sendDir == dMPM) {
        distributeData(rdata, index, minX1, maxX2, minX3);
        distributeData(rdata, index, minX1p1, maxX2m1, minX3p1);
        distributeData(rdata, index, minX1, maxX2m1, minX3p1);
        distributeData(rdata, index, minX1p1, maxX2, minX3p1);
        distributeData(rdata, index, minX1p1, maxX2m1, minX3);
        distributeData(rdata, index, minX1, maxX2, minX3p1);
        distributeData(rdata, index, minX1, maxX2m1, minX3);
        distributeData(rdata, index, minX1p1, maxX2, minX3);
    } else if (sendDir == dPMM) {
        distributeData(rdata, index, maxX1, minX2, minX3);
        distributeData(rdata, index, maxX1m1, minX2p1, minX3p1);
        distributeData(rdata, index, maxX1, minX2p1, minX3p1);
        distributeData(rdata, index, maxX1m1, minX2, minX3p1);
        distributeData(rdata, index, maxX1m1, minX2p1, minX3);
        distributeData(rdata, index, maxX1, minX2, minX3p1);
        distributeData(rdata, index, maxX1, minX2p1, minX3);
        distributeData(rdata, index, maxX1m1, minX2, minX3);
    } else if (sendDir == dMMM) {
        distributeData(rdata, index, minX1, minX2, minX3);
        distributeData(rdata, index, minX1p1, minX2p1, minX3p1);
        distributeData(rdata, index, minX1, minX2p1, minX3p1);
        distributeData(rdata, index, minX1p1, minX2, minX3p1);
        distributeData(rdata, index, minX1p1, minX2p1, minX3);
        distributeData(rdata, index, minX1, minX2, minX3p1);
        distributeData(rdata, index, minX1, minX2p1, minX3);
        distributeData(rdata, index, minX1p1, minX2, minX3);
    } else if (sendDir == dPPP) {
        distributeData(rdata, index, maxX1, maxX2, maxX3);
        distributeData(rdata, index, maxX1m1, maxX2m1, maxX3m1);
        distributeData(rdata, index, maxX1, maxX2m1, maxX3m1);
        distributeData(rdata, index, maxX1m1, maxX2, maxX3m1);
        distributeData(rdata, index, maxX1m1, maxX2m1, maxX3);
        distributeData(rdata, index, maxX1, maxX2, maxX3m1);
        distributeData(rdata, index, maxX1, maxX2m1, maxX3);
        distributeData(rdata, index, maxX1m1, maxX2, maxX3);
    } else if (sendDir == dMPP) {
        distributeData(rdata, index, minX1, maxX2, maxX3);
        distributeData(rdata, index, minX1p1, maxX2m1, maxX3m1);
        distributeData(rdata, index, minX1, maxX2m1, maxX3m1);
        distributeData(rdata, index, minX1p1, maxX2, maxX3m1);
        distributeData(rdata, index, minX1p1, maxX2m1, maxX3);
        distributeData(rdata, index, minX1, maxX2, maxX3m1);
        distributeData(rdata, index, minX1, maxX2m1, maxX3);
        distributeData(rdata, index, minX1p1, maxX2, maxX3);
    } else if (sendDir == dPMP) {
        distributeData(rdata, index, maxX1, minX2, maxX3);
        distributeData(rdata, index, maxX1m1, minX2p1, maxX3m1);
        distributeData(rdata, index, maxX1, minX2p1, maxX3m1);
        distributeData(rdata, index, maxX1m1, minX2, maxX3m1);
        distributeData(rdata, index, maxX1m1, minX2p1, maxX3);
        distributeData(rdata, index, maxX1, minX2, maxX3m1);
        distributeData(rdata, index, maxX1, minX2p1, maxX3);
        distributeData(rdata, index, maxX1m1, minX2, maxX3);
    } else if (sendDir == dMMP) {
        distributeData(rdata, index, minX1, minX2, maxX3);
        distributeData(rdata, index, minX1p1, minX2p1, maxX3m1);
        distributeData(rdata, index, minX1, minX2p1, maxX3m1);
        distributeData(rdata, index, minX1p1, minX2, maxX3m1);
        distributeData(rdata, index, minX1p1, minX2p1, maxX3);
        distributeData(rdata, index, minX1, minX2, maxX3m1);
        distributeData(rdata, index, minX1, minX2p1, maxX3);
        distributeData(rdata, index, minX1p1, minX2, maxX3);
    } else
        UB_THROW(UbException(UB_EXARGS, "unknown dir"));

}
//////////////////////////////////////////////////////////////////////////



//! \}
