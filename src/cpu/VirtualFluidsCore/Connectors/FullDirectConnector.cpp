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
//! \file FullDirectConnector.cpp
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#include "FullDirectConnector.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
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
    // EAST
    if (sendDir == D3Q27System::E) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                exchangeData(maxX1 - 1, x2, x3, 0, x2, x3);
            }
        }
    }
    // WEST
    else if (sendDir == D3Q27System::W) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x2 = 1; x2 < maxX2; x2++) {
                exchangeData(1, x2, x3, maxX1, x2, x3);
            }
        }
    }
    // NORTH
    else if (sendDir == D3Q27System::N) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                exchangeData(x1, maxX2 - 1, x3, x1, 0, x3);
            }
        }
    }
    // SOUTH
    else if (sendDir == D3Q27System::S) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                exchangeData(x1, 1, x3, x1, maxX2, x3);
            }
        }
    }

    // TOP
    else if (sendDir == D3Q27System::T) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                exchangeData(x1, x2, maxX3 - 1, x1, x2, 0);
            }
        }
    }
    // BOTTOM
    else if (sendDir == D3Q27System::B) {
        for (int x2 = 1; x2 < maxX2; x2++) {
            for (int x1 = 1; x1 < maxX1; x1++) {
                exchangeData(x1, x2, 1, x1, x2, maxX3);
            }
        }
    }
    // NORTHEAST
    else if (sendDir == D3Q27System::NE) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            exchangeData(maxX1 - 1, maxX2 - 1, x3, 0, 0, x3);
        }
    }
    // NORTHWEST
    else if (sendDir == D3Q27System::NW) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            exchangeData(1, maxX2 - 1, x3, maxX1, 0, x3);
        }
    }
    // SOUTHWEST
    else if (sendDir == D3Q27System::SW) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            exchangeData(1, 1, x3, maxX1, maxX2, x3);
        }
    }
    // SOUTHEAST
    else if (sendDir == D3Q27System::SE) {
        for (int x3 = 1; x3 < maxX3; x3++) {
            exchangeData(maxX1 - 1, 1, x3, 0, maxX2, x3);
        }
    } else if (sendDir == D3Q27System::TE)
        for (int x2 = 1; x2 < maxX2; x2++) {
            exchangeData(maxX1 - 1, x2, maxX3 - 1, 0, x2, 0);
        }
    else if (sendDir == D3Q27System::BW)
        for (int x2 = 1; x2 < maxX2; x2++) {
            exchangeData(1, x2, 1, maxX1, x2, maxX3);
        }
    else if (sendDir == D3Q27System::BE)
        for (int x2 = 1; x2 < maxX2; x2++) {
            exchangeData(maxX1 - 1, x2, 1, 0, x2, maxX3);
        }
    else if (sendDir == D3Q27System::TW)
        for (int x2 = 1; x2 < maxX2; x2++) {
            exchangeData(1, x2, maxX3 - 1, maxX1, x2, 0);
        }
    else if (sendDir == D3Q27System::TN)
        for (int x1 = 1; x1 < maxX1; x1++) {
            exchangeData(x1, maxX2 - 1, maxX3 - 1, x1, 0, 0);
        }
    else if (sendDir == D3Q27System::BS)
        for (int x1 = 1; x1 < maxX1; x1++) {
            exchangeData(x1, 1, 1, x1, maxX2, maxX3);
        }
    else if (sendDir == D3Q27System::BN)
        for (int x1 = 1; x1 < maxX1; x1++) {
            exchangeData(x1, maxX2 - 1, 1, x1, 0, maxX3);
        }

    else if (sendDir == D3Q27System::TS)
        for (int x1 = 1; x1 < maxX1; x1++) {
            exchangeData(x1, 1, maxX3 - 1, x1, maxX2, 0);
        }

    else if (sendDir == D3Q27System::TSW) {
        exchangeData(1, 1, maxX3 - 1, maxX1, maxX2, 0);
    } else if (sendDir == D3Q27System::TSE) {
        exchangeData(maxX1 - 1, 1, maxX3 - 1, 0, maxX2, 0);
    } else if (sendDir == D3Q27System::TNW) {
        exchangeData(1, maxX2 - 1, maxX3 - 1, maxX1, 0, 0);
    } else if (sendDir == D3Q27System::TNE) {
        exchangeData(maxX1 - 1, maxX2 - 1, maxX3 - 1, 0, 0, 0);
    } else if (sendDir == D3Q27System::BSW) {
        exchangeData(1, 1, 1, maxX1, maxX2, maxX3);
    } else if (sendDir == D3Q27System::BSE) {
        exchangeData(maxX1 - 1, 1, 1, 0, maxX2, maxX3);
    } else if (sendDir == D3Q27System::BNW) {
        exchangeData(1, maxX2 - 1, 1, maxX1, 0, maxX3);
    } else if (sendDir == D3Q27System::BNE) {
        exchangeData(maxX1 - 1, maxX2 - 1, 1, 0, 0, maxX3);
    } else
        UB_THROW(UbException(UB_EXARGS, "unknown dir"));
}
