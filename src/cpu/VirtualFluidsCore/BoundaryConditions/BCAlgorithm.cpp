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
//! \file BCAlgorithm.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================

#include "BCAlgorithm.h"

#include "BCArray3D.h"
#include "BoundaryConditions.h"
#include "EsoTwist3D.h"
#include "BCArray3D.h"


BCAlgorithm::BCAlgorithm() : compressible(false), thixotropy(false)
{
}
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::setNodeIndex(int x1, int x2, int x3)
{
    this->x1 = x1;
    this->x2 = x2;
    this->x3 = x3;
}
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::setBcPointer(SPtr<BoundaryConditions> bcPtr) { this->bcPtr = bcPtr; }
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::setCompressible(bool c)
{
    compressible = c;

    if (this->compressible) {
        calcFeqsForDirFct  = &D3Q27System::getCompFeqForDirection;
        calcMacrosFct      = &D3Q27System::calcCompMacroscopicValues;
        calcFeqFct         = &D3Q27System::calcCompFeq;
        compressibleFactor = 1.0;
    } else {
        calcFeqsForDirFct  = &D3Q27System::getIncompFeqForDirection;
        calcMacrosFct      = &D3Q27System::calcIncompMacroscopicValues;
        calcFeqFct         = &D3Q27System::calcIncompFeq;
        compressibleFactor = 0.0;
    }
}
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::setCollFactor(LBMReal cf) { collFactor = cf; }

//////////////////////////////////////////////////////////////////////////
char BCAlgorithm::getType() { return type; }
//////////////////////////////////////////////////////////////////////////
bool BCAlgorithm::isPreCollision() { return preCollision; }
//////////////////////////////////////////////////////////////////////////
SPtr<BCArray3D> BCAlgorithm::getBcArray() { return bcArray; }
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::setBcArray(SPtr<BCArray3D> bcarray) { bcArray = bcarray; }
