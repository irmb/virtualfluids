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
//! \file BCStrategy.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================

#include "BCStrategy.h"

#include "BCArray3D.h"
#include "BoundaryConditions.h"
#include "EsoTwist3D.h"
#include "Block3D.h"

void BCStrategy::setBlock(SPtr<Block3D> block) 
{ 
    this->block = block; 
}
//////////////////////////////////////////////////////////////////////////
void BCStrategy::setNodeIndex(int x1, int x2, int x3)
{
    this->x1 = x1;
    this->x2 = x2;
    this->x3 = x3;
}
//////////////////////////////////////////////////////////////////////////
void BCStrategy::setBcPointer(SPtr<BoundaryConditions> bcPtr) { this->bcPtr = bcPtr; }
//////////////////////////////////////////////////////////////////////////
void BCStrategy::setCompressible(bool c)
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
void BCStrategy::setCollFactor(real cf) { collFactor = cf; }
//////////////////////////////////////////////////////////////////////////
void BCStrategy::setCollFactorL(real cf) { collFactorL = cf; }
//////////////////////////////////////////////////////////////////////////
void BCStrategy::setCollFactorG(real cf) { collFactorG = cf; }
//////////////////////////////////////////////////////////////////////////
void BCStrategy::setCollFactorPh(real cf) { collFactorPh = cf; }
//////////////////////////////////////////////////////////////////////////
void BCStrategy::setDensityRatio(real dr) { densityRatio = dr; }
//////////////////////////////////////////////////////////////////////////
void BCStrategy::setPhiBound(real phiL, real phiH)
{
    this->phiL = phiL;
    this->phiH = phiH;
}
//////////////////////////////////////////////////////////////////////////
char BCStrategy::getType() { return type; }
//////////////////////////////////////////////////////////////////////////
bool BCStrategy::isPreCollision() { return preCollision; }
//////////////////////////////////////////////////////////////////////////
SPtr<BCArray3D> BCStrategy::getBcArray() { return bcArray; }
//////////////////////////////////////////////////////////////////////////
void BCStrategy::setBcArray(SPtr<BCArray3D> bcarray) { bcArray = bcarray; }
