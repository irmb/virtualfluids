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
//! \file ThreeDistributionsFullVectorConnector.cpp
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#include "ThreeDistributionsFullVectorConnector.h"
#include "Block3D.h"
#include "LBMKernel.h"
#include "EsoTwist3D.h"
#include "DataSet3D.h"

//////////////////////////////////////////////////////////////////////////
ThreeDistributionsFullVectorConnector::ThreeDistributionsFullVectorConnector(SPtr<Block3D> block,
                                                                         VectorTransmitterPtr sender,
                                                                         VectorTransmitterPtr receiver, int sendDir)
    : FullVectorConnector(block, sender, receiver, sendDir)
{
   if (!block || !sender || !receiver)
      UB_THROW(UbException(UB_EXARGS, "sender or receiver == NULL!!"));

}
//////////////////////////////////////////////////////////////////////////
void ThreeDistributionsFullVectorConnector::init()
{
   FullVectorConnector::init();

   fDis = dynamicPointerCast<EsoTwist3D>(block.lock()->getKernel()->getDataSet()->getFdistributions());
   hDis = dynamicPointerCast<EsoTwist3D>(block.lock()->getKernel()->getDataSet()->getHdistributions());
   h2Dis = dynamicPointerCast<EsoTwist3D>(block.lock()->getKernel()->getDataSet()->getH2distributions());

   int anz = 3*27;
   switch (sendDir)
   {
   case D3Q27System::REST: UB_THROW(UbException(UB_EXARGS, "ZERO not allowed")); break;
   case D3Q27System::E:
   case D3Q27System::W: sender->getData().resize(maxX2*maxX3*anz, 0.0);   break;
   case D3Q27System::N:
   case D3Q27System::S: sender->getData().resize(maxX1*maxX3*anz, 0.0);   break;
   case D3Q27System::T:
   case D3Q27System::B: sender->getData().resize(maxX1*maxX2*anz, 0.0);   break;

   case D3Q27System::NE:
   case D3Q27System::SW:
   case D3Q27System::SE:
   case D3Q27System::NW:  sender->getData().resize(maxX3*anz, 0.0);   break;

   case D3Q27System::TE:
   case D3Q27System::BW:
   case D3Q27System::BE:
   case D3Q27System::TW:  sender->getData().resize(maxX2*anz, 0.0);   break;

   case D3Q27System::TN:
   case D3Q27System::BS:
   case D3Q27System::BN:
   case D3Q27System::TS:  sender->getData().resize(maxX1*anz, 0.0);   break;

   case D3Q27System::TNE:
   case D3Q27System::BSW:
   case D3Q27System::BNE:
   case D3Q27System::TSW:
   case D3Q27System::TSE:
   case D3Q27System::BNW:
   case D3Q27System::BSE:
   case D3Q27System::TNW:  sender->getData().resize(anz, 0.0);   break;

   default: UB_THROW(UbException(UB_EXARGS, "unknown sendDir"));
   }
}



