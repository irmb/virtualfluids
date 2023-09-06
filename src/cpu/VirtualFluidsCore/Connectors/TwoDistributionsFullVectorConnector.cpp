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
//! \file TwoDistributionsFullVectorConnector.cpp
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#include "TwoDistributionsFullVectorConnector.h"
#include "Block3D.h"
#include "LBMKernel.h"
#include "EsoTwist3D.h"
#include "DataSet3D.h"

//////////////////////////////////////////////////////////////////////////
TwoDistributionsFullVectorConnector::TwoDistributionsFullVectorConnector(SPtr<Block3D> block,
                                                                         VectorTransmitterPtr sender,
                                                                         VectorTransmitterPtr receiver, int sendDir)
    : FullVectorConnector(block, sender, receiver, sendDir)
{
   if (!block || !sender || !receiver)
      UB_THROW(UbException(UB_EXARGS, "sender or receiver == NULL!!"));

}
//////////////////////////////////////////////////////////////////////////
void TwoDistributionsFullVectorConnector::init()
{
   using namespace vf::lbm::dir;
   using namespace vf::basics::constant;
  
   FullVectorConnector::init();

   fDis = dynamicPointerCast<EsoTwist3D>(block.lock()->getKernel()->getDataSet()->getFdistributions());
   hDis = dynamicPointerCast<EsoTwist3D>(block.lock()->getKernel()->getDataSet()->getHdistributions());

   int anz = 2*27;
   switch (sendDir)
   {
   case DIR_000: UB_THROW(UbException(UB_EXARGS, "ZERO not allowed")); break;
   case DIR_P00:
   case DIR_M00: sender->getData().resize(maxX2*maxX3*anz, c0o1);   break;
   case DIR_0P0:
   case DIR_0M0: sender->getData().resize(maxX1*maxX3*anz, c0o1);   break;
   case DIR_00P:
   case DIR_00M: sender->getData().resize(maxX1*maxX2*anz, c0o1);   break;

   case DIR_PP0:
   case DIR_MM0:
   case DIR_PM0:
   case DIR_MP0:  sender->getData().resize(maxX3*anz, c0o1);   break;

   case DIR_P0P:
   case DIR_M0M:
   case DIR_P0M:
   case DIR_M0P:  sender->getData().resize(maxX2*anz, c0o1);   break;

   case DIR_0PP:
   case DIR_0MM:
   case DIR_0PM:
   case DIR_0MP:  sender->getData().resize(maxX1*anz, c0o1);   break;

   case DIR_PPP:
   case DIR_MMM:
   case DIR_PPM:
   case DIR_MMP:
   case DIR_PMP:
   case DIR_MPM:
   case DIR_PMM:
   case DIR_MPP:  sender->getData().resize(anz, c0o1);   break;

   default: UB_THROW(UbException(UB_EXARGS, "unknown sendDir"));
   }
}



