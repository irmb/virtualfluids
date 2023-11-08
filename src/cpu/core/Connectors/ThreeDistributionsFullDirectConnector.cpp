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
//! \file ThreeDistributionsFullDirectConnector.cpp
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#include "ThreeDistributionsFullDirectConnector.h"
#include "LBMKernel.h"
#include "DataSet3D.h"

ThreeDistributionsFullDirectConnector::ThreeDistributionsFullDirectConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir)
    : FullDirectConnector(from, to, sendDir)
{

}
//////////////////////////////////////////////////////////////////////////
void ThreeDistributionsFullDirectConnector::init()
{
    FullDirectConnector::init();

	fFrom =dynamicPointerCast<EsoTwist3D>(from.lock()->getKernel()->getDataSet()->getFdistributions());
	fTo = dynamicPointerCast<EsoTwist3D>(to.lock()->getKernel()->getDataSet()->getFdistributions());
	hFrom = dynamicPointerCast<EsoTwist3D>(from.lock()->getKernel()->getDataSet()->getHdistributions());
	hTo = dynamicPointerCast<EsoTwist3D>(to.lock()->getKernel()->getDataSet()->getHdistributions());
    hFrom2 = dynamicPointerCast<EsoTwist3D>(from.lock()->getKernel()->getDataSet()->getH2distributions());
    hTo2  = dynamicPointerCast<EsoTwist3D>(to.lock()->getKernel()->getDataSet()->getH2distributions());
}
