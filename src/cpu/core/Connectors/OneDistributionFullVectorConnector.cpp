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
//! \file OneDistributionFullVectorConnector.cpp
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#include "OneDistributionFullVectorConnector.h"
#include "DataSet3D.h"
#include "LBMKernel.h"
//////////////////////////////////////////////////////////////////////////
OneDistributionFullVectorConnector::OneDistributionFullVectorConnector(SPtr<Block3D> block, VectorTransmitterPtr sender,
                                                       VectorTransmitterPtr receiver, int sendDir)
    : FullVectorConnector(block, sender, receiver, sendDir)
{
    if (!block || !sender || !receiver)
        UB_THROW(UbException(UB_EXARGS, "sender or receiver == NULL!!"));
}
//////////////////////////////////////////////////////////////////////////
void OneDistributionFullVectorConnector::init()
{
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;

    FullVectorConnector::init();
    
    fDis = dynamicPointerCast<EsoTwist3D>(block.lock()->getKernel()->getDataSet()->getFdistributions());

    int anz = 27;
    switch (sendDir) {
        case d000:
            UB_THROW(UbException(UB_EXARGS, "ZERO not allowed"));
            break;
        case dP00:
        case dM00:
            sender->getData().resize(maxX2 * maxX3 * anz, c0o1);
            break;
        case d0P0:
        case d0M0:
            sender->getData().resize(maxX1 * maxX3 * anz, c0o1);
            break;
        case d00P:
        case d00M:
            sender->getData().resize(maxX1 * maxX2 * anz, c0o1);
            break;

        case dPP0:
        case dMM0:
        case dPM0:
        case dMP0:
            sender->getData().resize(maxX3 * anz, c0o1);
            break;

        case dP0P:
        case dM0M:
        case dP0M:
        case dM0P:
            sender->getData().resize(maxX2 * anz, c0o1);
            break;

        case d0PP:
        case d0MM:
        case d0PM:
        case d0MP:
            sender->getData().resize(maxX1 * anz, c0o1);
            break;

        case dPPP:
        case dMMM:
        case dPPM:
        case dMMP:
        case dPMP:
        case dMPM:
        case dPMM:
        case dMPP:
            sender->getData().resize(anz, c0o1);
            break;

        default:
            UB_THROW(UbException(UB_EXARGS, "unknown sendDir"));
    }
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
