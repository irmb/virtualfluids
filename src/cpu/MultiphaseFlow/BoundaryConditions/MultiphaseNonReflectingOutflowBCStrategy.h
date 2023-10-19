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
//! \file MultiphaseNonReflectingOutflowBCStrategy.h
//! \ingroup BoundarConditions
//! \author Hesameddin Safari
//=======================================================================================

#ifndef MultiphaseNonReflectingOutflowBCStrategy_h__
#define MultiphaseNonReflectingOutflowBCStrategy_h__

#include "MultiphaseBCStrategy.h"
//! A class implements non reflecting outflow boundary condition for multiphase simulations
class MultiphaseNonReflectingOutflowBCStrategy : public MultiphaseBCStrategy
{
public:
    MultiphaseNonReflectingOutflowBCStrategy();
    ~MultiphaseNonReflectingOutflowBCStrategy();
    SPtr<BCStrategy> clone();
    void addDistributions(SPtr<DistributionArray3D> distributions);
    void addDistributionsH(SPtr<DistributionArray3D> distributionsH);
    void addDistributionsH2(SPtr<DistributionArray3D> distributionsH2);
    void applyBC();
};
#endif // MultiphaseNonReflectingOutflowBCStrategy_h__
