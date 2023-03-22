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
//! \file RheologyVelocityBCStrategy.h
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef RheologyVelocityBCStrategy_h__
#define RheologyVelocityBCStrategy_h__

#include "BCStrategy.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class RheologyVelocityBCStrategy : public BCStrategy
{
public:
   RheologyVelocityBCStrategy();
   ~RheologyVelocityBCStrategy();
   virtual SPtr<BCStrategy> clone() override { UB_THROW(UbException("real clone() - belongs in the derived class")); }
   void addDistributions(SPtr<DistributionArray3D> distributions) override;
   void applyBC() override;
protected:
   virtual real getRheologyCollFactor(real omegaInf, real shearRate, real drho) const = 0; // { UB_THROW(UbException("real getRheologyCollFactor() - belongs in the derived class")); }
};

#endif // RheologyVelocityBCStrategy_h__

