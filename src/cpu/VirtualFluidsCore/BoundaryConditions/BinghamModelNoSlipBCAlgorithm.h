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
//! \file BinghamModelNoSlipBCAlgorithm.h
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef BinghamModelNoSlipBCAlgorithm_h__
#define BinghamModelNoSlipBCAlgorithm_h__

#include "ThixotropyNoSlipBCAlgorithm.h"
#include "Thixotropy.h"

class BinghamModelNoSlipBCAlgorithm : public ThixotropyNoSlipBCAlgorithm
{
public:
   BinghamModelNoSlipBCAlgorithm()
   {
      BCAlgorithm::type = BCAlgorithm::BinghamModelNoSlipBCAlgorithm;
      BCAlgorithm::preCollision = true;
   }
   ~BinghamModelNoSlipBCAlgorithm() {}
   SPtr<BCAlgorithm> clone() override
   {
      SPtr<BCAlgorithm> bc(new BinghamModelNoSlipBCAlgorithm());
      return bc;
   }
protected:
   LBMReal getThyxotropyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho) const override 
   { 
      return Thixotropy::getBinghamCollFactor(omegaInf, shearRate, drho);
   }
};
#endif // BinghamModelNoSlipBCAlgorithm_h__
