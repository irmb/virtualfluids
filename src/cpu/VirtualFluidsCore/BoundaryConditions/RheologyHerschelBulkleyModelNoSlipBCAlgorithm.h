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
//! \file RheologyHerschelBulkleyModelNoSlipBCAlgorithm.h
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef RheologyHerschelBulkleyModelNoSlipBCAlgorithm_h__
#define RheologyHerschelBulkleyModelNoSlipBCAlgorithm_h__

#include "RheologyNoSlipBCAlgorithm.h"
#include "Rheology.h"

class RheologyHerschelBulkleyModelNoSlipBCAlgorithm : public RheologyNoSlipBCAlgorithm
{
public:
   RheologyHerschelBulkleyModelNoSlipBCAlgorithm() 
   {
      BCAlgorithm::type = BCAlgorithm::RheologyHerschelBulkleyModelNoSlipBCAlgorithm;
      BCAlgorithm::preCollision = true;
   }
   ~RheologyHerschelBulkleyModelNoSlipBCAlgorithm() {}
   SPtr<BCAlgorithm> clone() override
   {
      SPtr<BCAlgorithm> bc(new RheologyHerschelBulkleyModelNoSlipBCAlgorithm());
      return bc;
   }
protected:
   LBMReal getRheologyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho) const override
   {
      return Rheology::getHerschelBulkleyCollFactor(omegaInf, shearRate, drho);
   }
};
#endif // RheologyHerschelBulkleyModelNoSlipBCAlgorithm_h__