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
//! \file RheologyBinghamModelVelocityBCAlgorithm.h
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef BinghamModelVelocityBCAlgorithm_h__
#define BinghamModelVelocityBCAlgorithm_h__

#include "RheologyVelocityBCAlgorithm.h"
#include "Rheology.h"

class RheologyBinghamModelVelocityBCAlgorithm : public RheologyVelocityBCAlgorithm
{
public:
   RheologyBinghamModelVelocityBCAlgorithm()
   {
      BCAlgorithm::type = BCAlgorithm::RheologyBinghamModelVelocityBCAlgorithm;
      BCAlgorithm::preCollision = true;
   }
   ~RheologyBinghamModelVelocityBCAlgorithm() {}
   SPtr<BCAlgorithm> clone() override
   {
      SPtr<BCAlgorithm> bc(new RheologyBinghamModelVelocityBCAlgorithm());
      return bc;
   }
protected:
   real getRheologyCollFactor(real omegaInf, real shearRate, real drho) const override 
   { 
      return Rheology::getBinghamCollFactor(omegaInf, shearRate, drho);
   }
};
#endif // BinghamModelVelocityBCAlgorithm_h__
