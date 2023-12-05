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
//! \file ChangeRandomQs.hpp
//! \ingroup Utilities
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef ChangeRandomQs_h__
#define ChangeRandomQs_h__

#include "LBMKernel.h"
#include "IntegrateValuesHelper.h"
#include "BoundaryConditions.h"
#include "BCArray3D.h"
#include "BCSet.h"

namespace Utilities
{
   void ChangeRandomQs(SPtr<IntegrateValuesHelper> integrateValues)
   {
      std::vector<IntegrateValuesHelper::CalcNodes> cnodes = integrateValues->getCNodes();
      
      for(IntegrateValuesHelper::CalcNodes cn : cnodes)
      {
         SPtr<ILBMKernel> kernel = cn.block->getKernel();
         SPtr<BCArray3D> bcArray = kernel->getBCSet()->getBCArray();
         for(UbTupleInt3 node : cn.nodes)
         {
            SPtr<BoundaryConditions> bc = bcArray->getBC(val<1>(node), val<2>(node), val<3>(node));
            if (bc)
            {
                for (int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
                {
                  if (bc->hasNoSlipBoundaryFlag(fdir))
                  {
                     const int invDir = D3Q27System::INVDIR[fdir];
                     float q = bc->getQ(invDir);
                     //double r = (double)UbRandom::rand(-50, 50);
                     float r = (float)UbRandom::rand(-10, 10);
                     float q_temp = q + q/r;
                     if (q_temp < 0.0)
                     {
                        q_temp = 0.0001f;
                     }
                     else if (q_temp > 1.0)
                     {
                        q_temp = 0.9999f;
                     }
                     bc->setQ(q_temp, fdir);
                  }
                }
            }
         }
      }
   }

}
#endif // ChangeRandomQs_h__
