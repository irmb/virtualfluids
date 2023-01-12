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
//! \file SlipBCAdapter.cpp
//! \ingroup BoundarConditions
//! \author Sören Freudiger
//=======================================================================================
#include "SlipBCAdapter.h"
#include "D3Q27Interactor.h"
#include "D3Q27System.h"
#include "geometry3d/GbCuboid3D.h"

//*==========================================================*/
// ObObject* D3Q27SlipBCAdapterCreator::createObObject()
//{
//   return new D3Q27SlipBCAdapter;
//}
//*==========================================================*/
// ObObjectCreator* D3Q27SlipBCAdapter::getCreator()
//{
//   return D3Q27SlipBCAdapterCreator::getInstance();
//}
//*==========================================================*/
void SlipBCAdapter::adaptBC(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const double & /*worldX1*/,
                            const double & /*worldX2*/, const double & /*worldX3*/, const double & /*time*/)
{
    using namespace vf::lbm::dir;

    //////////////////////////////////////////////////////////////////////////
    //>>> nur workaround! -> Hendrick nach normalen berechnung aus qs fragen

    GbCuboid3DPtr geo = dynamicPointerCast<GbCuboid3D>(interactor.getGbObject3D());
    if (!geo)
        throw UbException(UB_EXARGS, "derzeit nur fuer Cubes valide");

    if (bc->hasSlipBoundaryFlag(DIR_P00))
        bc->setNormalVector(1.0, 0.0, 0.0);
    else if (bc->hasSlipBoundaryFlag(DIR_M00))
        bc->setNormalVector(-1.0, 0.0, 0.0);
    else if (bc->hasSlipBoundaryFlag(DIR_0P0))
        bc->setNormalVector(0.0, 1.0, 0.0);
    else if (bc->hasSlipBoundaryFlag(DIR_0M0))
        bc->setNormalVector(0.0, -1.0, 0.0);
    else if (bc->hasSlipBoundaryFlag(DIR_00P))
        bc->setNormalVector(0.0, 0.0, 1.0);
    else if (bc->hasSlipBoundaryFlag(DIR_00M))
        bc->setNormalVector(0.0, 0.0, -1.0);

    bc->setBcAlgorithmType(algorithmType);
}
