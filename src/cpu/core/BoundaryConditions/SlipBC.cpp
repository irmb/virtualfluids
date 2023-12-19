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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_BoundaryConditions BoundaryConditions
//! \ingroup cpu_core core
//! \{
//! \author Soeren Freudiger
//=======================================================================================
#include "SlipBC.h"
#include "D3Q27Interactor.h"
#include "D3Q27System.h"
#include "geometry3d/GbCuboid3D.h"

//*==========================================================*/
void SlipBC::adaptBC(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const real & /*worldX1*/,
                            const real & /*worldX2*/, const real & /*worldX3*/, const real & /*time*/)
{
    using namespace vf::lbm::dir;

    //////////////////////////////////////////////////////////////////////////
    //>>> nur workaround! -> Hendrick nach normalen berechnung aus qs fragen

    GbCuboid3DPtr geo = dynamicPointerCast<GbCuboid3D>(interactor.getGbObject3D());
    if (!geo)
        throw UbException(UB_EXARGS, "derzeit nur fuer Cubes valide");

    if (bc->hasSlipBoundaryFlag(dP00))
        bc->setNormalVector(vf::basics::constant::c1o1, vf::basics::constant::c0o1, vf::basics::constant::c0o1);
    else if (bc->hasSlipBoundaryFlag(dM00))
        bc->setNormalVector(-vf::basics::constant::c1o1, vf::basics::constant::c0o1, vf::basics::constant::c0o1);
    else if (bc->hasSlipBoundaryFlag(d0P0))
        bc->setNormalVector(vf::basics::constant::c0o1, vf::basics::constant::c1o1, vf::basics::constant::c0o1);
    else if (bc->hasSlipBoundaryFlag(d0M0))
        bc->setNormalVector(vf::basics::constant::c0o1, -vf::basics::constant::c1o1, vf::basics::constant::c0o1);
    else if (bc->hasSlipBoundaryFlag(d00P))
        bc->setNormalVector(vf::basics::constant::c0o1, vf::basics::constant::c0o1, vf::basics::constant::c1o1);
    else if (bc->hasSlipBoundaryFlag(d00M))
        bc->setNormalVector(vf::basics::constant::c0o1, vf::basics::constant::c0o1, -vf::basics::constant::c1o1);

    bc->setBCStrategyKey(bcStrategyKey);
}

//! \}
