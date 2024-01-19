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
//! \author Konstantin Kutscher
//=======================================================================================

#include "BCSet.h"
#include "BCStrategy.h"
#include "BCArray3D.h"
#include "EsoSplit.h"
#include "DataSet3D.h"
#include "ILBMKernel.h"

BCSet::BCSet() = default;
//////////////////////////////////////////////////////////////////////////
BCSet::BCSet(SPtr<ILBMKernel> kernel)
{
    SPtr<DistributionArray3D> distributions =
        std::dynamic_pointer_cast<EsoTwist3D>(kernel->getDataSet()->getFdistributions());
    bcArray = std::make_shared<BCArray3D>(distributions->getNX1(), distributions->getNX2(), distributions->getNX3(),
                                          BCArray3D::FLUID);
}
//////////////////////////////////////////////////////////////////////////
BCSet::~BCSet() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCSet> BCSet::clone(SPtr<ILBMKernel> kernel)
{
    SPtr<BCSet> bcSet(new BCSet(kernel));
    return bcSet;
}
//////////////////////////////////////////////////////////////////////////
SPtr<BCArray3D> BCSet::getBCArray() { return bcArray; }
//////////////////////////////////////////////////////////////////////////
void BCSet::setBCArray(SPtr<BCArray3D> bcarray) { bcArray = bcarray; }
//////////////////////////////////////////////////////////////////////////
void BCSet::addBC(SPtr<BCStrategy> bc)
{
    if (bc->isPreCollision()) {
        preBC.push_back(bc);
    } else {
        postBC.push_back(bc);
    }
}
//////////////////////////////////////////////////////////////////////////
void BCSet::applyPreCollisionBC()
{
    for (SPtr<BCStrategy> bc : preBC)
        bc->applyBC();
}
//////////////////////////////////////////////////////////////////////////
void BCSet::applyPostCollisionBC()
{
    for (SPtr<BCStrategy> bc : postBC)
        bc->applyBC();
}
//////////////////////////////////////////////////////////////////////////
void BCSet::clearBC()
{
    preBC.clear();
    postBC.clear();
}

//! \}
