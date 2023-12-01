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
//! \file NoSlipBC.h
//! \ingroup BoundarConditions
//! \author Sören Freudiger
//=======================================================================================

#ifndef NoSlipBC_H
#define NoSlipBC_H

#include "BC.h"

//! A class provides an interface for no-slip boundary condition in grid generator
class NoSlipBC : public BC
{
public:
    NoSlipBC() : BC() {}
    NoSlipBC(const short &secondaryBcOption) : BC(secondaryBcOption) {}

    void init(const D3Q27Interactor *const &interactor, const real &time = 0) override {}
    void update(const D3Q27Interactor *const &interactor, const real &time = 0) override {}

    void adaptBCForDirection(const D3Q27Interactor & /*interactor*/, SPtr<BoundaryConditions> bc,
                             const real & /*worldX1*/, const real & /*worldX2*/, const real & /*worldX3*/,
                             const real &q, const int &fdirection, const real & /*time*/ = 0) override
    {
        bc->setNoSlipBoundaryFlag(D3Q27System::INVDIR[fdirection], secondaryBcOption);
        bc->setQ((real)q, fdirection);
    }
    void adaptBC(const D3Q27Interactor & /*interactor*/, SPtr<BoundaryConditions> bc, const real & /*worldX1*/,
                 const real & /*worldX2*/, const real & /*worldX3*/, const real & /*time*/ = 0) override
    {
        bc->setBCStrategyKey(bcStrategyKey);
    }

private:
};
#endif // NoSlipBC_H
