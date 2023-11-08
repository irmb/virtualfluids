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
//! \file BC.h
//! \ingroup BoundarConditions
//! \author Sören Freudiger
//=======================================================================================
#ifndef BC_H
#define BC_H

#include <PointerDefinitions.h>

#include "BCStrategy.h"
#include "BoundaryConditions.h"

class D3Q27Interactor;

//! \brief Abstract class of baundary conditions adapter
//! \details  BC supports the definition of boundary conditions in grid generation
class BC
{
public:
    BC() = default;

    //! \param secondaryBcOption additional option of boundary conditions
    BC(const short &secondaryBcOption) : secondaryBcOption(secondaryBcOption) {}
    virtual ~BC() = default;

    // methods
    bool isTimeDependent() { return ((this->type & TIMEDEPENDENT) == TIMEDEPENDENT); }

    virtual short getSecondaryBcOption() { return this->secondaryBcOption; }
    virtual void setSecondaryBcOption(const short &val) { this->secondaryBcOption = val; }

    virtual void init(const D3Q27Interactor *const &interactor, const real &time = 0)   = 0;
    virtual void update(const D3Q27Interactor *const &interactor, const real &time = 0) = 0;

    virtual void adaptBC(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const real &worldX1,
                         const real &worldX2, const real &worldX3, const real &time = 0)       = 0;
    virtual void adaptBCForDirection(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc,
                                     const real &worldX1, const real &worldX2, const real &worldX3,
                                     const real &q, const int &fdirection, const real &time = 0) = 0;

    void setBCStrategy(SPtr<BCStrategy> alg)
    {
        algorithmType = alg->getType();
        algorithm     = alg;
    }
    SPtr<BCStrategy> getAlgorithm() { return algorithm; }
    char getBCStrategyType() { return algorithmType; }

protected:
    short secondaryBcOption{ 0 };

    char type{ 0 };

    SPtr<BCStrategy> algorithm;
    char algorithmType{ -1 };

    static const char TIMEDEPENDENT = 1 << 0; //'1';
    static const char TIMEPERIODIC  = 1 << 1; //'2';
};

#endif // D3Q27BOUNDARYCONDITIONADAPTER_H
