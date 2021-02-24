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
//! \file DensityBCAdapter.h
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef DensityBCAdapter_H
#define DensityBCAdapter_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "basics/utilities/UbMath.h"
#include "basics/utilities/UbTuple.h"

#include "BCAdapter.h"
#include "BCFunction.h"

//*  DensityBCAdapter                                                            */
//*                                                                         */
//**
//<BR><BR>
//@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
//@version 1.0 - 06.09.06
//*/
//
//*
// usage: ...
//*/

class DensityBCAdapter : public BCAdapter
{
public:
    // constructors
    DensityBCAdapter() { this->init(); }
    DensityBCAdapter(const double &dens, const double &startTime = 0.0, const double &endTime = BCFunction::INFCONST);
    DensityBCAdapter(const BCFunction &densBC);
    DensityBCAdapter(const std::vector<BCFunction> &densBCs);
    DensityBCAdapter(const mu::Parser &function, const double &startTime = 0.0,
                     const double &endTime = BCFunction::INFCONST);

    //------------- implements D3Q27BoundaryConditionAdapter ----- start
    std::string toString();

    void init(const D3Q27Interactor *const &interactor, const double &time = 0) override;
    void update(const D3Q27Interactor *const &interactor, const double &time = 0) override;

    void adaptBCForDirection(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const double &worldX1,
                             const double &worldX2, const double &worldX3, const double &q, const int &fdirection,
                             const double &time = 0) override;
    void adaptBC(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const double &worldX1,
                 const double &worldX2, const double &worldX3, const double &time = 0) override;

    double getDensity(const double &x1, const double &x2, const double &x3, const double &timeStep);

    //------------- implements D3Q27BoundaryConditionAdapter ----- end

protected:
    void init();

    // time dependency wird automatisch ueber D3Q27BCFunction Intervalle ermittelt!
    void setTimeDependent() { (this->type |= TIMEDEPENDENT); }
    void unsetTimeDependent() { (this->type &= ~TIMEDEPENDENT); }

    void clear() { densBCs.clear(); }
    void setNodeDensity(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const double &worldX1,
                        const double &worldX2, const double &worldX3, const double &timestep);

private:
    mu::value_type x1, x2, x3; // brauch man nicht serialisieren!
    mu::value_type timeStep;   // brauch man nicht serialisieren!

    mu::Parser *tmpDensityFunction; // brauch man nicht serialisieren!

    std::vector<BCFunction> densBCs;

private:
};

#endif
