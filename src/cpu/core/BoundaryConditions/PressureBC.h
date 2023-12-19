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

#ifndef PressureBC_H
#define PressureBC_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

//#include "basics/utilities/UbMath.h"
#include "basics/utilities/UbTuple.h"

#include "BC.h"
#include "BCFunction.h"

class PressureBC : public BC
{
public:
    // constructors
    PressureBC() { this->init(); }
    PressureBC(const real &dens, const real &startTime = 0.0, const real &endTime = BCFunction::INFCONST);
    PressureBC(const BCFunction &densBC);
    PressureBC(const std::vector<BCFunction> &densBCs);
    PressureBC(const mu::Parser &function, const real &startTime = 0.0,
                     const real &endTime = BCFunction::INFCONST);

    //------------- implements D3Q27BoundaryConditionAdapter ----- start
    std::string toString();

    void init(const D3Q27Interactor *const &interactor, const real &time = 0) override;
    void update(const D3Q27Interactor *const &interactor, const real &time = 0) override;

    void adaptBCForDirection(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const real &worldX1,
                             const real &worldX2, const real &worldX3, const real &q, const int &fdirection,
                             const real &time = 0) override;
    void adaptBC(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const real &worldX1,
                 const real &worldX2, const real &worldX3, const real &time = 0) override;

    real getDensity(const real &x1, const real &x2, const real &x3, const real &timeStep);

    //------------- implements D3Q27BoundaryConditionAdapter ----- end

protected:
    void init();

    // time dependency wird automatisch ueber D3Q27BCFunction Intervalle ermittelt!
    void setTimeDependent() { (this->type |= TIMEDEPENDENT); }
    void unsetTimeDependent() { (this->type &= ~TIMEDEPENDENT); }

    void clear() { densBCs.clear(); }
    void setNodeDensity(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const real &worldX1,
                        const real &worldX2, const real &worldX3, const real &timestep);

private:
    mu::value_type x1, x2, x3; // brauch man nicht serialisieren!
    mu::value_type timeStep;   // brauch man nicht serialisieren!

    mu::Parser *tmpDensityFunction; // brauch man nicht serialisieren!

    std::vector<BCFunction> densBCs;

private:
};

#endif

//! \}
