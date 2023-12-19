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
//! \author Sören Freudiger
//=======================================================================================
#ifndef VelocityBC_H
#define VelocityBC_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <basics/utilities/UbInfinity.h>

#include <BC.h>
#include <BCFunction.h>

//! \brief A class provides an interface for velocity boundary condition in grid generator.

//! \details
//! Example:
//! \code{.cpp}  vector<BCFunction> vx1BCs,vx2BCs,vx3BCs;
//!        vx1BCs.push_back(BCFunction(0.01 , 0  , 100) );   //t=[0  ..100[ -> vx1 = 0.01
//!        vx1BCs.push_back(BCFunction(0.004, 100, 200) );   //t=[100..200[ -> vx1 = 0.004
//!        vx1BCs.push_back(BCFunction(0.03 , 200, 400) );   //t=[200..400] -> vx1 = 0.03
//!
//!        vx2BCs.push_back(BCFunction(0.02 , 0  , 200) );   //t=[0  ..200[ -> vx2 = 0.02
//!        vx2BCs.push_back(BCFunction(0.002, 200, 300) );   //t=[200..300[ -> vx2 = 0.002
//!        vx2BCs.push_back(BCFunction(0.043, 300, 600) );   //t=[300..600] -> vx2 = 0.043
//!
//!        VelocityBC BC(vx1BCs,vx2BCs,vx3BCs);
//!        BC.setTimePeriodic(); //->  t=[0  ..100[ -> vx1 = 0.01
//!                                     //    t=[100..200[ -> vx1 = 0.004
//!                                     //    t=[200..400[ -> vx1 = 0.03
//!                                     //    t=[400..500[ -> vx1 = 0.01
//!                                     //    t=[500..600[ -> vx1 = 0.004
//!                                     //    t=[600..800[ -> vx1 = 0.03  ...
//!                                     //    t=[0  ..200[ -> vx2 = 0.02
//!                                     //    t=[200..300[ -> vx2 = 0.002
//!                                     //    t=[300..600] -> vx2 = 0.043
//!                                     //    t=[600..800[ -> vx2 = 0.02
//!                                     //    t=[800..900[ -> vx2 = 0.002
//!                                     //    t=[900..1200]-> vx2 = 0.043  ...
//! \endcode
//! Example of parabolic inflow:
//! \code{.cpp}
//!    mu::Parser fct;
//!    fct.SetExpr("max(vmax*(1.0-4.0*((x2-x2_vmax)^2+(x3-x3_vmax)^2)/H^2),0.0)"); //paraboloid (with vmax for
//!    (0/x2_vmax/x3_vmax) fct.DefineConst("x2Vmax", 0.0            ); //x2-Pos für vmax fct.DefineConst("x3Vmax", 0.0
//!    ); //x3-Pos für vmax fct.DefineConst("H"     , diameterOfPipe); fct.DefineConst("vmax"  , vmax           );
//!    VelocityBC velBC(true, false ,false ,fct, 0, BCFunction::INFCONST);
//! \endcode

class VelocityBC : public BC
{
public:
    // constructors
    VelocityBC() { this->init(); }

    VelocityBC(const bool &vx1, const bool &vx2, const bool &vx3, const BCFunction &velVxBC);

    VelocityBC(const bool &vx1, const bool &vx2, const bool &vx3, const mu::Parser &function,
                      const real &startTime, const real &endTime);

    VelocityBC(const bool &vx1, const bool &vx2, const bool &vx3, const mu::Parser &function1,
                      const mu::Parser &function2, const mu::Parser &function3, const real &startTime,
                      const real &endTime);

    VelocityBC(const bool &vx1, const bool &vx2, const bool &vx3, const std::string &functionstring,
                      const real &startTime, const real &endTime);

    VelocityBC(const BCFunction &velBC, bool x1Dir, bool x2Dir, bool x3Dir);

    VelocityBC(const BCFunction &velVx1BC, const BCFunction &velVx2BC, const BCFunction &velVx3BC);

    VelocityBC(const std::vector<BCFunction> &velVx1BCs, const std::vector<BCFunction> &velVx2BCs,
                      const std::vector<BCFunction> &velVx3BCs);

    VelocityBC(const real &vx1, const real &vx1StartTime, const real &vx1EndTime, const real &vx2,
                      const real &vx2StartTime, const real &vx2EndTime, const real &vx3,
                      const real &vx3StartTime, const real &vx3EndTime);

    VelocityBC(const std::string &vx1Function, const real &vx1StartTime, const real &vx1EndTime,
                      const std::string &vx2Function, const real &vx2StartTime, const real &vx2EndTime,
                      const std::string &vx3Function, const real &vx3StartTime, const real &vx3EndTime);

    // methods
    void setTimePeriodic() { (this->type |= TIMEPERIODIC); }
    void unsetTimePeriodic() { (this->type &= ~TIMEPERIODIC); }
    bool isTimePeriodic() { return ((this->type & TIMEPERIODIC) == TIMEPERIODIC); }

    // The following is meant for moving objects...
    void setNewVelocities(const real &vx1, const real &vx1StartTime, const real &vx1EndTime, const real &vx2,
                          const real &vx2StartTime, const real &vx2EndTime, const real &vx3,
                          const real &vx3StartTime, const real &vx3EndTime);

    //------------- implements BC ----- start
    std::string toString();

    void init(const D3Q27Interactor *const &interactor, const real &time = 0) override;
    void update(const D3Q27Interactor *const &interactor, const real &time = 0) override;

    void adaptBCForDirection(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const real &worldX1,
                             const real &worldX2, const real &worldX3, const real &q, const int &fdirection,
                             const real &time = 0) override;
    void adaptBC(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const real &worldX1,
                 const real &worldX2, const real &worldX3, const real &time = 0) override;

    //------------- implements BC ----- end

    UbTupleDouble3 getVelocity(const real &x1, const real &x2, const real &x3, const real &timeStep) const;

protected:
    void init();
    void init(std::vector<BCFunction> &vxBCs);

    // time dependency is determined automatically via BCFunction intervals!
    void setTimeDependent() { (this->type |= TIMEDEPENDENT); }
    void unsetTimeDependent() { (this->type &= ~TIMEDEPENDENT); }

    void clear()
    {
        vx1BCs.clear();
        vx2BCs.clear();
        vx3BCs.clear();
        this->init();
    }
    void setNodeVelocity(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const real &worldX1,
                         const real &worldX2, const real &worldX3, const real &timestep);

private:
    mutable mu::value_type x1, x2, x3;
    mutable mu::value_type timeStep;

    mu::Parser *tmpVx1Function;
    mu::Parser *tmpVx2Function;
    mu::Parser *tmpVx3Function;

    std::vector<BCFunction> vx1BCs;
    std::vector<BCFunction> vx2BCs;
    std::vector<BCFunction> vx3BCs;
};

#endif

//! \}
