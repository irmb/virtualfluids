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
//! \addtogroup NumericalTests
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#ifndef INITIAL_CONDITION_IMP_H
#define INITIAL_CONDITION_IMP_H

#include "InitialCondition.h"

#include <vector>
#include <memory>

class Parameter;

class InitialConditionImp : public InitialCondition
{
public:
    void setParameter(std::shared_ptr<Parameter> para);
    void init(const int level);
    virtual real getInitVX(int i, int level) = 0;
    virtual real getInitVY(int i, int level) = 0;
    virtual real getInitVZ(int i, int level) = 0;
    virtual real getInitROH(int i, int level) = 0;
    virtual real getInitPRESS(int i, int level) = 0;

protected:
    InitialConditionImp() {};
    real getXCoord(int i, int level);
    real getYCoord(int i, int level);
    real getZCoord(int i, int level);

    std::shared_ptr<Parameter> para;
    real XCoordStopNode, YCoordStopNode, ZCoordStopNode;

};
#endif
//! \}
