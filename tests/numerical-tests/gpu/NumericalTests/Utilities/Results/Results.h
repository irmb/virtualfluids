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
#ifndef RESULTS_H
#define RESULTS_H

#include <vector>

class Results
{
public:
    virtual ~Results() = default;
    virtual int getNumberOfTimeSteps() = 0;
    virtual std::vector<std::vector<double> > getVx() = 0;
    virtual std::vector<std::vector<double> > getVy() = 0;
    virtual std::vector<std::vector<double> > getVz() = 0;
    virtual int getNumberOfXNodes() = 0;
    virtual int getNumberOfYNodes() = 0;
    virtual int getNumberOfZNodes() = 0;
    virtual std::vector<std::vector<double> > getXNodes() = 0;
    virtual std::vector<std::vector<double> > getYNodes() = 0;
    virtual std::vector<std::vector<double> > getZNodes() = 0;
    virtual int getTimeStepLength() = 0;
    virtual std::vector<unsigned int> getTimeSteps() = 0;
    virtual std::vector < std::vector<unsigned int> > getLevels() = 0;

    virtual bool checkYourData() = 0;

    virtual int getL0() = 0;

};
#endif
//! \}
