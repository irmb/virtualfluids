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
#ifndef RESULTS_IMP_H
#define RESULTS_IMP_H

#include "Results.h"

class ResultsImp : public Results
{
public:
    int getNumberOfTimeSteps();
    int getTimeStepLength();
    std::vector<unsigned int> getTimeSteps();
    std::vector<int> getTime();
    std::vector<std::vector<double> > getVx();
    std::vector<std::vector<double> > getVy();
    std::vector<std::vector<double> > getVz();
    int getNumberOfXNodes();
    int getNumberOfYNodes();
    int getNumberOfZNodes();
    std::vector<std::vector<double> > getXNodes();
    std::vector<std::vector<double> > getYNodes();
    std::vector<std::vector<double> > getZNodes();
    std::vector<std::vector<unsigned int> > getLevels();
    std::vector<std::vector<double> > getPress();
    std::vector<std::vector<double> > getRho();
    int getL0();

    bool checkYourData();

protected:
    ResultsImp(int l0);
    ResultsImp() = default;

    unsigned int numberOfTimeSteps;
    unsigned int timeStepLength;
    unsigned int xNodes, yNodes, zNodes;
    unsigned int numberOfNodes;

    std::vector<unsigned int> timeStep;
    std::vector<int> time;
    std::vector<std::vector<double> > x, y, z;
    std::vector<std::vector<double> > vx, vy, vz;
    std::vector<std::vector<double> > press;
    std::vector<std::vector<double> > rho;
    std::vector<std::vector<unsigned int> > level;

    int l0;

private:
};
#endif
//! \}
