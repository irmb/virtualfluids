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
#ifndef FFTCALCULATOR_H
#define FFTCALCULATOR_H

#include <memory>
#include <vector>
#include <fftw3.h>

class SimulationResults;
class EvaluationParameter;
class TestResults;
class PhiAndNuTest;

class FFTCalculator
{
public:
    static std::shared_ptr<FFTCalculator> getInstance();
    
    double calcNy(std::vector<std::vector<double> > data, bool transposeData, int lx, int lz, int timeStepLength);
    double calcPhiDiff(std::vector<std::vector<double> > data, bool transposeData, int lx, int lz, int timeStepLength);

    double calcAmplitudeForTimeStep(std::vector<double> data, bool transposeData, int lx, int lz);

private:
    FFTCalculator();
    void init();
    double calcNy();
    double calcPhiDiff();
    std::vector<double> calcPhiForAllSteps();
    std::vector<double> calcLinReg(std::vector<double> y);
    std::vector<double> calcLogAmplitudeForAllSteps();
    std::vector<double> calcAmplitudeForAllSteps();
    void calcFFT2D(unsigned int step);
    std::vector<std::vector<double> > transpose(std::vector<std::vector<double> >);
    void initDataForFFT(fftw_complex* input, unsigned int step);
    void setFFTResults(fftw_complex* result, unsigned int step);

    std::vector<std::vector<double> > data;
    std::vector<std::vector<double> > fftResultsIm;
    std::vector<std::vector<double> > fftResultsRe;
    
    bool fftCalculated;
    bool transposeData;
    double lx, lz;
    double timeStepLength;
};
#endif

//! \}
