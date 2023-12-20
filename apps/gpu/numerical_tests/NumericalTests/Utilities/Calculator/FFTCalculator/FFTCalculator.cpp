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
#include "FFTCalculator.h"

#include "Utilities/Results/SimulationResults/SimulationResults.h"

#define _USE_MATH_DEFINES
#include <fstream>
#include <math.h>

std::shared_ptr<FFTCalculator> FFTCalculator::getInstance()
{
    static std::shared_ptr<FFTCalculator> uniqueInstance;
    if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<FFTCalculator>(new FFTCalculator());
    return uniqueInstance;
}

double FFTCalculator::calcNy(std::vector<std::vector<double>> data, bool transposeData, int lx, int lz,
                             int timeStepLength)
{
    this->lx = (double)lx;
    this->lz = (double)lz;
    this->timeStepLength = (double)timeStepLength;
    this->transposeData = transposeData;
    if (!transposeData)
        this->data = data;
    else
        this->data = transpose(data);

    init();

    double ny = calcNy();
    return ny;
}

double FFTCalculator::calcPhiDiff(std::vector<std::vector<double>> data, bool transposeData, int lx, int lz,
                                  int timeStepLength)
{
    this->lx = (double)lx;
    this->lz = (double)lz;
    this->timeStepLength = (double)timeStepLength;
    this->transposeData = transposeData;
    if (!transposeData)
        this->data = data;
    else
        this->data = transpose(data);

    init();

    double phidiff = calcPhiDiff();

    return abs(phidiff);
}

FFTCalculator::FFTCalculator()
{
}

double FFTCalculator::calcAmplitudeForTimeStep(std::vector<double> data, bool transposeData, int lx, int lz)
{
    this->lx = (double)lx;
    this->lz = (double)lz;
    init();
    this->transposeData = transposeData;
    this->data.resize(0);
    this->data.push_back(data);
    std::vector<double> amplitude = calcAmplitudeForAllSteps();
    return amplitude.at(0);
}

void FFTCalculator::init()
{
    fftResultsIm.clear();
    fftResultsRe.clear();
    fftCalculated = false;
}

double FFTCalculator::calcNy()
{
    std::vector<double> logAmplitude = calcLogAmplitudeForAllSteps();
    std::vector<double> linReg = calcLinReg(logAmplitude);
    double nu =
        -(1.0 / (((2.0 * M_PI / lz) * (2.0 * M_PI / lz) + (2.0 * M_PI / lx) * (2.0 * M_PI / lx)) * timeStepLength)) *
        linReg.at(0);
    return nu;
}

double FFTCalculator::calcPhiDiff()
{
    std::vector<double> phi = calcPhiForAllSteps();
    std::vector<double> linReg = calcLinReg(phi);

    return linReg.at(0);
}

std::vector<double> FFTCalculator::calcLinReg(std::vector<double> y)
{
    std::vector<double> result;
    std::vector<double> x(y.size());
    double sumX = 0.0;
    double sumY = 0.0;

    for (int i = 0; i < y.size(); i++) {
        sumY += y.at(i);
        x.at(i) = i;
        sumX += i;
    }
    double avgX = sumX / y.size();
    double avgY = sumY / y.size();
    double zaehler = 0.0;
    double nenner = 0.0;
    for (int i = 0; i < y.size(); i++) {
        zaehler += (x.at(i) - avgX) * (y.at(i) - avgY);
        nenner += (x.at(i) - avgX) * (x.at(i) - avgX);
    }
    double a1 = zaehler / nenner;
    result.push_back(a1);
    double a0 = avgY - a1 * avgX;
    result.push_back(a0);

    double ess = 0;
    double tss = 0;
    for (int i = 0; i < y.size(); i++) {
        ess += ((a0 + a1 * x.at(i)) - avgY) * ((a0 + a1 * x.at(i)) - avgY);
        tss += (y.at(i) - avgY) * (y.at(i) - avgY);
    }
    double r2 = ess / tss;
    result.push_back(r2);
    return result;
}

std::vector<double> FFTCalculator::calcLogAmplitudeForAllSteps()
{
    std::vector<double> amplitude = calcAmplitudeForAllSteps();
    std::vector<double> logAmplitude;
    for (int i = 0; i < amplitude.size(); i++)
        logAmplitude.push_back(log(amplitude.at(i)));

    return logAmplitude;
}

std::vector<double> FFTCalculator::calcAmplitudeForAllSteps()
{
    std::vector<double> amplitude;
    if (fftCalculated == false) {
        for (int step = 0; step < data.size(); step++)
            calcFFT2D(step);
        fftCalculated = true;
    }
    int pos;
    if (!transposeData)
        pos = 2 + (lx - 1);
    else
        pos = 2 + (lz - 1);

    for (int step = 0; step < data.size(); step++)
        amplitude.push_back(4.0 / (lx * lz) *
                            sqrt(fftResultsRe.at(step).at(pos) * fftResultsRe.at(step).at(pos) +
                                 fftResultsIm.at(step).at(pos) * fftResultsIm.at(step).at(pos)));

    return amplitude;
}

std::vector<double> FFTCalculator::calcPhiForAllSteps()
{
    std::vector<double> phi;
    if (fftCalculated == false) {
        for (int step = 0; step < data.size(); step++)
            calcFFT2D(step);
        fftCalculated = true;
    }
    int pos;
    if (!transposeData)
        pos = 2 + (lx - 1);
    else
        pos = 2 + (lz - 1);

    for (int step = 0; step < data.size(); step++) {
        phi.push_back(atan(fftResultsIm.at(step).at(pos) / fftResultsRe.at(step).at(pos)));
    }

    return phi;
}

void FFTCalculator::calcFFT2D(unsigned int timeStep)
{
    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * lx * lz);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * lx * lz);

    initDataForFFT(in, timeStep);

    fftw_plan p;
    if (!transposeData)
        p = fftw_plan_dft_2d(lz, lx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    else
        p = fftw_plan_dft_2d(lx, lz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    setFFTResults(out, timeStep);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

std::vector<std::vector<double>> FFTCalculator::transpose(std::vector<std::vector<double>> dataToTranspose)
{
    std::vector<std::vector<std::vector<double>>> dataInLx;
    dataInLx.resize(dataToTranspose.size());
    for (int i = 0; i < dataInLx.size(); i++) {
        dataInLx.at(i).resize(lz);
        for (int j = 0; j < dataInLx.at(i).size(); j++)
            dataInLx.at(i).at(j).resize(lx);
    }
    for (int timeStep = 0; timeStep < dataInLx.size(); timeStep++) {
        for (int posInLZ = 0; posInLZ < lz; posInLZ++)
            for (int posInLX = 0; posInLX < lx; posInLX++)
                dataInLx.at(timeStep).at(posInLZ).at(posInLX) = dataToTranspose.at(timeStep).at(posInLX + posInLZ * lx);
    }

    std::vector<std::vector<std::vector<double>>> dataInLz;
    dataInLz.resize(dataToTranspose.size());
    for (int i = 0; i < dataInLx.size(); i++) {
        dataInLz.at(i).resize(lx);
        for (int j = 0; j < dataInLz.at(i).size(); j++)
            dataInLz.at(i).at(j).resize(lz);
    }

    for (int timeStep = 0; timeStep < dataInLz.size(); timeStep++) {
        for (int posInLX = 0; posInLX < lx; posInLX++)
            for (int posInLZ = 0; posInLZ < lz; posInLZ++)
                dataInLz.at(timeStep).at(posInLX).at(posInLZ) = dataInLx.at(timeStep).at(posInLZ).at(posInLX);
    }

    std::vector<std::vector<double>> result;
    result.resize(dataToTranspose.size());

    for (int timeStep = 0; timeStep < dataInLz.size(); timeStep++) {
        result.at(timeStep).resize(0);
        for (int posInLX = 0; posInLX < lx; posInLX++)
            for (int posInLZ = 0; posInLZ < lz; posInLZ++)
                result.at(timeStep).push_back(dataInLz.at(timeStep).at(posInLX).at(posInLZ));
    }
    return result;
}

void FFTCalculator::initDataForFFT(fftw_complex *input, unsigned int step)
{
    for (int i = 0; i < data.at(step).size(); i++) {
        input[i][0] = data.at(step).at(i);
        input[i][1] = 0;
    }
}

void FFTCalculator::setFFTResults(fftw_complex *result, unsigned int step)
{
    std::vector<double> fftRe, fftIm;
    fftRe.resize(data.at(step).size());
    fftIm.resize(data.at(step).size());

    for (int i = 0; i < data.at(step).size(); i++) {
        fftRe.at(i) = result[i][0];
        fftIm.at(i) = result[i][1];
    }
    fftResultsIm.push_back(fftIm);
    fftResultsRe.push_back(fftRe);
}
//! \}
