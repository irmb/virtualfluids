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
//! \addtogroup NumericalTestPostProcessing
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#include "LogFileDataImp.h"

std::shared_ptr<LogFileDataImp> LogFileDataImp::getNewInstance()
{
    return std::shared_ptr<LogFileDataImp>(new LogFileDataImp());
}

void LogFileDataImp::setBasisTimeStepLength(int basisTimeStepLength)
{
    this->basisTimeStepLength = basisTimeStepLength;
}

void LogFileDataImp::setAnalyticalVTKWritingTime(std::vector<int> analyticalVTKWritingTime)
{
    this->analyticalVTKWritingTime = analyticalVTKWritingTime;
}

void LogFileDataImp::setBasicGridLengths(std::vector<double> basicGridLenghts)
{
    this->basicGridLenghts = basicGridLenghts;
}

void LogFileDataImp::setDate(std::string date)
{
    this->date = date;
}

void LogFileDataImp::setGpuDevices(std::vector<std::string> gpuDevices)
{
    this->gpuDevices = gpuDevices;
}

void LogFileDataImp::setSimName(std::string simName)
{
    this->simName = simName;
}

void LogFileDataImp::setSimTime(std::vector<int> simTime)
{
    this->simTime = simTime;
}

void LogFileDataImp::setTestTime(std::vector<double> testTime)
{
    this->testTime = testTime;
}

void LogFileDataImp::setTime(std::string time)
{
    this->time = time;
}

void LogFileDataImp::setKernel(std::string kernelName)
{
    this->kernelName = kernelName;
}

void LogFileDataImp::setL2NormTestRun(bool l2NormTestRun)
{
    this->l2NormTestRun = l2NormTestRun;
}

void LogFileDataImp::setL2NormTestBetweenKernelRun(bool l2NormTestBetweenKernelRun)
{
    this->l2NormTestBetweenKernelRun = l2NormTestBetweenKernelRun;
}

void LogFileDataImp::setNumberOfTimeSteps(int numberOfTimeSteps)
{
    this->numberOfTimeSteps = numberOfTimeSteps;
}

void LogFileDataImp::setPhiTestRun(bool phiTestRun)
{
    this->phiTestRun = phiTestRun;
}

void LogFileDataImp::setNyTestRun(bool nyTestRun)
{
    this->nyTestRun = nyTestRun;
}

void LogFileDataImp::setResultCheckTime(std::vector<double> resultsCheckTime)
{
    this->resultsCheckTime = resultsCheckTime;
}

void LogFileDataImp::setViscosity(double viscosity)
{
    this->viscosity = viscosity;
}

void LogFileDataImp::setVTKFileWriting(bool vtkFileWriting)
{
    this->vtkFileWriting = vtkFileWriting;
}

void LogFileDataImp::setFilePath(std::string filePath)
{
    this->filePath = filePath;
}

void LogFileDataImp::setSimulationSigniture(std::string simulationSigniture)
{
    this->simulationSigniture = simulationSigniture;
}

void LogFileDataImp::setL2NormLogFileData(std::vector<std::shared_ptr<L2NormLogFileData>> data)
{
    this->l2NormLogFileData = data;
}

void LogFileDataImp::setL2NormBetweenKernelsLogFileData(std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData>> data)
{
    this->l2NormBetweenKernelsDataLogFileData = data;
}

void LogFileDataImp::setNyLogFileData(std::vector<std::shared_ptr<NyLogFileData>> data)
{
    this->nyLogFileData = data;
}

void LogFileDataImp::setPhiLogFileData(std::vector<std::shared_ptr<PhiLogFileData>> data)
{
    this->phiLogFileData = data;
}

void LogFileDataImp::setTaylorGreenVortexUxLogFileData(std::shared_ptr<TaylorGreenVortexUxLogFileData> data)
{
    this->tgvUxLogFileData = data;
}

void LogFileDataImp::setTaylorGreenVortexUzLogFileData(std::shared_ptr<TaylorGreenVortexUzLogFileData> data)
{
    this->tgvUzLogFileData = data;
}

void LogFileDataImp::setShearWaveLogFileData(std::shared_ptr<ShearWaveLogFileData> data)
{
    this->shearWaveLogFileData = data;
}

void LogFileDataImp::setBasicSimulation(BasicSimulation sim)
{
    this->sim = sim;
}

LogFileDataImp::~LogFileDataImp()
{
}

int LogFileDataImp::getBasisTimeStepLength()
{
    return basisTimeStepLength;
}

std::vector<int> LogFileDataImp::getAnalyticalVTKWritingTime()
{
    return analyticalVTKWritingTime;
}

std::vector<double> LogFileDataImp::getBasicGridLengths()
{
    return basicGridLenghts;
}

std::string LogFileDataImp::getDate()
{
    return date;
}

std::vector<std::string> LogFileDataImp::getGpuDevices()
{
    return gpuDevices;
}

std::string LogFileDataImp::getKernel()
{
    return kernelName;
}

bool LogFileDataImp::getL2NormTestRun()
{
    return l2NormTestRun;
}

bool LogFileDataImp::getL2NormTestBetweenKernelRun()
{
    return l2NormTestBetweenKernelRun;
}

int LogFileDataImp::getNumberOfTimeSteps()
{
    return numberOfTimeSteps;
}

bool LogFileDataImp::getNyTestRun()
{
    return nyTestRun;
}

bool LogFileDataImp::getPhiTestRun()
{
    return phiTestRun;
}

std::vector<double> LogFileDataImp::getResultCheckTime()
{
    return resultsCheckTime;
}

std::string LogFileDataImp::getSimName()
{
    return simName;
}

std::vector<int> LogFileDataImp::getSimTime()
{
    return simTime;
}

std::vector<double> LogFileDataImp::getTestTime()
{
    return testTime;
}

std::string LogFileDataImp::getTime()
{
    return time;
}

double LogFileDataImp::getViscosity()
{
    return viscosity;
}

bool LogFileDataImp::getVTKFileWriting()
{
    return vtkFileWriting;
}

std::string LogFileDataImp::getFilePath()
{
    return filePath;
}

std::string LogFileDataImp::getSimulationSigniture()
{
    return simulationSigniture;
}

std::vector<std::shared_ptr<L2NormLogFileData>> LogFileDataImp::getL2NormLogFileData()
{
    return l2NormLogFileData;
}

std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData>> LogFileDataImp::getL2NormBetweenKernelsLogFileData()
{
    return l2NormBetweenKernelsDataLogFileData;
}

std::vector<std::shared_ptr<PhiLogFileData>> LogFileDataImp::getPhiLogFileData()
{
    return phiLogFileData;
}

std::vector<std::shared_ptr<NyLogFileData>> LogFileDataImp::getNyLogFileData()
{
    return nyLogFileData;
}

std::shared_ptr<TaylorGreenVortexUxLogFileData> LogFileDataImp::getTaylorGreenVortexUxLogFileData()
{
    return tgvUxLogFileData;
}

std::shared_ptr<TaylorGreenVortexUzLogFileData> LogFileDataImp::getTaylorGreenVortexUzLogFileData()
{
    return tgvUzLogFileData;
}

std::shared_ptr<ShearWaveLogFileData> LogFileDataImp::getShearWaveLogFileData()
{
    return shearWaveLogFileData;
}

BasicSimulation LogFileDataImp::getBasicSimulation()
{
    return sim;
}

LogFileDataImp::LogFileDataImp()
{
}

//! \}
