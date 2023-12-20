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
#ifndef LOGFILE_DATA_IMP_H
#define LOGFILE_DATA_IMP_H

#include "LogFileData.h"

#include <memory>

class LogFileDataImp : public LogFileData
{
public:
    static std::shared_ptr<LogFileDataImp> getNewInstance();

    std::vector<int> getAnalyticalVTKWritingTime();
    std::vector<double> getBasicGridLengths();
    int getBasisTimeStepLength();
    std::string getDate();
    std::vector<std::string> getGpuDevices();
    std::string getKernel();
    bool getL2NormTestRun();
    bool getL2NormTestBetweenKernelRun();
    int getNumberOfTimeSteps();
    bool getNyTestRun();
    bool getPhiTestRun();
    std::vector<double> getResultCheckTime();
    std::string getSimName();
    std::vector<int> getSimTime();
    std::vector<double> getTestTime();
    std::string getTime();
    double getViscosity();
    bool getVTKFileWriting();
    std::string getFilePath();
    std::string getSimulationSigniture();

    std::vector<std::shared_ptr<L2NormLogFileData> > getL2NormLogFileData();
    std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > getL2NormBetweenKernelsLogFileData();
    std::vector<std::shared_ptr<PhiLogFileData> > getPhiLogFileData();
    std::vector<std::shared_ptr<NyLogFileData> > getNyLogFileData();
    
    std::shared_ptr<TaylorGreenVortexUxLogFileData> getTaylorGreenVortexUxLogFileData();
    std::shared_ptr<TaylorGreenVortexUzLogFileData> getTaylorGreenVortexUzLogFileData();
    std::shared_ptr<ShearWaveLogFileData> getShearWaveLogFileData();
    BasicSimulation getBasicSimulation();

    void setAnalyticalVTKWritingTime(std::vector<int> analyticalVTKWritingTime);
    void setBasicGridLengths(std::vector<double> basicGridLenghts);
    void setBasisTimeStepLength(int basisTimeStepLength);
    void setDate(std::string date);
    void setGpuDevices(std::vector<std::string> gpuDevices);
    void setKernel(std::string kernelName);
    void setL2NormTestRun(bool l2NormTestRun);
    void setL2NormTestBetweenKernelRun(bool l2NormTestBetweenKernelRun);
    void setNumberOfTimeSteps(int numberOfTimeSteps);
    void setPhiTestRun(bool phiTestRun);
    void setNyTestRun(bool nyTestRun);
    void setResultCheckTime(std::vector<double> resultsCheckTime);
    void setSimName(std::string simName);
    void setSimTime(std::vector<int> simTime);
    void setTestTime(std::vector<double> testTime);
    void setTime(std::string time);
    void setViscosity(double viscosity);
    void setVTKFileWriting(bool vtkFileWriting);
    void setFilePath(std::string filePath);
    void setSimulationSigniture(std::string simulationSigniture);

    void setL2NormLogFileData(std::vector<std::shared_ptr<L2NormLogFileData> > data);
    void setL2NormBetweenKernelsLogFileData(std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > data);
    void setNyLogFileData(std::vector<std::shared_ptr<NyLogFileData> > data);
    void setPhiLogFileData(std::vector<std::shared_ptr<PhiLogFileData> > data);
    
    void setTaylorGreenVortexUxLogFileData(std::shared_ptr<TaylorGreenVortexUxLogFileData> data);
    void setTaylorGreenVortexUzLogFileData(std::shared_ptr<TaylorGreenVortexUzLogFileData> data);
    void setShearWaveLogFileData(std::shared_ptr<ShearWaveLogFileData> data);
    void setBasicSimulation(BasicSimulation sim);

    ~LogFileDataImp();

private:
    LogFileDataImp();

    std::vector<int> analyticalVTKWritingTime;
    std::vector<double> basicGridLenghts;
    int basisTimeStepLength;
    std::string date;
    std::vector<std::string> gpuDevices;
    std::string kernelName;
    bool l2NormTestRun;
    bool l2NormTestBetweenKernelRun;
    int numberOfTimeSteps;
    bool phiTestRun;
    bool nyTestRun;
    std::vector<double> resultsCheckTime;
    std::string simName;
    std::vector<int> simTime;
    std::vector<double> testTime;
    std::string time;
    double viscosity;
    bool vtkFileWriting;
    std::string filePath;
    std::string simulationSigniture;

    std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > l2NormBetweenKernelsDataLogFileData;
    std::vector<std::shared_ptr<PhiLogFileData> > phiLogFileData;
    std::vector<std::shared_ptr<NyLogFileData> > nyLogFileData;
    std::vector<std::shared_ptr<L2NormLogFileData> > l2NormLogFileData;

    std::shared_ptr<ShearWaveLogFileData> shearWaveLogFileData;
    std::shared_ptr<TaylorGreenVortexUxLogFileData> tgvUxLogFileData;
    std::shared_ptr<TaylorGreenVortexUzLogFileData> tgvUzLogFileData;

    BasicSimulation sim;
};
#endif
//! \}
