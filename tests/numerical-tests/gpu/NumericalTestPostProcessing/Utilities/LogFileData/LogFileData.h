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
#ifndef LOGFILE_DATA_H
#define LOGFILE_DATA_H

#include "Simulation/BasicSimulation.h"

#include <memory>
#include <string>
#include <vector>

class L2NormLogFileData;
class L2NormBetweenKernelsLogFileData;
class NyLogFileData;
class PhiLogFileData;
class ShearWaveLogFileData;
class TaylorGreenVortexUxLogFileData;
class TaylorGreenVortexUzLogFileData;

class LogFileData
{
public:
    virtual std::vector<int> getAnalyticalVTKWritingTime() = 0;
    virtual std::vector<double> getBasicGridLengths() = 0;
    virtual int getBasisTimeStepLength() = 0;
    virtual std::string getDate() = 0;
    virtual std::vector<std::string> getGpuDevices() = 0;
    virtual std::string getKernel() = 0;
    virtual bool getL2NormTestRun() = 0;
    virtual bool getL2NormTestBetweenKernelRun() = 0;
    virtual int getNumberOfTimeSteps() = 0;
    virtual bool getNyTestRun() = 0;
    virtual bool getPhiTestRun() = 0;
    virtual std::vector<double> getResultCheckTime() = 0;
    virtual std::string getSimName() = 0;
    virtual std::vector<int> getSimTime() = 0;

    virtual std::vector<double> getTestTime() = 0;
    virtual std::string getTime() = 0;
    virtual double getViscosity() = 0;
    virtual bool getVTKFileWriting() = 0;

    virtual std::string getFilePath() = 0;
    virtual std::string getSimulationSigniture() = 0;

    virtual std::vector<std::shared_ptr<L2NormLogFileData> > getL2NormLogFileData() = 0;
    virtual std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > getL2NormBetweenKernelsLogFileData() = 0;
    virtual std::vector<std::shared_ptr<NyLogFileData> > getNyLogFileData() = 0;
    virtual std::vector<std::shared_ptr<PhiLogFileData> > getPhiLogFileData() = 0;

    virtual std::shared_ptr<TaylorGreenVortexUxLogFileData> getTaylorGreenVortexUxLogFileData() = 0;
    virtual std::shared_ptr<TaylorGreenVortexUzLogFileData> getTaylorGreenVortexUzLogFileData() = 0;
    virtual std::shared_ptr<ShearWaveLogFileData> getShearWaveLogFileData() = 0;

    virtual BasicSimulation getBasicSimulation() = 0;


    
};
#endif
//! \}
