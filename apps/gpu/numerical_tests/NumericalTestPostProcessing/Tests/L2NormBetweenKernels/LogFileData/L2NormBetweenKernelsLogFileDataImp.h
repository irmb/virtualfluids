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
#ifndef L2NORM_BETWEEN_KERNELS_LOG_FILE_DATA_IMP_H
#define L2NORM_BETWEEN_KERNELS_LOG_FILE_DATA_IMP_H

#include "L2NormBetweenKernelsLogFileData.h"

#include <memory>

class L2NormBetweenKernelsLogFileDataImp : public L2NormBetweenKernelsLogFileData
{
public:
    static std::shared_ptr<L2NormBetweenKernelsLogFileDataImp> getNewInstance();

    std::vector<double> getBasicGridLengths();
    std::string getBasicKernel();
    std::string getDivergentKernel();
    std::string getDataToCalculate();
    int getTimeStep();
    std::vector<double> getL2NormForBasicKernel();
    std::vector<double> getL2NormForDivergentKernel();
    std::vector<double> getL2NormBetweenKernels();
    std::string getNormalizeData();

    void setBasicGridLengths(std::vector<double> basicGridLengths);
    void setBasicKernel(std::string basicKernel);
    void setDivergentKernel(std::string divergentKernel);
    void setDataToCalculate(std::string dataToCalc);
    void setTimeStep(int timeStep);
    void setL2NormForBasicKernel(std::vector<double> l2Norm);
    void setL2NormForDivergentKernel(std::vector<double> l2Norm);
    void setL2NormBetweenKernels(std::vector<double> l2Norm);
    void setNormalizeData(std::string normalizeData);

    ~L2NormBetweenKernelsLogFileDataImp();

private:
    L2NormBetweenKernelsLogFileDataImp();

    std::vector<double> basicGridLengths;
    std::string basicKernel;
    std::string divergentKernel;
    std::string dataToCalc;
    int timeStep;
    std::vector<double> l2NormForBasicKernel;
    std::vector<double> l2NormForDivergentKernel;
    std::vector<double> l2NormBetweenKernels;
    std::string normalizeData;
};
#endif
//! \}
