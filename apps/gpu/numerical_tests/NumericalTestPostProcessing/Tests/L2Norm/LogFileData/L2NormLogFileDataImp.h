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
#ifndef L2_NORM_LOG_FILE_DATA_IMP_H
#define L2_NORM_LOG_FILE_DATA_IMP_H

#include "L2NormLogFileData.h"

#include <memory>

class L2NormLogFileDataImp : public L2NormLogFileData
{
public:
    static std::shared_ptr<L2NormLogFileDataImp> getNewInstance();

    std::vector<double> getBasicGridLengths();
    std::string getDataToCalc();
    std::string getNormalizeData();
    int getBasicTimeStep();
    int getDivergentTimeStep();
    std::vector<double> getL2NormForBasicTimeStep();
    std::vector<double> getL2NormForDivergentTimeStep();
    std::vector<double> getL2NormDiff();

    void setBasicGridLengths(std::vector<double> basicGridLengths);
    void setDataToCalc(std::string dataToCalc);
    void setNormalizeData(std::string normalizeData);
    void setBasicTimeStep(int basicTimeStep);
    void setDivergentTimeStep(int divergentTimeStep);
    void setL2NormForBasicTimeStep(std::vector<double> l2NormForBasicTimeStep);
    void setL2NormForDivergentTimeStep(std::vector<double> l2NormForDivergentTimeStep);
    void setL2NormDiff(std::vector<double> l2NormDiff);

    ~L2NormLogFileDataImp();

private:
    L2NormLogFileDataImp();

    std::vector<double> basicGridLengths;
    std::string dataToCalc;
    std::string normalizeData;
    int basicTimeStep;
    int divergentTimeStep;
    std::vector<double> l2NormForBasicTimeStep;
    std::vector<double> l2NormForDivergentTimeStep;
    std::vector<double> l2NormDiff;
};
#endif
//! \}
