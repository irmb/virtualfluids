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
#ifndef PHI_LOG_FILE_DATA_IMP_H
#define PHI_LOG_FILE_DATA_IMP_H

#include "PhiLogFileData.h"

#include <memory>

class PhiLogFileDataImp : public PhiLogFileData
{
public:
    static std::shared_ptr<PhiLogFileDataImp> getNewInstance();

    std::vector<double> getBasicGridLengths();
    int getStartTimeStepCalculation();
    int getEndTimeStepCalculation();
    std::string getDataToCalc();
    std::vector<double> getPhiDiff();
    std::vector<std::vector<double> > getOrderOfAccuracy();

    void setBasicGridLengths(std::vector<double> basicGridLengths);
    void setStartTimeStepCalculation(int startTimeStepCalculation);
    void setEndTimeStepCalculation(int endTimeStepCalculation);
    void setDataToCalc(std::string dataToCalcPhiAndNu);
    void setPhiDiff(std::vector<double> phiDiff);
    void setOrderOfAccuracy(std::vector<std::vector<double> > orderOfAccuracy);

private:
    PhiLogFileDataImp();

    std::vector<double> basicGridLengths;
    int startTimeStepCalculation;
    int endTimeStepCalculation;
    std::string dataToCalc;
    std::vector<double> phiDiff;
    std::vector<std::vector<double> > orderOfAccuracy;
};
#endif
//! \}
