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
#ifndef POST_PROCESSING_CONFIG_DATA_IMP_H
#define POST_PROCESSING_CONFIG_DATA_IMP_H

#include "PostProcessingConfigData.h"

class PostProcessingConfigDataImp : public PostProcessingConfigData
{
public:
    static std::shared_ptr<PostProcessingConfigDataImp> getNewInstance();

    std::vector<BasicSimulation> getSimulations();
    std::vector<Assistant> getAssistants();
    std::vector<DataCombination> getDataCombinations();
    std::string getMathematicaFilePath();
    std::string getLogFilesPath();

    void setSimulations(std::vector<BasicSimulation> simulations);
    void setAssistants(std::vector<Assistant> assis);
    void setDataCombinations(std::vector<DataCombination> dataComb);
    void setMathematicaFilePath(std::string mathematicaFilePath);
    void setLogFilesPath(std::string logFilesPath);

private:
    PostProcessingConfigDataImp();

    std::vector<BasicSimulation> simulations;
    std::vector<Assistant> assistants;
    std::vector<DataCombination> dataCombinations;
    std::string mathematicaFilePath;
    std::string logFilesPath;
};
#endif
//! \}
