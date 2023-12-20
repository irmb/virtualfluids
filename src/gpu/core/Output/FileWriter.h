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
//! \addtogroup gpu_Output Output
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "DataWriter.h"

class Parameter;
class CudaMemoryManager;

class FileWriter : public DataWriter
{
public:
    void writeInit(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaMemoryManager) override;
    void writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep) override;

private:
    void writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep, int level) override;
    std::vector<std::string> writeUnstructuredGridLT(std::shared_ptr<Parameter> para, int level,
                                                     std::vector<std::string>& fname);
    std::vector<std::string> writeUnstructuredGridMeanLT(std::shared_ptr<Parameter> para, int level,
                                                         std::vector<std::string>& fname);

    std::string writeCollectionFile(std::shared_ptr<Parameter> para, unsigned int timestep);

    std::string writeCollectionFileMean(std::shared_ptr<Parameter> para, unsigned int timestep);

    std::string writePvdCollectionFileForTimeSeries(const Parameter& para);

    std::vector<std::string> getNodeDataNames(std::shared_ptr<Parameter> para);
    std::vector<std::string> getMeanNodeDataNames(std::shared_ptr<Parameter> para);

    std::vector<std::string> fileNamesForCollectionFile;
    std::vector<std::string> fileNamesForCollectionFileMean;

    std::map<uint, std::vector<std::string>>
        fileNamesForCollectionFileTimeSeries; // key: timeStep, value: fileNames for timeStep
};
#endif

//! \}
