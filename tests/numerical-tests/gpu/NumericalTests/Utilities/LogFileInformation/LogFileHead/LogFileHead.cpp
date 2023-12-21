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
#include "LogFileHead.h"

#include <iomanip>
#include <ctime>

#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

std::shared_ptr<LogFileHead> LogFileHead::getNewInstance(std::vector<int> devices)
{
    return std::shared_ptr<LogFileHead>(new LogFileHead(devices));
}

std::string LogFileHead::getOutput()
{
    calcDateAndTime();

    makeCenterHead("LogFile Information");
    oss << "Date=" << std::setw(2) << std::setfill('0') << nowLocal.tm_mday << "." << std::setw(2) << nowLocal.tm_mon + 1 << "." << nowLocal.tm_year + 1900 << std::endl;
    oss << "Time=" << std::setw(2) << std::setfill('0') << nowLocal.tm_hour << ":" << std::setw(2) << nowLocal.tm_min << ":" << std::setw(2) << nowLocal.tm_sec << std::endl;
    oss << std::endl;

    oss << "GPU_Devices=\"";
    for (int i = 0; i < devices.size(); i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, devices.at(i));
        std::string deviceName = prop.name;
        deviceName.assign(deviceName.begin(), remove_if(deviceName.begin(), deviceName.end(), &isspace));
        oss << deviceName;
        if (i < devices.size() - 1)
            oss << " ";
        else
            oss << "\"" << std::endl;
    }
    oss << std::endl;

    return oss.str();
}

void LogFileHead::calcDateAndTime()
{
    now = time(NULL);
    nowLocal = *localtime(&now);
}

LogFileHead::LogFileHead(std::vector<int> devices) : devices(devices)
{

}

//! \}
