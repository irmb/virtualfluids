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
//! \addtogroup gpu_Samplers Utilities
//! \ingroup gpu_core core
//! \{

#ifndef SAMPLE_UTILITIES_H
#define SAMPLE_UTILITIES_H
#include <basics/StringUtilities/StringUtil.h>
#include <filesystem>
#include <fstream>
#include <ios>
#include <stdexcept>
#include <string>
#include <vector>

template <typename T>
inline std::string nameComponent(const std::string& name, T value)
{
    return "_" + name + "_" + StringUtil::toString<T>(value);
}

inline std::string makeParallelFileName(const std::string& probeName, int id, int t)
{
    return probeName + "_bin" + nameComponent("ID", id) + nameComponent("t", t) + ".vtk";
}

inline std::string makeGridFileName(const std::string& probeName, int level, int id, int t, uint part)
{
    return probeName + "_bin" + nameComponent("lev", level) + nameComponent("ID", id) + nameComponent<int>("Part", part) +
           nameComponent("t", t) + ".vtk";
}

inline std::string makeTimeseriesFileName(const std::string& probeName, int level, int id)
{
    return probeName + "_timeseries" + nameComponent("lev", level) + nameComponent("ID", id) + ".txt";
}

template <typename T>
__host__ __device__ inline T computeNewTimeAverage(T oldAverage, T newValue, real inverseNumberOfTimesteps)
{
    return oldAverage + (newValue - oldAverage) * inverseNumberOfTimesteps;
}

inline void writeTimeseriesFileHeader(const std::string& fileName, int numberOfPoints,
                                      std::vector<std::string>& variableNames, const real* coordsX, const real* coordsY,
                                      const real* coordsZ)
{
    std::filesystem::create_directories(std::filesystem::path(fileName).parent_path());
    std::ofstream out(fileName.c_str(), std::ios::out | std::ios::binary);

    if (!out.is_open())
        throw std::runtime_error("Could not open timeseries file " + fileName + "!");

    out << "TimeseriesOutput \n";
    out << "Quantities: ";
    for (const std::string& name : variableNames)
        out << name << ", ";
    out << "\n";
    out << "Number of points in this file: \n";
    out << numberOfPoints << "\n";
    out << "Positions: x, y, z\n";
    for (int i = 0; i < numberOfPoints; i++)
        out << coordsX[i] << ", " << coordsY[i] << ", " << coordsZ[i] << "\n";

    out.close();
}

//! \brief Write data to timeseries file, that can be read by TimeseriesFileReader.
//! Layout of the file is:
//! TimeseriesOutput
//! Quantities: Quant1 Quant2 Quant3
//! Positions:
//! point1.x, point1.y, point1.z
//! point2.x, point2.y, point2.z
//! ...
//! t0 point1.quant1 point2.quant1 ... point1.quant2 point2.quant2 ...
//! t1 point1.quant1 point2.quant1 ... point1.quant2 point2.quant2 ...
inline void appendDataToTimeseriesFile(const std::string& fileName, std::vector<std::vector<real>>& data)
{
    std::ofstream out(fileName.c_str(), std::ios::app | std::ios::binary);
    for (auto& timestepData : data) {
        out.write((char*)timestepData.data(), sizeof(real) * timestepData.size());
    }
    out.close();
}

#endif

//! \}