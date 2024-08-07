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
//! \addtogroup gpu_Restart Restart
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "RestartObject.h"

#include <fstream>
#include <filesystem>

#include "Parameter/Parameter.h"
#include "basics/utilities/UbMath.h"
#include "logger/Logger.h"

void RestartObject::deserialize(const std::string &filename, std::shared_ptr<Parameter>& para)
{
    deserialize_internal(filename);

    for (int index1 = para->getCoarse(); index1 <= para->getFine(); index1++) {
        std::vector<real> vec;
        fs.push_back(vec);

        for (size_t index2 = 0; index2 < (para->getD3Qxx() * para->getParH(index1)->numberOfNodes); index2++) {
            para->getParH(index1)->distributions.f[0][index2] = fs[index1][index2];
        }
    }
}

void RestartObject::serialize(const std::string &filename, const std::shared_ptr<Parameter>& para)
{
    if (!fs.empty()) {
        clear(para);
    }
    for (int index1 = para->getCoarse(); index1 <= para->getFine(); index1++) {
        std::vector<real> vec;
        fs.push_back(vec);

        for (size_t index2 = 0; index2 < (para->getD3Qxx() * para->getParH(index1)->numberOfNodes); index2++) {
            if (ub_math::isNaN(para->getParH(index1)->distributions.f[0][index2])) {
                fs[index1].push_back((real)0.0);
            } else {
                fs[index1].push_back(para->getParH(index1)->distributions.f[0][index2]);
            }
        }
    }

    serialize_internal(filename);
}

void RestartObject::clear(const std::shared_ptr<Parameter>& para)
{
    for (int j = para->getCoarse(); j <= para->getFine(); j++) {
        fs[j].resize(0);
    }
    fs.resize(0);
}

void RestartObject::delete_restart_file(const std::string& filename)
{
    const std::string file = filename + getFileExtension();
    const bool wasDeleted = std::filesystem::remove(file.data());
    if (!wasDeleted)
        VF_LOG_WARNING("Could not delete the restart object file \"{}\".", file);
}

//////////////////////////////////////////////////////////////////////////

std::string ASCIIRestartObject::getFileExtension() const
{
    return ".txt";
}

void ASCIIRestartObject::serialize_internal(const std::string &filename)
{
    std::ofstream stream(filename + getFileExtension());

    if (!stream) {
        stream.clear();
        std::string path = ub_system::getPathFromString(filename);
        if (!path.empty()) {
            ub_system::makeDirectory(path);
            stream.open(filename.c_str());
        }

        if (!stream)
            throw UbException(UB_EXARGS, "couldn't open file " + filename);
    }

    for (std::vector<std::vector<real>>::size_type i = 0; i < fs.size(); i++) {
        for (std::vector<real>::size_type j = 0; j < fs[i].size(); j++) {
            stream << fs[i][j] << " ";
        }
        stream << "\n";
    }

    stream.close();
}

void ASCIIRestartObject::deserialize_internal(const std::string &filename)
{
    std::ifstream stream(filename + getFileExtension(), std::ios_base::in);

    if (!stream.is_open())
        throw UbException(UB_EXARGS, "couldn't open check point file " + filename);

    std::string line;
    while (std::getline(stream, line)) {
        std::vector<real> lineData;
        std::stringstream lineStream(line);

        real value;
        while (lineStream >> value) {
            lineData.push_back(value);
        }
        fs.push_back(lineData);
    }

    stream.close();
}

std::string BinaryRestartObject::getFileExtension() const
{
    return ".bin";
}

void BinaryRestartObject::serialize_internal(const std::string &filename)
{
    std::ofstream stream(filename + getFileExtension(), std::ios_base::binary);

    if (!stream) {
        stream.clear();
        std::string path = ub_system::getPathFromString(filename);
        if (!path.empty()) {
            ub_system::makeDirectory(path);
            stream.open(filename.c_str());
        }

        if (!stream)
            throw UbException(UB_EXARGS, "couldn't open file " + filename);
    }

    // Store size of the outer vector
    int s1 = (int)fs.size();
    stream.write(reinterpret_cast<const char *>(&s1), sizeof(s1));

    // Now write each vector one by one
    for (auto &v : fs) {
        // Store its size
        int size = (int)v.size();
        stream.write(reinterpret_cast<const char *>(&size), sizeof(size));

        // Store its contents
        stream.write(reinterpret_cast<const char *>(&v[0]), v.size() * sizeof(real));
    }

    stream.close();
}

void BinaryRestartObject::deserialize_internal(const std::string &filename)
{
    std::ifstream stream(filename + getFileExtension(), std::ios_base::in | std::ifstream::binary);

    if (!stream.is_open())
        throw UbException(UB_EXARGS, "couldn't open check point file " + filename);

    int size = 0;
    stream.read(reinterpret_cast<char *>(&size), sizeof(size));
    fs.resize(size);
    for (int n = 0; n < size; ++n) {
        int size2 = 0;
        stream.read(reinterpret_cast<char *>(&size2), sizeof(size2));
        real f;
        for (int k = 0; k < size2; ++k) {
            stream.read(reinterpret_cast<char *>(&f), sizeof(f));
            fs[n].push_back(f);
        }
    }

    stream.close();
}

//! \}
