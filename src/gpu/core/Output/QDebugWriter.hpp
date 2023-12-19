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
#ifndef QDEBUGWRITER_HPP
#define QDEBUGWRITER_HPP

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

#include <lbm/constants/D3Q27.h>

#include <basics/StringUtilities/StringUtil.h>
#include <basics/utilities/UbSystem.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"

namespace QDebugWriter
{

void writeQValues(QforBoundaryConditions& Q, int* k, int kq, const std::string& name)
{
    std::vector<std::vector<real>> qs;
    for (int j = 0; j < kq; j++) {
        uint32_t qKey = 0;
        std::vector<real> qNode;

        for (int i = 26; i >= 0; i--) {
            real q = Q.q27[i][j];
            if (q > 0) {
                qKey += (uint32_t)pow(2, 26 - i);
                qNode.push_back(q);
            }
        }
        if (qKey > 0) {
            float transportKey = *((float*)&qKey);
            qNode.push_back((real)transportKey);
            qNode.push_back((real)k[j]);
            qs.push_back(qNode);
        }
        qNode.clear();
    }

    SPtr<std::ofstream> outQ(new std::ofstream);
    outQ->open(name.c_str(), std::ios::out | std::ios::binary);

    for (std::size_t index = 0; index < qs.size(); index++) {
        std::vector<real> bcs = qs[index];
        uint32_t key = *((uint32_t*)&bcs[bcs.size() - 2]);
        int qIndex = (int)bcs[bcs.size() - 1];

        *outQ << qIndex << " " << key;

        for (std::size_t i = 0; i < bcs.size() - 2; i++) {
            *outQ << " " << std::fixed << std::setprecision(16) << bcs[i];
        }

        *outQ << "\n";
    }

    outQ->close();
}

} // namespace QDebugWriter

#endif

//! \}
