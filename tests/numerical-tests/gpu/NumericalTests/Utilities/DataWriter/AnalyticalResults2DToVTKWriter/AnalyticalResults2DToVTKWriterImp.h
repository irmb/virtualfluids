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
#ifndef ANALYTICAL_RESULTS_2D_TO_VTK_WRITER_IMP_H
#define ANALYTICAL_RESULTS_2D_TO_VTK_WRITER_IMP_H

#include "AnalyticalResults2DToVTKWriter.h"

#include <vector>


class AnalyticalResults2DToVTKWriterImp : public AnalyticalResults2DToVTKWriter
{
public:
    static std::shared_ptr<AnalyticalResults2DToVTKWriterImp> getInstance(bool writeAnalyticalResults);

    void writeAnalyticalResult(std::shared_ptr<Parameter> para, std::shared_ptr<AnalyticalResults> analyticalResult);

private:
    AnalyticalResults2DToVTKWriterImp() {};
    AnalyticalResults2DToVTKWriterImp(bool writeAnalyticalResults);

    void writeTimeStep(std::shared_ptr<Parameter> para, std::shared_ptr<AnalyticalResults> analyticalResult, int level, std::vector<std::string> & fname, int timeStep);
    int CoordResults2DTo1D(int x, int z);

    std::shared_ptr<Parameter> para;
    int maxX, maxY, maxZ;
    bool writeAnalyticalResults;
};
#endif
//! \}
