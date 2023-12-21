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
#ifndef Y_2D_Slice_To_RESULTS_H
#define Y_2D_Slice_To_RESULTS_H

#include "../ToVectorWriter.h"

#include <vector>
#include <memory>

class Parameter;
class SimulationResults;
struct VectorWriterInformationStruct;

class Y2dSliceToResults : public ToVectorWriter
{
public:
    static std::shared_ptr<Y2dSliceToResults> getNewInstance(std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo, unsigned int timeStepLength, std::shared_ptr<SimulationResults> simResults, unsigned int ySliceForCalculation);

    
private:
    Y2dSliceToResults();
    Y2dSliceToResults(std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo, unsigned int timeStepLength, std::shared_ptr<SimulationResults> simResults, unsigned int ySliceForCalculation);
    
    void writeTimestep(std::shared_ptr<Parameter> para, unsigned int t, int level);
    int CoordPara3DTo1D(int x, int y, int z);
    int CoordResults2DTo1D(int x, int z);

    std::shared_ptr<SimulationResults> simResults;
    unsigned int ySliceForCalculation;
    unsigned int maxX, maxY, maxZ;
};
#endif
//! \}
