
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
//! \author Anna Wellmann
//=======================================================================================
#ifndef NEIGHBORDEBUGWRITER_H
#define NEIGHBORDEBUGWRITER_H

#include <string>

class Parameter;
struct LBMSimulationParameter;
class WbWriter;
struct QforDirectionalBoundaryCondition;
struct QforBoundaryConditions;

namespace neighbor_debug_writer
{

//! \brief Write the links to the neighbors as lines for all 27 directions.
void writeNeighborLinkLines(Parameter* para);

//! \brief Write the links to the neighbors as lines for the specified direction.
void writeNeighborLinkLinesForDirection(LBMSimulationParameter* parH, int direction, const std::string& filePath,
                                        WbWriter* writer);

void writeBoundaryConditionNeighbors(QforDirectionalBoundaryCondition* boundaryCondition, LBMSimulationParameter* parH,
                                     std::string& filePathBase);

void writeBoundaryConditionNeighbors(QforBoundaryConditions* boundaryCondition, LBMSimulationParameter* parH,
                                     std::string& filePathBase);

} // namespace neighbor_debug_writer

#endif

//! \}
