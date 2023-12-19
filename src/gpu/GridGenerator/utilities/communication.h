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
//! \addtogroup gpu_utilities utilities
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef Communication_H
#define Communication_H

#include <basics/geometry3d/Axis.h>

#include "grid/BoundaryConditions/Side.h"

// has to have the same order as SideType in Side.h
namespace CommunicationDirections
{
enum CommunicationDirection {
    MX = static_cast<int>(SideType::MX),
    PX = static_cast<int>(SideType::PX),
    MY = static_cast<int>(SideType::MY),
    PY = static_cast<int>(SideType::PY),
    MZ = static_cast<int>(SideType::MZ),
    PZ = static_cast<int>(SideType::PZ)
};

bool isNegative(CommunicationDirection direction);
bool isPositive(CommunicationDirection direction);

const std::map<CommunicationDirection, Axis> communicationDirectionToAxes { { MX, Axis::x }, { PX, Axis::x },
                                                                            { MY, Axis::y }, { PY, Axis::y },
                                                                            { MZ, Axis::z }, { PZ, Axis::z } };

CommunicationDirection getNegativeDirectionAlongAxis(Axis axis);
CommunicationDirection getPositiveDirectionAlongAxis(Axis axis);

} // namespace CommunicationDirections

#endif // Communication_H
//! \}
