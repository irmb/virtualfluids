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
#include "communication.h"

using namespace communication_directions;

bool communication_directions::isNegative(CommunicationDirection direction)
{
    return direction == CommunicationDirection::MX || direction == CommunicationDirection::MY ||
           direction == CommunicationDirection::MZ;
}

bool communication_directions::isPositive(CommunicationDirection direction)
{
    return direction == CommunicationDirection::PX || direction == CommunicationDirection::PY ||
           direction == CommunicationDirection::PZ;
}

CommunicationDirection communication_directions::getNegativeDirectionAlongAxis(Axis axis)
{
    switch (axis) {
        case Axis::x:
            return MX;
            break;
        case Axis::y:
            return MY;
            break;
        case Axis::z:
            return MZ;
            break;
        default:
            throw std::runtime_error("Unknown coordinate direction" + axis::to_string(axis));
    }
}

CommunicationDirection communication_directions::getPositiveDirectionAlongAxis(Axis axis)
{
    switch (axis) {
        case Axis::x:
            return PX;
            break;
        case Axis::y:
            return PY;
            break;
        case Axis::z:
            return PZ;
            break;
        default:
            throw std::runtime_error("Unknown coordinate direction" + axis::to_string(axis));
    }
}
//! \}
