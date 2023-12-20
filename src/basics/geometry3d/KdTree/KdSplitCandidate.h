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
//! \addtogroup geometry3d
//! \ingroup basics
//! \{
//! \author Soeren Textor, Sebastian Bindick
//=======================================================================================
#ifndef KDSPLITCANDIDATE_H
#define KDSPLITCANDIDATE_H

#include <basics/utilities/UbMath.h>

namespace Kd
{
template <typename T>
class SplitCandidate
{
public:
    SplitCandidate() = default;
    /* ======================================================================================= */
    SplitCandidate(const int &axis, const T &position, const int &starting, const int &ending, const int &insidePlane)
        : axis(axis), position(position), starting(starting), ending(ending), np(insidePlane),
          isValid(true) // FIXME: isValid default false is correct?
    {
    }
    /* ======================================================================================= */
    bool operator!() const { return isValid; }
    /* ======================================================================================= */
    friend inline bool operator<(const SplitCandidate &lhs, const SplitCandidate &rhs)
    {
        return lhs.position < rhs.position;
    }
    /* ======================================================================================= */

public:
    int axis{ 0 };
    T Cn{ 0.0 };
    T position{ 0.0 };
    int nl{ 0 };
    int nr{ 0 };
    int np;
    int starting{ 0 };
    int ending{ 0 };
    bool np_left{ false };
    bool np_right{ false };
    bool isValid{ false };
};
} // namespace Kd

#endif // KDSPLITCANDIDATE_H

//! \}
