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
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#include <GbObject3D.h>
#include <GbPoint3D.h>
#include <basics/utilities/UbMath.h>

using namespace std;

/*======================================================================*/
bool GbObject3D::isPointInGbObject3D(GbPoint3D *p)
{
    return this->isPointInGbObject3D(p->getX1Centroid(), p->getX2Coordinate(), p->getX3Coordinate());
}
/*======================================================================*/
bool GbObject3D::isCellInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                        const double &x2b, const double &x3b)
{

    if (this->isPointInGbObject3D(x1a, x2a, x3a) && this->isPointInGbObject3D(x1b, x2a, x3a) &&
        this->isPointInGbObject3D(x1b, x2b, x3a) && this->isPointInGbObject3D(x1a, x2b, x3a) &&
        this->isPointInGbObject3D(x1a, x2a, x3b) && this->isPointInGbObject3D(x1b, x2a, x3b) &&
        this->isPointInGbObject3D(x1b, x2b, x3b) && this->isPointInGbObject3D(x1a, x2b, x3b)) {
        return true;
    }

    return false;
}
/*======================================================================*/
bool GbObject3D::isCellCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b)
{
    if (this->isPointInGbObject3D(x1a, x2a, x3a) || this->isPointInGbObject3D(x1b, x2a, x3a) ||
        this->isPointInGbObject3D(x1b, x2b, x3a) || this->isPointInGbObject3D(x1a, x2b, x3a) ||
        this->isPointInGbObject3D(x1a, x2a, x3b) || this->isPointInGbObject3D(x1b, x2a, x3b) ||
        this->isPointInGbObject3D(x1b, x2b, x3b) || this->isPointInGbObject3D(x1a, x2b, x3b)) {
        if (!this->isPointInGbObject3D(x1a, x2a, x3a) || !this->isPointInGbObject3D(x1b, x2a, x3a) ||
            !this->isPointInGbObject3D(x1b, x2b, x3a) || !this->isPointInGbObject3D(x1a, x2b, x3a) ||
            !this->isPointInGbObject3D(x1a, x2a, x3b) || !this->isPointInGbObject3D(x1b, x2a, x3b) ||
            !this->isPointInGbObject3D(x1b, x2b, x3b) || !this->isPointInGbObject3D(x1a, x2b, x3b))
            return true;
    }
    return false;
}
/*======================================================================*/
bool GbObject3D::isCellInsideOrCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a,
                                                 const double &x1b, const double &x2b, const double &x3b)
{
    if (this->isPointInGbObject3D(x1a, x2a, x3a) || this->isPointInGbObject3D(x1b, x2a, x3a) ||
        this->isPointInGbObject3D(x1b, x2b, x3a) || this->isPointInGbObject3D(x1a, x2b, x3a) ||
        this->isPointInGbObject3D(x1a, x2a, x3b) || this->isPointInGbObject3D(x1b, x2a, x3b) ||
        this->isPointInGbObject3D(x1b, x2b, x3b) || this->isPointInGbObject3D(x1a, x2b, x3b)) {
        return true;
    }

    return false;
}
/*=======================================================*/
bool GbObject3D::isInsideCell(const double &minX1, const double &minX2, const double &minX3, const double &maxX1,
                              const double &maxX2, const double &maxX3)
{
    if (UbMath::greaterEqual(this->getX1Minimum(), minX1) && UbMath::greaterEqual(this->getX2Minimum(), minX2) &&
        UbMath::greaterEqual(this->getX3Minimum(), minX3) && UbMath::lessEqual(this->getX1Maximum(), maxX1) &&
        UbMath::lessEqual(this->getX2Maximum(), maxX2) && UbMath::lessEqual(this->getX2Maximum(), maxX3))
        return true;

    return false;
}

//! \}
