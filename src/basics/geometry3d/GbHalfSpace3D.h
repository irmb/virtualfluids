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
//! \author Konstantin Kutscher, Soeren Textor, Sebastian Geller
//=======================================================================================
#ifndef GBHALFSPACE3D_H
#define GBHALFSPACE3D_H

#include <iostream>
#include <sstream>

#include <basics/utilities/UbMath.h>

#include <geometry3d/GbPoint3D.h>
#include <geometry3d/GbTriangle3D.h>
#include <geometry3d/GbVector3D.h>

#include <PointerDefinitions.h>

/*=========================================================================*/
/* GbHalfSpace3D                                                             */
/*                                                                         */
/**
 * This Class helps in performing some operations on a halfspace defined by 2 or 3 points
 */

class GbHalfSpace3D
{
public:
    GbHalfSpace3D(GbTriangle3D *triangle);

    GbHalfSpace3D(GbPoint3D *PointA, GbPoint3D *PointB, GbPoint3D *PointC);

    GbHalfSpace3D(GbPoint3D *PointA, GbPoint3D *PointB);

    GbHalfSpace3D(const double &p1x, const double &p1y, const double &p1z, const double &p2x, const double &p2y,
                  const double &p2z, const double &p3x, const double &p3y, const double &p3z);
    GbHalfSpace3D(const double &p1x, const double &p1y, const double &p1z, const double &nx, const double &ny,
                  const double &nz);

    /*=======================================================*/
    std::string getTypeID() { return "GbHalfSpace3D"; }
    /*=============================================*/
    bool ptInside(const double &x, const double &y, const double &z)
    {
        return UbMath::greaterEqual(normalX * x + normalY * y + normalZ * z, this->d);
    }
    /*=============================================*/
    bool ptInside(GbPoint3D *pointX)
    {
        // GbVector3D X(PointX->x1, PointX->x2, PointX->x3 );
        // return UbMath::greaterEqual(this->Normal.Dot(X), this->d);
        return UbMath::greaterEqual(normalX * pointX->x1 + normalY * pointX->x2 + normalZ * pointX->x3, this->d);
    }
    /*=============================================*/
    bool ptInside(GbVector3D &x)
    {
        // return UbMath::greaterEqual(this->Normal.Dot(X), this->d);
        return UbMath::greaterEqual(normalX * x[0] + normalY * x[1] + normalZ * x[2], this->d);
    }
    /*=============================================*/
    double getDistance(const double &x1p, const double &x2p, const double &x3p)
    {
        return (normalX * x1p + normalY * x2p + normalZ * x3p) - this->d;
    }

    const double &getNormalX() { return this->normalX; }
    const double &getNormalY() { return this->normalY; }
    const double &getNormalZ() { return this->normalZ; }
    const double &getD() { return this->d; }

private:
    // GbVector3D Normal;
    double normalX;
    double normalY;
    double normalZ;
    double d;
};
/*=========================================================================*/

#endif // GBHALFSPACE3D_H

//! \}
