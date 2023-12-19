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
#ifndef KDUTILIES_H
#define KDUTILIES_H

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbTuple.h>

#include <geometry3d/GbTriFaceMesh3D.h>

#include <algorithm>
#include <vector>

namespace Kd
{
struct Axis {
    static const int X; // = 0;
    static const int Y; // = 1;
    static const int Z; // = 2;
};
/* ======================================================================================= */
struct Intersection {
    static const int ON_BOUNDARY;     // = -2;
    static const int INTERSECT_EDGE;  // = -1;
    static const int INTERSECTION;    // = 1;
    static const int NO_INTERSECTION; // = 0;
};
/* ======================================================================================= */
template <typename T>
inline void project2Axis(GbTriFaceMesh3D::TriFace &triFace, std::vector<GbTriFaceMesh3D::Vertex> &nodes,
                         const int &axis, std::vector<T> &projection)
{
    projection.resize(3);

    if (axis == Axis::X) {
        projection[0] = triFace.getV1x(nodes);
        projection[1] = triFace.getV2x(nodes);
        projection[2] = triFace.getV3x(nodes);
    } else if (axis == Axis::Y) {
        projection[0] = triFace.getV1y(nodes);
        projection[1] = triFace.getV2y(nodes);
        projection[2] = triFace.getV3y(nodes);
    } else if (axis == Axis::Z) {
        projection[0] = triFace.getV1z(nodes);
        projection[1] = triFace.getV2z(nodes);
        projection[2] = triFace.getV3z(nodes);
    } else
        throw UbException(UB_EXARGS, "unknown axis");

    std::sort(projection.begin(), projection.end(), std::less<>());
}
/* ======================================================================================= */
template <typename T>
inline bool isPointOnPlane(const T &px, const T &py, const T &pz, const T &precision,
                           GbTriFaceMesh3D::Vertex &pointOfTriFace, GbTriFaceMesh3D::TriFace &triFace)
{
    return std::fabs((px - pointOfTriFace.x) * triFace.nx + (py - pointOfTriFace.y) * triFace.ny +
                     (pz - pointOfTriFace.z) * triFace.nz) < precision;
}
/* ======================================================================================= */
template <typename T>
inline bool isPointOnTriangle(const T &px, const T &py, const T &pz, const T &precision, GbTriFaceMesh3D::Vertex &p1,
                              GbTriFaceMesh3D::Vertex &p2, GbTriFaceMesh3D::Vertex &p3,
                              GbTriFaceMesh3D::TriFace &triFace)
{
    if (Kd::isPointOnPlane(px, py, pz, precision, p1, triFace)) {
        T a_x = p1.x - px;
        T a_y = p1.y - py;
        T a_z = p1.z - pz;
        T b_x = p2.x - px;
        T b_y = p2.y - py;
        T b_z = p2.z - pz;
        T c_x = p3.x - px;
        T c_y = p3.y - py;
        T c_z = p3.z - pz;

        const T factor = 0.5;
        T Q1_x         = (a_y * b_z - a_z * b_y) * factor;
        T Q1_y         = (a_z * b_x - a_x * b_z) * factor;
        T Q1_z         = (a_x * b_y - a_y * b_x) * factor;

        T Q2_x = (b_y * c_z - b_z * c_y) * factor;
        T Q2_y = (b_z * c_x - b_x * c_z) * factor;
        T Q2_z = (b_x * c_y - b_y * c_x) * factor;

        T Q3_x = (c_y * a_z - c_z * a_y) * factor;
        T Q3_y = (c_z * a_x - c_x * a_z) * factor;
        T Q3_z = (c_x * a_y - c_y * a_x) * factor;

        T Q_x = Q1_x + Q2_x + Q3_x;
        T Q_y = Q1_y + Q2_y + Q3_y;
        T Q_z = Q1_z + Q2_z + Q3_z;

        if (UbMath::zero(Q_x * Q1_x + Q_y * Q1_y + Q_z * Q1_z))
            return true;
        else if (UbMath::zero(Q_x * Q2_x + Q_y * Q2_y + Q_z * Q2_z))
            return true;
        else if (UbMath::zero(Q_x * Q3_x + Q_y * Q3_y + Q_z * Q3_z))
            return true;
        else if (UbMath::less(Q_x * Q1_x + Q_y * Q1_y + Q_z * Q1_z, T(0.0)))
            return false;
        else if (UbMath::less(Q_x * Q2_x + Q_y * Q2_y + Q_z * Q2_z, T(0.0)))
            return false;
        else if (UbMath::less(Q_x * Q3_x + Q_y * Q3_y + Q_z * Q3_z, T(0.0)))
            return false;

        return true;
    }

    return false;
}
/* ======================================================================================= */
template <typename T>
inline bool intersectLine(const UbTuple<T, T, T> &n1, const UbTuple<T, T, T> &n2, GbTriFaceMesh3D::TriFace &triFace,
                          std::vector<GbTriFaceMesh3D::Vertex> &nodes)
{
    GbTriFaceMesh3D::Vertex &p0 = triFace.getNode(0, nodes);

    const T &n1X = val<1>(n1);
    const T &n1Y = val<2>(n1);
    const T &n1Z = val<3>(n1);

    const T &n2X = val<1>(n2);
    const T &n2Y = val<2>(n2);
    const T &n2Z = val<3>(n2);

    // if(   Kd::isPointOnPlane(n1X, n1Y, n1Z, T(1.0E-6), p0, triFace)
    //   && Kd::isPointOnPlane(n2X, n2Y, n2Z, T(1.0E-6), p0, triFace))
    //{
    //   return true;
    //}

    T denom = (n2X - n1X) * triFace.nx + (n2Y - n1Y) * triFace.ny + (n2Z - n1Z) * triFace.nz;

    if (UbMath::zero(denom)) // line does not intersect the plane of the triangle !
    {
        return false;
    } else {
        T d  = -triFace.nx * p0.x - triFace.ny * p0.y - triFace.nz * p0.z;
        T mu = T(-1.0 * (d + n1X * triFace.nx + n1Y * triFace.ny + n1Z * triFace.nz) / denom);

        if (!UbMath::inClosedInterval(mu, T(0.0),
                                      T(1.0))) // Point of intersection of line and plane does not lie on the triangle
        {
            return false;
        } else {
            // intersection with plane

            // Test whether Point lies inside the triangle or not
            GbTriFaceMesh3D::Vertex &p1 = triFace.getNode(1, nodes);
            GbTriFaceMesh3D::Vertex &p2 = triFace.getNode(2, nodes);

            return Kd::isPointOnTriangle(n1X + ((n2X - n1X) * mu) // intersectionPointX
                                         ,
                                         n1Y + ((n2Y - n1Y) * mu) // intersectionPointY
                                         ,
                                         n1Z + ((n2Z - n1Z) * mu) // intersectionPointZ
                                         ,
                                         T(0.001), p0, p1, p2, triFace);
        }
    }
}
} // namespace Kd

#endif // KDUTILIES_H

//! \}
