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
#ifndef GBVECTOR3D_H
#define GBVECTOR3D_H

#include <cassert>
#include <cfloat>
#include <string>

#include <PointerDefinitions.h>

class GbPoint3D;

//! \brief This Class provides basic 3D vector objects.
class GbVector3D
{
public:
    // construction
    GbVector3D();
    GbVector3D(const double &fX1, const double &fX2, const double &fX3);
    GbVector3D(const GbVector3D &rkV);
    GbVector3D(const GbPoint3D &rkV);

    std::string toString();

    // coordinate access
    operator const double *() const;
    operator double *();
    double operator[](int i) const;
    double &operator[](int i);
    double X1() const;
    double &X1();
    double X2() const;
    double &X2();
    double X3() const;
    double &X3();

    // assignment
    GbVector3D &operator=(const GbVector3D &rkV);

    // comparison
    bool operator==(const GbVector3D &rkV) const;
    bool operator!=(const GbVector3D &rkV) const;
    bool operator<(const GbVector3D &rkV) const;
    bool operator<=(const GbVector3D &rkV) const;
    bool operator>(const GbVector3D &rkV) const;
    bool operator>=(const GbVector3D &rkV) const;

    // arithmetic operations
    GbVector3D operator+(const GbVector3D &rkV) const;
    GbVector3D operator-(const GbVector3D &rkV) const;
    GbVector3D operator*(const double &fScalar) const;
    GbVector3D operator/(const double &fScalar) const;
    GbVector3D operator-() const;

    // arithmetic updates
    GbVector3D &operator+=(const GbVector3D &rkV);
    GbVector3D &operator-=(const GbVector3D &rkV);
    GbVector3D &operator*=(const double &fScalar);
    GbVector3D &operator/=(const double &fScalar);

    GbVector3D Add(GbVector3D &vector);
    GbVector3D Subtract(GbVector3D &vector);
    GbVector3D Scale(const double &x);

    // vector operations
    double Length() const;
    double SquaredLength() const;
    double Dot(const GbVector3D &rkV) const;
    double Normalize();

    // The cross products are computed using the right-handed rule.  Be aware
    // that some graphics APIs use a left-handed rule.  If you have to compute
    // a cross product with these functions and send the result to the API
    // that expects left-handed, you will need to change sign on the vector
    // (replace each component value c by -c).
    GbVector3D Cross(const GbVector3D &rkV) const;
    GbVector3D UnitCross(const GbVector3D &rkV) const;

    // Compute the barycentric coordinates of the point with respect to the
    // tetrahedron <V0,V1,V2,V3>, P = b0*V0 + b1*V1 + b2*V2 + b3*V3, where
    // b0 + b1 + b2 + b3 = 1.
    void GetBarycentrics(const GbVector3D &rkV0, const GbVector3D &rkV1, const GbVector3D &rkV2, const GbVector3D &rkV3,
                         double afBary[4]) const;

    // Gram-Schmidt orthonormalization.  Take linearly independent vectors
    // U, V, and W and compute an orthonormal set (unit length, mutually
    // perpendicular).
    static void Orthonormalize(GbVector3D &rkU, GbVector3D &rkV, GbVector3D &rkW);
    static void Orthonormalize(GbVector3D *akV);

    // Input W must be initialized to a nonzero vector, output is {U,V,W},
    // an orthonormal basis.  A hint is provided about whether or not W
    // is already unit length.
    static void GenerateOrthonormalBasis(GbVector3D &rkU, GbVector3D &rkV, GbVector3D &rkW, bool bUnitLengthW);

    // special vectors
    static const GbVector3D ZERO;
    static const GbVector3D UNIT_X1;
    static const GbVector3D UNIT_X2;
    static const GbVector3D UNIT_X3;

private:
    // support for comparisons
    int CompareArrays(const GbVector3D &rkV) const;

    double m_afTuple[3];
};

GbVector3D operator*(const double &fScalar, const GbVector3D &rkV);

#endif // GBVECTOR3D_H

//! \}
