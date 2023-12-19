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
#ifndef GBPOINT3D_H
#define GBPOINT3D_H

#include <cmath>
#include <sstream>
#include <string>

#include <GbObject3D.h>

#include <PointerDefinitions.h>

class GbTriangle3D;

//! \brief This Class provides basic 3D point objects.
class GbPoint3D : public GbObject3D
{
public:
    GbPoint3D();
    GbPoint3D(const double &x1, const double &x2, const double &x3);
    GbPoint3D(GbPoint3D *point);
    ~GbPoint3D() override = default;

    GbPoint3D *clone() override { return new GbPoint3D(this); }
    void finalize() override {}

    void setCoordinates(const double &x1, const double &x2, const double &x3)
    {
        this->x1 = x1;
        this->x2 = x2;
        this->x3 = x3;
        this->notifyObserversObjectChanged();
    }
    void setX1(const double &x1)
    {
        this->x1 = x1;
        this->notifyObserversObjectChanged();
    }
    void setX2(const double &x2)
    {
        this->x2 = x2;
        this->notifyObserversObjectChanged();
    }
    void setX3(const double &x3)
    {
        this->x3 = x3;
        this->notifyObserversObjectChanged();
    }

    double getX1Coordinate() const { return this->x1; }
    double getX2Coordinate() const { return this->x2; }
    double getX3Coordinate() const { return this->x3; }

    void transform(const double matrix[4][4]);

    double getX1Centroid() override { return this->x1; }
    double getX1Minimum() override { return this->x1; }
    double getX1Maximum() override { return this->x1; }
    double getX2Centroid() override { return this->x2; }
    double getX2Minimum() override { return this->x2; }
    double getX2Maximum() override { return this->x2; }
    double getX3Centroid() override { return this->x3; }
    double getX3Minimum() override { return this->x3; }
    double getX3Maximum() override { return this->x3; }

    void translate(const double &x1, const double &x2, const double &x3) override;
    void rotate(const double &rx1, const double &rx2, const double &rx3) override;
    void scale(const double &sx1, const double &sx2, const double &sx3) override;

    double getDistance(GbPoint3D *p);
    bool equals(const GbPoint3D *point) const;
    bool isPointInGbObject3D(const double &x1, const double &x2, const double &x3, bool &pointIsOnBoundary) override;
    bool isPointInGbObject3D(const double &x1, const double &x2, const double &x3) override;
    bool isCellInsideGbObject3D(const double & /*x11*/, const double & /*x21*/, const double & /*x31*/,
                                const double & /*x12*/, const double & /*x22*/, const double & /*x23*/) override
    {
        return false;
    }

    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override;
    GbLine3D *createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2) override;
    std::string toString() override;

    using GbObject3D::isPointInGbObject3D; // Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht
                                           // ausprogrammieren , welche sonst hier "ueberdeckt" waere,da es dieselbe
                                           //methode mit anderen args gibt!

    // member
    double x1;
    double x2;
    double x3;
};

#endif

//! \}
