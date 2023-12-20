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
#ifndef GBLINE3D_H
#define GBLINE3D_H

#include <cmath>
#include <sstream>

#include <basics/utilities/UbObserver.h>

#include <GbObject3D.h>
#include <GbPoint3D.h>

class GbCuboid3D;

#include <PointerDefinitions.h>

//////////////////////////////////////////////////////////////////////////
//!
//!  \class GbLine3D
//!
//! \brief This Class provides basic 3D line objects.
//! \details The describing points are observed by 3D line objects.
//!
//////////////////////////////////////////////////////////////////////////

class GbLine3D : public GbObject3D, public UbObserver
{
public:
    GbLine3D();
    GbLine3D(GbPoint3D *point1, GbPoint3D *point2);
    GbLine3D(GbLine3D *line);
    ~GbLine3D() override;

    GbLine3D *clone() override { return new GbLine3D(this); }
    void finalize() override;

    void setPoint1(GbPoint3D *point1);
    void setPoint2(GbPoint3D *point2);
    void setPoints(GbPoint3D *point1, GbPoint3D *point2);

    void deletePoint1()
    {
        if (this->p1) {
            this->p1->removeObserver(this);
            delete this->p1;
            this->p1 = NULL;
        }
    }
    void deletePoint2()
    {
        if (this->p2) {
            this->p2->removeObserver(this);
            delete this->p2;
            this->p2 = NULL;
        }
    }
    void deletePoints()
    {
        this->deletePoint1();
        this->deletePoint2();
    }

    GbPoint3D *getPoint1() { return this->p1; }
    GbPoint3D *getPoint2() { return this->p2; }

    double getLength() { return (this->length); }

    double getX1Centroid() override { return ((this->p1->x1 + this->p2->x1) * 0.5); }
    double getX2Centroid() override { return ((this->p1->x2 + this->p2->x2) * 0.5); };
    double getX3Centroid() override { return ((this->p1->x3 + this->p2->x3) * 0.5); }

    double getX1Minimum() override { return (this->p1->x1 < this->p2->x1 ? this->p1->x1 : this->p2->x1); }
    double getX2Minimum() override { return (this->p1->x2 < this->p2->x2 ? this->p1->x2 : this->p2->x2); }
    double getX3Minimum() override { return (this->p1->x3 < this->p2->x3 ? this->p1->x3 : this->p2->x3); }

    double getX1Maximum() override { return (this->p1->x1 > this->p2->x1 ? this->p1->x1 : this->p2->x1); }
    double getX2Maximum() override { return (this->p1->x2 > this->p2->x2 ? this->p1->x2 : this->p2->x2); }
    double getX3Maximum() override { return (this->p1->x3 > this->p2->x3 ? this->p1->x3 : this->p2->x3); }

    void scale(const double &sx1, const double &sx2, const double &sx3) override;
    void translate(const double &tx1, const double &tx2, const double &tx3) override;

    GbPoint3D *calculateIntersectionPoint3D(GbLine3D *line);
    GbLine3D *createClippedLine3D(GbCuboid3D *cuboid);
    GbLine3D *createClippedLine3D(GbPoint3D *pA, GbPoint3D *pE);

    double getDistance(const GbPoint3D &point);
    double getDistance(const double &x1, const double &x2, const double &x3);

    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override;
    bool isPointInGbObject3D(const double & /*x1*/, const double & /*x2*/, const double & /*x3*/) override
    {
        throw UbException(UB_EXARGS, "not implemented");
    }
    bool isPointInGbObject3D(const double & /*x1*/, const double & /*x2*/, const double & /*x3*/,
                             bool & /*pointIsOnBoundary*/) override
    {
        throw UbException(UB_EXARGS, "not implemented");
    }
    bool isCellInsideGbObject3D(const double & /*x11*/, const double & /*x21*/, const double & /*x31*/,
                                const double & /*x12*/, const double & /*x22*/, const double & /*x32*/) override
    {
        return false;
    }

    GbLine3D *createClippedLine3D(GbPoint3D & /*point1*/, GbPoint3D & /*point2*/) override
    {
        throw UbException(UB_EXARGS, "not implemented");
    }

    // virtuelle Methoden von UbObserver
    void objectChanged(UbObservable *changedObject) override;
    void objectWillBeDeleted(UbObservable *objectForDeletion) override;

    std::string toString() override;

    using GbObject3D::isPointInGbObject3D; // Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht
                                           // ausprogrammieren, welche sonst hier "ueberdeckt" waere
protected:
    GbPoint3D *p1;
    GbPoint3D *p2;
    double length;

private:
    void calculateValues();
};

#endif

//! \}
