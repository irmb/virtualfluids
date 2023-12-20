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
#ifndef GBCUBOID3D_H
#define GBCUBOID3D_H

#include <cmath>
#include <vector>

#include <GbPoint3D.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbObserver.h>

class GbLine3D;
class GbObject3DCreator;

#include <PointerDefinitions.h>
class GbCuboid3D;
using GbCuboid3DPtr = SPtr<GbCuboid3D>;

//! \brief This Class provides basic 3D box objects.
class GbCuboid3D : public GbObject3D, public UbObserver
{
public:
    GbCuboid3D();
    GbCuboid3D(const double &minX1, const double &minX2, const double &minX3, const double &maxX1, const double &maxX2,
               const double &maxX3);
    GbCuboid3D(GbPoint3D *p1, GbPoint3D *p2);
    GbCuboid3D(GbCuboid3D *cuboid);
    ~GbCuboid3D() override;

    GbCuboid3D *clone() override { return new GbCuboid3D(this); }
    void finalize() override;

    GbPoint3D *getPoint1() { return this->p1; }
    GbPoint3D *getPoint2() { return this->p2; }

    void setPoint1(GbPoint3D *point1);
    void setPoint2(GbPoint3D *point2);
    void setPoints(GbPoint3D *point1, GbPoint3D *point2);

    double getX1Centroid() override;
    double getX1Minimum() override;
    double getX1Maximum() override;
    double getX2Centroid() override;
    double getX2Minimum() override;
    double getX2Maximum() override;
    double getX3Centroid() override;
    double getX3Minimum() override;
    double getX3Maximum() override;
    void setCenterCoordinates(const UbTupleDouble3 & /*position*/) override
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    void setCenterCoordinates(const double &x1, const double &x2, const double &x3) override;

    void translate(const double &x1, const double &x2, const double &x3) override;
    void rotate(const double &rx1, const double &rx2, const double &rx3) override { (void)rx1; (void)rx2; (void)rx3; }
    void scale(const double &sx1, const double &sx2, const double &sx3) override;

    double getLengthX1();
    double getLengthX2();
    double getLengthX3();

    bool isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p, bool &pointIsOnBoundary) override;
    bool isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p) override;
    bool isCellInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                const double &x2b, const double &x3b) override;
    bool isCellCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                 const double &x2b, const double &x3b) override;
    bool isCellInsideOrCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b) override;
    double getCellVolumeInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b) override;

    GbPoint3D *calculateInterSectionPoint3D(GbPoint3D &point1, GbPoint3D &point2);
    // GbCuboid3D* createClippedRectangle3D(const double& x1a,const double& x2a,const double& x3a,const double&
    // x1b,const double& x2b,const double& x3b);
    GbLine3D *createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2) override;

    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override;
    void addSurfaceTriangleSet(std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt3> &triangles) override;

    bool hasRaytracing() override { return true; }
    /*|r| must be 1! einheitsvector!!*/
    double getIntersectionRaytraceFactor(const double &x1, const double &x2, const double &x3, const double &rx1,
                                         const double &rx2, const double &rx3) override;

    double getDistance(GbPoint3D *p)
    {
        return this->getDistance(p->getX1Coordinate(), p->getX2Coordinate(), p->getX3Coordinate());
    }
    double getDistance(const double &x1p, const double &x2p, const double &x3p)
    {
        throw UbException(UB_EXARGS, "not implemented");

        // falls punkt innerhalt ist: minimalen abstand ausrechnen
        if (this->isPointInGbObject3D(x1p, x2p, x3p)) {
            double x1Dist = UbMath::min(std::abs(x1p - this->getX1Minimum()), std::abs(x1p - this->getX1Maximum()));
            double x2Dist = UbMath::min(std::abs(x2p - this->getX2Minimum()), std::abs(x2p - this->getX2Maximum()));
            double x3Dist = UbMath::min(std::abs(x3p - this->getX3Minimum()), std::abs(x3p - this->getX3Maximum()));

            return UbMath::min(x1Dist, x2Dist, x3Dist);
        } else {
        }
    }

    std::string toString() override;

    // virtuelle Methoden von UbObserver
    void objectChanged(UbObservable *changedObject) override;
    void objectWillBeDeleted(UbObservable *objectForDeletion) override;

    using GbObject3D::isPointInGbObject3D; // Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht
                                           // ausprogrammieren, welche sonst hier "ueberdeckt" waere

protected:
    GbPoint3D *p1;
    GbPoint3D *p2;
};

#endif

//! \}
