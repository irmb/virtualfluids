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
#ifndef GBCYLINDER3D_H
#define GBCYLINDER3D_H

#include <cmath>
#include <vector>

#include <basics/utilities/UbObserver.h>
#include <geometry3d/GbLine3D.h>
#include <geometry3d/GbObject3D.h>

class GbPoint3D;
class GbLine3D;
class GbTriangle3D;

class GbObject3DCreator;

#include <PointerDefinitions.h>
class GbCylinder3D;
using GbCylinder3DPtr = SPtr<GbCylinder3D>;

class GbCylinder3D : public GbObject3D, public UbObserver
{
public:
    GbCylinder3D();
    GbCylinder3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b, const double &x2b,
                 const double &x3b, const double &radius);
    GbCylinder3D(GbPoint3D *p1, GbPoint3D *p2, const double &radius);
    GbCylinder3D(GbLine3D *line, const double &rad);
    GbCylinder3D(GbCylinder3D *cylinder);
    ~GbCylinder3D() override;

    GbCylinder3D *clone() override { return new GbCylinder3D(this); }
    void finalize() override;

    double getRadius() { return this->mRad; };
    GbLine3D *getLine() { return mLine; }
    GbPoint3D *getPoint1();
    GbPoint3D *getPoint2();

    void setRadius(const double &radius);
    void setLine(GbLine3D *line);
    void setPoint1(const double &x1, const double &x2, const double &x3);
    void setPoint2(const double &x1, const double &x2, const double &x3);

    bool isParallelToX1Axis() { return ((this->cylinderType & X1PARALLEL) == X1PARALLEL); }
    bool isParallelToX2Axis() { return ((this->cylinderType & X2PARALLEL) == X2PARALLEL); }
    bool isParallelToX3Axis() { return ((this->cylinderType & X3PARALLEL) == X3PARALLEL); }
    bool isNotParallelToAxis() { return ((this->cylinderType & NOTPARALLELTOAXIS) == NOTPARALLELTOAXIS); }

    double getHeight();

    void scale(const double &sx1, const double &sx2, const double &sx3) override;

    void translate(const double &x1, const double &x2, const double &x3) override
    {
        this->mLine->translate(x1, x2, x3);
        this->calculateValues();
        // this->notifyObserversObjectChanged();
    }

    double getX1Centroid() override { return centerX1; }
    double getX1Minimum() override { return minX1; }
    double getX1Maximum() override { return maxX1; }
    double getX2Centroid() override { return centerX2; }
    double getX2Minimum() override { return minX2; }
    double getX2Maximum() override { return maxX2; }
    double getX3Centroid() override { return centerX3; }
    double getX3Minimum() override { return minX3; }
    double getX3Maximum() override { return maxX3; }

    bool isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p) override;
    bool isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p, bool &pointIsOnBoundary) override;
    bool isCellInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                const double &x2b, const double &x3b) override;
    bool isCellCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                 const double &x2b, const double &x3b) override;
    bool isCellInsideOrCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b) override;

    GbLine3D *createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2) override;

    // SG ausdokumentieren, da der nur unendlcihe Zylinder macht ...
    // bool hasRaytracing() { return true; }
    bool hasRaytracing() override { return false; }
    bool raytracingSupportsPointsInside() override { return true; }

    /*|r| must be 1! einheitsvector!!*/
    double getIntersectionRaytraceFactor(const double &x1, const double &x2, const double &x3, const double &rx1,
                                         const double &rx2, const double &rx3) override;

    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override;
    void addSurfaceTriangleSet(std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt3> &triangles) override;
    void addSurfaceTriangleSetSegments(std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt3> &triangles,
                                       int segmentsRound, int segmentsHeight);

    std::string toString() override;

    // virtuelle Methoden von UbObserver
    void objectChanged(UbObservable *changedObject) override;
    void objectWillBeDeleted(UbObservable *objectForDeletion) override;

    using GbObject3D::isPointInGbObject3D; // Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht
                                           // ausprogrammieren, welche sonst hier "ueberdeckt" waere

protected:
    void calculateValues();

    GbLine3D *mLine;
    double mRad;

    double minX1{ 0.0 }, minX2{ 0.0 }, minX3{ 0.0 };
    double maxX1{ 0.0 }, maxX2{ 0.0 }, maxX3{ 0.0 };
    double centerX1{ 0.0 }, centerX2{ 0.0 }, centerX3{ 0.0 };

    int cylinderType;

    // void berechneQuerschnittsWerte();
    static const int NOTPARALLELTOAXIS = (1 << 0); // 1
    static const int X1PARALLEL        = (1 << 1); // 2
    static const int X2PARALLEL        = (1 << 2); // 4
    static const int X3PARALLEL        = (1 << 3); // 8
};

#endif

//! \}
