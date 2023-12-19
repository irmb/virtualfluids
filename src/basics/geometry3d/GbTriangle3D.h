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
#ifndef GBTRIANGLE3D_H
#define GBTRIANGLE3D_H

#include <sstream>

#include <GbObject3D.h>
#include <GbPoint3D.h>
#include <GbVector3D.h>

#include <PointerDefinitions.h>

class GbCuboid3D;
class GbPolygon3D;
class GbObject3DCreator;

//////////////////////////////////////////////////////////////////////////
//!
//! \class GbTriangle3D
//!
//! \brief This Class provides basic 3D triangle objects.
//! \details The describing points are observed by 2D triangle objects.
//!
//////////////////////////////////////////////////////////////////////////

class GbTriangle3D : public GbObject3D, public UbObserver
{
public:
    /*======================================================================*/
    /*  Konstruktoren                                                       */
    /*                                                                      */
    GbTriangle3D();
    GbTriangle3D(GbPoint3D *point1, GbPoint3D *point2, GbPoint3D *point3);
    GbTriangle3D(GbTriangle3D *triangle);
    ~GbTriangle3D() override;
    /*======================================================================*/
    /*  Methoden                                                            */
    /*                                                                      */
    GbTriangle3D *clone() override;
    void finalize() override { this->deletePoints(); }

    GbPoint3D *getPoint1() { return this->points[0]; }
    GbPoint3D *getPoint2() { return this->points[1]; }
    GbPoint3D *getPoint3() { return this->points[2]; }

    GbVector3D getNormal();
    void calculateNormal();

    void deletePoints();

    int contains(GbPoint3D *point);
    int containsEqual(GbPoint3D *point);
    GbPoint3D *getPoint(const int &index);
    std::vector<GbPoint3D> getPoints();
    double getArea();
    double getX1Centroid() override;
    double getX1Minimum() override;
    double getX1Maximum() override;
    double getX2Centroid() override;
    double getX2Minimum() override;
    double getX2Maximum() override;
    double getX3Centroid() override;
    double getX3Minimum() override;
    double getX3Maximum() override;

    void setInconsistent() { this->consistent = false; }

    void setPoint(GbPoint3D *point, int index);

    // bool equals(GbObject3D *object)
    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override;
    bool isPointInGbObject3D(const double & /*x1*/, const double & /*x2*/, const double & /*x3*/) override
    {
        // der einfachheit halber ...
        return false;
        // throw UbException(__FILE__, __LINE__, "GbTriangle3D::isPointInObject3D- not implemented");
    }
    bool isPointInGbObject3D(const double & /*x1*/, const double & /*x2*/, const double & /*x3*/,
                             bool &pointIsOnBoundary) override
    {
        // der einfachheit halber ...
        pointIsOnBoundary = false;
        return false;
        // throw UbException(__FILE__, __LINE__, "GbTriangle3D::isPointInObject3D- not implemented");
    }
    bool isCellInsideGbObject3D(const double & /*x11*/, const double & /*x21*/, const double & /*x31*/,
                                const double & /*x12*/, const double & /*x22*/, const double & /*x23*/) override
    {
        return false;
    }

    // get distance from a point to the triangle
    // todo CHANGE...
    double getDistanceFromPoint(GbVector3D punct);

    std::string toString() override;

    /*======================================================================*/
    /*  Calculation                                                         */
    /*                                                                      */
    //   std::vector<GbPoint3D> calculateIntersectionPoints3D(GbLine3D *line);
    bool hasRaytracing() override { return true; }
    /*|r| must be 1! einheitsvector!!*/
    double getIntersectionRaytraceFactor(const double &x1, const double &x2, const double &x3, const double &rx1,
                                         const double &rx2, const double &rx3) override;
    //   bool isPointOnEdge(GbVector3D& q);

    GbPoint3D *calculateIntersectionPoints3D(GbLine3D *line);
    GbPoint3D *calculateIntersectionPoints3D(GbPoint3D *linePoint1, GbPoint3D *linePoint2);
    double calculateDistanceToPoint3D(GbPoint3D *point);
    double calculateDistanceToPoint3D(const double &x1, const double &x2, const double &x3);
    double calculateNormalizedDistanceToPoint3D(const double &x1, const double &y1, const double &z1, const double &x2,
                                                const double &y2, const double &z2);

    bool enclosesPoint2D(double x1, double x2);
    GbPolygon3D *createClippedPolygon3D(GbCuboid3D *cube);
    GbLine3D *createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2) override;
    // public GbPolygon2D createClippedPolygon2D(GbPoint2D p1, GbPoint2D p2);
    GbPolygon3D *createClippedPolygon3D(const double &p1x1, const double &p1x2, const double &p1x3, const double &p2x1,
                                        const double &p2x2, const double &p2x3);
    // bool enclosesRectangle2D(GbRectangle2D *rectangle);
    // public boolean enclosesRectangle2D(GbPoint2D p1, GbPoint2D p2);
    // public boolean enclosesRectangle2D(double p1x1, double p1x2, double p2x1, double p2x2);
    // public boolean crossesRectangle2D(GbRectangle2D rectangle);
    // public boolean crossesRectangle2D(GbPoint2D p1, GbPoint2D p2);
    // public boolean crossesRectangle2D(double p1x1, double p1x2, double p2x1, double p2x2);
    // public boolean enclosesOrCrossesRectangle2D(GbRectangle2D rectangle);
    // public boolean enclosesOrCrossesRectangle2D(GbPoint2D p1, GbPoint2D p2);
    // public boolean enclosesOrCrossesRectangle2D(double p1x1, double p1x2, double p2x1, double p2x2);
    /*======================================================================*/
    /*======================================================================*/
    /*  Private Methoden                                                    */
    /*                                                                      */
    virtual void calculateValues();

    // virtuelle Methoden von UbObserver
    //!! quick and dirty von sirann !!
    void objectChanged(UbObservable *changedObject) override
    {
        GbPoint3D *point = dynamic_cast<GbPoint3D *>(changedObject);
        if (!point || (this->points[0] != point && this->points[1] != point && this->points[2] != point))
            return;

        this->consistent = false;
    }
    void objectWillBeDeleted(UbObservable *objectForDeletion) override
    {
        if (this->points[0]) {
            UbObservable *observedObj = dynamic_cast<UbObservable *>(this->points[0]);
            if (objectForDeletion == observedObj) {
                this->points[0] = NULL;
            }
        }
        if (this->points[1]) {
            UbObservable *observedObj = dynamic_cast<UbObservable *>(this->points[1]);
            if (objectForDeletion == observedObj) {
                this->points[1] = NULL;
            }
        }
        if (this->points[2]) {
            UbObservable *observedObj = dynamic_cast<UbObservable *>(this->points[2]);
            if (objectForDeletion == observedObj) {
                this->points[2] = NULL;
            }
        }
        // ACHTUNG: eigentlich muessten in allen methoden von GbLine if abfragen fuer NULL pointer hin... toDo
    }
    using GbObject3D::isPointInGbObject3D; // Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht
                                           // ausprogrammieren, welche sonst hier "ueberdeckt" waere

protected:
    bool consistent;
    double x1s;
    double x2s;
    double x3s;
    double x1min;
    double x1max;
    double x2min;
    double x2max;
    double x3min;
    double x3max;
    double area;

    GbVector3D normal;
    std::vector<GbPoint3D *> points;

private:
    void init();
};
/*=========================================================================*/

#endif

//! \}
