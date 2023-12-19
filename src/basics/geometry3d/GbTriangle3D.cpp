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
#include <GbCuboid3D.h>
#include <GbLine3D.h>
#include <GbSystem3D.h>
#include <GbTriangle3D.h>

#include <basics/utilities/UbMath.h>

using namespace std;

GbTriangle3D::GbTriangle3D()
{
    this->init();
    this->consistent = false;
}
/*======================================================================*/
/*
 * Creates an empty 2D triangle with the specified points.
 * @param point1 the 1st point
 * @param point2 the 2nd point
 * @param point3 the 3nd point
 */
GbTriangle3D::GbTriangle3D(GbPoint3D *point1, GbPoint3D *point2, GbPoint3D *point3)
{
    this->init();
    this->points[0] = point1;
    this->points[1] = point2;
    this->points[2] = point3;

    this->calculateNormal();
    this->consistent = false;

    this->points[0]->addObserver(this);
    this->points[1]->addObserver(this);
    this->points[2]->addObserver(this);

    // this.po        = new PointObserver(this);
    // this.points[0].addObserver(this.po);
    // this.points[1].addObserver(this.po);
    // this.points[2].addObserver(this.po);
}
/*======================================================================*/
/*
 * Creates a 3D triangle as clone of the specified 2D triangle.
 * @param triangle the 3D triangle to be cloned
 */
GbTriangle3D::GbTriangle3D(GbTriangle3D *triangle)
{
    this->init();
    this->points[0] = triangle->points[0]->clone();
    this->points[1] = triangle->points[1]->clone();
    this->points[2] = triangle->points[2]->clone();

    this->consistent = false;
    this->calculateNormal();
    this->calculateValues();
}
/*======================================================================*/
GbTriangle3D::~GbTriangle3D()
{
    if (this->points[0])
        this->points[0]->removeObserver(this);
    if (this->points[1])
        this->points[1]->removeObserver(this);
    if (this->points[2])
        this->points[2]->removeObserver(this);
}
/*======================================================================*/
void GbTriangle3D::deletePoints()
{
    if (points[0]) {
        delete points[0];
        points[0] = NULL;
    }
    if (points[1]) {
        delete points[1];
        points[1] = NULL;
    }
    if (points[2]) {
        delete points[2];
        points[2] = NULL;
    }
}

/*======================================================================*/
/*  Methoden                                                            */
/*                                                                      */
/*
 * Creates a 3D triangle as clone of this 3D triangle.
 */
GbTriangle3D *GbTriangle3D::clone() { return (new GbTriangle3D(this)); }
/*======================================================================*/
/*
 * Returns the number of times this 2D triangle contains the specified point.
 * @param point the point
 * @return the number of times this 2D triangle contains the specified point
 */
int GbTriangle3D::contains(GbPoint3D *point)
{
    int n = 0;
    for (int i = 0; i < 3; i++)
        if (this->points[i]->equals(point))
            n++;
    return (n);
}
/*======================================================================*/
/*
 * Returns the number of times this 2D triangle contains a point equal to the specified point.
 * @param point the point
 * @return the number of times this 2D triangle contains a point equal to the specified point
 */
int GbTriangle3D::containsEqual(GbPoint3D *point)
{
    int n = 0;
    for (int i = 0; i < 3; i++)
        if (this->points[i]->equals(point))
            n++;
    return (n);
}
/*======================================================================*/
/*
 * Returns the specified point.
 * @param index the index (must be 0, 1, or 2)
 * @return the specified point
 * @exception ArrayIndexOutOfBoundsException if the specified index is not valid
 */
GbPoint3D *GbTriangle3D::getPoint(const int &index)
{
    if (index < 0 || index > 2)
        throw UbException(UB_EXARGS, "invalid index specified: ");
    return ((this->points[index]));
}
/*======================================================================*/
vector<GbPoint3D> GbTriangle3D::getPoints()
{
    vector<GbPoint3D> p(3);
    p[0] = *(points[0]);
    p[1] = *(points[1]);
    p[2] = *(points[2]);
    return p;
    //
    // vector<GbPoint3D> p(3);// = new vector<GbPoint3D*>;
    // p.resize(3);//, NULL);
    // p[0] = this->points[0];
    // p[1] = this->points[1];
    // p[2] = this->points[2];
    // return(p);
}
/*======================================================================*/
/*
 * Returns the area of this triangle.
 * The area is positive for positive ordered points, otherwise negative.
 * @return the area of this triangle
 */
double GbTriangle3D::getArea()
{
    if (!this->consistent)
        this->calculateValues();
    // throw UbException(UB_EXARGS,"not correct calculated ...");
    return (this->area);
}
/*
 * Returns the centroid x1 coordinate of this triangle.
 * @return the centroid x1 coordinate of this triangle
 */
double GbTriangle3D::getX1Centroid()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x1s);
}
/*
 * Returns the minimum x1 coordinate of this triangle.
 * @return the minimum x1 coordinate of this triangle
 */
double GbTriangle3D::getX1Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x1min);
}
/*
 * Returns the maximum x1 coordinate of this triangle.
 * @return the maximum x1 coordinate of this triangle
 */
double GbTriangle3D::getX1Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x1max);
}
/*
 * Returns the centroid x2 coordinate of this triangle.
 * @return the centroid x2 coordinate of this triangle
 */
double GbTriangle3D::getX2Centroid()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x2s);
}
/*
 * Returns the minimum x2 coordinate of this triangle.
 * @return the minimum x2 coordinate of this triangle
 */
double GbTriangle3D::getX2Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x2min);
}
/*
 * Returns the maximum x2 coordinate of this triangle.
 * @return the maximum x2 coordinate of this triangle
 */
double GbTriangle3D::getX2Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x2max);
}
double GbTriangle3D::getX3Centroid()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x3s);
}
double GbTriangle3D::getX3Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x3min);
}
double GbTriangle3D::getX3Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x3max);
}

/*
 * Sets the specified point.
 * @param point the point
 * @param index the index (must be 0, 1, or 2)
 * @exception ArrayIndexOutOfBoundsException if the specified index is not valid
 */
void GbTriangle3D::setPoint(GbPoint3D *point, int index)
{
    if (index < 0 || index > 2)
        throw UbException(UB_EXARGS, "invalid index specified: ");
    this->points[index] = point;
    this->consistent    = false;
    this->calculateNormal();
}

/*
 * Returns the surface triangle set with new nodes !!!
 * @returns the surface triangle set with new nodes !!!
 */
vector<GbTriangle3D *> GbTriangle3D::getSurfaceTriangleSet()
{
    vector<GbTriangle3D *> triangles;

    triangles.push_back(
        new GbTriangle3D(new GbPoint3D(getPoint1()), new GbPoint3D(getPoint2()), new GbPoint3D(getPoint3())));

    return triangles;
}

/*
 * Returns the string representation of the triangle
 * @returns the string representation of the triangle
 */

string GbTriangle3D::toString()
{
    stringstream ss;
    ss << "GbTriangle3D[area=";
    ss << this->getArea();

    ss << ", x1s=" << this->x1s;
    ss << ", x2s=" << this->x2s;
    ss << ", x3s=" << this->x3s;
    ss << ", x1min=" << this->x1min;
    ss << ", x1max=" << this->x1max;
    ss << ", x2min=" << this->x2min;
    ss << ", x2max=" << this->x2max;
    ss << ", x3min=" << this->x3min;
    ss << ", x3max=" << this->x3max;
    ss << ", points1=" << this->points[0]->toString();
    ss << ", points2=" << this->points[1]->toString();
    ss << ", points3=" << this->points[2]->toString();
    ss << "]";
    return ((ss.str()).c_str());
}
/*======================================================================*/
double GbTriangle3D::getIntersectionRaytraceFactor(const double &x1, const double &x2, const double &x3,
                                                   const double &rx1, const double &rx2, const double &rx3)
{
    // e1 = v1 - v0
    double e1x1 = this->points[1]->x1 - this->points[0]->x1;
    double e1x2 = this->points[1]->x2 - this->points[0]->x2;
    double e1x3 = this->points[1]->x3 - this->points[0]->x3;

    // e2 = v2 - v0
    double e2x1 = this->points[2]->x1 - this->points[0]->x1;
    double e2x2 = this->points[2]->x2 - this->points[0]->x2;
    double e2x3 = this->points[2]->x3 - this->points[0]->x3;

    // p = d x e2
    double px1 = rx2 * e2x3 - rx3 * e2x2;
    double px2 = rx3 * e2x1 - rx1 * e2x3;
    double px3 = rx1 * e2x2 - rx2 * e2x1;

    // a = e1 dot p
    double a = e1x1 * px1 + e1x2 * px2 + e1x3 * px3;
    if (fabs(a) < 1.E-10)
        return -1.0;
    double f = 1.0 / a;

    // s = o - v0
    double sx1 = x1 - this->points[0]->x1;
    double sx2 = x2 - this->points[0]->x2;
    double sx3 = x3 - this->points[0]->x3;

    // u = f * ( s dot p)
    double u = f * (sx1 * px1 + sx2 * px2 + sx3 * px3);
    if (u < -1.E-10 || u > 1.0 + 1.E-10)
        return -1.0;

    // q = s x e1
    double qx1 = sx2 * e1x3 - sx3 * e1x2;
    double qx2 = sx3 * e1x1 - sx1 * e1x3;
    double qx3 = sx1 * e1x2 - sx2 * e1x1;

    // v = f*(e2 dot q)
    double v = f * (rx1 * qx1 + rx2 * qx2 + rx3 * qx3);
    if (v < -1.E-10 || (u + v) > 1.0 + 1.E-10)
        return -1.0;

    // t = f * (e2 dot q)
    return f * (e2x1 * qx1 + e2x2 * qx2 + e2x3 * qx3);
}

/*===========================================================*/

GbLine3D *GbTriangle3D::createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2)
{
    GbPoint3D *result = this->calculateIntersectionPoints3D(&point1, &point2);
    if (!result)
        return NULL;

    return new GbLine3D(result, new GbPoint3D(point2));

    // return GbSystem::createClipLine3D(point1, point2,
    // p1->getX1Coordinate(),p1->getX2Coordinate(),p1->getX3Coordinate(),
    // p2->getX1Coordinate(),p2->getX2Coordinate(),p2->getX3Coordinate() );
}

// von Navodit ...
/*===========================================================*/
GbPoint3D *GbTriangle3D::calculateIntersectionPoints3D(GbLine3D *line)
{
    return this->calculateIntersectionPoints3D(line->getPoint1(), line->getPoint2());
}
/*===========================================================*/
GbPoint3D *GbTriangle3D::calculateIntersectionPoints3D(GbPoint3D *linePoint1, GbPoint3D *linePoint2)
{
    GbVector3D Point1(linePoint1->x1, linePoint1->x2, linePoint1->x3);
    GbVector3D Point2(linePoint2->x1, linePoint2->x2, linePoint2->x3);
    GbVector3D direction = Point2 - Point1;
    GbVector3D GbPoint3D1(this->getPoint1()->x1, this->getPoint1()->x2, this->getPoint1()->x3);
    GbVector3D GbPoint3D2(this->getPoint2()->x1, this->getPoint2()->x2, this->getPoint2()->x3);
    GbVector3D GbPoint3D3(this->getPoint3()->x1, this->getPoint3()->x2, this->getPoint3()->x3);
    GbVector3D V2V1     = GbPoint3D2 - GbPoint3D1;
    GbVector3D V3V1     = GbPoint3D3 - GbPoint3D1;
    GbVector3D V2V1V3V1 = V2V1.Cross(V3V1);
    V2V1V3V1.Normalize();
    GbVector3D Normal = V2V1V3V1;

    double d     = -Normal.Dot(GbPoint3D1);
    double denom = Normal.Dot(direction);

    if (UbMath::zero(denom))
        return NULL; // line does not intersect the plane of the triangle !
    else {
        double mu = -1. * (d + Point1.Dot(Normal)) / denom; // mu = -(d+ Normal.Point1)/denom

        //   GbVector3D p1 = Point2-Point1;
        //   GbVector3D p2 = p1*mu;
        //   GbVector3D p3 = Point1+p2;
        GbVector3D point = Point1 + mu * (Point2 - Point1);

        if (mu < 0.0 || mu > 1.0)
            return NULL; // Point of intersection of line and plane does not lie on the triangle
        else {
            // Test whether Point lies inside the triangle or not
            bool test       = true;
            GbVector3D a    = GbPoint3D1 - point;
            GbVector3D b    = GbPoint3D2 - point;
            GbVector3D c    = GbPoint3D3 - point;
            GbVector3D ab   = a.Cross(b);
            GbVector3D bc   = b.Cross(c);
            GbVector3D ca   = c.Cross(a);
            GbVector3D Q1   = ab * 0.5;
            GbVector3D Q2   = bc * 0.5;
            GbVector3D Q3   = ca * 0.5;
            GbVector3D Q1Q2 = Q1 + Q2;
            GbVector3D Q    = Q1Q2 + Q3;

            if (UbMath::less(Q.Dot(Q1), 0.0))
                test = false;
            if (UbMath::less(Q.Dot(Q2), 0.0))
                test = false;
            if (UbMath::less(Q.Dot(Q3), 0.0))
                test = false;

            if (test == true)
                return (new GbPoint3D(point.X1(), point.X2(), point.X3()));
            else
                return NULL;
        }
    }
}

/**
 * Returns the distance between the 3D triangle and the specified 3D Point
 * @param point the 3D point from whom the distance is to be calculated
 * @return the distance of the specified point from the triangle
 */
double GbTriangle3D::calculateDistanceToPoint3D(GbPoint3D *point)
{
    return this->calculateDistanceToPoint3D(point->x1, point->x2, point->x3);
}
/*=======================================================================*/
double GbTriangle3D::calculateDistanceToPoint3D(const double &x1, const double &x2, const double &x3)
{
    //
    // throw UbException(UB_EXARGS,"Ich glaub GbTriangle3D::calculateDistanceToPoint3D(...) kann man so nicht
    // nehmen,jedenfalls nicht fuer die q's");
    cout << "??? ch glaub GbTriangle3D::calculateDistanceToPoint3D(...) kann man so nicht nehmen,jedenfalls nicht fuer "
            "die q's"
         << endl;
    GbVector3D P0(x1, x2, x3);
    GbVector3D P1(this->points[0]->x1, this->points[0]->x2, this->points[0]->x3);
    GbVector3D P2(this->points[1]->x1, this->points[1]->x2, this->points[1]->x3);
    GbVector3D P3(this->points[2]->x1, this->points[2]->x2, this->points[2]->x3);

    // Determine normal to triangle
    GbVector3D Normal = (P1 - P2).Cross(P1 - P3);
    double alpha      = UbMath::ACos((P1 - P0).Dot(Normal) / ((P1 - P0).Length() * Normal.Length()));

    double P0P0dash = (P0 - P1).Length() * cos(alpha);
    Normal.Normalize();
    GbVector3D Projection = Normal * (-P0P0dash);

    GbVector3D P0dash = P0 + Projection;

    // Check if point P0dash lies within the triangle P1P2P3.
    bool test = false;
    if (((P1 - P0).Cross(P2 - P0)).Dot(Normal) > 0)
        test = true;
    if (((P2 - P0).Cross(P3 - P0)).Dot(Normal) > 0)
        test = true;
    if (((P3 - P0).Cross(P1 - P0)).Dot(Normal) > 0)
        test = true;

    if (test == true)
        return (P0 - P0dash).Length();
    else
    // Determine the distance of point P0 from all edges and vertices and return the minimum distance
    {
        double dP0P1 = (P0 - P1).Length(); // Distance of Point P0 from Point P1
        double dP0P2 = (P0 - P2).Length(); // Distance of Point P0 from Point P2
        double dP0P3 = (P0 - P3).Length(); // Distance of Point P0 from Point P3

        GbVector3D MP1P2 = P2 - P1; // Direction vector for line P1P2
        GbVector3D MP2P3 = P3 - P2; // Direction vector for line P2P3
        GbVector3D MP3P1 = P1 - P3; // Direction vector for line P3P1

        double tP1P2 = MP1P2.Dot(P0 - P1) / MP1P2.Dot(MP1P2);
        double tP2P3 = MP2P3.Dot(P0 - P2) / MP2P3.Dot(MP2P3);
        double tP3P1 = MP3P1.Dot(P0 - P3) / MP3P1.Dot(MP3P1);

        double dP1P2 = (P0 - (P1 + (MP1P2 * tP1P2))).Length(); // Distance of Point P0 from line P1P2
        double dP2P3 = (P0 - (P2 + (MP2P3 * tP2P3))).Length(); // Distance of Point P0 from line P2P3
        double dP3P1 = (P0 - (P3 + (MP3P1 * tP3P1))).Length(); // Distance of Point P0 from line P3P1

        double distanceP0[6]; // Array to store all the distances from Point P0
        distanceP0[0] = dP0P1;
        distanceP0[1] = dP0P2;
        distanceP0[2] = dP0P3;
        distanceP0[3] = dP1P2;
        distanceP0[4] = dP2P3;
        distanceP0[5] = dP3P1;

        double d = 0.0;
        // Find the minimum distance from Point P0
        for (int i = 0; i < 6; i++) {
            if (distanceP0[i] < d)
                d = distanceP0[i];
        }
        return d;
    }
}
/**
 * Returns the normalized distance between the 3D triangle and the specified 3D Point
 * copied from Benjamin A.
 * @param point the 3D point from whom the distance is to be calculated
 * @return the distance of the specified point from the triangle
 */
double GbTriangle3D::calculateNormalizedDistanceToPoint3D(const double &x1, const double &y1, const double &z1,
                                                          const double &x2, const double &y2, const double &z2)
{
    // face* pf
    double xa, xb, xc, ya, yb, yc, za, zb, zc;
    // double xp, yp, zp;
    double tt = 0, xi = 0, eta = 0;
    double zaehler, nenner;
    double wurzel3 = sqrt(3.);

    // Weltkoordinaten der Dreiecke
    xa = this->points[0]->x1;
    xb = this->points[1]->x1;
    xc = this->points[2]->x1;

    ya = this->points[0]->x2;
    yb = this->points[1]->x2;
    yc = this->points[2]->x2;

    za = this->points[0]->x3;
    zb = this->points[1]->x3;
    zc = this->points[2]->x3;

    // Shape-Funktionen zum Berechnen der Schnittpunkte
    zaehler = static_cast<double>(((-1.0 * zc + zb) * ya + (yc - 1.0 * yb) * za + zc * yb - 1.0 * zb * yc) * x1 +
                                  ((-1.0 * zb + zc) * xa + (xb - 1.0 * xc) * za - 1.0 * xb * zc + xc * zb) * y1 +
                                  ((-1.0 * yc + yb) * xa + (-1.0 * xb + xc) * ya - 1.0 * xc * yb + xb * yc) * z1 +
                                  ((-1.0 * zc + zb) * ya + (yc - 1.0 * yb) * za + zc * yb - 1.0 * zb * yc) * x2 +
                                  ((-1.0 * zb + zc) * xa + (xb - 1.0 * xc) * za - 1.0 * xb * zc + xc * zb) * y2 +
                                  ((-1.0 * yc + yb) * xa + (-1.0 * xb + xc) * ya - 1.0 * xc * yb + xb * yc) * z2 +
                                  (2.0 * zb * yc - 2.0 * zc * yb) * xa + (2.0 * xb * zc - 2.0 * xc * zb) * ya +
                                  (-2.0 * xb * yc + 2.0 * xc * yb) * za);
    nenner  = static_cast<double>((((-1.0 * zc + zb) * ya + (yc - 1.0 * yb) * za + zc * yb - 1.0 * zb * yc) * x1 +
                                  ((-1.0 * zb + zc) * xa + (xb - 1.0 * xc) * za - 1.0 * xb * zc + xc * zb) * y1 +
                                  ((-1.0 * yc + yb) * xa + (-1.0 * xb + xc) * ya - 1.0 * xc * yb + xb * yc) * z1 +
                                  ((-1.0 * zb + zc) * ya + (-1.0 * yc + yb) * za - 1.0 * zc * yb + zb * yc) * x2 +
                                  ((-1.0 * zc + zb) * xa + (-1.0 * xb + xc) * za + xb * zc - 1.0 * xc * zb) * y2 +
                                  ((yc - 1.0 * yb) * xa + (xb - 1.0 * xc) * ya + xc * yb - 1.0 * xb * yc) * z2));
    if (UbMath::greater(nenner, 0.0))
        tt = zaehler / nenner;
    else
        tt = -999.;

    zaehler = static_cast<double>(((-2.0 * zc + za + zb) * y2 + (-1.0 * yb - 1.0 * ya + 2.0 * yc) * z2 + zc * ya -
                                   1.0 * zb * yc + zc * yb - 1.0 * za * yc) *
                                      x1 +
                                  ((-1.0 * za + 2.0 * zc - 1.0 * zb) * x2 + (xa - 2.0 * xc + xb) * z2 - 1.0 * xa * zc -
                                   1.0 * xb * zc + xc * za + xc * zb) *
                                      y1 +
                                  ((-2.0 * yc + ya + yb) * x2 + (-1.0 * xa - 1.0 * xb + 2.0 * xc) * y2 - 1.0 * xc * yb +
                                   xa * yc + xb * yc - 1.0 * xc * ya) *
                                      z1 +
                                  (zb * yc - 1.0 * zc * ya - 1.0 * zc * yb + za * yc) * x2 +
                                  (-1.0 * xc * za + xb * zc + xa * zc - 1.0 * xc * zb) * y2 +
                                  (xc * yb - 1.0 * xa * yc - 1.0 * xb * yc + xc * ya) * z2);
    nenner  = static_cast<double>((((zc - 1.0 * zb) * ya + (yb - 1.0 * yc) * za + zb * yc - 1.0 * zc * yb) * x1 +
                                  ((zb - 1.0 * zc) * xa + (xc - 1.0 * xb) * za - 1.0 * xc * zb + xb * zc) * y1 +
                                  ((-1.0 * yb + yc) * xa + (xb - 1.0 * xc) * ya - 1.0 * xb * yc + xc * yb) * z1 +
                                  ((zb - 1.0 * zc) * ya + (-1.0 * yb + yc) * za + zc * yb - 1.0 * zb * yc) * x2 +
                                  ((zc - 1.0 * zb) * xa + (xb - 1.0 * xc) * za - 1.0 * xb * zc + xc * zb) * y2 +
                                  ((yb - 1.0 * yc) * xa + (xc - 1.0 * xb) * ya + xb * yc - 1.0 * xc * yb) * z2));
    if (UbMath::greater(nenner, 0.0))
        xi = zaehler / nenner;
    else
        xi = -999.;

    zaehler = static_cast<double>(((za - 1.0 * zb) * y2 + (-1.0 * ya + yb) * z2 - 1.0 * za * yb + zb * ya) * x1 +
                                  ((-1.0 * za + zb) * x2 + (xa - 1.0 * xb) * z2 - 1.0 * xa * zb + xb * za) * y1 +
                                  ((ya - 1.0 * yb) * x2 + (xb - 1.0 * xa) * y2 + xa * yb - 1.0 * xb * ya) * z1 +
                                  (-1.0 * zb * ya + za * yb) * x2 + (-1.0 * xb * za + xa * zb) * y2 +
                                  (-1.0 * xa * yb + xb * ya) * z2);
    nenner  = static_cast<double>((((zc - 1.0 * zb) * ya + (yb - 1.0 * yc) * za + zb * yc - 1.0 * zc * yb) * x1 +
                                  ((zb - 1.0 * zc) * xa + (xc - 1.0 * xb) * za - 1.0 * xc * zb + xb * zc) * y1 +
                                  ((-1.0 * yb + yc) * xa + (xb - 1.0 * xc) * ya - 1.0 * xb * yc + xc * yb) * z1 +
                                  ((zb - 1.0 * zc) * ya + (-1.0 * yb + yc) * za + zc * yb - 1.0 * zb * yc) * x2 +
                                  ((zc - 1.0 * zb) * xa + (xb - 1.0 * xc) * za - 1.0 * xb * zc + xc * zb) * y2 +
                                  ((yb - 1.0 * yc) * xa + (xc - 1.0 * xb) * ya + xb * yc - 1.0 * xc * yb) * z2));
    if (UbMath::greater(nenner, 0.0))
        eta = static_cast<double>((zaehler / nenner) * wurzel3 * -1.);
    else
        eta = -999.;

    if (tt >= -1.0 - UbMath::Epsilon<double>::val() && tt <= 1.0) {
        if (xi >= -1.0 + eta / wurzel3 - UbMath::Epsilon<double>::val() &&
            xi <= 1.0 - eta / wurzel3 + UbMath::Epsilon<double>::val()) {
            if (eta >= 0 - UbMath::Epsilon<double>::val() && eta <= wurzel3 + UbMath::Epsilon<double>::val()) {
                /*xp = x1*(0.5-tt/2)+x2*(0.5+tt/2);
                yp = y1*(0.5-tt/2)+y2*(0.5+tt/2);
                zp = z1*(0.5-tt/2)+z2*(0.5+tt/2);*/
                return static_cast<double>((sqrt(pow((x1 * (0.5 - tt / 2) + x2 * (0.5 + tt / 2)) - x1, 2) +
                                                 pow((y1 * (0.5 - tt / 2) + y2 * (0.5 + tt / 2)) - y1, 2) +
                                                 pow((z1 * (0.5 - tt / 2) + z2 * (0.5 + tt / 2)) - z1, 2))));
            }
        }
    }
    return (-999.);
}
/*
 * Returns true if the specified 2D point lies within (or on the border of) this 2D triangle.
 * @param point the 2D point to check
 * @return true if the specified 2D point lies within (or on the border of) this 2D triangle
 */
bool GbTriangle3D::enclosesPoint2D(double x1, double x2)
{
    int i = 0;
    // Punkt(x1,x2) liegt auf einem der Eckpunkte
    if (x1 == this->getPoint(0)->getX1Coordinate() && x2 == this->getPoint(0)->getX2Coordinate())
        return true;
    if (x1 == this->getPoint(1)->getX1Coordinate() && x2 == this->getPoint(1)->getX2Coordinate())
        return true;
    if (x1 == this->getPoint(2)->getX1Coordinate() && x2 == this->getPoint(2)->getX2Coordinate())
        return true;

    // Erste Grade aus dem zu pruefenden Punkt(x,y) und einem zweiten Punkt(x+0.333,y+2.333)
    GbPoint3D p1;
    p1.setX1(x1);
    p1.setX2(x2);
    p1.setX3(0.0);
    GbPoint3D p2;
    p2.setX1(x1 + 0.333);
    p2.setX2(x2 + 3.333);
    p2.setX3(0.0);
    // Punkte des Dreiecks auf 2D reduziert
    GbPoint3D dp1;
    dp1.setX1(this->getPoint(0)->getX1Coordinate());
    dp1.setX2(this->getPoint(0)->getX2Coordinate());
    dp1.setX3(0.0);
    GbPoint3D dp2;
    dp2.setX1(this->getPoint(1)->getX1Coordinate());
    dp2.setX2(this->getPoint(1)->getX2Coordinate());
    dp2.setX3(0.0);
    GbPoint3D dp3;
    dp3.setX1(this->getPoint(2)->getX1Coordinate());
    dp3.setX2(this->getPoint(2)->getX2Coordinate());
    dp3.setX3(0.0);
    // ueberpruefen, ob der Punkt(x,y) innerhalt der Boundingbox des Dreiecks liegt
    if (x1 < this->getX1Maximum() && x1 > getX1Minimum() && x2 < this->getX2Maximum() && x2 > getX2Minimum()) {
        GbPoint3D *dummy = NULL;
        // ueberpruefen, ob der Punkt innerhalb des Dreiecks liegt
        dummy = GbSystem3D::calculateIntersectionPoint3D(p1, p2, dp1, dp2);
        if (dummy != NULL) {
            if (dummy->getX1Coordinate() == p1.getX1Coordinate() && dummy->getX2Coordinate() == p1.getX2Coordinate()) {
                delete dummy;
                return true;
            } else if (dummy->getX1Coordinate() > p1.getX1Coordinate()) {
                i++;
            } else {
                i--;
            }
        }
        if (dummy)
            delete dummy;

        dummy = GbSystem3D::calculateIntersectionPoint3D(p1, p2, dp2, dp3);
        if (dummy != NULL) {
            if (dummy->getX1Coordinate() == p1.getX1Coordinate() && dummy->getX2Coordinate() == p1.getX2Coordinate()) {
                if (dummy)
                    delete dummy;
                return true;
            } else if (dummy->getX1Coordinate() > p1.getX1Coordinate()) {
                i++;
            } else {
                i--;
            }
        }
        if (dummy)
            delete dummy;

        dummy = GbSystem3D::calculateIntersectionPoint3D(p1, p2, dp3, dp1);
        if (dummy != NULL) {
            if (dummy->getX1Coordinate() == p1.getX1Coordinate() && dummy->getX2Coordinate() == p1.getX2Coordinate()) {
                if (dummy)
                    delete dummy;
                return true;
            } else if (dummy->getX1Coordinate() > p1.getX1Coordinate()) {
                i++;
            } else {
                i--;
            }
        }
        if (dummy)
            delete dummy;
    }
    if (i == -1)
        return true;
    if (i == 1)
        return true;

    return false;
}

///*
//* Returns a new 2D polygon clipped by the specified 2D rectangle (result may be null!).
//* @param rectangle the 2D rectangle
//* @return a new 2D polygon clipped by the specified 2D rectangle
//*/
GbPolygon3D *GbTriangle3D::createClippedPolygon3D(GbCuboid3D *cube)
{
    return (GbSystem3D::clipPolygon3D(this->getPoints(), cube->getPoint1()->getX1Coordinate(),
                                      cube->getPoint1()->getX2Coordinate(), cube->getPoint1()->getX3Coordinate(),
                                      cube->getPoint2()->getX1Coordinate(), cube->getPoint2()->getX2Coordinate(),
                                      cube->getPoint2()->getX3Coordinate()));
}
///*
//* Returns a new 2D polygon clipped by the specified 2D rectangle (result may be null!).
//* @param p1 the 1st point of the rectangle
//* @param p2 the 2nd point of the rectangle
//* @return a new 2D polygon clipped by the specified 2D rectangle
//*/
// public GbPolygon2D createClippedPolygon2D(GbPoint2D p1, GbPoint2D p2)
//{
//   return(GbSystem.clipPolygon2D(this.points, p1.x1, p1.x2, p2.x1, p2.x2));
//}
/*
 * Returns a new 2D polygon clipped by the specified 2D rectangle (result may be null!).
 * @param p1x1 the 1st x1 coordinate of the rectangle
 * @param p1x2 the 1st x2 coordinate of the rectangle
 * @param p2x1 the 2nd x1 coordinate of the rectangle
 * @param p2x2 the 2nd x2 coordinate of the rectangle
 * @return a new 2D polygon clipped by the specified 2D rectangle
 */
GbPolygon3D *GbTriangle3D::createClippedPolygon3D(const double &p1x1, const double &p1x2, const double &p1x3,
                                                  const double &p2x1, const double &p2x2, const double &p2x3)
{
    return (GbSystem3D::clipPolygon3D(this->getPoints(), p1x1, p1x2, p1x3, p2x1, p2x2, p2x3));
}

/*
 * Returns true if the specified 2D rectangle lies completely within this 2D triangle.
 * @param rectangle the 2D rectangle to check
 * @return true if the specified 2D rectangle lies completely within this 2D triangle
 */
// bool enclosesRectangle2D(GbRectangle2D *rectangle)
//{
//   GbPolygon2D p = GbSystem.clipPolygon2D(this.points, rectangle.p1.x1, rectangle.p1.x2, rectangle.p2.x1,
//   rectangle.p2.x2); return(p!=null && GbSystem.equal(Math.abs(p.getArea()), rectangle.getArea()));
//}
/*
 * Returns true if the specified 2D rectangle lies completely within this 2D triangle.
 * @param p1 the 1st point of the rectangle to check
 * @param p2 the 2nd point of the rectangle to check                         triangle
 * @return true if the specified 2D rectangle lies completely within this 2D
 */
// public boolean enclosesRectangle2D(GbPoint2D p1, GbPoint2D p2)
//{
//   GbPolygon2D p = GbSystem.clipPolygon2D(this.points, p1.x1, p1.x2, p2.x1, p2.x2);
//   return(p!=null && GbSystem.equal(Math.abs(p.getArea()), Math.abs((p1.x1-p2.x1)*(p1.x2-p2.x2))));
//}
/*
 * Returns true if the specified 2D rectangle lies completely within this 2D triangle.
 * @param p1x1 the 1st x1 coordinate of the rectangle to check
 * @param p1x2 the 1st x2 coordinate of the rectangle to check
 * @param p2x1 the 2nd x1 coordinate of the rectangle to check
 * @param p2x2 the 2nd x2 coordinate of the rectangle to check
 * @return true if the specified 2D rectangle lies completely within this 2D triangle
 */
// public boolean enclosesRectangle2D(double p1x1, double p1x2, double p2x1, double p2x2)
//{
//   GbPolygon2D p = GbSystem.clipPolygon2D(this.points, p1x1, p1x2, p2x1, p2x2);
//   return(p!=null && GbSystem.equal(Math.abs(p.getArea()), Math.abs((p1x1-p2x1)*(p1x2-p2x2))));
//}

/*
 * Returns true if the specified 2D rectangle is crossed by this 2D triangle.
 * @param rectangle the 2D rectangle to check
 * @return true if the specified 2D rectangle is crossed by this 2D triangle
 */
// public boolean crossesRectangle2D(GbRectangle2D rectangle)
//{
//   GbPolygon2D p = GbSystem.clipPolygon2D(this.points, rectangle.p1.x1, rectangle.p1.x2, rectangle.p2.x1,
//   rectangle.p2.x2); return(p!=null && GbSystem.inOpenInterval(Math.abs(p.getArea()), 0.0, rectangle.getArea()));
//}
/*
 * Returns true if the specified 2D rectangle is crossed by this 2D triangle.
 * @param p1 the 1st point of the rectangle to check
 * @param p2 the 2nd point of the rectangle to check
 * @return true if the specified 2D rectangle is crossed by this 2D triangle
 */
// public boolean crossesRectangle2D(GbPoint2D p1, GbPoint2D p2)
//{
//   GbPolygon2D p = GbSystem.clipPolygon2D(this.points, p1.x1, p1.x2, p2.x1, p2.x2);
//   return(p!=null && GbSystem.inOpenInterval(Math.abs(p.getArea()), 0.0, Math.abs((p1.x1-p2.x1)*(p1.x2-p2.x2))));
//}
/*
 * Returns true if the specified 2D rectangle is crossed by this 2D triangle.
 * @param p1x1 the 1st x1 coordinate of the rectangle to check
 * @param p1x2 the 1st x2 coordinate of the rectangle to check
 * @param p2x1 the 2nd x1 coordinate of the rectangle to check
 * @param p2x2 the 2nd x2 coordinate of the rectangle to check
 * @return true if the specified 2D rectangle is crossed by this 2D triangle
 */
// public boolean crossesRectangle2D(double p1x1, double p1x2, double p2x1, double p2x2)
//{
//   GbPolygon2D p = GbSystem.clipPolygon2D(this.points, p1x1, p1x2, p2x1, p2x2);
//   return(p!=null && GbSystem.inOpenInterval(Math.abs(p.getArea()), 0.0, Math.abs((p1x1-p2x1)*(p1x2-p2x2))));
//}

/*
 * Returns true if the specified 2D rectangle lies (at least partly) within this 2D triangle.
 * @param rectangle the 2D rectangle to check
 * @return true if the specified 2D rectangle lies (at least partly) within this 2D triangle
 */
// public boolean enclosesOrCrossesRectangle2D(GbRectangle2D rectangle)
//{
//   GbPolygon2D p = GbSystem.clipPolygon2D(this.points, rectangle.p1.x1, rectangle.p1.x2, rectangle.p2.x1,
//   rectangle.p2.x2); return(p!=null && GbSystem.greater(Math.abs(p.getArea()), 0.0));
//}
/*
 * Returns true if the specified 2D rectangle lies (at least partly) within this 2D triangle.
 * @param p1 the 1st point of the rectangle to check
 * @param p2 the 2nd point of the rectangle to check
 * @return true if the specified 2D rectangle lies (at least partly) within this 2D triangle
 */
// public boolean enclosesOrCrossesRectangle2D(GbPoint2D p1, GbPoint2D p2)
//{
//   GbPolygon2D p = GbSystem.clipPolygon2D(this.points, p1.x1, p1.x2, p2.x1, p2.x2);
//   return(p!=null && GbSystem.greater(Math.abs(p.getArea()), 0.0));
//}
/*
 * Returns true if the specified 2D rectangle lies (at least partly) within this 2D triangle.
 * @param p1x1 the 1st x1 coordinate of the rectangle to check
 * @param p1x2 the 1st x2 coordinate of the rectangle to check
 * @param p2x1 the 2nd x1 coordinate of the rectangle to check
 * @param p2x2 the 2nd x2 coordinate of the rectangle to check
 * @return true if the specified 2D rectangle lies (at least partly) within this 2D triangle
 */
// public boolean enclosesOrCrossesRectangle2D(double p1x1, double p1x2, double p2x1, double p2x2)
//{
//   GbPolygon2D p = GbSystem.clipPolygon2D(this.points, p1x1, p1x2, p2x1, p2x2);
//   return(p!=null && GbSystem.greater(Math.abs(p.getArea()), 0.0));
//}
/*======================================================================*/

/*======================================================================*/
/*  Private Methoden                                                    */
/*                                                                      */
void GbTriangle3D::calculateValues()
{
    this->x1min = this->points[0]->x1;
    this->x1max = this->points[0]->x1;
    this->x2min = this->points[0]->x2;
    this->x2max = this->points[0]->x2;
    this->x3min = this->points[0]->x3;
    this->x3max = this->points[0]->x3;

    if (this->points[1]->x1 < this->x1min)
        this->x1min = this->points[1]->x1;
    if (this->points[1]->x1 > this->x1max)
        this->x1max = this->points[1]->x1;
    if (this->points[1]->x2 < this->x2min)
        this->x2min = this->points[1]->x2;
    if (this->points[1]->x2 > this->x2max)
        this->x2max = this->points[1]->x2;
    if (this->points[1]->x3 < this->x3min)
        this->x3min = this->points[1]->x3;
    if (this->points[1]->x3 > this->x3max)
        this->x3max = this->points[1]->x3;

    if (this->points[2]->x1 < this->x1min)
        this->x1min = this->points[2]->x1;
    if (this->points[2]->x1 > this->x1max)
        this->x1max = this->points[2]->x1;
    if (this->points[2]->x2 < this->x2min)
        this->x2min = this->points[2]->x2;
    if (this->points[2]->x2 > this->x2max)
        this->x2max = this->points[2]->x2;
    if (this->points[2]->x3 < this->x3min)
        this->x3min = this->points[2]->x3;
    if (this->points[2]->x3 > this->x3max)
        this->x3max = this->points[2]->x3;

    this->x1s = (this->points[0]->x1 + this->points[1]->x1 + this->points[2]->x1) / 3.0;
    this->x2s = (this->points[0]->x2 + this->points[1]->x2 + this->points[2]->x2) / 3.0;
    this->x3s = (this->points[0]->x3 + this->points[1]->x3 + this->points[2]->x3) / 3.0;

    GbVector3D A(points[0]->x1, points[0]->x2, points[0]->x3);
    GbVector3D B(points[1]->x1, points[1]->x2, points[1]->x3);
    GbVector3D C(points[2]->x1, points[2]->x2, points[2]->x3);
    GbVector3D AB    = B - A;
    GbVector3D AC    = C - A;
    GbVector3D N     = AB.Cross(AC);
    this->area       = 0.5 * N.Length();
    this->consistent = true;
}
/*======================================================================*/

/*======================================================================*/
GbVector3D GbTriangle3D::getNormal()
{
    this->calculateNormal();
    return normal;
}
/*======================================================================*/
void GbTriangle3D::init()
{
    x1s        = 0.0;
    x2s        = 0.0;
    x3s        = 0.0;
    x1min      = 0.0;
    x1max      = 0.0;
    x2min      = 0.0;
    x2max      = 0.0;
    area       = 0.0;
    consistent = false;
    points.resize(3, NULL);
}
/*=======================================================*/
void GbTriangle3D::calculateNormal()
{
    GbPoint3D *&a = points[0];
    GbPoint3D *&b = points[1];
    GbPoint3D *&c = points[2];
    normal[0]     = (c->getX3Coordinate() - a->getX3Coordinate()) * (b->getX2Coordinate() - a->getX2Coordinate()) -
                (b->getX3Coordinate() - a->getX3Coordinate()) * (c->getX2Coordinate() - a->getX2Coordinate());
    normal[1] = (b->getX3Coordinate() - a->getX3Coordinate()) * (c->getX1Coordinate() - a->getX1Coordinate()) -
                (b->getX1Coordinate() - a->getX1Coordinate()) * (c->getX3Coordinate() - a->getX3Coordinate());
    normal[2] = (b->getX1Coordinate() - a->getX1Coordinate()) * (c->getX2Coordinate() - a->getX2Coordinate()) -
                (b->getX2Coordinate() - a->getX2Coordinate()) * (c->getX1Coordinate() - a->getX1Coordinate());
    normal.Normalize();
}
/*=======================================================*/
// toDo:
double GbTriangle3D::getDistanceFromPoint(GbVector3D punct)
{
    GbVector3D Point1(this->getPoint1()->getX1Coordinate(), this->getPoint1()->getX2Coordinate(),
                      this->getPoint1()->getX3Coordinate());
    GbVector3D Point2(this->getPoint2()->getX1Coordinate(), this->getPoint2()->getX2Coordinate(),
                      this->getPoint2()->getX3Coordinate());
    GbVector3D Point3(this->getPoint3()->getX1Coordinate(), this->getPoint3()->getX2Coordinate(),
                      this->getPoint3()->getX3Coordinate());

    GbVector3D kDiff  = Point1 - punct;
    GbVector3D kEdge0 = Point2 - Point1;
    GbVector3D kEdge1 = Point3 - Point1;
    double fA00       = kEdge0.SquaredLength();
    double fA01       = kEdge0.Dot(kEdge1);
    double fA11       = kEdge1.SquaredLength();
    double fB0        = kDiff.Dot(kEdge0);
    double fB1        = kDiff.Dot(kEdge1);
    double fC         = kDiff.SquaredLength();
    double fDet       = fabs(fA00 * fA11 - fA01 * fA01);
    double fS         = fA01 * fB1 - fA11 * fB0;
    double fT         = fA01 * fB0 - fA00 * fB1;
    double fSqrDistance;

    if (fS + fT <= fDet) {
        if (fS < (double)0.0) {
            if (fT < (double)0.0) // region 4
            {
                if (fB0 < (double)0.0) {
                    fT = (double)0.0;
                    if (-fB0 >= fA00) {
                        fS           = (double)1.0;
                        fSqrDistance = fA00 + ((double)2.0) * fB0 + fC;
                    } else {
                        fS           = -fB0 / fA00;
                        fSqrDistance = fB0 * fS + fC;
                    }
                } else {
                    fS = (double)0.0;
                    if (fB1 >= (double)0.0) {
                        fT           = (double)0.0;
                        fSqrDistance = fC;
                    } else if (-fB1 >= fA11) {
                        fT           = (double)1.0;
                        fSqrDistance = fA11 + ((double)2.0) * fB1 + fC;
                    } else {
                        fT           = -fB1 / fA11;
                        fSqrDistance = fB1 * fT + fC;
                    }
                }
            } else // region 3
            {
                fS = (double)0.0;
                if (fB1 >= (double)0.0) {
                    fT           = (double)0.0;
                    fSqrDistance = fC;
                } else if (-fB1 >= fA11) {
                    fT           = (double)1.0;
                    fSqrDistance = fA11 + ((double)2.0) * fB1 + fC;
                } else {
                    fT           = -fB1 / fA11;
                    fSqrDistance = fB1 * fT + fC;
                }
            }
        } else if (fT < (double)0.0) // region 5
        {
            fT = (double)0.0;
            if (fB0 >= (double)0.0) {
                fS           = (double)0.0;
                fSqrDistance = fC;
            } else if (-fB0 >= fA00) {
                fS           = (double)1.0;
                fSqrDistance = fA00 + ((double)2.0) * fB0 + fC;
            } else {
                fS           = -fB0 / fA00;
                fSqrDistance = fB0 * fS + fC;
            }
        } else // region 0
        {
            // minimum at interior point
            double fInvDet = ((double)1.0) / fDet;
            fS *= fInvDet;
            fT *= fInvDet;
            fSqrDistance = fS * (fA00 * fS + fA01 * fT + ((double)2.0) * fB0) +
                           fT * (fA01 * fS + fA11 * fT + ((double)2.0) * fB1) + fC;
        }
    } else {
        double fTmp0, fTmp1, fNumer, fDenom;

        if (fS < (double)0.0) // region 2
        {
            fTmp0 = fA01 + fB0;
            fTmp1 = fA11 + fB1;
            if (fTmp1 > fTmp0) {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00 - 2.0f * fA01 + fA11;
                if (fNumer >= fDenom) {
                    fS           = (double)1.0;
                    fT           = (double)0.0;
                    fSqrDistance = fA00 + ((double)2.0) * fB0 + fC;
                } else {
                    fS           = fNumer / fDenom;
                    fT           = (double)1.0 - fS;
                    fSqrDistance = fS * (fA00 * fS + fA01 * fT + 2.0f * fB0) +
                                   fT * (fA01 * fS + fA11 * fT + ((double)2.0) * fB1) + fC;
                }
            } else {
                fS = (double)0.0;
                if (fTmp1 <= (double)0.0) {
                    fT           = (double)1.0;
                    fSqrDistance = fA11 + ((double)2.0) * fB1 + fC;
                } else if (fB1 >= (double)0.0) {
                    fT           = (double)0.0;
                    fSqrDistance = fC;
                } else {
                    fT           = -fB1 / fA11;
                    fSqrDistance = fB1 * fT + fC;
                }
            }
        } else if (fT < (double)0.0) // region 6
        {
            fTmp0 = fA01 + fB1;
            fTmp1 = fA00 + fB0;
            if (fTmp1 > fTmp0) {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00 - ((double)2.0) * fA01 + fA11;
                if (fNumer >= fDenom) {
                    fT           = (double)1.0;
                    fS           = (double)0.0;
                    fSqrDistance = fA11 + ((double)2.0) * fB1 + fC;
                } else {
                    fT           = fNumer / fDenom;
                    fS           = (double)1.0 - fT;
                    fSqrDistance = fS * (fA00 * fS + fA01 * fT + ((double)2.0) * fB0) +
                                   fT * (fA01 * fS + fA11 * fT + ((double)2.0) * fB1) + fC;
                }
            } else {
                fT = (double)0.0;
                if (fTmp1 <= (double)0.0) {
                    fS           = (double)1.0;
                    fSqrDistance = fA00 + ((double)2.0) * fB0 + fC;
                } else if (fB0 >= (double)0.0) {
                    fS           = (double)0.0;
                    fSqrDistance = fC;
                } else {
                    fS           = -fB0 / fA00;
                    fSqrDistance = fB0 * fS + fC;
                }
            }
        } else // region 1
        {
            fNumer = fA11 + fB1 - fA01 - fB0;
            if (fNumer <= (double)0.0) {
                fS           = (double)0.0;
                fT           = (double)1.0;
                fSqrDistance = fA11 + ((double)2.0) * fB1 + fC;
            } else {
                fDenom = fA00 - 2.0f * fA01 + fA11;
                if (fNumer >= fDenom) {
                    fS           = (double)1.0;
                    fT           = (double)0.0;
                    fSqrDistance = fA00 + ((double)2.0) * fB0 + fC;
                } else {
                    fS           = fNumer / fDenom;
                    fT           = (double)1.0 - fS;
                    fSqrDistance = fS * (fA00 * fS + fA01 * fT + ((double)2.0) * fB0) +
                                   fT * (fA01 * fS + fA11 * fT + ((double)2.0) * fB1) + fC;
                }
            }
        }
    }

    // account for numerical round-off error
    if (fSqrDistance < (double)0.0) {
        fSqrDistance = (double)0.0;
    }
    /*
        m_kClosestPoint0 = punct;
        m_kClosestPoint1 = m_rkTriangle.V[0] + fS*kEdge0 + fT*kEdge1;
        m_afTriangleBary[1] = fS;
        m_afTriangleBary[2] = fT;
        m_afTriangleBary[0] = (double)1.0 - fS - fT;
    */
    return sqrt(fSqrDistance);
}

//! \}
