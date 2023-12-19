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
#include <geometry3d/GbLine3D.h>
#include <geometry3d/GbPoint3D.h>
#include <geometry3d/GbSphere3D.h>
#include <geometry3d/GbSystem3D.h>
#include <geometry3d/GbTriangle3D.h>

using namespace std;

/*=====================================================*/
GbSphere3D::GbSphere3D() : GbObject3D(), UbObserver()
{
    this->setName("sphere");
    radius   = 0;
    midPoint = new GbPoint3D(0, 0, 0);
}
/*=====================================================*/
GbSphere3D::GbSphere3D(const double &x1, const double &x2, const double &x3, const double &radius)
    : GbObject3D(), UbObserver()
{
    this->setName("sphere");
    midPoint = new GbPoint3D(x1, x2, x3);
    midPoint->addObserver(this);

    this->radius      = radius;
    triangulationMode = RAYPROJECTION;
    // triangulationMode = CUBOIDPROJECTION;
}
/*=====================================================*/
GbSphere3D::GbSphere3D(const GbSphere3D &sphere) : GbObject3D(sphere), UbObserver()
{
    this->setName("sphere");

    this->midPoint    = sphere.midPoint->clone();
    this->radius      = sphere.radius;
    triangulationMode = RAYPROJECTION;

    this->midPoint->addObserver(this);
}
/*=====================================================*/
GbSphere3D::GbSphere3D(GbSphere3D *sphere) : GbObject3D(), UbObserver()
{
    this->setName(sphere->getName());
    midPoint = sphere->midPoint->clone();
    midPoint->addObserver(this);

    this->radius      = sphere->getRadius();
    triangulationMode = RAYPROJECTION;
}
/*=====================================================*/
GbSphere3D::~GbSphere3D()
{
    if (this->midPoint)
        this->midPoint->removeObserver(this);
}
/*=====================================================*/
void GbSphere3D::finalize()
{
    if (this->midPoint) {
        this->midPoint->removeObserver(this);
        this->midPoint->finalize();
        delete this->midPoint;
        this->midPoint = NULL;
    }

    if (this->midPoint)
        this->midPoint->removeObserver(this);
}
/*=====================================================*/
bool GbSphere3D::intersects(SPtr<GbSphere3D> sphere)
{
    return this->getDistance(sphere->midPoint) < radius + sphere->radius;
}
/*=======================================================*/
void GbSphere3D::setCenterCoordinates(const double &x1, const double &x2, const double &x3)
{
    this->translate(x1 - getX1Centroid(), x2 - getX2Centroid(), x3 - getX3Centroid());
}

void GbSphere3D::setCenterCoordinates(const UbTupleDouble3 &position)
{
    this->setCenterCoordinates(val<1>(position), val<2>(position), val<3>(position));
}

/*=====================================================*/
double GbSphere3D::getDistance(GbPoint3D *p)
{
    return this->getDistance(p->getX1Centroid(), p->getX2Coordinate(), p->getX3Coordinate());
}
/*=====================================================*/
void GbSphere3D::setCenterX1Coordinate(const double &value)
{
    if (this->midPoint)
        this->midPoint->setX1(value);
    else
        throw UbException(UB_EXARGS, "Sphere has no midPoint");
    // kein notifyObserver(), da der knoten notifyObserver() ausfuehrt und die GbSphere dieses event
    // abfaengt und dann selbst notifyObservers ausfuehrt ;-)
}
/*=====================================================*/
void GbSphere3D::setCenterX2Coordinate(const double &value)
{
    if (this->midPoint)
        this->midPoint->setX2(value);
    else
        throw UbException(UB_EXARGS, "Sphere has no midPoint");
    // kein notifyObserver(), da der knoten notifyObserver() ausfuehrt und die GbSphere dieses event
    // abfaengt und dann selbst notifyObservers ausfuehrt ;-)
}
/*=====================================================*/
void GbSphere3D::setCenterX3Coordinate(const double &value)
{
    if (this->midPoint)
        this->midPoint->setX3(value);
    else
        throw UbException(UB_EXARGS, "sphere has no midPoint");
    // kein notifyObserver(), da der knoten notifyObserver() ausfuehrt und die GbSphere dieses event
    // abfaengt und dann selbst notifyObservers ausfuehrt ;-)
}
/*=====================================================*/
void GbSphere3D::setRadius(const double &radius)
{
    if (radius != this->radius) {
        this->radius = radius;
        this->notifyObserversObjectChanged();
    }
}
/*=====================================================*/
double GbSphere3D::getDistance(const double &x1p, const double &x2p, const double &x3p)
{
    double deltaX1 = x1p - midPoint->getX1Coordinate();
    double deltaX2 = x2p - midPoint->getX2Coordinate();
    double deltaX3 = x3p - midPoint->getX3Coordinate();
    return sqrt(deltaX1 * deltaX1 + deltaX2 * deltaX2 + deltaX3 * deltaX3) - this->radius;
}
/*=====================================================*/
// true, wenn 'in Object' oder 'auf Boundary'!
bool GbSphere3D::isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p)
{
    double deltaX1 = x1p - midPoint->getX1Coordinate();
    double deltaX2 = x2p - midPoint->getX2Coordinate();
    double deltaX3 = x3p - midPoint->getX3Coordinate();

    return UbMath::lessEqual(deltaX1 * deltaX1 + deltaX2 * deltaX2 + deltaX3 * deltaX3, this->radius * this->radius);
}
/*=====================================================*/
// true, wenn 'in Object' oder 'auf Boundary'!
bool GbSphere3D::isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p, bool &pointIsOnBoundary)
{
    double deltaX1 = x1p - midPoint->getX1Coordinate();
    double deltaX2 = x2p - midPoint->getX2Coordinate();
    double deltaX3 = x3p - midPoint->getX3Coordinate();

    double distanceSquare = deltaX1 * deltaX1 + deltaX2 * deltaX2 + deltaX3 * deltaX3;
    double radiusSquare   = this->radius * this->radius;

    pointIsOnBoundary = UbMath::equal(distanceSquare, radiusSquare);

    return UbMath::lessEqual(distanceSquare, radiusSquare);
}
/*=====================================================*/
// bool GbSphere3D::crossCellCrossSection(double x11,double x21,double x12,double x22, double ra)
//{
//   if(this->isPointInCrossection(x11, x12) || this->isPointInCrossection(x21, x22) || this->isPointInCrossection(x11,
//   x22) || this->isPointInCrossection(x21, x12))
//   {
//        if(!this->isPointInCrossection(x11, x12) || !this->isPointInCrossection(x21, x22) ||
//!this->isPointInCrossection(x11, x22) || !this->isPointInCrossection(x21, x12)) return true;
//   }
//   return false;
//}
//
///*=====================================================*/
// bool GbSphere3D::cellCrossAndInsideCrossSection(double x11,double x21,double x12,double x22, double ra)
//{
//   if(this->isPointInCrossection(x11, x12) || this->isPointInCrossection(x21, x22) || this->isPointInCrossection(x11,
//   x22) || this->isPointInCrossection(x21, x12))  return true; return false;
//}
/*=====================================================*/
string GbSphere3D::toString()
{
    stringstream ss;
    ss << "GbSphere3D[";
    ss << "mid=" << midPoint->toString() << ", r=" << radius << "]";
    return ss.str();
}
/*=====================================================*/
GbLine3D *GbSphere3D::createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2)
{
    double factor = 100.0; // um rundungsfehler beim wurzelterm zu minimieren
    double xa = factor * point1.getX1Coordinate();
    double ya = factor * point1.getX2Coordinate();
    double za = factor * point1.getX3Coordinate();
    double xb = factor * point2.getX1Coordinate();
    double yb = factor * point2.getX2Coordinate();
    double zb = factor * point2.getX3Coordinate();
    double xm = factor * this->midPoint->getX1Coordinate();
    double ym = factor * this->midPoint->getX2Coordinate();
    double zm = factor * this->midPoint->getX3Coordinate();
    double r  = factor * this->radius;

    double xa2 = xa * xa;
    double ya2 = ya * ya;
    double za2 = za * za;
    double xb2 = xb * xb;
    double yb2 = yb * yb;
    double zb2 = zb * zb;
    double xm2 = xm * xm;
    double ym2 = ym * ym;
    double zm2 = zm * zm;
    double r2  = r * r;

    double wurzel =
        2.0 * xa * xb * ym2 - 2.0 * ya * yb * r2 + 2.0 * ya * ym * xb2 + 2.0 * yb * ym * za2 + 2.0 * ya * ym * zb2 +
        2.0 * xb * xm * za2 + 2.0 * za * zb * ym2 + 2.0 * xb * xm * ya2 + 2.0 * xa * xm * yb2 + 2.0 * yb * ym * xa2 +
        2.0 * zb * zm * ya2 + 2.0 * xa * xm * zb2 + 2.0 * za * zm * xb2 + 2.0 * za * zm * yb2 + 2.0 * xa * xb * zm2 -
        2.0 * xa * xb * r2 - 2.0 * za * zb * r2 + 2.0 * za * zb * xm2 - 2.0 * ya * yb * xa * xm +
        2.0 * ya * yb * xa * xb + 2.0 * zb * zm * xa2 - 2.0 * ya * yb * xb * xm + 2.0 * ya * yb * xm2 -
        2.0 * ya * yb * zb * zm + 2.0 * ya * yb * zm2 + 2.0 * zb * zm * yb * ym - 2.0 * zb * zm * ya * ym +
        2.0 * zb * zm * xb * xm - 2.0 * xa * xm * yb * ym + 2.0 * xa * xm * za * zm + 2.0 * xa * xm * ya * ym -
        2.0 * yb * ym * za * zm + 2.0 * yb * ym * xb * xm + 2.0 * za * zm * ya * ym - 2.0 * za * zm * xb * xm -
        2.0 * ya * ym * xb * xm + 2.0 * za * zb * xa * xb - 2.0 * za * zb * xa * xm - 2.0 * za * zb * xb * xm +
        2.0 * za * zb * ya * yb - 2.0 * za * zb * ya * ym - 2.0 * za * zb * yb * ym - 2.0 * ya * yb * za * zm -
        xa2 * zb2 - xa2 * yb2 - zb2 * ya2 - za2 * xb2 - za2 * yb2 - xb2 * ya2 - 2.0 * zb * zm * xa * xm -
        2.0 * xa * xb * za * zm - 2.0 * xa * xb * zb * zm - 2.0 * xa * xb * ya * ym - 2.0 * xa * xb * yb * ym +
        za2 * r2 - za2 * xm2 - za2 * ym2 + zb2 * r2 - zb2 * xm2 - zb2 * ym2 + xa2 * r2 - xa2 * zm2 - xa2 * ym2 +
        xb2 * r2 - xb2 * zm2 - xb2 * ym2 + ya2 * r2 - ya2 * zm2 - ya2 * xm2 + yb2 * r2 - yb2 * zm2 - yb2 * xm2;
    double nenner  = -2.0 * za * zb - 2.0 * ya * yb - 2.0 * xa * xb + za2 + zb2 + xa2 + xb2 + ya2 + yb2;
    double zaehler = 2.0 * zb * zm - 2.0 * xa * xm + 2.0 * yb * ym - 2.0 * za * zm + xa2 - 2.0 * ya * ym +
                     2.0 * xb * xm - zb2 + za2 - xb2 + ya2 - yb2;

    vector<GbPoint3D *> schnittpunkte;

    if (fabs(nenner) > 1.E-13 && UbMath::greaterEqual(wurzel, 0.0)) {
        double t1 = (zaehler + 2.0 * sqrt(wurzel)) / nenner;
        double t2 = (zaehler - 2.0 * sqrt(wurzel)) / nenner;

        if (UbMath::inClosedInterval(t1, -1.0, 1.0)) {
            double x = (xa * (0.5 - 0.5 * t1) + xb * (0.5 + 0.5 * t1)) / factor;
            double y = (ya * (0.5 - 0.5 * t1) + yb * (0.5 + 0.5 * t1)) / factor;
            double z = (za * (0.5 - 0.5 * t1) + zb * (0.5 + 0.5 * t1)) / factor;

            schnittpunkte.push_back(new GbPoint3D(x, y, z));
        }
        if (fabs(t2 - t1) > 1.E-13 && UbMath::inClosedInterval(t2, -1.0, 1.0)) {
            double x = (xa * (0.5 - 0.5 * t2) + xb * (0.5 + 0.5 * t2)) / factor;
            double y = (ya * (0.5 - 0.5 * t2) + yb * (0.5 + 0.5 * t2)) / factor;
            double z = (za * (0.5 - 0.5 * t2) + zb * (0.5 + 0.5 * t2)) / factor;

            schnittpunkte.push_back(new GbPoint3D(x, y, z));
        }
    }

    int nofSchnittpunkte = (int)schnittpunkte.size();
    if (nofSchnittpunkte == 1) {
        if (this->isPointInGbObject3D(&point1))
            return new GbLine3D(schnittpunkte[0], new GbPoint3D(point1));
        else if (this->isPointInGbObject3D(&point2))
            return new GbLine3D(schnittpunkte[0], new GbPoint3D(point2));
        else // line beruehrt kugel! -> clippedLine reduziert sich zu einem Punkt!
        {
            if (std::fabs(this->getDistance(schnittpunkte[0]) - this->radius) < 1.E-13)
                throw UbException(UB_EXARGS, "Beide LinenPunkte ausserhalb des Kreises, der berechnete Punkt ist "
                                             "jedoch KEIN Beruhrungspunkt der Sphere...");
            return new GbLine3D(schnittpunkte[0], new GbPoint3D(*(schnittpunkte[0])));
        }
    } else if (nofSchnittpunkte == 2)
        return new GbLine3D(schnittpunkte[0], schnittpunkte[1]);

    return NULL;
}
/*=========================================================================*/
vector<GbTriangle3D *> GbSphere3D::getSurfaceTriangleSet()
{
    if (triangulationMode == RAYPROJECTION) {
        double x1m = midPoint->getX1Coordinate();
        double x2m = midPoint->getX2Coordinate();
        double x3m = midPoint->getX3Coordinate();

        vector<GbTriangle3D *> triangles;

        int segments    = 30;
        double deltaPhi = UbMath::PI / (double)segments;
        double phiX1a, phiX1b, phiX3a, phiX3b;
        double x1a, x2a, x3a, x1b, x2b, x3b, x1c, x2c, x3c, x1d, x2d, x3d;

        for (phiX3a = 0.5 * UbMath::PI; phiX3a > -1.5 * UbMath::PI; phiX3a -= deltaPhi) {
            for (phiX1a = 0.0; phiX1a < UbMath::PI; phiX1a += deltaPhi) {
                phiX1b = phiX1a + deltaPhi;
                phiX3b = phiX3a + deltaPhi;

                x1a = x1m + radius * cos(phiX3a) * std::cos(phiX1a);
                x2a = x2m + radius * cos(phiX3a) * std::sin(phiX1a);
                x3a = x3m + radius * sin(phiX3a);
                x1b = x1m + radius * cos(phiX3a) * std::cos(phiX1b);
                x2b = x2m + radius * cos(phiX3a) * std::sin(phiX1b);
                x3b = x3m + radius * sin(phiX3a);
                x1c = x1m + radius * cos(phiX3b) * std::cos(phiX1b);
                x2c = x2m + radius * cos(phiX3b) * std::sin(phiX1b);
                x3c = x3m + radius * sin(phiX3b);
                x1d = x1m + radius * cos(phiX3b) * std::cos(phiX1a);
                x2d = x2m + radius * cos(phiX3b) * std::sin(phiX1a);
                x3d = x3m + radius * sin(phiX3b);

                if (UbMath::greater(phiX3b, -0.5 * UbMath::PI) && UbMath::less(phiX3a, 0.5 * UbMath::PI)) {
                    triangles.push_back(new GbTriangle3D(new GbPoint3D(x1a, x2a, x3a), new GbPoint3D(x1b, x2b, x3b),
                                                         new GbPoint3D(x1c, x2c, x3c)));
                    triangles.push_back(new GbTriangle3D(new GbPoint3D(x1a, x2a, x3a), new GbPoint3D(x1c, x2c, x3c),
                                                         new GbPoint3D(x1d, x2d, x3d)));
                } else {
                    triangles.push_back(new GbTriangle3D(new GbPoint3D(x1d, x2d, x3d), new GbPoint3D(x1c, x2c, x3c),
                                                         new GbPoint3D(x1a, x2a, x3a)));
                    triangles.push_back(new GbTriangle3D(new GbPoint3D(x1c, x2c, x3c), new GbPoint3D(x1b, x2b, x3b),
                                                         new GbPoint3D(x1a, x2a, x3a)));
                }
            }
        }
        return triangles;
    } else if (triangulationMode == CUBOIDPROJECTION) {
        vector<GbTriangle3D *> triangles;
        vector<GbPoint3D *> points;
        double x1min = this->getX1Minimum();
        double x2min = this->getX2Minimum();
        double x3min = this->getX3Minimum();
        double x1max = this->getX1Maximum();
        double x2max = this->getX2Maximum();
        double x3max = this->getX3Maximum();
        double ax1   = x1min;
        double bx2   = x2min;
        double cx1   = x1min;
        double ax2   = x2min;
        double bx3   = x3min;
        double cx3   = x3min;

        int anzahl = 20;
        double dx1 = (x1max - x1min) / (double)(anzahl - 1);
        double dx2 = (x2max - x2min) / (double)(anzahl - 1);
        double dx3 = (x3max - x3min) / (double)(anzahl - 1);

        for (int u = 0; u < anzahl; u++) {
            ax2 = x2min;
            bx2 = x2min;
            cx3 = x3min;
            for (int v = 0; v < anzahl; v++) {
                GbPoint3D p1 = GbPoint3D(ax1, ax2, x3max);
                GbPoint3D p2 = GbPoint3D(ax1, ax2, x3min);
                GbPoint3D p3 = GbPoint3D(cx1, x2min, cx3);
                GbPoint3D p4 = GbPoint3D(cx1, x2max, cx3);
                GbPoint3D p5 = GbPoint3D(x1min, bx2, bx3);
                GbPoint3D p6 = GbPoint3D(x1max, bx2, bx3);

                GbLine3D *clippedline1 = this->createClippedLine3D(*this->midPoint, p1);
                GbLine3D *clippedline2 = this->createClippedLine3D(*this->midPoint, p2);
                GbLine3D *clippedline3 = this->createClippedLine3D(*this->midPoint, p3);
                GbLine3D *clippedline4 = this->createClippedLine3D(*this->midPoint, p4);
                GbLine3D *clippedline5 = this->createClippedLine3D(*this->midPoint, p5);
                GbLine3D *clippedline6 = this->createClippedLine3D(*this->midPoint, p6);
                points.push_back(new GbPoint3D(clippedline1->getPoint1()));
                points.push_back(new GbPoint3D(clippedline2->getPoint1()));
                points.push_back(new GbPoint3D(clippedline3->getPoint1()));
                points.push_back(new GbPoint3D(clippedline4->getPoint1()));
                points.push_back(new GbPoint3D(clippedline5->getPoint1()));
                points.push_back(new GbPoint3D(clippedline6->getPoint1()));
                clippedline1->deletePoints();
                delete clippedline1;
                clippedline2->deletePoints();
                delete clippedline2;
                clippedline3->deletePoints();
                delete clippedline3;
                clippedline4->deletePoints();
                delete clippedline4;
                clippedline5->deletePoints();
                delete clippedline5;
                clippedline6->deletePoints();
                delete clippedline6;
                ax2 += dx2;
                cx3 += dx3;
                bx2 += dx2;
            }
            ax1 += dx1;
            cx1 += dx1;
            bx3 += dx3;
        }

        int anz           = anzahl * anzahl * 6;
        GbPoint3D *point1 = NULL;
        GbPoint3D *point2 = NULL;
        GbPoint3D *point3 = NULL;
        int anzahl2       = anzahl * 6;
        int anzahl3       = anzahl2 + 6;
        for (int u = 0; u < anz - anzahl3; u++) {
            point1 = new GbPoint3D(points[u + 6]);
            point2 = new GbPoint3D(points[u]);
            point3 = new GbPoint3D(points[u + anzahl2]);
            if (u % 2 == 0)
                triangles.push_back(new GbTriangle3D(point1, point2, point3));
            else
                triangles.push_back(new GbTriangle3D(point2, point1, point3));

            point1 = new GbPoint3D(points[u + 6]);
            point2 = new GbPoint3D(points[u + anzahl2]);
            point3 = new GbPoint3D(points[u + anzahl3]);
            if (u % 2 == 0)
                triangles.push_back(new GbTriangle3D(point1, point2, point3));
            else
                triangles.push_back(new GbTriangle3D(point2, point1, point3));
        }
        for (int u = 0; u < anz; u++)
            delete points[u];

        return triangles;
    } else
        throw UbException(UB_EXARGS, "undefined triangulationmode");
}
/*=======================================================*/
void GbSphere3D::addSurfaceTriangleSet(vector<UbTupleFloat3> &nodes, vector<UbTupleInt3> &triangles)
{
    // wenn ich viele Kugeln bei der PE rausschreibe sollten die vektoren nicht geresized werden
    // nodes.resize(0);
    // triangles.resize(0);

    if (triangulationMode == RAYPROJECTION) {
        float x1m = (float)midPoint->getX1Coordinate();
        float x2m = (float)midPoint->getX2Coordinate();
        float x3m = (float)midPoint->getX3Coordinate();

        int segments   = 30;
        float deltaPhi = (float)UbMath::PI / (float)segments;
        float phiX1a, phiX1b, phiX3a, phiX3b;
        float x1a, x2a, x3a, x1b, x2b, x3b, x1c, x2c, x3c, x1d, x2d, x3d;
        int nodeNr = int(nodes.size());
        for (phiX3a = (float)(0.5 * UbMath::PI); phiX3a > (float)(-1.5 * UbMath::PI); phiX3a -= deltaPhi) {
            for (phiX1a = 0.0; phiX1a < UbMath::PI; phiX1a += deltaPhi) {
                phiX1b = phiX1a + deltaPhi;
                phiX3b = phiX3a + deltaPhi;

                x1a = x1m + (float)(radius * cos(phiX3a) * std::cos(phiX1a));
                x2a = x2m + (float)(radius * cos(phiX3a) * std::sin(phiX1a));
                x3a = x3m + (float)(radius * sin(phiX3a));
                x1b = x1m + (float)(radius * cos(phiX3a) * std::cos(phiX1b));
                x2b = x2m + (float)(radius * cos(phiX3a) * std::sin(phiX1b));
                x3b = x3m + (float)(radius * sin(phiX3a));
                x1c = x1m + (float)(radius * cos(phiX3b) * std::cos(phiX1b));
                x2c = x2m + (float)(radius * cos(phiX3b) * std::sin(phiX1b));
                x3c = x3m + (float)(radius * sin(phiX3b));
                x1d = x1m + (float)(radius * cos(phiX3b) * std::cos(phiX1a));
                x2d = x2m + (float)(radius * cos(phiX3b) * std::sin(phiX1a));
                x3d = x3m + (float)(radius * sin(phiX3b));

                if (UbMath::greater(phiX3b, -0.5 * UbMath::PI) && UbMath::less(phiX3a, 0.5 * UbMath::PI)) {
                    nodes.push_back(makeUbTuple(x1a, x2a, x3a));
                    nodes.push_back(makeUbTuple(x1b, x2b, x3b));
                    nodes.push_back(makeUbTuple(x1c, x2c, x3c));

                    nodes.push_back(makeUbTuple(x1a, x2a, x3a));
                    nodes.push_back(makeUbTuple(x1c, x2c, x3c));
                    nodes.push_back(makeUbTuple(x1d, x2d, x3d));
                } else {
                    nodes.push_back(makeUbTuple(x1d, x2d, x3d));
                    nodes.push_back(makeUbTuple(x1c, x2c, x3c));
                    nodes.push_back(makeUbTuple(x1a, x2a, x3a));

                    nodes.push_back(makeUbTuple(x1c, x2c, x3c));
                    nodes.push_back(makeUbTuple(x1b, x2b, x3b));
                    nodes.push_back(makeUbTuple(x1a, x2a, x3a));
                }
                triangles.push_back(makeUbTuple(nodeNr, nodeNr + 1, nodeNr + 2));
                triangles.push_back(makeUbTuple(nodeNr + 3, nodeNr + 4, nodeNr + 5));
                nodeNr += 6;
            }
        }
    } else if (triangulationMode == CUBOIDPROJECTION) {
        vector<GbPoint3D *> points;
        double x1min = this->getX1Minimum();
        double x2min = this->getX2Minimum();
        double x3min = this->getX3Minimum();
        double x1max = this->getX1Maximum();
        double x2max = this->getX2Maximum();
        double x3max = this->getX3Maximum();
        double ax1   = x1min;
        double bx2   = x2min;
        double cx1   = x1min;
        double ax2   = x2min;
        double bx3   = x3min;
        double cx3   = x3min;

        int anzahl = 20;
        double dx1 = (x1max - x1min) / (double)(anzahl - 1);
        double dx2 = (x2max - x2min) / (double)(anzahl - 1);
        double dx3 = (x3max - x3min) / (double)(anzahl - 1);

        for (int u = 0; u < anzahl; u++) {
            ax2 = x2min;
            bx2 = x2min;
            cx3 = x3min;
            for (int v = 0; v < anzahl; v++) {
                GbPoint3D p1 = GbPoint3D(ax1, ax2, x3max);
                GbPoint3D p2 = GbPoint3D(ax1, ax2, x3min);
                GbPoint3D p3 = GbPoint3D(cx1, x2min, cx3);
                GbPoint3D p4 = GbPoint3D(cx1, x2max, cx3);
                GbPoint3D p5 = GbPoint3D(x1min, bx2, bx3);
                GbPoint3D p6 = GbPoint3D(x1max, bx2, bx3);

                GbLine3D *clippedline1 = this->createClippedLine3D(*this->midPoint, p1);
                GbLine3D *clippedline2 = this->createClippedLine3D(*this->midPoint, p2);
                GbLine3D *clippedline3 = this->createClippedLine3D(*this->midPoint, p3);
                GbLine3D *clippedline4 = this->createClippedLine3D(*this->midPoint, p4);
                GbLine3D *clippedline5 = this->createClippedLine3D(*this->midPoint, p5);
                GbLine3D *clippedline6 = this->createClippedLine3D(*this->midPoint, p6);
                points.push_back(new GbPoint3D(clippedline1->getPoint1()));
                points.push_back(new GbPoint3D(clippedline2->getPoint1()));
                points.push_back(new GbPoint3D(clippedline3->getPoint1()));
                points.push_back(new GbPoint3D(clippedline4->getPoint1()));
                points.push_back(new GbPoint3D(clippedline5->getPoint1()));
                points.push_back(new GbPoint3D(clippedline6->getPoint1()));
                clippedline1->deletePoints();
                delete clippedline1;
                clippedline2->deletePoints();
                delete clippedline2;
                clippedline3->deletePoints();
                delete clippedline3;
                clippedline4->deletePoints();
                delete clippedline4;
                clippedline5->deletePoints();
                delete clippedline5;
                clippedline6->deletePoints();
                delete clippedline6;
                ax2 += dx2;
                cx3 += dx3;
                bx2 += dx2;
            }
            ax1 += dx1;
            cx1 += dx1;
            bx3 += dx3;
        }

        int anz     = anzahl * anzahl * 6;
        int anzahl2 = anzahl * 6;
        int anzahl3 = anzahl2 + 6;
        int nodeNr  = 0;
        for (int u = 0; u < anz - anzahl3; u++) {
            nodes.push_back(makeUbTuple((float)points[u + 6]->x1, (float)points[u + 6]->x2, (float)points[u + 6]->x3));
            nodes.push_back(makeUbTuple((float)points[u]->x1, (float)points[u]->x2, (float)points[u]->x3));
            nodes.push_back(makeUbTuple((float)points[u + anzahl2]->x1, (float)points[u + anzahl2]->x2,
                                        (float)points[u + anzahl2]->x3));

            if (u % 2 == 0)
                triangles.push_back(makeUbTuple(nodeNr, nodeNr + 1, nodeNr + 2));
            else
                triangles.push_back(makeUbTuple(nodeNr, nodeNr + 1, nodeNr + 2));

            nodes.push_back(makeUbTuple((float)points[u + 6]->x1, (float)points[u + 6]->x2, (float)points[u + 6]->x3));
            nodes.push_back(makeUbTuple((float)points[u + anzahl2]->x1, (float)points[u + anzahl2]->x2,
                                        (float)points[u + anzahl2]->x3));
            nodes.push_back(makeUbTuple((float)points[u + anzahl3]->x1, (float)points[u + anzahl3]->x2,
                                        (float)points[u + anzahl3]->x3));
            if (u % 2 == 0)
                triangles.push_back(makeUbTuple(nodeNr + 3, nodeNr + 4, nodeNr + 5));
            else
                triangles.push_back(makeUbTuple(nodeNr + 3, nodeNr + 4, nodeNr + 5));

            nodeNr += 6;
        }
        for (int u = 0; u < anz; u++)
            delete points[u];
    } else
        throw UbException(UB_EXARGS, "undefined triangulationmode");
}
/*=======================================================*/
void GbSphere3D::transform(const double matrix[4][4])
{
    midPoint->transform(matrix);
    this->setRadius(this->getRadius() * matrix[0][0]);
    this->notifyObserversObjectChanged();
}
/*=======================================================*/
bool GbSphere3D::hasIntersectionWithDirectedLine(GbPoint3D origin, GbPoint3D direction)
{
    GbVector3D vecOrigin(origin.getX1Coordinate(), origin.getX2Coordinate(), origin.getX3Coordinate());
    GbVector3D vecDirection(direction.getX1Coordinate(), direction.getX2Coordinate(), direction.getX3Coordinate());
    GbVector3D vecSfereCenter(getX1Centroid(), getX2Centroid(), getX3Centroid());
    GbVector3D diff = vecOrigin - vecSfereCenter;
    float a         = (float)(vecDirection.Dot(vecDirection));
    float b         = (float)(2.0 * vecDirection.Dot(diff));
    float c         = (float)(diff.Dot(diff) - this->getRadius() * this->getRadius());

    // use 'abc'-formula for finding root t_1,2 = (-b +/- sqrt(b^2-4ac))/(2a)
    float inRoot = (float)(b * b - 4.0 * a * c);
    if (inRoot < 0)
        return false;
    float root = sqrt(inRoot);

    float dist = (float)((-b - root) / (2.0 * a));

    double infinity = DBL_MAX;
    double eps      = 1E-4;

    if (dist > infinity)
        return false;

    if (dist < eps) {
        dist = (float)((-b + root) / (2.0 * a));
        if (dist < eps || dist > infinity)
            return false;
    }
    return true;
}
/*=======================================================*/
bool GbSphere3D::isCellCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b)
// Merksatz: cell oder deren Volumen schneidet oder beinhaltet komplette oder Teile der SphereUmrandung
// returns true:
//  - cell cuts  sphere3D
//  - cell boxes sphere3D
// returns false:
//  - cell completely inside sphere3D ( = sphere3D boxes cell)
//  - cell und sphere3D haben kein gemeinsames Volumen
{
    double midX[] = { this->getX1Centroid(), this->getX2Centroid(), this->getX3Centroid() };

    double Bmin[] = { UbMath::min(x1a, x1b), UbMath::min(x2a, x2b), UbMath::min(x3a, x3b) };

    double Bmax[] = { UbMath::max(x1a, x1b), UbMath::max(x2a, x2b), UbMath::max(x3a, x3b) };

    /* Solid Box - Hollow Sphere */
    double dmin = 0.0;
    double dmax = 0.0;
    double r2   = radius * radius;

    for (int i = 0; i < 3; i++) {
        double a = pow(midX[i] - Bmin[i], 2.0);
        double b = pow(midX[i] - Bmax[i], 2.0);
        dmax += UbMath::max(a, b);
        if (UbMath::less(midX[i], Bmin[i]))
            dmin += a;
        else if (UbMath::greater(midX[i], Bmax[i]))
            dmin += b;
    }
    if (UbMath::lessEqual(dmin, r2) && UbMath::lessEqual(r2, dmax)) {
        return true;
    }
    return false;
}
/*=======================================================*/
bool GbSphere3D::isCellInsideOrCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a,
                                                 const double &x1b, const double &x2b, const double &x3b)
// returns true:
//  - cell completely inside sphere3D ( = sphere3D boxes cell)
//  - cell cuts  sphere3D
//  - cell boxes sphere3D
// returns false:
//  - cell und sphere3D haben kein gemeinsames Volumen
{
    // URL: http://tog.acm.org/GraphicsGems/gems/BoxSphere.c (mode=4, beides solids!!!)
    // solid - solid
    // this routine tests for intersection between an 3-dimensional
    // axis-aligned box and an 3-dimensional sphere.

    // true:
    //  - wenn Schnitt
    //  - Cell komplett innerhalb GbSphere3D
    //  - Cell umhuellt GbSphere3D

    double midX1 = this->getX1Centroid();
    double midX2 = this->getX2Centroid();
    double midX3 = this->getX3Centroid();

    double dmin = 0.0;

    if (UbMath::less(midX1, x1a))
        dmin += std::pow(midX1 - x1a, 2.0);
    else if (UbMath::greater(midX1, x1b))
        dmin += std::pow(midX1 - x1b, 2.0);

    if (UbMath::less(midX2, x2a))
        dmin += std::pow(midX2 - x2a, 2.0);
    else if (UbMath::greater(midX2, x2b))
        dmin += std::pow(midX2 - x2b, 2.0);

    if (UbMath::less(midX3, x3a))
        dmin += std::pow(midX3 - x3a, 2.0);
    else if (UbMath::greater(midX3, x3b))
        dmin += std::pow(midX3 - x3b, 2.0);

    if (UbMath::lessEqual(dmin, radius * radius)) {
        return true;
    }

    return false;
}
/*==========================================================*/
double GbSphere3D::getCellVolumeInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a,
                                                 const double &x1b, const double &x2b, const double &x3b)
{
    double deltaX1 = (x1b - x1a);
    double deltaX2 = (x2b - x2a);
    double deltaX3 = (x3b - x3a);

    if (this->isCellInsideGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b))
        return 1.0 * deltaX1 * deltaX2 * deltaX3;
    if (!(this->isCellCuttingGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b)))
        return 0.0;

    double tempResult = 0.0;

    int iMax = 10;
    int jMax = 10;
    int kMax = 10;

    for (int i = 0; i < iMax; i++) {
        for (int j = 0; j < jMax; j++) {
            for (int k = 0; k < kMax; k++) {

                tempResult += getCellVolumeInsideGbObject3DHelperFunction(
                    x1a + ((double)i) * deltaX1 / ((double)iMax), x2a + ((double)j) * deltaX2 / ((double)jMax),
                    x3a + ((double)k) * deltaX3 / ((double)kMax), x1a + ((double)(i + 1)) * deltaX1 / ((double)iMax),
                    x2a + ((double)(j + 1)) * deltaX2 / ((double)jMax),
                    x3a + ((double)(k + 1)) * deltaX3 / ((double)kMax));
            }
        }
    }

    // double resultWithOneCell = getCellVolumeInsideGbObject3DHelperFunction( x1a, x2a, x3a, x1b, x2b, x3b );
    // cout << tempResult << " vs. " << resultWithOneCell << endl;

    return tempResult;
}
/*==========================================================*/
double GbSphere3D::getCellVolumeInsideGbObject3DHelperFunction(const double &x1a, const double &x2a, const double &x3a,
                                                               const double &x1b, const double &x2b, const double &x3b)
{

    double deltaX1 = x1b - x1a;
    double deltaX2 = x2b - x2a;
    double deltaX3 = x3b - x3a;

    if (this->isCellInsideGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b))
        return 1.0 * deltaX1 * deltaX2 * deltaX3;
    if (!(this->isCellCuttingGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b)))
        return 0.0;

    double alpha = 0.0;
    double internX1, internX2, internX3;

    for (int x1vers = 0; x1vers < 2; x1vers++) {
        for (int x2vers = 0; x2vers < 2; x2vers++) {
            for (int x3vers = 0; x3vers < 2; x3vers++) {
                internX1 = x1a + (x1b - x1a) * x1vers;
                internX2 = x2a + (x2b - x2a) * x2vers;
                internX3 = x3a + (x3b - x3a) * x3vers;

                if (UbMath::lessEqual(this->getDistance(internX1, internX2, internX3), alpha))
                    alpha = this->getDistance(internX1, internX2, internX3);
                // cout<<zelltyp<<" "<<kugel->getDistance(internX1,internX2,internX3)<<" "<<alpha<<endl;
            } // end first for
        }     // end second for
    }         // end third for

    alpha = (-1) * alpha;

    double n[3];
    n[0] = 0.5 * (x1b + x1a) - this->getX1Centroid();
    n[1] = 0.5 * (x2b + x2a) - this->getX2Centroid();
    n[2] = 0.5 * (x3b + x3a) - this->getX3Centroid();

    // cout << "Koordinaten:  "<<x1<<" "<<x2<<" "<<x3<<endl;
    // cout << "Deltas:       "<<deltaX1<<" "<<deltaX2<<" "<<deltaX3<<endl;
    // cout << "Halbe Zelle:  "<<halfcelldelta<<endl;

    // cout<<"Centroid:  "<<kugel->getX1Centroid()<<" "<<kugel->getX2Centroid()<<" "<<kugel->getX3Centroid()<<endl;

    // cout<<"Normals: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;

    double normLength;
    normLength = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] /= normLength;
    n[1] /= normLength;
    n[2] /= normLength;

    if (UbMath::less(n[0], 0.0))
        n[0] = -n[0];
    if (UbMath::less(n[1], 0.0))
        n[1] = -n[1];
    if (UbMath::less(n[2], 0.0))
        n[2] = -n[2];

    // cout<<"Normals: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;

    double dummy;
    if (UbMath::greater(n[0], n[1])) {
        dummy = n[1];
        n[1]  = n[0];
        n[0]  = dummy;
    }
    if (UbMath::greater(n[1], n[2])) {
        dummy = n[2];
        n[2]  = n[1];
        n[1]  = dummy;
    }
    if (UbMath::greater(n[0], n[1])) {
        dummy = n[1];
        n[1]  = n[0];
        n[0]  = dummy;
    }

    // cout<<"Normals: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;

    double n1, n2, n3;
    n1 = n[0];
    n2 = n[1];
    n3 = n[2];

    double maxVol = deltaX1 * deltaX2 * deltaX3;

    double result = 0.0, preresult = 0.0;

    if (UbMath::lessEqual(maxVol, 0.000001))
        return 0.0;

    // 1D Check
    if (UbMath::lessEqual(n1, 0.001) && UbMath::lessEqual(n2, 0.001)) {
        result = alpha * deltaX1 * deltaX2;
    }
    // 2D Check
    else if (UbMath::lessEqual(n1, 0.001)) {
        preresult = (2 * n2 * n3);
        result    = (alpha * alpha) / preresult;

        if (UbMath::greater(alpha, n2 * deltaX2)) {
            result += -(alpha - n2 * deltaX2) * (alpha - n2 * deltaX2) / preresult;
        }
        if (UbMath::greater(alpha, n3 * deltaX3)) {
            result += -(alpha - n3 * deltaX3) * (alpha - n3 * deltaX3) / preresult;
        }
        if (UbMath::greater(alpha, n2 * deltaX2 + n3 * deltaX3)) {
            result += (alpha - n2 * deltaX2 - n3 * deltaX3) * (alpha - n2 * deltaX2 - n3 * deltaX3) / preresult;
        }

        // tiefenrichtung mit einmultiplizieren...
        result *= deltaX1;
    }
    // 3D Check
    else {
        preresult = 6 * n1 * n2 * n3;

        result = alpha * alpha * alpha;

        if (UbMath::greaterEqual(alpha, n1 * deltaX1)) {
            result += -((alpha - n1 * deltaX1) * (alpha - n1 * deltaX1) * (alpha - n1 * deltaX1));
        }
        if (UbMath::greaterEqual(alpha, n2 * deltaX2)) {
            result += -((alpha - n2 * deltaX2) * (alpha - n2 * deltaX2) * (alpha - n2 * deltaX2));
        }
        if (UbMath::greaterEqual(alpha, n3 * deltaX3)) {
            result += -((alpha - n3 * deltaX3) * (alpha - n3 * deltaX3) * (alpha - n3 * deltaX3));
        }
        if (UbMath::greaterEqual(alpha, (n1 * deltaX1 + n2 * deltaX2))) {
            result += ((alpha - (n1 * deltaX1 + n2 * deltaX2)) * (alpha - (n1 * deltaX1 + n2 * deltaX2)) *
                       (alpha - (n1 * deltaX1 + n2 * deltaX2)));
        }
        if (UbMath::greaterEqual(alpha, (n1 * deltaX1 + n3 * deltaX3))) {
            result += ((alpha - (n1 * deltaX1 + n3 * deltaX3)) * (alpha - (n1 * deltaX1 + n3 * deltaX3)) *
                       (alpha - (n1 * deltaX1 + n3 * deltaX3)));
        }
        if (UbMath::greaterEqual(alpha, (n2 * deltaX2 + n3 * deltaX3))) {
            result += ((alpha - (n2 * deltaX2 + n3 * deltaX3)) * (alpha - (n2 * deltaX2 + n3 * deltaX3)) *
                       (alpha - (n2 * deltaX2 + n3 * deltaX3)));
        }

        // NEW
        if (UbMath::greaterEqual(alpha, (n1 * deltaX1 + n2 * deltaX2 + n3 * deltaX3))) {
            result += -((alpha - (n1 * deltaX1 + n2 * deltaX2 + n3 * deltaX3)) *
                        (alpha - (n1 * deltaX1 + n2 * deltaX2 + n3 * deltaX3)) *
                        (alpha - (n1 * deltaX1 + n2 * deltaX2 + n3 * deltaX3)));
        }

        result = result / preresult;
    }
    return (result);

    // cout << "alpha ist " << alpha << endl;
    // cout << "fillLevel ist " << eps << endl;
}
/*==========================================================*/
double GbSphere3D::getIntersectionRaytraceFactor(const double &x1, const double &x2, const double &x3,
                                                 const double &rx1, const double &rx2, const double &rx3)
{
    double lx1  = midPoint->x1 - x1;
    double lx2  = midPoint->x2 - x2;
    double lx3  = midPoint->x3 - x3;
    double l_sq = lx1 * lx1 + lx2 * lx2 + lx3 * lx3; // l = abstand Punkt(x1,x2,x3)<->kreismittelpunkt

    double s    = lx1 * rx1 + lx2 * rx2 + lx3 * rx3; // s= l*ray_dir)
    double r_sq = this->radius * this->radius;       // r� =r*r
    // if (d<0 (fuer die Richtung falls sie gegen das Kreis dann haben wir ein negativer Zahl)
    //     && l� > r� (point outside ))
    // wenn s<0->Punkt liegt rechts vom mittelpunkt, wenn nun punkt ausserhalb des kreises liegt, kann es keinen SP mehr
    // geben
    if (s < -1.E-10 && l_sq > r_sq + 1.E-10)
        return -1.0;
    // Pythagor on Triangle Rectangle (point, center of the cercle, intersection of the direction on point and m)
    // l� = m� + d�
    double m_sq = l_sq - s * s;
    // if (m� > r� (dann gibt es kein schnittpunt zwischen direction und circle))
    if (m_sq > r_sq + 1.E-10)
        return -1.0;
    // Pythagoras on Triangle Rectangle in cercle (direction , m, r)
    // r� = m� + h�

    // patch: rundungsfehler bei kleinen delta!!!
    //-> wenn wurzel minimal null->
    double wurzelTerm = r_sq - m_sq;
    if (wurzelTerm < 0.0) {
        if (wurzelTerm < -1E-10)
            return -1.0; // definitiv kein SP
        else
            return s; // im rundungsfehler-bereich. SP liegt dierkt auf sphere umrandung
    }

    // if point outside of the circle
    if (l_sq > r_sq)
        return s - sqrt(wurzelTerm);

    return s + sqrt(wurzelTerm);
}
/*=======================================================*/

//! \}
