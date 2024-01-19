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

using namespace std;

/*=======================================================*/
GbLine3D::GbLine3D()
{
    p1     = NULL;
    p2     = NULL;
    length = 0.0;
}
/*=======================================================*/
GbLine3D::GbLine3D(GbPoint3D *point1, GbPoint3D *point2)
{
    this->p1 = point1;
    this->p2 = point2;
    this->p1->addObserver(this);
    this->p2->addObserver(this);
    this->calculateValues();
}
/*=======================================================*/
GbLine3D::GbLine3D(GbLine3D *line)
{
    this->p1 = line->p1->clone();
    this->p2 = line->p2->clone();
    this->p1->addObserver(this);
    this->p2->addObserver(this);
    this->calculateValues();
}
/*=======================================================*/
GbLine3D::~GbLine3D()
{
    if (this->p1)
        this->p1->removeObserver(this);
    if (this->p2)
        this->p2->removeObserver(this);
}
/*=======================================================*/
void GbLine3D::finalize()
{
    if (this->p1) {
        this->p1->removeObserver(this);
        this->p1->finalize();
        delete this->p1;
        this->p1 = NULL;
    }
    if (this->p2) {
        this->p2->removeObserver(this);
        this->p2->finalize();
        delete this->p2;
        this->p2 = NULL;
    }
}
/*=======================================================*/
vector<GbTriangle3D *> GbLine3D::getSurfaceTriangleSet()
{
    vector<GbTriangle3D *> triangles;
    GbPoint3D p1(getX1Minimum(), getX2Minimum(), getX3Minimum());
    GbPoint3D p2(getX1Centroid(), getX2Centroid(), getX3Centroid());
    GbPoint3D p3(getX1Maximum(), getX2Maximum(), getX3Maximum());

    triangles.push_back(new GbTriangle3D(new GbPoint3D(p1), new GbPoint3D(p2), new GbPoint3D(p3)));

    return triangles;
}
/*=======================================================*/
void GbLine3D::setPoint1(GbPoint3D *point1)
{
    if (this->p1)
        this->p1->removeObserver(this);
    this->p1 = point1;
    this->p1->addObserver(this);

    if (this->p1 && this->p2)
        this->calculateValues();
}
/*=======================================================*/
void GbLine3D::setPoint2(GbPoint3D *point2)
{
    if (this->p2)
        this->p2->removeObserver(this);
    this->p2 = point2;
    this->p2->addObserver(this);

    if (this->p1 && this->p2)
        this->calculateValues();
}
/*=======================================================*/
void GbLine3D::setPoints(GbPoint3D *point1, GbPoint3D *point2)
{
    if (this->p1)
        this->p1->removeObserver(this);
    if (this->p2)
        this->p2->removeObserver(this);

    this->p1 = point1;
    this->p2 = point2;

    this->p1->addObserver(this);
    this->p2->addObserver(this);

    this->calculateValues();
}
/*=======================================================*/
string GbLine3D::toString()
{
    stringstream ss;
    ss << "GbLine3D[p1=";
    ss << this->p1->toString() << ",p2=" << this->p2->toString() << ",l=" << this->getLength() << "]";
    return (ss.str());
}
/*=======================================================*/
GbPoint3D *GbLine3D::calculateIntersectionPoint3D(GbLine3D * /*line*/)
{
    throw UbException(UB_EXARGS, " not implemented");
    // return(GbSystem::calculateIntersectionPoint3D(*this->p1, *this->p2, *line->p1, *line->p2));
}
/*======================================================================*/
GbLine3D *GbLine3D::createClippedLine3D(GbCuboid3D *cuboid)
{
    return GbSystem3D::createClipLine3D(*this->p1, *this->p2, cuboid->getPoint1()->x1, cuboid->getPoint1()->x2,
                                        cuboid->getPoint1()->x3, cuboid->getPoint2()->x1, cuboid->getPoint2()->x2,
                                        cuboid->getPoint2()->x3);
}
/*======================================================================*/
GbLine3D *GbLine3D::createClippedLine3D(GbPoint3D *pA, GbPoint3D *pE)
{
    return GbSystem3D::createClipLine3D(*this->p1, *this->p2, pA->x1, pA->x2, pA->x3, pE->x1, pE->x2, pE->x3);
}
/*======================================================================*/
double GbLine3D::getDistance(const GbPoint3D &point) { return this->getDistance(point.x1, point.x2, point.x3); }
/*======================================================================*/
double GbLine3D::getDistance(const double &x1, const double &x2, const double &x3)
{
    double dx1 = this->p2->x1 - this->p1->x1;
    double dx2 = this->p2->x2 - this->p1->x2;
    double dx3 = this->p2->x3 - this->p1->x3;

    // double vec[3];
    double a0 = x1 - p1->x1;
    double a1 = x2 - p1->x2;
    double a2 = x3 - p1->x3;

    double kreuzProd0 = a1 * dx3 - a2 * dx2;
    double kreuzProd1 = a2 * dx1 - a0 * dx3;
    double kreuzProd2 = a0 * dx2 - a1 * dx1;

    return (std::sqrt(kreuzProd0 * kreuzProd0 + kreuzProd1 * kreuzProd1 + kreuzProd2 * kreuzProd2)) / length;
}
/*=======================================================*/
void GbLine3D::calculateValues()
{
    double dx1   = this->p2->x1 - this->p1->x1;
    double dx2   = this->p2->x2 - this->p1->x2;
    double dx3   = this->p2->x3 - this->p1->x3;
    this->length = std::sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
}
/*==========================================================*/
void GbLine3D::objectChanged(UbObservable *changedObject)
{
    GbPoint3D *point = dynamic_cast<GbPoint3D *>(changedObject);
    if (!point || (this->p1 != point && this->p2 != point))
        return;

    this->calculateValues();
}
/*==========================================================*/
void GbLine3D::objectWillBeDeleted(UbObservable *objectForDeletion)
{
    if (this->p1) {
        UbObservable *observedObj = dynamic_cast<UbObservable *>(this->p1);
        if (objectForDeletion == observedObj) {
            this->p1 = NULL;
            length   = 0.0;
        }
    }
    if (this->p2) {
        UbObservable *observedObj = dynamic_cast<UbObservable *>(this->p2);
        if (objectForDeletion == observedObj) {
            this->p2 = NULL;
            length   = 0.0;
        }
    }
    // ACHTUNG: eigentlich muessten in allen methoden von GbLine if abfragen fuer NULL pointer hin... toDo
}
/*==========================================================*/
void GbLine3D::scale(const double &sx1, const double &sx2, const double &sx3)
{
    double p1X1 = this->p1->getX1Coordinate();
    double p1X2 = this->p1->getX2Coordinate();
    double p1X3 = this->p1->getX3Coordinate();

    double p2X1 = this->p2->getX1Coordinate();
    double p2X2 = this->p2->getX2Coordinate();
    double p2X3 = this->p2->getX3Coordinate();

    double lenX1 = fabs(p1X1 - p2X1);
    double lenX2 = fabs(p1X2 - p2X2);
    double lenX3 = fabs(p1X3 - p2X3);

    double deltaX1 = lenX1 * sx1 - lenX1;
    double deltaX2 = lenX2 * sx2 - lenX2;
    double deltaX3 = lenX3 * sx3 - lenX3;

    if (p1X1 < p2X1) {
        p1X1 -= 0.5 * deltaX1;
        p2X1 += 0.5 * deltaX1;
    } else {
        p1X1 += 0.5 * deltaX1;
        p2X1 -= 0.5 * deltaX1;
    }
    if (p1X2 < p2X2) {
        p1X2 -= 0.5 * deltaX2;
        p2X2 += 0.5 * deltaX2;
    } else {
        p1X2 += 0.5 * deltaX2;
        p2X2 -= 0.5 * deltaX2;
    }
    if (p1X3 < p2X3) {
        p1X3 -= 0.5 * deltaX3;
        p2X3 += 0.5 * deltaX3;
    } else {
        p1X3 += 0.5 * deltaX3;
        p2X3 -= 0.5 * deltaX3;
    }

    this->p1->setCoordinates(p1X1, p1X2, p1X3);
    this->p2->setCoordinates(p2X1, p2X2, p2X3);
}
/*=======================================================*/
void GbLine3D::translate(const double &tx1, const double &tx2, const double &tx3)
{
    this->p1->translate(tx1, tx2, tx3);
    this->p2->translate(tx1, tx2, tx3);
    // this->notifyObserversObjectChanged();
}

//! \}
