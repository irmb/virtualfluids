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
#include <GbPolygon3D.h>

using namespace std;

int GbPolygon3D::counter = 0;

GbPolygon3D::GbPolygon3D()
{
    init();
    counter++;
    this->ps = new GbSystem3D::PointSet3(0);
}
void GbPolygon3D::init()
{
    x1s   = 0.0;
    x2s   = 0.0;
    x3s   = 0.0;
    x1min = 0.0;
    x1max = 0.0;
    x2min = 0.0;
    x2max = 0.0;
    x3min = 0.0;
    x3max = 0.0;
    //        points   = NULL;
    consistent = false;
    ps         = NULL;
}

/*!
 * Creates an empty 3D polygon with the specified capacity.
 * @param capacity the initial capacity
 */
GbPolygon3D::GbPolygon3D(int capacity)
{
    init();
    counter++;
    this->ps = new GbSystem3D::PointSet3(capacity);
    //     this.po = new PointObserver(this);
}
/**
 * Creates a 3D polygon with the specified points.
 * @param points the initial points of the polygon
 */
GbPolygon3D::GbPolygon3D(vector<GbPoint3D> &points)
{
    init();
    counter++;
    this->ps = new GbSystem3D::PointSet3((int)points.size());
    this->addPoints(points);
}
/**
 * Creates a 3D polygon as clone of the specified 3D polygon.
 * @param polygon the 3D polygon to be cloned
 */
GbPolygon3D::GbPolygon3D(GbPolygon3D *polygon)
{
    this->init();
    counter++;
    this->ps               = new GbSystem3D::PointSet3((int)polygon->size());
    vector<GbPoint3D> temp = polygon->getPoints();
    this->addPoints(temp);
}

GbPolygon3D::~GbPolygon3D()
{
    counter--;
    // if(points)
    // for(unsigned u=0; u<points->size(); u++)
    //{
    //    delete (*points)[u];
    //}
    //        delete this->points;
    delete this->ps;
}

/*======================================================================*/
/**
 * Returns the number of points.
 * @return the number of points
 */
int GbPolygon3D::size() { return (this->ps->size()); }
/**
 * Returns the number of times this 3D polygon contains the specified point.
 * @param point the point
 * @return the number of times this 3D polygon contains the specified point
 */
int GbPolygon3D::contains(GbPoint3D *point) { return (this->ps->contains(point)); }
/**
 * Returns the number of times this 3D polygon contains a point equal to the specified point.
 * @param point the point
 * @return the number of times this 3D polygon contains a point equal to the specified point
 */
int GbPolygon3D::containsEqual(GbPoint3D *point) { return (this->ps->containsEqual(point)); }
/**
 * Returns true, if this 3D polygon contains the specified line.
 * @param point1 the first point
 * @param point2 the second point
 * @return true, if this 3D polygon contains the specified line
 */
bool GbPolygon3D::containsLine(GbPoint3D *point1, GbPoint3D *point2)
{
    return (this->ps->containsLine(point1, point2));
}
/**
 * Returns true, if this 3D polygon contains the specified line.
 * @param line the line
 * @return true, if this 3D polygon contains the specified line
 */
bool GbPolygon3D::containsLine(GbLine3D *line)
{
    return (this->ps->containsLine(line->getPoint1(), line->getPoint2()));
}
/**
 * Returns the first point.
 * @return the first point
 */
GbPoint3D *GbPolygon3D::getFirstPoint() { return (this->ps->getFirstPoint()); }
/**
 * Returns the last point.
 * @return the last point
 */
GbPoint3D *GbPolygon3D::getLastPoint() { return (this->ps->getLastPoint()); }
/**
 * Returns the specified point.
 * @param index the index
 * @return the specified point
 * @exception ArrayIndexOutOfBoundsException if the specified index is not valid
 */
GbPoint3D *GbPolygon3D::getPoint(const int &index)
{
    if (index < 0 || index > this->ps->size())
        throw UbException(UB_EXARGS, "ArrayIndexOutOfBoundsException-GbPolygon3D.getPoint()");
    return (this->ps->getPoint(index));
}
/**
 * Returns the points.
 * @return the points
 */
vector<GbPoint3D> GbPolygon3D::getPoints()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->points);
}
/**
 * Returns the points within the specified rectangle.
 * @param p1 the 1st point of the rectangle
 * @param p2 the 2nd point of the rectangle
 * @return the points within the specified rectangle
 */
vector<GbPoint3D> GbPolygon3D::getPoints(GbPoint3D *p1, GbPoint3D *p2)
{
    return (this->getPoints(p1->x1, p1->x2, p1->x3, p2->x1, p2->x2, p2->x3));
}
/**
 * Returns the points within the specified rectangle.
 * @param p1x1 the 1st x1 coordinate of the rectangle
 * @param p1x2 the 1st x2 coordinate of the rectangle
 * @param p1x3 the 1st x3 coordinate of the rectangle
 * @param p2x1 the 2nd x1 coordinate of the rectangle
 * @param p2x2 the 2nd x2 coordinate of the rectangle
 * @param p2x3 the 2nd x3 coordinate of the rectangle
 * @return the points within the specified rectangle
 */
vector<GbPoint3D> GbPolygon3D::getPoints(const double &p1x1, const double &p1x2, const double &p1x3, const double &p2x1,
                                         const double &p2x2, const double &p2x3)
{
    double x1min, x1max, x2min, x2max, x3min, x3max;

    if (UbMath::less(p1x1, p2x1)) {
        x1min = p1x1;
        x1max = p2x1;
    } else {
        x1min = p2x1;
        x1max = p1x1;
    }
    if (UbMath::less(p1x2, p2x2)) {
        x2min = p1x2;
        x2max = p2x2;
    } else {
        x2min = p2x2;
        x2max = p1x2;
    }
    if (UbMath::less(p1x3, p2x3)) {
        x3min = p1x3;
        x3max = p2x3;
    } else {
        x3min = p2x3;
        x3max = p1x3;
    }

    GbSystem3D::PointSet3 *pts = new GbSystem3D::PointSet3(1);

    if (!this->consistent)
        this->calculateValues();
    for (int i = this->size() - 1; i >= 0; i--) {
        if (UbMath::lessEqual(x1min, (this->points)[i].x1) && UbMath::greaterEqual(x1max, (this->points)[i].x1) &&
            UbMath::lessEqual(x2min, (this->points)[i].x2) && UbMath::greaterEqual(x2max, (this->points)[i].x2) &&
            UbMath::lessEqual(x3min, (this->points)[i].x3) && UbMath::greaterEqual(x3max, (this->points)[i].x3))
            pts->add((this->points)[i]);
    }
    return (pts->getPoints());
}
/**
 * Returns the area of this polygon.
 * The area is positive for positive ordered points, otherwise negative.
 * @return the area of this polygon
 */
// double getArea()
//{
//   if(!this.consistent) this.calculateValues();
//   return(this.area);
//}
double GbPolygon3D::getX1Centroid()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x1s);
}
double GbPolygon3D::getX1Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x1min);
}
double GbPolygon3D::getX1Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x1max);
}
double GbPolygon3D::getX2Centroid()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x2s);
}
double GbPolygon3D::getX2Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x2min);
}
double GbPolygon3D::getX2Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x2max);
}
double GbPolygon3D::getX3Centroid()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x3s);
}
double GbPolygon3D::getX3Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x3min);
}
double GbPolygon3D::getX3Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x3max);
}

/**
 * Adds a point to the end of this polygon. Notifies the observers of this 3D polygon.
 * @param point the point
 */
void GbPolygon3D::addPoint(GbPoint3D *point)
{
    // if((this instanceof GbPolygon3D) && !(point instanceof GbPoint3D)) throw new
    // IllegalArgumentException("GbPolygon3D.addPoint(): points of 3D polygons have to be 3D points!");

    this->ps->add(point);
    // point.addObserver(this.po);
    this->consistent = false;
    // super.notifyObservers();
}
/**
 * Adds a number of points to the end of this polygon. Notifies the observers of this 3D polygon.
 * @param points the points
 */
void GbPolygon3D::addPoints(vector<GbPoint3D> &points)
{
    // if((this instanceof GbPolygon3D) && (points.getClass().getComponentType() != GbPoint3D.class)) throw new
    // IllegalArgumentException("GbPolygon3D.addPoints(): points of 3D polygons have to be 3D points!");

    this->ps->add(points);
    // for(int i=0; i<points.length; i++) points[i].addObserver(this.po);
    this->consistent = false;
    // super.notifyObservers();
}
/**
 * Removes all points from this polygon. Notifies the observers of this 3D polygon.
 */
void GbPolygon3D::clear()
{
    //        delete this->points;
    this->ps->clearAndTrim();
    delete this->ps;

    // for(int i=points.length-1; i>=0; i--) points[i].removeObserver(this.po);
    this->consistent = false;
    // super.notifyObservers();
}
/**
 * Returns a string representation of this 3D polygon.
 * @return a string representation of this 3D polygon
 */
string GbPolygon3D::toString()
{
    stringstream ss;
    ss << "GbPolygon3D[";
    ss << this->size() << " points";
    ss << "]" << endl;
    for (int u = 0; u < this->size(); u++)
        ss << this->ps->getPoint(u)->toString() << endl;

    return (ss.str());
}
/*======================================================================*/

void GbPolygon3D::calculateValues()
{
    this->x1s   = 0.0;
    this->x2s   = 0.0;
    this->x3s   = 0.0;
    this->x1min = 0.0;
    this->x1max = 0.0;
    this->x2min = 0.0;
    this->x2max = 0.0;
    this->x3min = 0.0;
    this->x3max = 0.0;
    throw UbException(UB_EXARGS, "should be implemented");
}
/*======================================================================*/

//! \}
