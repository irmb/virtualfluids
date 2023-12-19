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
#ifndef GBPOLYGON3D_H
#define GBPOLYGON3D_H

#include <iostream>
#include <sstream>

#include <GbLine3D.h>
#include <GbObject3D.h>
#include <GbSystem3D.h>
#include <GbTriangle3D.h>

#include <PointerDefinitions.h>

/*=========================================================================*/
//! \class GbPolygon2D
/*                                                                         */
//! \brief This Class provides basic 3D polygon objects.

class GbPolygon3D : public GbObject3D
{
public:
    using GbObject3D::isPointInGbObject3D; // Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht
                                           // ausprogrammieren, welche sonst hier "ueberdeckt" waere
private:
    /*======================================================================*/
    double x1s;
    double x2s;
    double x3s;
    double x1min;
    double x1max;
    double x2min;
    double x2max;
    double x3min;
    double x3max;

    std::vector<GbPoint3D> points;
    bool consistent;

    GbSystem3D::PointSet3 *ps;
    // private PointObserver     po         = null;

    void init();

    /*======================================================================*/

    /*======================================================================*/
    /*  Construcrors                                                       */
    /*                                                                      */
    /*
     * Creates an empty 2D polygon.
     */
public:
    static int counter;
    GbPolygon3D();
    /*
     * Creates an empty 2D polygon with the specified capacity.
     * @param capacity the initial capacity
     */
    GbPolygon3D(int capacity);
    /*
     * Creates a 2D polygon with the specified points.
     * @param points the initial points of the polygon
     */
    GbPolygon3D(std::vector<GbPoint3D> &points);
    /*
     * Creates a 2D polygon as clone of the specified 2D polygon.
     * @param polygon the 2D polygon to be cloned
     */
    GbPolygon3D(GbPolygon3D *polygon);

    ~GbPolygon3D() override;

    /*======================================================================*/

    /*======================================================================*/
    /*  Methoden                                                            */
    /*                                                                      */
    /*
     * Creates a 2D polygon as clone of this 2D polygon.
     */
    GbPolygon3D *clone() override { return (new GbPolygon3D(this)); }
    void finalize() override { throw UbException(UB_EXARGS, "toDo"); }

    /*
     * Returns the number of points.
     * @return the number of points
     */
    int size();
    /*
     * Returns the number of times this 2D polygon contains the specified point.
     * @param point the point
     * @return the number of times this 2D polygon contains the specified point
     */
    int contains(GbPoint3D *point);
    /*
     * Returns the number of times this 2D polygon contains a point equal to the specified point.
     * @param point the point
     * @return the number of times this 2D polygon contains a point equal to the specified point
     */
    int containsEqual(GbPoint3D *point);
    /*
     * Returns true, if this 2D polygon contains the specified line.
     * @param point1 the first point
     * @param point2 the second point
     * @return true, if this 2D polygon contains the specified line
     */
    bool containsLine(GbPoint3D *point1, GbPoint3D *point2);
    /*
     * Returns true, if this 2D polygon contains the specified line.
     * @param line the line
     * @return true, if this 2D polygon contains the specified line
     */
    bool containsLine(GbLine3D *line);
    /*
     * Returns the first point.
     * @return the first point
     */
    GbPoint3D *getFirstPoint();
    /*
     * Returns the last point.
     * @return the last point
     */
    GbPoint3D *getLastPoint();
    /*
     * Returns the specified point.
     * @param index the index
     * @return the specified point
     * @exception ArrayIndexOutOfBoundsException if the specified index is not valid
     */
    GbPoint3D *getPoint(const int &index);
    /*
     * Returns the points.
     * @return the points
     */
    std::vector<GbPoint3D> getPoints();
    /*
     * Returns the points within the specified rectangle.
     * @param p1 the 1st point of the rectangle
     * @param p2 the 2nd point of the rectangle
     * @return the points within the specified rectangle
     */
    std::vector<GbPoint3D> getPoints(GbPoint3D *p1, GbPoint3D *p2);
    /*
     * Returns the points within the specified rectangle.
     * @param p1x1 the 1st x1 coordinate of the rectangle
     * @param p1x2 the 1st x2 coordinate of the rectangle
     * @param p2x1 the 2nd x1 coordinate of the rectangle
     * @param p2x2 the 2nd x2 coordinate of the rectangle
     * @return the points within the specified rectangle
     */
    std::vector<GbPoint3D> getPoints(const double &p1x1, const double &p1x2, const double &p1x3, const double &p2x1,
                                     const double &p2x2, const double &p2x3);
    /*
     * Returns the area of this polygon.
     * The area is positive for positive ordered points, otherwise negative.
     * @return the area of this polygon
     */
    // double getArea()
    //{
    //   if(!this.consistent) this.calculateValues();
    //   return(this.area);
    //}
    double getX1Centroid() override;
    double getX1Minimum() override;
    double getX1Maximum() override;
    double getX2Centroid() override;
    double getX2Minimum() override;
    double getX2Maximum() override;
    double getX3Centroid() override;
    double getX3Minimum() override;
    double getX3Maximum() override;

    /*
     * Adds a point to the end of this polygon. Notifies the observers of this 2D polygon.
     * @param point the point
     */
    void addPoint(GbPoint3D *point);
    /*
     * Adds a number of points to the end of this polygon. Notifies the observers of this 2D polygon.
     * @param points the points
     */
    void addPoints(std::vector<GbPoint3D> &points);
    /*
     * Removes all points from this polygon. Notifies the observers of this 2D polygon.
     */
    void clear();

    /*
     * Returns true if this 2D polygon equals the specified object.
     * Two polygon are equal, if their points are equal.
     * <BR>Note that the order of points is recognized!
     * @return true if this 2D polygon equals the specified object
     * @see GbPoint2D#equals(java.lang.Object)
     * @see GbPoint3D#equals(java.lang.Object)
     */
    // bool equals(Object object)
    // {
    //    try
    //    {
    //    GbPolygon2D polygon = (GbPolygon2D) object;
    // int         n       = this.size();

    // if(n != polygon.size()) return(false);
    // for(int i=0; i<n; i++) if(!this.getPoint(i).equals(polygon.getPoint(i))) return(false);
    // return(true);
    //    }
    //    catch(Exception e){ return(false); }
    // }
    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override
    {
        std::cout << "GbPolygon3D::getSurfaceTriangleSet() - not implemented\n";
        std::vector<GbTriangle3D *> tmp;
        return tmp;
    }
    bool isPointInGbObject3D(const double & /*x1*/, const double & /*x2*/, const double & /*x3*/) override
    {
        throw UbException(__FILE__, __LINE__, "GbPolygon3D::isPointInObject3D- not implemented");
    }
    bool isPointInGbObject3D(const double & /*x1*/, const double & /*x2*/, const double & /*x3*/,
                             bool & /*pointIsOnBoundary*/) override
    {
        throw UbException(__FILE__, __LINE__, "GbPolygon3D::isPointInObject3D- not implemented");
    }
    bool isCellInsideGbObject3D(const double & /*x1a*/, const double & /*x2a*/, const double & /*x3a*/,
                                const double & /*x1b*/, const double & /*x2b*/, const double & /*x3b*/) override
    {
        return false;
    }

    GbLine3D *createClippedLine3D(GbPoint3D & /*point1*/, GbPoint3D & /*point2*/) override
    {
        throw UbException(__FILE__, __LINE__, "GbPolygon3D::createClippedLine3D - not implemented");
    }

    /*
     * Returns a string representation of this 2D polygon.
     * @return a string representation of this 2D polygon
     */
    std::string toString() override;

    /*======================================================================*/
    /*  Private Methoden                                                    */
    /*                                                                      */
    void calculateValues();
    /*======================================================================*/
};
/*=========================================================================*/
#endif

//! \}
