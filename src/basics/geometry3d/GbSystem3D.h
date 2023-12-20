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
#ifndef GBSYSTEM3D_H
#define GBSYSTEM3D_H

#include <cmath>
#include <iostream>
#include <vector>

#include <GbObject3D.h>
#include <GbPoint3D.h>
#include <basics/utilities/UbMath.h>
#include <basics/writer/WbWriter.h>

class GbPolygon3D;
class GbCuboid3D;
class GbLine3D;

namespace GbSystem3D
{
extern double getDistance(const GbPoint3D &p11, const GbPoint3D &p12);
extern GbPoint3D *calculateIntersectionPoint3D(GbPoint3D &p11, GbPoint3D &p12, GbPoint3D &p21, GbPoint3D &p22);
extern GbPolygon3D *clipPolygon3D(std::vector<GbPoint3D> points, double x11, double x12, double x13, double x21,
                                  double x22, double x23);
extern GbLine3D *createClipLine3D(GbPoint3D &p1, GbPoint3D &p2, double x11, double x12, double x13, double x21,
                                  double x22, double x23);
extern bool hasIntersectionPoint3D(GbPoint3D &p11, GbPoint3D &p12, GbPoint3D &p21, GbPoint3D &p22);
extern bool isParallelIn3D(GbPoint3D &p11, GbPoint3D &p12, GbPoint3D &p21, GbPoint3D &p22);
extern GbCuboid3D *clipRectangle3D(GbPoint3D &p1, GbPoint3D &p2, double x11, double x12, double x13, double x21,
                                   double x22, double x23);

/*========================================================================================*/
inline static std::string writeGeoObject(GbObject3D *gbobject, const std::string &filename, WbWriter *writer)
{
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt3> triangles;
    gbobject->addSurfaceTriangleSet(nodes, triangles);

    std::string outFilename = writer->writeTriangles(filename, nodes, triangles);
    return outFilename;
}
// the same as before
/*========================================================================================*/
inline static std::string writeGeoObject(SPtr<GbObject3D> gbobject, const std::string &filename, WbWriter *writer)
{
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt3> triangles;
    gbobject->addSurfaceTriangleSet(nodes, triangles);

    std::string outFilename = writer->writeTriangles(filename, nodes, triangles);
    return outFilename;
}
/*========================================================================================*/
inline static std::vector<std::string> writeGeoObject(GbObject3D *gbobject, const std::string &filename,
                                                      std::vector<WbWriter *> writer)
{
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt3> triangles;
    gbobject->addSurfaceTriangleSet(nodes, triangles);

    std::vector<std::string> outFilenames;
    for (std::size_t i = 0; i < writer.size(); ++i)
        outFilenames.push_back(writer[i]->writeTriangles(filename, nodes, triangles));
    return outFilenames;
}
/*========================================================================================*/
inline static std::string writeGeoObjects(std::vector<GbObject3D *> gbobjects, const std::string &filename,
                                          WbWriter *writer)
{
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt3> triangles;

    for (std::size_t i = 0; i < gbobjects.size(); ++i) {
        // std::cout<<i<<", "<<gbobjects[i]<<std::endl;
        gbobjects[i]->addSurfaceTriangleSet(nodes, triangles);
    }

    std::string outFilename = writer->writeTriangles(filename, nodes, triangles);
    return outFilename;
}
/*========================================================================================*/

//////////////////////////////////////////////////////////////////////////
class PointSet3
{
public:
    PointSet3(int n);
    PointSet3(const std::vector<GbPoint3D> &points);
    ~PointSet3() = default;
    void add(const GbPoint3D &point);
    void addUnequal(const GbPoint3D &point);
    void add(const std::vector<GbPoint3D> &p);
    int contains(GbPoint3D *point);
    void insert(const GbPoint3D &point, int index);
    void clear();
    void clearAndTrim();
    int containsEqual(const GbPoint3D &point);
    bool containsLine(GbPoint3D *point1, GbPoint3D *point2);
    bool containsEqualLine(const GbPoint3D &point1, const GbPoint3D &point2);
    double getX1Minimum();
    double getX1Maximum();
    double getX2Minimum();
    double getX2Maximum();
    double getX3Minimum();
    double getX3Maximum();
    void calculateValues();
    int size();
    GbPoint3D *getPoint(int index);
    GbPoint3D *getFirstPoint();
    GbPoint3D *getLastPoint();
    std::vector<GbPoint3D> getPoints();

private:
    double x1min;
    double x1max;
    double x2min;
    double x2max;
    double x3min;
    double x3max;
    bool consistent;
    std::vector<GbPoint3D> points;

    void init()
    {
        consistent = false;
        x1min = x2min = x3min = 0.0;
        x1max = x2max = x3max = 0.0;
    }
};
/*=================================================================*/
class OldPointSet3
{
private:
    int sizet;
    double x1min;
    double x1max;
    double x2min;
    double x2max;
    double x3min;
    double x3max;
    bool consistent;
    std::vector<GbPoint3D> points;

    void init()
    {
        sizet      = 0;
        x1min      = 0.0;
        x1max      = 0.0;
        x2min      = 0.0;
        x2max      = 0.0;
        x3min      = 0.0;
        x3max      = 0.0;
        consistent = false;
        // points   = NULL;
    };

public:
    OldPointSet3(int n)
    {
        this->init();
        this->points.resize(n);
    }
    OldPointSet3(std::vector<GbPoint3D> &points)
    {
        this->init();
        this->points.resize(0); //, NULL);
        this->add(points);
    };
    ~OldPointSet3() = default;
    void add(GbPoint3D point)
    {
        if (this->sizet > 0 && point.equals(&(this->points)[this->sizet - 1]))
            return;
        if (this->sizet == (int)this->points.size()) {
            std::vector<GbPoint3D> a;
            a.resize(1 + (this->sizet << 1));
            for (int u = 0; u < this->sizet; u++) {
                (a)[u] = (points)[u];
            }
            this->points = a;
        }
        (this->points)[this->sizet] = point;
        this->consistent            = false;
        this->sizet++;
    }
    void addUnequal(GbPoint3D *point)
    {
        if (this->containsEqual(point) > 0)
            return;
        if (this->sizet == (int)this->points.size()) {
            std::vector<GbPoint3D> a;
            a.resize(1 + (this->sizet << 1));
            for (int u = 0; u < this->sizet; u++) {
                (a)[u] = (points)[u];
            }
            this->points = a;
        }
        (this->points)[this->sizet] = point;
        this->consistent            = false;
        this->sizet++;
    }
    void add(std::vector<GbPoint3D> &p)
    {
        if (this->sizet + p.size() >= this->points.size()) {
            std::vector<GbPoint3D> a;
            a.resize(this->sizet + p.size());
            for (int u = 0; u < (int)this->points.size(); u++) {
                (a)[u] = (this->points)[u];
            }
            this->points = a;
        }
        int u = this->sizet; // (int)this->points->size();
        for (int b = 0; b < (int)p.size(); b++)
            (this->points)[u++] = (p)[b];
        // u = this->sizet;
        // for(int b=0; b<(int)p->size(); b++)
        //    cout<<(this->points)[u++].toString()<<endl;
        this->consistent = false;
        this->sizet += (int)p.size();
    };
    //        void insert(GbPoint3D *point, int index)
    //      {
    //     if(this.size == this.points.length)
    //     {
    //        GbPoint3D a[] = new GbPoint3D[1+(this.size<<1)];
    //        System.arraycopy(this.points, 0, a, 0, this.size);
    //        this.points = a;
    //     }
    //     System.arraycopy(this.points, index, this.points, index+1, this.size-index);
    //     this.points[index] = point;
    //     this.consistent    = false;
    //     this.size++;
    //      }
    //      void delete(GbPoint3D point)
    //      {
    //     for(int i=this.size-1; i>=0; i--) if(this.points[i] == point) this.delete(i);
    //      }
    //      void delete(int index)
    //      {
    //     int j = this.size - index - 1;
    //     if(j > 0) System.arraycopy(this.points, index + 1, this.points, index, j);
    //     this.consistent = false;
    //     this.size--;
    //      }
    void clear()
    {
        this->sizet      = 0;
        this->consistent = false;
    }
    void clearAndTrim()
    {
        this->sizet = 0;
        this->points.resize(0);
        this->consistent = false;
    }

    double getX1Minimum()
    {
        if (!this->consistent)
            this->calculateValues();
        return (this->x1min);
    }
    double getX1Maximum()
    {
        if (!this->consistent)
            this->calculateValues();
        return (this->x1max);
    }
    double getX2Minimum()
    {
        if (!this->consistent)
            this->calculateValues();
        return (this->x2min);
    }
    double getX2Maximum()
    {
        if (!this->consistent)
            this->calculateValues();
        return (this->x2max);
    }
    double getX3Minimum()
    {
        if (!this->consistent)
            this->calculateValues();
        return (this->x3min);
    }
    double getX3Maximum()
    {
        if (!this->consistent)
            this->calculateValues();
        return (this->x3max);
    }
    void calculateValues()
    {
        this->x1min      = 0.0;
        this->x1max      = 0.0;
        this->x2min      = 0.0;
        this->x2max      = 0.0;
        this->x3min      = 0.0;
        this->x3max      = 0.0;
        this->consistent = true;

        if (this->sizet == 0)
            return;

        this->x1min = (this->points)[0].x1;
        this->x1max = (this->points)[0].x1;
        this->x2min = (this->points)[0].x2;
        this->x2max = (this->points)[0].x2;
        this->x3min = (this->points)[0].x3;
        this->x3max = (this->points)[0].x3;

        for (int i = this->sizet - 1; i > 0; i--) {
            if ((this->points)[i].x1 < this->x1min)
                this->x1min = (this->points)[i].x1;
            if ((this->points)[i].x1 > this->x1max)
                this->x1max = (this->points)[i].x1;
            if ((this->points)[i].x2 < this->x2min)
                this->x2min = (this->points)[i].x2;
            if ((this->points)[i].x2 > this->x2max)
                this->x2max = (this->points)[i].x2;
            if ((this->points)[i].x3 < this->x3min)
                this->x3min = (this->points)[i].x3;
            if ((this->points)[i].x3 > this->x3max)
                this->x3max = (this->points)[i].x3;
        }
    };

    int contains(GbPoint3D *point)
    {
        int n = 0;
        for (int i = this->sizet - 1; i >= 0; i--)
            if (&(this->points)[i] == point)
                n++;
        return (n);
    };
    int containsEqual(GbPoint3D *point)
    {
        int n = 0;
        for (int i = this->sizet - 1; i >= 0; i--)
            if ((this->points)[i].equals(point))
                n++;
        return (n);
    }
    bool containsLine(GbPoint3D *point1, GbPoint3D *point2)
    {
        for (int i = this->sizet - 1; i >= 0; i--)
            if (&(this->points)[i] == point1) {
                if (i == 0) {
                    if (&(this->points)[i + 1] == point2)
                        return (true);
                    if (&(this->points)[this->sizet - 1] == point2)
                        return (true);
                } else if (i == this->sizet - 1) {
                    if (&(this->points)[0] == point2)
                        return (true);
                    if (&(this->points)[i - 1] == point2)
                        return (true);
                } else {
                    if (&(this->points)[i + 1] == point2)
                        return (true);
                    if (&(this->points)[i - 1] == point2)
                        return (true);
                }
            }
        return (false);
    };
    //      boolean containsEqualLine(GbPoint2D point1, GbPoint2D point2)
    //      {
    //     for(int i=this.size-1; i>=0; i--) if(this.points[i].equals(point1))
    //     {
    //        if(i == 0)
    //        {
    //           if(this.points[i+1].equals(point2))         return(true);
    //           if(this.points[this.size-1].equals(point2)) return(true);
    //        }
    //        else if(i == this.size-1)
    //        {
    //           if(this.points[0].equals(point2))   return(true);
    //           if(this.points[i-1].equals(point2)) return(true);
    //        }
    //        else
    //        {
    //           if(this.points[i+1].equals(point2)) return(true);
    //           if(this.points[i-1].equals(point2)) return(true);
    //        }
    //     }
    //     return(false);
    //      }
    GbPoint3D *getPoint(int index) { return (&(this->points)[index]); }
    GbPoint3D *getFirstPoint() { return (&(this->points)[0]); }
    GbPoint3D *getLastPoint() { return (&(this->points)[this->sizet - 1]); }
    int size() { return (this->sizet); }
    std::vector<GbPoint3D> getPoints()
    {
        points.resize(sizet);
        return points;
        // int l = this->sizet;
        // if(l > 1 && (this->points)[0].equals(&(this->points)[l-1])) l--;

        // vector<GbPoint3D*> *a = new vector<GbPoint3D*>;
        // a->resize(l, NULL);
        // for(int u=0; u<l; u++)  { (*a)[u] = &((points)[u]);    }
        // return(a);
    }
};
/*=================================================================*/
} // namespace GbSystem3D

#endif // GBSYSTEM3D_H

//! \}
