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
#include <GbPolygon3D.h>
#include <GbSystem3D.h>

using namespace std;

double GbSystem3D::getDistance(const GbPoint3D &p11, const GbPoint3D &p12)
{
    double dx1 = p11.x1 - p12.x1;
    double dx2 = p11.x2 - p12.x2;
    double dx3 = p11.x3 - p12.x3;
    return std::sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
}

GbPoint3D *GbSystem3D::calculateIntersectionPoint3D(GbPoint3D &p11, GbPoint3D &p12, GbPoint3D &p21, GbPoint3D &p22)
{
    if (UbMath::less2(p11.x1, p12.x1, p21.x1, p22.x1))
        return NULL;
    if (UbMath::less2(p11.x2, p12.x2, p21.x2, p22.x2))
        return NULL;
    if (UbMath::less2(p11.x3, p12.x3, p21.x3, p22.x3))
        return NULL;
    if (UbMath::greater2(p11.x1, p12.x1, p21.x1, p22.x1))
        return NULL;
    if (UbMath::greater2(p11.x2, p12.x2, p21.x2, p22.x2))
        return NULL;
    if (UbMath::greater2(p11.x3, p12.x3, p21.x3, p22.x3))
        return NULL;

    double a11 = p12.x1 - p11.x1; //..HOW PARAMETERS ARE USED.........
    double a12 = p12.x2 - p11.x2; //
    double a13 = p12.x3 - p11.x3; //  p11 and p12 represent line 1
    double a21 = p21.x1 - p22.x1; //  p21 and p22 represent line 2
    double a22 = p21.x2 - p22.x2; //
    double a23 = p21.x3 - p22.x3; //..................................
    double b1  = p21.x1 - p11.x1;
    double b2  = p21.x2 - p11.x2;
    double b3  = p21.x3 - p11.x3;
    double d1  = a11 * a22 - a12 * a21;
    double d2  = a11 * a23 - a13 * a21;
    double d3  = a12 * a23 - a13 * a22;
    double t;

    if (UbMath::zero(d1) && UbMath::zero(d2) && UbMath::zero(d3))
        return NULL;
    if (UbMath::zero(d1)) {
        if (!UbMath::zero(d2))
            t = (a23 * b1 - a21 * b3) / d2;
        else
            t = (a23 * b2 - a22 * b3) / d3;
    } else if (UbMath::zero(d2)) {
        if (!UbMath::zero(d1))
            t = (a22 * b1 - a21 * b2) / d1;
        else
            t = (a23 * b2 - a22 * b3) / d3;
    } else if (UbMath::zero(d3)) {
        if (!UbMath::zero(d1))
            t = (a22 * b1 - a21 * b2) / d1;
        else
            t = (a23 * b1 - a21 * b3) / d2;
    } else
        return NULL;

    double x1 = p11.x1 + t * a11;
    double x2 = p11.x2 + t * a12;
    double x3 = p11.x3 + t * a13;

    if (UbMath::inClosedInterval(x1, p11.x1, p12.x1) && UbMath::inClosedInterval(x1, p21.x1, p22.x1) &&
        UbMath::inClosedInterval(x2, p11.x2, p12.x2) && UbMath::inClosedInterval(x2, p21.x2, p22.x2) &&
        UbMath::inClosedInterval(x3, p11.x3, p12.x3) && UbMath::inClosedInterval(x3, p21.x3, p22.x3))
        return new GbPoint3D(x1, x2, x3);

    return NULL;
}
/*=================================================================*/
// Line1: p11 -> p12 and Line2: p21 -> p22
bool GbSystem3D::hasIntersectionPoint3D(GbPoint3D &p11, GbPoint3D &p12, GbPoint3D &p21, GbPoint3D &p22)
{
    if (UbMath::less2(p11.x1, p12.x1, p21.x1, p22.x1))
        return false;
    if (UbMath::less2(p11.x2, p12.x2, p21.x2, p22.x2))
        return false;
    if (UbMath::less2(p11.x3, p12.x3, p21.x3, p22.x3))
        return false;
    if (UbMath::greater2(p11.x1, p12.x1, p21.x1, p22.x1))
        return false;
    if (UbMath::greater2(p11.x2, p12.x2, p21.x2, p22.x2))
        return false;
    if (UbMath::greater2(p11.x3, p12.x3, p21.x3, p22.x3))
        return false;

    double a11 = p12.x1 - p11.x1; //..HOW PARAMETERS ARE USED.........
    double a12 = p12.x2 - p11.x2; //
    double a13 = p12.x3 - p11.x3; //  p11 and p12 represent line 1
    double a21 = p21.x1 - p22.x1; //  p21 and p22 represent line 2
    double a22 = p21.x2 - p22.x2; //
    double a23 = p21.x3 - p22.x3; //..................................
    double b1  = p21.x1 - p11.x1;
    double b2  = p21.x2 - p11.x2;
    double b3  = p21.x3 - p11.x3;
    double d1  = a11 * a22 - a12 * a21;
    double d2  = a11 * a23 - a13 * a21;
    double d3  = a12 * a23 - a13 * a22;
    double t;

    if (UbMath::zero(d1) && UbMath::zero(d2) && UbMath::zero(d3))
        return false;
    if (UbMath::zero(d1)) {
        if (!UbMath::zero(d2))
            t = (a23 * b1 - a21 * b3) / d2;
        else
            t = (a23 * b2 - a22 * b3) / d3;
    } else if (UbMath::zero(d2)) {
        if (!UbMath::zero(d1))
            t = (a22 * b1 - a21 * b2) / d1;
        else
            t = (a23 * b2 - a22 * b3) / d3;
    } else if (UbMath::zero(d3)) {
        if (!UbMath::zero(d1))
            t = (a22 * b1 - a21 * b2) / d1;
        else
            t = (a23 * b1 - a21 * b3) / d2;
    } else
        return false;

    double x1 = p11.x1 + t * a11;
    double x2 = p11.x2 + t * a12;
    double x3 = p11.x3 + t * a13;

    if (UbMath::inClosedInterval(x1, p11.x1, p12.x1) && UbMath::inClosedInterval(x1, p21.x1, p22.x1) &&
        UbMath::inClosedInterval(x2, p11.x2, p12.x2) && UbMath::inClosedInterval(x2, p21.x2, p22.x2) &&
        UbMath::inClosedInterval(x3, p11.x3, p12.x3) && UbMath::inClosedInterval(x3, p21.x3, p22.x3))
        return true;
    return false;
}
/*======================================================================*/
//
//
//   /*======================================================================*/
//   /*  Private Methoden (Parallelism)                                      */
//   /*                                                                      */
bool GbSystem3D::isParallelIn3D(GbPoint3D &p11, GbPoint3D &p12, GbPoint3D &p21, GbPoint3D &p22)
{
    double a11 = p12.x1 - p11.x1; //..HOW PARAMETERS ARE USED.........
    double a12 = p12.x2 - p11.x2; //
    double a13 = p12.x3 - p11.x3; //  p11 and p12 represent line 1
    double a21 = p21.x1 - p22.x1; //  p21 and p22 represent line 2
    double a22 = p21.x2 - p22.x2; //
    double a23 = p21.x3 - p22.x3; //..................................

    return (UbMath::zero(a11 * a22 - a12 * a21) && UbMath::zero(a11 * a23 - a13 * a21) &&
            UbMath::zero(a12 * a23 - a13 * a22));
}
/*======================================================================*/

/*======================================================================*/
/*  General Clipping Methods                                            */
//......................................................................*/
//
//  Method       Parameters                                       Result      Remarks
//  ---------    ---------------------------------------------    ---------   -------------------
//  clip###2D   (2D objects to be clipped, 2+2 clipping values)   2D object   clipping x1, x2
//  clip###3D   (3D objects to be clipped, 2+2 clipping values)   3D object   clipping x1, x2
//  clip###3D   (3D objects to be clipped, 3+3 clipping values)   3D object   clipping x1, x2, x3
//  clip###3D   (3D objects to be clipped, 1+1 clipping values)   3D object   clipping x3
//
/*======================================================================*/
/*  Private Methoden (Clipping Lines)                                   */
/*                                                                      */
GbLine3D *GbSystem3D::createClipLine3D(GbPoint3D &pA, GbPoint3D &pB, double x1a, double x2a, double x3a, double x1b,
                                       double x2b, double x3b)
{
    GbPoint3D *p1 = new GbPoint3D(pA);
    GbPoint3D *p2 = new GbPoint3D(pB);

    if (UbMath::greater(x1a, x1b)) {
        double x1 = x1a;
        x1a       = x1b;
        x1b       = x1;
    }
    if (UbMath::greater(x2a, x2b)) {
        double x2 = x2a;
        x2a       = x2b;
        x2b       = x2;
    }
    if (UbMath::greater(x3a, x3b)) {
        double x3 = x3a;
        x3a       = x3b;
        x3b       = x3;
    }

    double f;

    /*-------------------------------------------------------------------*/
    /*  Schneiden an vorderer Kante                                      */
    /*                                                                   */
    if (UbMath::less(p1->x3, x3a)) {
        if (UbMath::less(p2->x3, x3a)) {
            delete p1;
            delete p2;
            return NULL;
        }

        f = (x3a - p1->x3) / (p1->x3 - p2->x3);
        p1->x1 += (p1->x1 - p2->x1) * f;
        p1->x2 += (p1->x2 - p2->x2) * f;
        p1->x3 = x3a;
    } else if (UbMath::less(p2->x3, x3a)) {
        f = (x3a - p2->x3) / (p2->x3 - p1->x3);
        p2->x1 += (p2->x1 - p1->x1) * f;
        p2->x2 += (p2->x2 - p1->x2) * f;
        p2->x3 = x3a;
    }
    /*-------------------------------------------------------------------*/
    /*  Schneiden an unterer Kante                                       */
    /*                                                                   */
    if (UbMath::less(p1->x2, x2a)) {
        if (UbMath::less(p2->x2, x2a)) {
            delete p1;
            delete p2;
            return NULL;
        }

        f = (x2a - p1->x2) / (p1->x2 - p2->x2);
        p1->x1 += (p1->x1 - p2->x1) * f;
        p1->x3 += (p1->x3 - p2->x3) * f;
        p1->x2 = x2a;
    } else if (UbMath::less(p2->x2, x2a)) {
        f = (x2a - p2->x2) / (p2->x2 - p1->x2);
        p2->x1 += (p2->x1 - p1->x1) * f;
        p2->x3 += (p2->x3 - p1->x3) * f;
        p2->x2 = x2a;
    }
    /*-------------------------------------------------------------------*/
    /*  Schneiden an rechter Kante                                       */
    /*                                                                   */
    if (UbMath::greater(p1->x1, x1b)) {
        if (UbMath::greater(p2->x1, x1b)) {
            delete p1;
            delete p2;
            return NULL;
        }

        f = (x1b - p1->x1) / (p1->x1 - p2->x1);
        p1->x2 += (p1->x2 - p2->x2) * f;
        p1->x3 += (p1->x3 - p2->x3) * f;
        p1->x1 = x1b;
    } else if (UbMath::greater(p2->x1, x1b)) {
        f = (x1b - p2->x1) / (p2->x1 - p1->x1);
        p2->x2 += (p2->x2 - p1->x2) * f;
        p2->x3 += (p2->x3 - p1->x3) * f;
        p2->x1 = x1b;
    }
    /*-------------------------------------------------------------------*/
    /*  Schneiden an hinterer Kante                                      */
    /*                                                                   */
    if (UbMath::greater(p1->x3, x3b)) {
        if (UbMath::greater(p2->x3, x3b)) {
            delete p1;
            delete p2;
            return NULL;
        }

        f = (x3b - p1->x3) / (p1->x3 - p2->x3);
        p1->x1 += (p1->x1 - p2->x1) * f;
        p1->x2 += (p1->x2 - p2->x2) * f;
        p1->x3 = x3b;
    } else if (UbMath::greater(p2->x3, x3b)) {
        f = (x3b - p2->x3) / (p2->x3 - p1->x3);
        p2->x1 += (p2->x1 - p1->x1) * f;
        p2->x2 += (p2->x2 - p1->x2) * f;
        p2->x3 = x3b;
    }
    /*-------------------------------------------------------------------*/
    /*  Schneiden an oberer Kante                                        */
    /*                                                                   */
    if (UbMath::greater(p1->x2, x2b)) {
        if (UbMath::greater(p2->x2, x2b)) {
            delete p1;
            delete p2;
            return NULL;
        }

        f = (x2b - p1->x2) / (p1->x2 - p2->x2);
        p1->x1 += (p1->x1 - p2->x1) * f;
        p1->x3 += (p1->x3 - p2->x3) * f;
        p1->x2 = x2b;
    } else if (UbMath::greater(p2->x2, x2b)) {
        f = (x2b - p2->x2) / (p2->x2 - p1->x2);
        p2->x1 += (p2->x1 - p1->x1) * f;
        p2->x3 += (p2->x3 - p1->x3) * f;
        p2->x2 = x2b;
    }
    /*-------------------------------------------------------------------*/
    /*  Schneiden an linker Kante                                        */
    /*                                                                   */
    if (UbMath::less(p1->x1, x1a)) {
        if (UbMath::less(p2->x1, x1a)) {
            delete p1;
            delete p2;
            return NULL;
        }

        f = (x1a - p1->x1) / (p1->x1 - p2->x1);
        p1->x2 += (p1->x2 - p2->x2) * f;
        p1->x3 += (p1->x3 - p2->x3) * f;
        p1->x1 = x1a;
    } else if (UbMath::less(p2->x1, x1a)) {
        f = (x1a - p2->x1) / (p2->x1 - p1->x1);
        p2->x2 += (p2->x2 - p1->x2) * f;
        p2->x3 += (p2->x3 - p1->x3) * f;
        p2->x1 = x1a;
    }
    /*-------------------------------------------------------------------*/
    return new GbLine3D(p1, p2);
}
//   /*======================================================================*/
//   /*  Private Methoden (Clipping Rectangles)                              */
//   /*                                                                      */
//   final static GbPolygon3D clipPolygon3D(GbPoint3D points[], double x11, double x12, double x21, double x22)
//   {
//      GbPoint3D last = null;
//      PointSet3 ps   = new PointSet3(points);
//      boolean   flag = false;
//      int       n    = points.length;
//      int       i;
//      double    f;
//
//      if(n == 0)              return(null);
//      if(greater(x11, x21)) { double ax = x11; x11 = x21; x21 = ax; }
//      if(greater(x12, x22)) { double ay = x12; x12 = x22; x22 = ay; }
//
//      /*-------------------------------------------------------------------*/
//      /*  Schneiden an unterer Kante                                       */
//      /*                                                                   */
//      if(less(ps.getX2Minimum(), x12))
//      {
//     ps.clear();
//     last = points[0];
//     if(less((*points)[0]->x2, x12)) flag = false;
//     else
//     {
//        ps.add(points[0]);
//        flag = true;
//     }
//     for(i=1; i<n; i++)
//     {
//        if(less((*points)[i]->x2, x12))
//        {
//           if(flag)
//           {
//              f = (x12-(*points)[i]->x2)/((*points)[i]->x2-last->x2);
//              ps.add(new GbPoint3D((*points)[i]->x1 + ((*points)[i]->x1-last->x1)*f, x12, (*points)[i]->x3 +
//((*points)[i]->x3-last->x3)*f));
//           }
//           flag = false;
//        }
//        else
//        {
//           if(!flag)
//           {
//              f = (x12-(*points)[i]->x2)/((*points)[i]->x2-last->x2);
//              ps.add(new GbPoint3D((*points)[i]->x1 + ((*points)[i]->x1-last->x1)*f, x12, (*points)[i]->x3 +
//((*points)[i]->x3-last->x3)*f));
//           }
//           ps.add((*points)[i]);
//           flag = true;
//        }
//        last = points[i];
//     }
//     if(!((less(points[0].x2, x12)) ^ flag))
//     {
//        f = (x12-points[0].x2)/(points[0].x2-last->x2);
//        ps.add(new GbPoint3D(points[0].x1 + (points[0].x1-last->x1)*f, x12, points[0].x3 + (points[0].x3-last->x3)*f));
//     }
//
//     points = ps.getPoints();
//     n      = points.length;
//
//     if(n == 0) return(null);
//      }
//      /*-------------------------------------------------------------------*/
//      /*  Schneiden an rechter Kante                                       */
//      /*                                                                   */
//      if(greater(ps.getX1Maximum(), x21))
//      {
//     ps.clear();
//     last = points[0];
//     if(greater(points[0].x1, x21)) flag = false;
//     else
//     {
//        ps.add(points[0]);
//        flag = true;
//     }
//     for(i=1; i<n; i++)
//     {
//        if(greater((*points)[i]->x1, x21))
//        {
//           if(flag)
//           {
//              f = (x21-(*points)[i]->x1)/((*points)[i]->x1-last->x1);
//              ps.add(new GbPoint3D(x21, (*points)[i]->x2 + ((*points)[i]->x2-last->x2)*f, (*points)[i]->x3 +
//((*points)[i]->x3-last->x3)*f));
//           }
//           flag = false;
//        }
//        else
//        {
//           if(!flag)
//           {
//              f = (x21-(*points)[i]->x1)/((*points)[i]->x1-last->x1);
//              ps.add(new GbPoint3D(x21, (*points)[i]->x2 + ((*points)[i]->x2-last->x2)*f, (*points)[i]->x3 +
//((*points)[i]->x3-last->x3)*f));
//           }
//           ps.add(points[i]);
//           flag = true;
//        }
//        last = points[i];
//     }
//     if(!((greater(points[0].x1, x21)) ^ flag))
//     {
//        f = (x21-points[0].x1)/(points[0].x1-last.x1);
//        ps.add(new GbPoint3D(x21, points[0].x2 + (points[0].x2-last.x2)*f, points[0].x3 + (points[0].x3-last.x3)*f));
//     }
//
//     points = ps.getPoints();
//     n      = points.length;
//
//     if(n == 0) return(null);
//      }
//      /*-------------------------------------------------------------------*/
//      /*  Schneiden an oberer Kante                                        */
//      /*                                                                   */
//      if(greater(ps.getX2Maximum(), x22))
//      {
//     ps.clear();
//     last = points[0];
//     if(greater(points[0].x2, x22)) flag = false;
//     else
//     {
//        ps.add(points[0]);
//        flag = true;
//     }
//     for(i=1; i<n; i++)
//     {
//        if(greater((*points)[i]->x2, x22))
//        {
//           if(flag)
//           {
//              f = (x22-(*points)[i]->x2)/(points[i].x2-last.x2);
//              ps.add(new GbPoint3D(points[i].x1 + (points[i].x1-last.x1)*f, x22, points[i].x3 +
//(points[i].x3-last.x3)*f));
//           }
//           flag = false;
//        }
//        else
//        {
//           if(!flag)
//           {
//              f = (x22-points[i].x2)/(points[i].x2-last.x2);
//              ps.add(new GbPoint3D(points[i].x1 + (points[i].x1-last.x1)*f, x22, points[i].x3 +
//(points[i].x3-last.x3)*f));
//           }
//           ps.add(points[i]);
//           flag = true;
//        }
//        last = points[i];
//     }
//     if(!((greater(points[0].x2, x22)) ^ flag))
//     {
//        f = (x22-points[0].x2)/(points[0].x2-last.x2);
//        ps.add(new GbPoint3D(points[0].x1 + (points[0].x1-last.x1)*f, x22, points[0].x3 + (points[0].x3-last.x3)*f));
//     }
//
//     points = ps.getPoints();
//     n      = points.length;
//
//     if(n == 0) return(null);
//      }
//      /*-------------------------------------------------------------------*/
//      /*  Schneiden an linker Kante                                        */
//      /*                                                                   */
//      if(less(ps.getX1Minimum(), x11))
//      {
//     ps.clear();
//     last = points[0];
//     if(less(points[0].x1, x11)) flag = false;
//     else
//     {
//        ps.add(points[0]);
//        flag = true;
//     }
//     for(i=1; i<n; i++)
//     {
//        if(less(points[i].x1, x11))
//        {
//           if(flag)
//           {
//              f = (x11-points[i].x1)/(points[i].x1-last.x1);
//              ps.add(new GbPoint3D(x11, points[i].x2 + (points[i].x2-last.x2)*f, points[i].x3 +
//(points[i].x3-last.x3)*f));
//           }
//           flag = false;
//        }
//        else
//        {
//           if(!flag)
//           {
//              f = (x11-points[i].x1)/(points[i].x1-last.x1);
//              ps.add(new GbPoint3D(x11, points[i].x2 + (points[i].x2-last.x2)*f, points[i].x3 +
//(points[i].x3-last.x3)*f));
//           }
//           ps.add(points[i]);
//           flag = true;
//        }
//        last = points[i];
//     }
//     if(!((less(points[0].x1, x11)) ^ flag))
//     {
//        f = (x11-points[0].x1)/(points[0].x1-last.x1);
//        ps.add(new GbPoint3D(x11, points[0].x2 + (points[0].x2-last.x2)*f, points[0].x3 + (points[0].x3-last.x3)*f));
//     }
//
//     points = ps.getPoints();
//     n      = points.length;
//
//     if(n == 0) return(null);
//      }
//      /*-------------------------------------------------------------------*/
//      GbPolygon3D polygon = new GbPolygon3D(points);
//
//      if(n > 2)
//      {
//     for(i=2; i<n; i++) if(zero(i_TA(points[i-2], points[i-1], points[i]))) polygon.deletePoint(points[i-1]);
//     if(zero(i_TA(points[n-2], points[n-1], points[0]))) polygon.deletePoint(points[n-1]);
//     if(zero(i_TA(points[n-1], points[0],   points[1]))) polygon.deletePoint(points[0]);
//      }
//      return(polygon);
//   }
//   final static GbPolygon3D clipPolygon3D(GbPoint3D points[], double x13, double x23)
//   {
//      GbPoint3D last = null;
//      PointSet3 ps   = new PointSet3(points);
//      boolean   flag = false;
//      int       n    = points.length;
//      int       i;
//      double    f;
//
//      if(n == 0)              return(null);
//      if(greater(x13, x23)) { double az = x13; x13 = x23; x23 = az; }
//
//      /*-------------------------------------------------------------------*/
//      /*  Schneiden an vorderer Kante                                      */
//      /*                                                                   */
//      if(less(ps.getX3Minimum(), x13))
//      {
//     ps.clear();
//     last = points[0];
//     if(less(points[0].x3, x13)) flag = false;
//     else
//     {
//        ps.add(points[0]);
//        flag = true;
//     }
//     for(i=1; i<n; i++)
//     {
//        if(less(points[i].x3, x13))
//        {
//           if(flag)
//           {
//              f = (x13-points[i].x3)/(points[i].x3-last.x3);
//              ps.add(new GbPoint3D(points[i].x1 + (points[i].x1-last.x1)*f, points[i].x2 + (points[i].x2-last.x2)*f,
//x13));
//           }
//           flag = false;
//        }
//        else
//        {
//           if(!flag)
//           {
//              f = (x13-points[i].x3)/(points[i].x3-last.x3);
//              ps.add(new GbPoint3D(points[i].x1 + (points[i].x1-last.x1)*f, points[i].x2 + (points[i].x2-last.x2)*f,
//x13));
//           }
//           ps.add(points[i]);
//           flag = true;
//        }
//        last = points[i];
//     }
//     if(!((less(points[0].x3, x13)) ^ flag))
//     {
//        f = (x13-points[0].x3)/(points[0].x3-last.x3);
//        ps.add(new GbPoint3D(points[0].x1 + (points[0].x1-last.x1)*f, points[0].x2 + (points[0].x2-last.x2)*f, x13));
//     }
//
//     points = ps.getPoints();
//     n      = points.length;
//
//     if(n == 0) return(null);
//      }
//      /*-------------------------------------------------------------------*/
//      /*  Schneiden an hinterer Kante                                      */
//      /*                                                                   */
//      if(greater(ps.getX3Maximum(), x23))
//      {
//     ps.clear();
//     last = points[0];
//     if(greater(points[0].x3, x23)) flag = false;
//     else
//     {
//        ps.add(points[0]);
//        flag = true;
//     }
//     for(i=1; i<n; i++)
//     {
//        if(greater(points[i].x3, x23))
//        {
//           if(flag)
//           {
//              f = (x23-points[i].x3)/(points[i].x3-last.x3);
//              ps.add(new GbPoint3D(points[i].x1 + (points[i].x1-last.x1)*f, points[i].x2 + (points[i].x2-last.x2)*f,
//x23));
//           }
//           flag = false;
//        }
//        else
//        {
//           if(!flag)
//           {
//              f = (x23-points[i].x3)/(points[i].x3-last.x3);
//              ps.add(new GbPoint3D(points[i].x1 + ((*points)[i]->x1-last.x1)*f, (*points)[i]->x2 +
//((*points)[i]->x2-last.x2)*f, x23));
//           }
//           ps.add(points[i]);
//           flag = true;
//        }
//        last = points[i];
//     }
//     if(!((greater(points[0].x3, x23)) ^ flag))
//     {
//        f = (x23-points[0].x3)/(points[0].x3-last.x3);
//        ps.add(new GbPoint3D(points[0].x1 + (points[0].x1-last.x1)*f, points[0].x2 + (points[0].x2-last.x2)*f, x23));
//     }
//
//     points = ps.getPoints();
//     n      = points.length;
//
//     if(n == 0) return(null);
//      }
//      /*-------------------------------------------------------------------*/
//      GbPolygon3D polygon = new GbPolygon3D(points);
//
//      return(polygon);
//   }
GbPolygon3D *GbSystem3D::clipPolygon3D(vector<GbPoint3D> points, double x11, double x12, double x13, double x21,
                                       double x22, double x23)
{
    GbPoint3D last;
    PointSet3 ps(points);
    bool flag = false;
    int n     = (int)points.size();
    int i;
    double f;

    if (n == 0)
        return NULL;
    if (UbMath::greater(x11, x21)) {
        double ax = x11;
        x11       = x21;
        x21       = ax;
    }
    if (UbMath::greater(x12, x22)) {
        double ay = x12;
        x12       = x22;
        x22       = ay;
    }
    if (UbMath::greater(x13, x23)) {
        double az = x13;
        x13       = x23;
        x23       = az;
    }

    /*-------------------------------------------------------------------*/
    /*  Schneiden an vorderer Kante                                      */
    /*                                                                   */
    if (UbMath::less(ps.getX3Minimum(), x13)) {
        ps.clear();
        last = (points)[0];
        if (UbMath::less((points)[0].x3, x13))
            flag = false;
        else {
            ps.add((points)[0]);
            flag = true;
        }
        for (i = 1; i < n; i++) {
            if (UbMath::less((points)[i].x3, x13)) {
                if (flag) {
                    f = (x13 - (points)[i].x3) / ((points)[i].x3 - last.x3);
                    ps.add(GbPoint3D((points)[i].x1 + ((points)[i].x1 - last.x1) * f,
                                     (points)[i].x2 + ((points)[i].x2 - last.x2) * f, x13));
                }
                flag = false;
            } else {
                if (!flag) {
                    f = (x13 - (points)[i].x3) / ((points)[i].x3 - last.x3);
                    ps.add(GbPoint3D((points)[i].x1 + ((points)[i].x1 - last.x1) * f,
                                     (points)[i].x2 + ((points)[i].x2 - last.x2) * f, x13));
                }
                ps.add((points)[i]);
                flag = true;
            }
            last = (points)[i];
        }
        if (!((UbMath::less((points)[0].x3, x13)) ^ flag)) {
            f = (x13 - (points)[0].x3) / ((points)[0].x3 - last.x3);
            ps.add(GbPoint3D((points)[0].x1 + ((points)[0].x1 - last.x1) * f,
                             (points)[0].x2 + ((points)[0].x2 - last.x2) * f, x13));
        }

        points = ps.getPoints();
        n      = (int)points.size();

        if (n == 0)
            return NULL;
    }

    /*-------------------------------------------------------------------*/
    /*  Schneiden an unterer Kante                                       */
    /*                                                                   */
    if (UbMath::less(ps.getX2Minimum(), x12)) {
        ps.clear();
        last = (points)[0];
        if (UbMath::less((points)[0].x2, x12))
            flag = false;
        else {
            ps.add((points)[0]);
            flag = true;
        }
        for (i = 1; i < n; i++) {
            if (UbMath::less((points)[i].x2, x12)) {
                if (flag) {
                    f = (x12 - (points)[i].x2) / ((points)[i].x2 - last.x2);
                    ps.add(GbPoint3D((points)[i].x1 + ((points)[i].x1 - last.x1) * f, x12,
                                     (points)[i].x3 + ((points)[i].x3 - last.x3) * f));
                }
                flag = false;
            } else {
                if (!flag) {
                    f = (x12 - (points)[i].x2) / ((points)[i].x2 - last.x2);
                    ps.add(GbPoint3D((points)[i].x1 + ((points)[i].x1 - last.x1) * f, x12,
                                     (points)[i].x3 + ((points)[i].x3 - last.x3) * f));
                }
                ps.add((points)[i]);
                flag = true;
            }
            last = (points)[i];
        }
        if (!((UbMath::less((points)[0].x2, x12)) ^ flag)) {
            f = (x12 - (points)[0].x2) / ((points)[0].x2 - last.x2);
            ps.add(GbPoint3D((points)[0].x1 + ((points)[0].x1 - last.x1) * f, x12,
                             (points)[0].x3 + ((points)[0].x3 - last.x3) * f));
        }

        points = ps.getPoints();
        n      = (int)points.size();

        if (n == 0)
            return NULL;
    }
    /*-------------------------------------------------------------------*/
    /*  Schneiden an rechter Kante                                       */
    /*                                                                   */

    if (UbMath::greater(ps.getX1Maximum(), x21)) {
        ps.clear();
        last = (points)[0];
        if (UbMath::greater((points)[0].x1, x21))
            flag = false;
        else {
            ps.add((points)[0]);
            flag = true;
        }
        for (i = 1; i < n; i++) {
            if (UbMath::greater((points)[i].x1, x21)) {
                if (flag) {
                    f = (x21 - (points)[i].x1) / ((points)[i].x1 - last.x1);
                    ps.add(GbPoint3D(x21, (points)[i].x2 + ((points)[i].x2 - last.x2) * f,
                                     (points)[i].x3 + ((points)[i].x3 - last.x3) * f));
                }
                flag = false;
            } else {
                if (!flag) {
                    f = (x21 - (points)[i].x1) / ((points)[i].x1 - last.x1);
                    ps.add(GbPoint3D(x21, (points)[i].x2 + ((points)[i].x2 - last.x2) * f,
                                     (points)[i].x3 + ((points)[i].x3 - last.x3) * f));
                }
                ps.add((points)[i]);
                flag = true;
            }
            last = (points)[i];
        }
        if (!((UbMath::greater((points)[0].x1, x21)) ^ flag)) {
            f = (x21 - (points)[0].x1) / ((points)[0].x1 - last.x1);
            ps.add(GbPoint3D(x21, (points)[0].x2 + ((points)[0].x2 - last.x2) * f,
                             (points)[0].x3 + ((points)[0].x3 - last.x3) * f));
        }

        points = ps.getPoints();
        n      = (int)points.size();

        if (n == 0)
            return NULL;
    }
    /*-------------------------------------------------------------------*/
    /*  Schneiden an hinterer Kante                                      */
    /*                                                                   */
    if (UbMath::greater(ps.getX3Maximum(), x23)) {
        ps.clear();
        last = (points)[0];
        if (UbMath::greater((points)[0].x3, x23))
            flag = false;
        else {
            ps.add((points)[0]);
            flag = true;
        }
        for (i = 1; i < n; i++) {
            if (UbMath::greater((points)[i].x3, x23)) {
                if (flag) {
                    f = (x23 - (points)[i].x3) / ((points)[i].x3 - last.x3);
                    ps.add(GbPoint3D((points)[i].x1 + ((points)[i].x1 - last.x1) * f,
                                     (points)[i].x2 + ((points)[i].x2 - last.x2) * f, x23));
                }
                flag = false;
            } else {
                if (!flag) {
                    f = (x23 - (points)[i].x3) / ((points)[i].x3 - last.x3);
                    ps.add(GbPoint3D((points)[i].x1 + ((points)[i].x1 - last.x1) * f,
                                     (points)[i].x2 + ((points)[i].x2 - last.x2) * f, x23));
                }
                ps.add((points)[i]);
                flag = true;
            }
            last = (points)[i];
        }
        if (!((UbMath::greater((points)[0].x3, x23)) ^ flag)) {
            f = (x23 - (points)[0].x3) / ((points)[0].x3 - last.x3);
            ps.add(GbPoint3D((points)[0].x1 + ((points)[0].x1 - last.x1) * f,
                             (points)[0].x2 + ((points)[0].x2 - last.x2) * f, x23));
        }

        points = ps.getPoints();
        n      = (int)points.size();

        if (n == 0)
            return NULL;
    }
    /*-------------------------------------------------------------------*/
    /*  Schneiden an oberer Kante                                        */
    /*                                                                   */

    if (UbMath::greater(ps.getX2Maximum(), x22)) {
        ps.clear();
        last = (points)[0];
        if (UbMath::greater((points)[0].x2, x22))
            flag = false;
        else {
            ps.add((points)[0]);
            flag = true;
        }
        for (i = 1; i < n; i++) {
            if (UbMath::greater((points)[i].x2, x22)) {
                if (flag) {
                    f = (x22 - (points)[i].x2) / ((points)[i].x2 - last.x2);
                    ps.add(GbPoint3D((points)[i].x1 + ((points)[i].x1 - last.x1) * f, x22,
                                     (points)[i].x3 + ((points)[i].x3 - last.x3) * f));
                }
                flag = false;
            } else {
                if (!flag) {
                    f = (x22 - (points)[i].x2) / ((points)[i].x2 - last.x2);
                    ps.add(GbPoint3D((points)[i].x1 + ((points)[i].x1 - last.x1) * f, x22,
                                     (points)[i].x3 + ((points)[i].x3 - last.x3) * f));
                }
                ps.add((points)[i]);
                flag = true;
            }
            last = (points)[i];
        }
        if (!((UbMath::greater((points)[0].x2, x22)) ^ flag)) {
            f = (x22 - (points)[0].x2) / ((points)[0].x2 - last.x2);
            ps.add(GbPoint3D((points)[0].x1 + ((points)[0].x1 - last.x1) * f, x22,
                             (points)[0].x3 + ((points)[0].x3 - last.x3) * f));
        }

        points = ps.getPoints();
        n      = (int)points.size();

        if (n == 0)
            return NULL;
    }
    /*-------------------------------------------------------------------*/
    /*  Schneiden an linker Kante                                        */
    /*                                                                   */
    if (UbMath::less(ps.getX1Minimum(), x11)) {
        ps.clear();
        last = (points)[0];
        if (UbMath::less((points)[0].x1, x11))
            flag = false;
        else {
            ps.add((points)[0]);
            flag = true;
        }
        for (i = 1; i < n; i++) {
            if (UbMath::less((points)[i].x1, x11)) {
                if (flag) {
                    f = (x11 - (points)[i].x1) / ((points)[i].x1 - last.x1);
                    ps.add(GbPoint3D(x11, (points)[i].x2 + ((points)[i].x2 - last.x2) * f,
                                     (points)[i].x3 + ((points)[i].x3 - last.x3) * f));
                }
                flag = false;
            } else {
                if (!flag) {
                    f = (x11 - (points)[i].x1) / ((points)[i].x1 - last.x1);
                    ps.add(GbPoint3D(x11, (points)[i].x2 + ((points)[i].x2 - last.x2) * f,
                                     (points)[i].x3 + ((points)[i].x3 - last.x3) * f));
                }
                ps.add((points)[i]);
                flag = true;
            }
            last = (points)[i];
        }
        if (!((UbMath::less((points)[0].x1, x11)) ^ flag)) {
            f = (x11 - (points)[0].x1) / ((points)[0].x1 - last.x1);
            ps.add(GbPoint3D(x11, (points)[0].x2 + ((points)[0].x2 - last.x2) * f,
                             (points)[0].x3 + ((points)[0].x3 - last.x3) * f));
        }

        points = ps.getPoints();
        n      = (int)points.size();

        if (n == 0)
            return NULL;
    }
    /*-------------------------------------------------------------------*/
    return new GbPolygon3D(points);
}
/*=========================================================================*/
GbCuboid3D *GbSystem3D::clipRectangle3D(GbPoint3D &p1, GbPoint3D &p2, double x11, double x12, double x13, double x21,
                                        double x22, double x23)
{
    double r11 = p1.x1;
    double r12 = p1.x2;
    double r13 = p1.x3;
    double r21 = p2.x1;
    double r22 = p2.x2;
    double r23 = p2.x3;

    if (UbMath::greater(x11, x21)) {
        double ax = x11;
        x11       = x21;
        x21       = ax;
    }
    if (UbMath::greater(x12, x22)) {
        double ay = x12;
        x12       = x22;
        x22       = ay;
    }
    if (UbMath::greater(x13, x23)) {
        double az = x13;
        x13       = x23;
        x23       = az;
    }
    if (UbMath::greater(r11, r21)) {
        double bx = r11;
        r11       = r21;
        r21       = bx;
    }
    if (UbMath::greater(r12, r22)) {
        double by = r12;
        r12       = r22;
        r22       = by;
    }
    if (UbMath::greater(r13, r23)) {
        double bz = r13;
        r13       = r23;
        r23       = bz;
    }

    double m11 = UbMath::greater(x11, r11) ? x11 : r11;
    double m12 = UbMath::greater(x12, r12) ? x12 : r12;
    double m13 = UbMath::greater(x13, r13) ? x13 : r13;
    double m21 = UbMath::greater(x21, r21) ? r21 : x21;
    double m22 = UbMath::greater(x22, r22) ? r22 : x22;
    double m23 = UbMath::greater(x23, r23) ? r23 : x23;

    if (UbMath::lessEqual(m11, m21) && UbMath::lessEqual(m12, m22) && UbMath::lessEqual(m13, m23))
        return (new GbCuboid3D(new GbPoint3D(m11, m12, m13), new GbPoint3D(m21, m22, m23)));
    else
        return (NULL);
}

/*=========================================================================*/
/*=========================================================================*/
/*=========================================================================*/

GbSystem3D::PointSet3::PointSet3(int n)
{
    this->init();
    this->points.reserve(n); // reserves n elements! but the size of the vector ist still "0"
}
/*=======================================================*/
GbSystem3D::PointSet3::PointSet3(const vector<GbPoint3D> &points)
{
    this->init();
    this->add(points);
}
/*=======================================================*/
void GbSystem3D::PointSet3::add(const GbPoint3D &point)
{
    // is point equal to last point in points then return
    if (!this->points.empty() && point.equals(&this->points.back()))
        return; // WHY???

    // push point to vector
    this->points.push_back(point);

    this->consistent = false;
}
/*=======================================================*/
void GbSystem3D::PointSet3::addUnequal(const GbPoint3D &point)
{
    if (this->containsEqual(point) > 0)
        return;

    this->points.push_back(point);
    this->consistent = false;
}
/*=======================================================*/
void GbSystem3D::PointSet3::add(const vector<GbPoint3D> &pointVector)
{
    for (int pos = 0; pos < (int)pointVector.size(); pos++)
        this->points.push_back(pointVector[pos]);

    this->consistent = false;
}
/*=======================================================*/
void GbSystem3D::PointSet3::insert(const GbPoint3D &point, int index)
{
    if (index < 0 || index >= (int)this->points.size())
        throw UbException(UB_EXARGS, "index out of range");

    // get iterator for index-position
    vector<GbPoint3D>::iterator pos = this->points.begin();
    for (int i = 1; i <= index; i++)
        ++pos;

    // insert point
    this->points.insert(pos, point);

    this->consistent = false;
}
/*=======================================================*/
// void delete(GbPoint3D point)
//{
//   for(int i=this.size-1; i>=0; i--) if(this.points[i] == point) this.delete(i);
//}
/*=======================================================*/
// void delete(int index)
//{
//   int j = this.size - index - 1;
//   if(j > 0) System.arraycopy(this.points, index + 1, this.points, index, j);
//   this.consistent = false;
//   this.size--;
//}
/*=======================================================*/
void GbSystem3D::PointSet3::clear()
{
    // clears points (size==0 but capacity is the old c1)
    this->points.clear();
    this->consistent = false;
}
/*=======================================================*/
void GbSystem3D::PointSet3::clearAndTrim()
{
    // clears points (size==0 AND capacity==0)
    this->points.resize(0);
    this->consistent = false;
}
/*=======================================================*/
double GbSystem3D::PointSet3::getX1Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return this->x1min;
}
/*=======================================================*/
double GbSystem3D::PointSet3::getX1Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return this->x1max;
}
/*=======================================================*/
double GbSystem3D::PointSet3::getX2Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return this->x2min;
}
/*=======================================================*/
double GbSystem3D::PointSet3::getX2Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return this->x2max;
}
/*=======================================================*/
double GbSystem3D::PointSet3::getX3Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return this->x3min;
}
/*=======================================================*/
double GbSystem3D::PointSet3::getX3Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return this->x3max;
}
/*=======================================================*/
int GbSystem3D::PointSet3::contains(GbPoint3D *point)
{
    // returns number of points which has the same adress (this should be 0 or 1!!!)
    int n = 0;

    for (int pos = (int)this->points.size() - 1; pos >= 0; pos--)
        if (&this->points[pos] == point)
            n++;

    return n;
}
/*=======================================================*/
int GbSystem3D::PointSet3::containsEqual(const GbPoint3D &point)
{
    // returns number of points which have the same coordinates with point (could be 0,1 or even more)
    int n = 0;

    for (int pos = (int)this->points.size() - 1; pos >= 0; pos--)
        if (this->points[pos].equals(&point))
            n++;

    return n;
}
/*=======================================================*/
bool GbSystem3D::PointSet3::containsLine(GbPoint3D *point1, GbPoint3D *point2)
{
    // returns true if pointset has c2 in "this->points"vector  neighboured points
    // wich have the same adress as point1 or point2
    vector<GbPoint3D>::iterator pos1 = this->points.begin();
    vector<GbPoint3D>::iterator pos2;

    for (pos2 = pos1++; pos2 != this->points.end(); ++pos2) {
        if (&(*pos1) == point1 && &(*pos2) == point2)
            return true;
        else if (&(*pos1) == point2 && &(*pos2) == point1)
            return true;

        pos1 = pos2;
    }

    return false;
}
/*=======================================================*/
bool GbSystem3D::PointSet3::containsEqualLine(const GbPoint3D &point1, const GbPoint3D &point2)
{
    // returns true if pointset has c2 in "this->points"vector  neighboured points
    // wich have the same coordinates as point1 or point2
    vector<GbPoint3D>::iterator pos1 = this->points.begin();
    vector<GbPoint3D>::iterator pos2;

    for (pos2 = pos1++; pos2 != this->points.end(); ++pos2) {
        if ((*pos1).equals(&point1) && (*pos2).equals(&point2))
            return true;
        else if ((*pos1).equals(&point2) && (*pos2).equals(&point1))
            return true;

        pos1 = pos2;
    }

    return false;
}
/*=======================================================*/
GbPoint3D *GbSystem3D::PointSet3::getPoint(int index)
{
    if (index < 0 || index >= (int)this->points.size())
        throw UbException(UB_EXARGS, "index out of range");
    return &(this->points)[index];
}
/*=======================================================*/
GbPoint3D *GbSystem3D::PointSet3::getFirstPoint() { return &(this->points.front()); }
/*=======================================================*/
GbPoint3D *GbSystem3D::PointSet3::getLastPoint() { return &(this->points.back()); }
/*=======================================================*/
int GbSystem3D::PointSet3::size() { return (int)this->points.size(); }
/*=======================================================*/
vector<GbPoint3D> GbSystem3D::PointSet3::getPoints()
{
    // is this right? it's another effect as at GbPoint3D* getNode(index)!!!
    // or do we want to have the next uncommented getPoints() funktion
    return this->points;
}
///*=======================================================*/
// vector<GbPoint3D*> GbSystem3D::PointSet3::getPoints()
//{
//   vector<GbPoint3D*> tmp;
//   for(int pos=0; pos<(int)this->points.size();pos++) tmp.push_back(&this->points[pos]);
//
//   return tmp;
//}
/*=======================================================*/
void GbSystem3D::PointSet3::calculateValues()
{
    if (this->points.empty()) {
        this->x1min = this->x2min = this->x3min = 0.0;
        this->x1max = this->x2max = this->x3max = 0.0;
    } else {
        this->x1min = (this->points)[0].x1;
        this->x1max = (this->points)[0].x1;
        this->x2min = (this->points)[0].x2;
        this->x2max = (this->points)[0].x2;
        this->x3min = (this->points)[0].x3;
        this->x3max = (this->points)[0].x3;

        for (int i = (int)this->points.size() - 1; i > 0; --i) {
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
    }
    this->consistent = true;
}

//! \}
