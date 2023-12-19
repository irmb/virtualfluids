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
#ifndef GbHalfSpaceKrischan3D_H
#define GbHalfSpaceKrischan3D_H

#include <iostream>
#include <sstream>

#include <basics/utilities/UbMath.h>

#include <geometry3d/GbLine3D.h>
#include <geometry3d/GbPoint3D.h>
#include <geometry3d/GbTriangle3D.h>
#include <geometry3d/GbVector3D.h>

/*=========================================================================*/
/* GbHalfSpaceKrischan3D                                                             */
/*                                                                         */
/**
 * This Class helps in performing some operations on a halfspace defined by 2 or 3 points
 */

class GbHalfSpaceKrischan3D : public GbObject3D, public UbObserver
{
public:
    GbHalfSpaceKrischan3D(GbTriangle3D *triangle);

    GbHalfSpaceKrischan3D(GbPoint3D *PointA, GbPoint3D *PointB, GbPoint3D *PointC);

    GbHalfSpaceKrischan3D(double nx, double ny, double nz, double dist);

    GbHalfSpaceKrischan3D(GbPoint3D *PointA, GbPoint3D *PointB);

    GbHalfSpaceKrischan3D(const double &p1x, const double &p1y, const double &p1z, const double &p2x, const double &p2y,
                          const double &p2z, const double &p3x, const double &p3y, const double &p3z);

    /*=======================================================*/
    ~GbHalfSpaceKrischan3D() override = default;
    /*=======================================================*/
    std::string getTypeID() { return "GbHalfSpaceKrischan3D"; }
    /*=============================================*/
    bool ptInside(const double &x, const double &y, const double &z)
    {
        return UbMath::lessEqual(Normal[0] * x + Normal[1] * y + Normal[2] * z, this->d);
    }
    /*=============================================*/
    bool ptInside(GbPoint3D *PointX)
    {
        GbVector3D X(PointX->x1, PointX->x2, PointX->x3);
        return UbMath::lessEqual(this->Normal.Dot(X), this->d);
    }
    /*=============================================*/
    bool ptInside(GbVector3D &X) { return UbMath::lessEqual(this->Normal.Dot(X), this->d); }

    /*=====================================================*/
    // true, wenn 'in Object' oder 'auf Boundary'!
    bool isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p) override
    {
        return (ptInside(x1p, x2p, x3p));
    }
    /*=====================================================*/
    // true, wenn 'in Object' oder 'auf Boundary'!
    bool isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p,
                             bool & /*pointIsOnBoundary*/) override
    {
        return (ptInside(x1p, x2p, x3p));
    }

    void finalize() override {}

    double getX1Centroid() override { return 0.0; }
    double getX1Minimum() override { return -99999.0; }
    double getX1Maximum() override { return 99999.0; }
    double getX2Centroid() override { return 0.0; }
    double getX2Minimum() override { return -99999.0; }
    double getX2Maximum() override { return 99999.0; }
    double getX3Centroid() override { return 0.0; }
    double getX3Minimum() override { return -99999.0; }
    double getX3Maximum() override { return 99999.0; }

    GbLine3D *createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2) override
    {
        GbPoint3D *p1 = new GbPoint3D(point1);
        GbPoint3D *p2 = new GbPoint3D(point2);

        GbVector3D p1p2(p2->x1 - p1->x1, p2->x2 - p1->x2, p2->x3 - p1->x3);

        double dist1 = getDistance(p1->x1, p1->x2, p1->x3);
        double dist2 = getDistance(p2->x1, p2->x2, p2->x3);

        double totalDist = std::abs(dist1) + std::abs(dist2);

        // Falls erster Punkt nicht drinliegt
        if (!ptInside(p1)) {
            if (!ptInside(p2))
                return NULL;

            // distance ausrechnen (groesser null)
            if (UbMath::less(dist1, 0.0))
                throw UbException(UB_EXARGS, "Punkt ausserhalb, aber Distanz kleiner null???");

            p1->x1 = p1->x1 + dist1 / totalDist * p1p2[0];
            p1->x2 = p1->x2 + dist1 / totalDist * p1p2[1];
            p1->x3 = p1->x3 + dist1 / totalDist * p1p2[2];
        }
        // Falls zweiter Punkt nicht drinliegt
        if (!ptInside(p2)) {
            if (!ptInside(p1))
                return NULL;

            // distance ausrechnen (groesser null)
            if (UbMath::less(dist2, 0.0))
                throw UbException(UB_EXARGS, "Punkt ausserhalb, aber Distanz kleiner null???");

            p2->x1 = p2->x1 - dist2 / totalDist * p1p2[0];
            p2->x2 = p2->x2 - dist2 / totalDist * p1p2[1];
            p2->x3 = p2->x3 - dist2 / totalDist * p1p2[2];
        }

        return new GbLine3D(p1, p2);
    }

    double getCellVolumeInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b) override;

    double getDistance(const double &x1p, const double &x2p, const double &x3p)
    {
        return (Normal[0] * x1p + Normal[1] * x2p + Normal[2] * x3p) - this->d;
    }

    void getNormal(double &n1, double &n2, double &n3)
    {
        n1 = this->Normal[0];
        n2 = this->Normal[1];
        n3 = this->Normal[2];
    }

    void addSurfaceTriangleSet(std::vector<UbTupleFloat3> & /*nodes*/,
                               std::vector<UbTupleInt3> & /*triangles*/) override
    {
        std::cout << " addSurfaceTriangleSet(): TO BE DONE AND CHECKED ... " << std::endl;
    }

    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override
    {
        std::vector<GbTriangle3D *> triangles;
        GbPoint3D p1(0.0, 0.0, 0.0);
        GbPoint3D p2(1.0, 0.0, 0.0);
        GbPoint3D p3(0.0, 1.0, 0.0);

        triangles.push_back(new GbTriangle3D(new GbPoint3D(p1), new GbPoint3D(p2), new GbPoint3D(p3)));

        return triangles;
    }

    void objectChanged(UbObservable * /*changedObject*/) override
    {
        return;

        // GbLine3D* line = dynamic_cast<GbLine3D*>(changedObject);
        // if(!line || this->mLine!=line) return;
        // this->notifyObserversObjectChanged();
    }
    /*==========================================================*/
    void objectWillBeDeleted(UbObservable * /*objectForDeletion*/) override
    {
        return;
        // if(this->mLine)
        //{
        //   UbObservable* observedObj = dynamic_cast<UbObservable*>(this->mLine);
        //   if(objectForDeletion == observedObj) { this->mLine = NULL; }
        //}
    }

    ObObject *clone() override { return NULL; };

    std::string toString() override
    {
        std::stringstream temp;

        temp << "GbHalfSpaceKrischan3D:   ";
        temp << " Distance   " << this->d;
        temp << " Norm vec   " << this->Normal[0];
        temp << " " << this->Normal[1];
        temp << " " << this->Normal[2];

        return temp.str();
    };

private:
    GbVector3D Normal;
    double d;
};
/*=========================================================================*/

#endif // GbHalfSpaceKrischan3D_H

//! \}
