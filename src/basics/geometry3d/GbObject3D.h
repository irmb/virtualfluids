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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file GbObject3D.h
//! \ingroup geometry3d
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef GBOBJECT3D_H
#define GBOBJECT3D_H

#include <string>
#include <vector>

#include <basics/objects/ObObject.h>
#include <basics/utilities/UbException.h>
#include <basics/utilities/UbObservable.h>
#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbTuple.h>

class GbPoint3D;
class GbLine3D;
class GbTriangle3D;
class GbObject3DCreator;

#include <PointerDefinitions.h>



//////////////////////////////////////////////////////////////////////////
//!
//! \class GbObject3D
//!
//! \brief This Interface provides basic 3D geometry objects methods.
//!
//////////////////////////////////////////////////////////////////////////

class GbObject3D : public ObObject
{
public:
    virtual ~GbObject3D() = default;
    // abstract Methods
    virtual void finalize() = 0; // destroys also all dynamic objects (e.g. GbPoints in GbLine)
    /**
     * Returns the centroid x1 coordinate of this 3D object.
     * @return the centroid x1 coordinate of this 3D object
     */
    virtual double getX1Centroid() = 0;
    /**
     * Returns the minimum x1 coordinate of this 3D object.
     * @return the minimum x1 coordinate of this 3D object
     */
    virtual double getX1Minimum() = 0;
    /**
     * Returns the maximum x1 coordinate of this 3D object.
     * @return the maximum x1 coordinate of this 3D object
     */
    virtual double getX1Maximum() = 0;
    /**
     * Returns the centroid x2 coordinate of this 3D object.
     * @return the centroid x2 coordinate of this 3D object
     */
    virtual double getX2Centroid() = 0;
    /**
     * Returns the minimum x2 coordinate of this 3D object.
     * @return the minimum x2 coordinate of this 3D object
     */
    virtual double getX2Minimum() = 0;
    /**
     * Returns the maximum x2 coordinate of this 3D object.
     * @return the maximum x2 coordinate of this 3D object
     */
    virtual double getX2Maximum() = 0;

    virtual double getX3Centroid() = 0;
    /**
     * Returns the minimum x2 coordinate of this 3D object.
     * @return the minimum x2 coordinate of this 3D object
     */
    virtual double getX3Minimum() = 0;
    /**
     * Returns the maximum x2 coordinate of this 3D object.
     * @return the maximum x2 coordinate of this 3D object
     */
    virtual double getX3Maximum() = 0;

    /*=======================================================*/
    double getLengthX1() { return (getX1Maximum() - getX1Minimum()); }
    double getLengthX2() { return (getX2Maximum() - getX2Minimum()); }
    double getLengthX3() { return (getX3Maximum() - getX3Minimum()); }

    virtual void setCenterX1Coordinate(const double & /*value*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual void setCenterX2Coordinate(const double & /*value*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual void setCenterX3Coordinate(const double & /*value*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual void setCenterCoordinates(const double & /*x1*/, const double & /*x2*/, const double & /*x3*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual void setCenterCoordinates(const UbTupleDouble3 & /*position*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }

    // Rotates the Point in relation to the origen.
    // Parameters must be radian measure.
    virtual void rotate(const double & /*rx1*/, const double & /*rx2*/, const double & /*rx3*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual void translate(const double & /*x1*/, const double & /*x2*/, const double & /*x3*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual void scale(const double & /*sx1*/, const double & /*sx2*/, const double & /*sx3*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }

    virtual bool isPointInGbObject3D(GbPoint3D *p);
    virtual bool isPointInGbObject3D(const double &x1, const double &x2, const double &x3, bool &pointIsOnBoundary) = 0;
    virtual bool isPointInGbObject3D(const double &x1, const double &x2, const double &x3)                          = 0;

    virtual bool isCellInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                        const double &x2b, const double &x3b);
    virtual bool isCellCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b);
    virtual bool isCellInsideOrCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a,
                                                 const double &x1b, const double &x2b, const double &x3b);
    virtual double getCellVolumeInsideGbObject3D(const double & /*x1a*/, const double & /*x2a*/, const double & /*x3a*/,
                                                 const double & /*x1b*/, const double & /*x2b*/, const double & /*x3b*/)
    {
        return -1.0;
    };

    virtual bool isInsideCell(const double &minX1, const double &minX2, const double &minX3, const double &maxX1,
                              const double &maxX2, const double &maxX3);

    virtual GbLine3D *createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2) = 0;
    virtual std::vector<GbTriangle3D *> getSurfaceTriangleSet()                 = 0;

    virtual void addSurfaceTriangleSet(std::vector<UbTupleFloat3> & /*nodes*/, std::vector<UbTupleInt3> & /*triangles*/)
    {
        throw UbException("GbObject3D::addSurfaceTriangleSet - not implemented for " +
                          (std::string) typeid(*this).name());
    }

    virtual bool hasRaytracing() { return false; }
    virtual bool raytracingSupportsPointsInside() { return false; }
    //|r| must be 1! einheitsvector!!
    // return negativ value oder zero if no intersection
    virtual double getIntersectionRaytraceFactor(const double & /*x1*/, const double & /*x2*/, const double & /*x3*/,
                                                 const double & /*rx1*/, const double & /*rx2*/, const double & /*rx3*/)
    {
        throw UbException("GbObject3D::getIntersectionRaytraceFactor - not implemented");
    }
};
/*=========================================================================*/

#endif
