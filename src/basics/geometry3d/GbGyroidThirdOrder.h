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
//! \file GbGyroidThirdOrder.h
//! \ingroup geometry3d
//! \author Hussein Alihussein
//=======================================================================================
#ifndef GbGyroidThirdOrder_H
#define GbGyroidThirdOrder_H

#ifdef BUILD_USE_BOOST

#include <vector>

#include <GbPoint3D.h>
#include <basics/utilities/UbObserver.h>
#include <basics/utilities/UbMath.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h> 

class GbLine3D;
class GbObject3DCreator;

#include <PointerDefinitions.h>
class GbGyroidThirdOrder;
using GbGyroidThirdOrderPtr = SPtr<GbGyroidThirdOrder>;


class GbGyroidThirdOrder : public GbObject3D, public UbObserver
{
public:
    GbGyroidThirdOrder();
    GbGyroidThirdOrder(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b, const double& edgeLength, const double& dx, const double& thickness=0);

    GbGyroidThirdOrder(const double & x1a, const double & x2a, const double & x3a, const double & x1b, const double & x2b, const double & x3b, const double & x1c, const double & x2c, const double & x3c, const double & x1d, const double & x2d, const double & x3d, const double & edgeLength, const double & dx);
    GbGyroidThirdOrder(GbGyroidThirdOrder *imp);
    ~GbGyroidThirdOrder();

    GbGyroidThirdOrder* clone() override { return new GbGyroidThirdOrder(this); }
    void finalize() override { throw UbException(UB_EXARGS, "finalize() - not implemented"); }

    double getX1Centroid() override;
    double getX1Minimum() override;
    double getX1Maximum() override;
    double getX2Centroid()override;
    double getX2Minimum() override;
    double getX2Maximum() override;
    double getX3Centroid()override;
    double getX3Minimum() override;
    double getX3Maximum() override;
    void setCenterCoordinates(const double &x1, const double &x2, const double &x3) override {throw UbException(UB_EXARGS, "finalize() - not implemented");
    }

    void translate(const double& x1, const double& x2, const double& x3) override { throw UbException(UB_EXARGS, "finalize() - not implemented"); }
    void rotate(const double& rx1, const double& rx2, const double& rx3) override{ throw UbException(UB_EXARGS, "finalize() - not implemented"); }
    void scale(const double& sx1, const double& sx2, const double& sx3) override { throw UbException(UB_EXARGS, "finalize() - not implemented"); }

    double getLengthX1();
    double getLengthX2();
    double getLengthX3();
    
    bool isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p, bool &pointinboundary) override { throw UbException(UB_EXARGS, "finalize() - not implemented"); }
    bool isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p) override;
    bool isCellInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                const double &x2b, const double &x3b) override;
    bool isCellCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                 const double &x2b, const double &x3b) override;
    bool isCellInsideOrCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b) override;
    double getCellVolumeInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b) override { throw UbException(UB_EXARGS, "finalize() - not implemented"); }
    GbPoint3D *calculateInterSectionPoint3D(GbPoint3D &point1, GbPoint3D &point2);
    GbLine3D *createClippedLine3D (GbPoint3D &point1, GbPoint3D &point2) override { throw UbException(UB_EXARGS, "finalize() - not implemented"); }
    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override { throw UbException(UB_EXARGS, "finalize() - not implemented"); }
    void addSurfaceTriangleSet(std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt3> &triangles) override { throw UbException(UB_EXARGS, "finalize() - not implemented"); }
    bool hasRaytracing() override { return true;  }
    double getIntersectionRaytraceFactor (const double& x1, const double& x2, const double& x3, const double& rx1, const double& rx2, const double& rx3) override;
    double evaluateImplicitFunction(const double & x1, const double & x2, const double & x3, const double & position);
    double getDistance(const double& x1p, const double& x2p, const double& x3p) { throw UbException(UB_EXARGS, "finalize() - not implemented"); }    
    std::string toString() override { throw UbException(UB_EXARGS, "finalize() - not implemented"); }
    void objectChanged(UbObservable *changedObject) override;
    void objectWillBeDeleted(UbObservable *objectForDeletion) override;
    using GbObject3D::isPointInGbObject3D; 

protected:
    GbPoint3D* p1;
    GbPoint3D* p2;
    GbPoint3D* p3;
    GbPoint3D* p4;
    double edgeLength;
    double dx;
    double thickness;
private:
};



#endif   
#endif
