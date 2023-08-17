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
//! \file Cylinder.h
//! \ingroup geometries
//! \author Anna Wellmann
//=======================================================================================


#ifndef CYLINDER_H
#define CYLINDER_H

#include "geometries/Object.h"
#include <basics/geometry3d/GbCylinder3D.h>
#include <map>

class GRIDGENERATOR_EXPORT Cylinder : public Object
{
public:
    enum PrincipalAxis {
        x = 0,
        y = 1,
        z = 2,
    };

    Cylinder(double centerX, double centerY, double centerZ, double radius, double height, PrincipalAxis axis);
    Cylinder(std::array<double, 3> center, double radius, double height, PrincipalAxis axis);

    SPtr<Object> clone() const override;

    double getX1Centroid() override;
    double getX1Minimum() override;
    double getX1Maximum() override;
    double getX2Centroid() override;
    double getX2Minimum() override;
    double getX2Maximum() override;
    double getX3Centroid() override;
    double getX3Minimum() override;
    double getX3Maximum() override;

    double getRadius() const;
    double getHeight() const;

    bool isPointInObject(const double &x1, const double &x2, const double &x3, const double &minOffset,
                         const double &maxOffset) override;
    void changeSizeByDelta(double delta) override;

private:
    double getCentroidCoordinate(PrincipalAxis coordinateDirection) const;
    double getMinimunCoordinate(PrincipalAxis coordinateDirection) const;
    double getMaximumCoordinate(PrincipalAxis coordinateDirection) const;

    bool isInCircle(double delta1, double delta2, double offset) const;

    const std::map<PrincipalAxis, std::array<double, 3>> unitVectors{ { x, { 1, 0, 0 } },
                                                                      { y, { 0, 1, 0 } },
                                                                      { z, { 0, 0, 1 } } };

    PrincipalAxis principalAxis;
    const std::array<double, 3> center;

    double radius;
    double height;
};

#endif