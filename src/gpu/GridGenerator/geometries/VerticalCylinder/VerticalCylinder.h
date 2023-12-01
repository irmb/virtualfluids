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
//! \file VerticalCylinder.h
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef VERTICAL_CYLINDER_H
#define VERTICAL_CYLINDER_H

#include "global.h"
#include "geometries/Object.h"

class VerticalCylinder : public Object
{
public:
    VerticalCylinder(const double& centerX, const double& centerY, const double& centerZ, const double& radius, const double& height);
    // the axis is in the z-direction

    static SPtr<VerticalCylinder> makeShared(double centerX, double centerY, double centerZ, double radius, double height);

    SPtr<Object> clone() const override;

    double getX1Centroid() const override;
    double getX1Minimum() const override;
    double getX1Maximum() const override;
    double getX2Centroid() const override;
    double getX2Minimum() const override;
    double getX2Maximum() const override;
    double getX3Centroid() const override;
    double getX3Minimum() const override;
    double getX3Maximum() const override;

    double getRadius() const;
    double getHeight() const;

    bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) override;


    void changeSizeByDelta(double delta) override;

protected:
    double centerX;
    double centerY;
    double centerZ;

    double radius;
    double height;
};



#endif   
