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
//! \file Sphere.h
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef SPHERE_H
#define SPHERE_H

#include "PointerDefinitions.h"
#include "global.h"
#include "geometries/Object.h"

class GRIDGENERATOR_EXPORT Sphere : public Object
{
public:
    Sphere(const double& centerX, const double& centerY, const double& centerZ, const double& radius);

    static SPtr<Sphere> makeShared(double centerX, double centerY, double centerZ, double radius);

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

    bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) override;


    void scale(double delta) override;
    
    int getIntersection(const Vertex &P, const Vertex &direction, Vertex &pointOnObject, real &qVal) override;


protected:
    double centerX;
    double centerY;
    double centerZ;

    double radius;
};



#endif   
