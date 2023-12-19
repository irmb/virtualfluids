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
//! \addtogroup gpu_geometries geometries
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef CUBOID_H
#define CUBOID_H

#include "global.h"

#include "geometries/Object.h"

class Cuboid : public Object
{
public:              
    Cuboid(const double& minX1, const double& minX2, const double& minX3, const double& maxX1,const double& maxX2, const double& maxX3);

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

    void changeSizeByDelta(double delta) override;

    bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) override;

private:
    static double getCenter(double x1, double x2);
    static double getMinimum(double x1, double x2);
    static double getMaximum(double x1, double x2);
    static bool isOn(const real& coord, const real& plane1, const real& plane2);
    static bool isBetween(const real& coord, const real& start, const real& end);

protected:
    double minX1;
    double minX2;
    double minX3;
    double maxX1;
    double maxX2;
    double maxX3;
};



#endif   

//! \}
