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
//! \file VerticalCylinder.cpp
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "VerticalCylinder.h"

VerticalCylinder::VerticalCylinder(const double& centerX, const double& centerY, const double& centerZ, const double& radius, const double& height)
    : centerX(centerX), centerY(centerY), centerZ(centerZ), radius(radius), height(height)
{

}

VerticalCylinder::~VerticalCylinder()
{
}

SPtr<VerticalCylinder> VerticalCylinder::makeShared(double centerX, double centerY, double centerZ, double radius, double height)
{
    return SPtr<VerticalCylinder>(new VerticalCylinder(centerX, centerY, centerZ, radius, height));
}

Object* VerticalCylinder::clone() const
{
    return new VerticalCylinder(centerX, centerY, centerZ, radius, height);
}

double VerticalCylinder::getX1Centroid()
{
    return centerX;
}

double VerticalCylinder::getX1Minimum()
{
    return centerX - radius;
}

double VerticalCylinder::getX1Maximum()
{
    return centerX + radius;
}

double VerticalCylinder::getX2Centroid()
{
    return centerY;
}

double VerticalCylinder::getX2Minimum()
{
    return centerY - radius;
}

double VerticalCylinder::getX2Maximum()
{
    return centerY + radius;
}

double VerticalCylinder::getX3Centroid()
{
    return centerZ;
}

double VerticalCylinder::getX3Minimum()
{
    return centerZ - 0.5 * height;
}

double VerticalCylinder::getX3Maximum()
{
    return centerZ + 0.5 * height;
}

bool VerticalCylinder::isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset)
{
    double offset = maxOffset;
    if (x1 < centerX || x2 < centerY || x3 < centerZ)
        offset = minOffset;
        

    const double deltaX1 = x1 - centerX;
    const double deltaX2 = x2 - centerY;
    const double deltaX3 = x3 - centerZ;

    if( deltaX3 > 0.5 * height || deltaX3 < - 0.5 * height )
        return false;

    return (deltaX1*deltaX1 + deltaX2*deltaX2) < ((this->radius - offset) * (this->radius - offset));
}


void VerticalCylinder::scale(double delta)
{
    this->radius += delta;
    this->height += delta;
}
