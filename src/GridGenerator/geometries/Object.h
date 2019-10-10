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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \file Object.h
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef OBJECT_H
#define OBJECT_H

#include <VirtualFluidsDefinitions.h>
#include "grid/Cell.h"
#include "global.h"

class GridImp;

class VF_PUBLIC Object
{
public:
    virtual ~Object() {}
    virtual Object* clone() const = 0;

    virtual double getX1Centroid() = 0;
    virtual double getX1Minimum()  = 0;
    virtual double getX1Maximum()  = 0;

    virtual double getX2Centroid() = 0;
    virtual double getX2Minimum()  = 0;
    virtual double getX2Maximum()  = 0;

    virtual double getX3Centroid() = 0;
    virtual double getX3Minimum()  = 0;
    virtual double getX3Maximum()  = 0;


    virtual void scale(double delta) = 0;


    virtual bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) = 0;

    virtual bool isCellInObject(const Cell& cell) {
        for (const auto point : cell)
        {
            const bool isInObject = isPointInObject(point.x, point.y, point.z, 0.0, 0.0);
            if (!isInObject)
                return false;
        }
        return true;
    }

    virtual void findInnerNodes(SPtr<GridImp> grid);
};


#endif
