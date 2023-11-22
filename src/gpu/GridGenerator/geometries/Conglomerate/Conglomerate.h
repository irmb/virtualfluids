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
//! \file Conglomerate.h
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef CONGLOMERATE_H
#define CONGLOMERATE_H

#include <array>

#include "global.h"

#include "geometries/Object.h"
#include "basics/PointerDefinitions.h"

#define MAX_NUMBER_OF_OBJECTS 20

class GRIDGENERATOR_EXPORT Conglomerate : public Object
{
public:
    static SPtr<Conglomerate> makeShared();

    void add(SPtr<Object> object);
    void subtract(SPtr<Object> objectStub);


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

    void findInnerNodes(SPtr<GridImp> grid) override;

protected:
    static double getMinimum(double val1, double val2);
    static double getMaximum(double val1, double val2);


    std::array<SPtr<Object>, MAX_NUMBER_OF_OBJECTS> addObjects;
    std::array<SPtr<Object>, MAX_NUMBER_OF_OBJECTS> subtractObjects;
    uint numberOfAddObjects = 0;
    uint numberOfSubtractObjects = 0;
};



#endif   
