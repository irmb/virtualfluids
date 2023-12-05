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
//! \file DistributionDebugInspector.h
//! \author Henrik Asmuth
//! \date 13/012/2022
//! \brief Basic debugging class to print out f's in a certain area of the domain
//!
//! Basic debugging class. Needs to be directly added in UpdateGrid (could potentially also be added as a proper Probe in the
//! future) How to use: Define a part of the domain via min/max x, y, and z. The DistributionDebugInspector will print out
//! all f's in that area.
//!
//=======================================================================================

#ifndef DISTRIBUTION_INSPECTOR_H
#define DISTRIBUTION_INSPECTOR_H

#include <basics/DataTypes.h>

class Parameter;

class DistributionDebugInspector
{
public:
    DistributionDebugInspector(uint _inspectionLevel, real _minX, real _maxX, real _minY, real _maxY, real _minZ, real _maxZ,
                               std::string _tag)
        : inspectionLevel(_inspectionLevel), minX(_minX), maxX(_maxX), minY(_minY), maxY(_maxY), minZ(_minZ), maxZ(_maxZ),
          tag(_tag) {};

    void inspect(Parameter* para, uint level, uint t);

private:
    uint inspectionLevel;
    real minX;
    real maxX;
    real minY;
    real maxY;
    real minZ;
    real maxZ;
    std::string tag;
};

#endif
