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
#include "BoundingBox.h"

#include "../Triangle/Triangle.h"
#include "../Vertex/Vertex.h"
#include <GridGenerator/utilities/math/Math.h>

#include <limits>



 BoundingBox::BoundingBox(real minX, real maxX, real minY, real maxY, real minZ, real maxZ) : minX(minX), maxX(maxX), minY(minY), maxY(maxY), minZ(minZ), maxZ(maxZ) {}

 BoundingBox BoundingBox::makeInvalidMinMaxBox()
 {
     BoundingBox box = BoundingBox(std::numeric_limits<real>::max(),
         std::numeric_limits<real>::lowest(),
         std::numeric_limits<real>::max(),
         std::numeric_limits<real>::lowest(),
         std::numeric_limits<real>::max(),
         std::numeric_limits<real>::lowest());
     return box;
 }

 void BoundingBox::setMinMax(const Triangle &t)
 {
     real minX, maxX, minY, maxY, minZ, maxZ;
     t.setMinMax(minX, maxX, minY, maxY, minZ, maxZ);
     if (minX < this->minX)
         this->minX = minX;
     if (minY < this->minY)
         this->minY = minY;
     if (minZ < this->minZ)
         this->minZ = minZ;

     if (maxX > this->maxX)
         this->maxX = maxX;
     if (maxY > this->maxY)
         this->maxY = maxY;
     if (maxZ > this->maxZ)
         this->maxZ = maxZ;
 }

 bool BoundingBox::intersect(const Triangle &t) const
 {
     if (isInside(t.v1) || isInside(t.v2) || isInside(t.v3))
         return true;
     return false;
 }

 bool BoundingBox::isInside(const Triangle &t) const
 {
     if (isInside(t.v1) && isInside(t.v2) && isInside(t.v3))
         return true;
     return false;
 }

 bool BoundingBox::isInside(const real x, const real y, const real z) const
 {
     return this->isInside(Vertex(x,y,z));
 }

 bool BoundingBox::isInside(const Vertex &v) const
 {
     if (v.isXbetween(minX, maxX) && v.isYbetween(minY, maxY) && v.isZbetween(minZ, maxZ))
         return true;
     return false;
 }

 std::vector<std::vector<Vertex> > BoundingBox::getIntersectionPoints(const BoundingBox &b) const
 {
     std::vector<std::vector<Vertex> > intersectionBox;
     intersectionBox.resize(6);

     int intersects = 0;
     if (b.minX < maxX && b.maxX > maxX) { //maxX is intersect
         intersectionBox[intersects].push_back(Vertex((real)maxX, (real)minY, (real)minZ));
         intersectionBox[intersects].push_back(Vertex((real)maxX, (real)maxY, (real)minZ));
         intersectionBox[intersects].push_back(Vertex((real)maxX, (real)minY, (real)maxZ));
         intersects++;
     }
     if (b.minX < minX && b.maxX > minX) { //minX is intersect
         intersectionBox[intersects].push_back(Vertex((real)minX, (real)minY, (real)minZ));
         intersectionBox[intersects].push_back(Vertex((real)minX, (real)maxY, (real)minZ));
         intersectionBox[intersects].push_back(Vertex((real)minX, (real)minY, (real)maxZ));
         intersects++;
     }
     if (b.minY < minY && b.maxY > minY) { //minY is intersect
         intersectionBox[intersects].push_back(Vertex((real)minX, (real)minY, (real)minZ));
         intersectionBox[intersects].push_back(Vertex((real)maxX, (real)minY, (real)minZ));
         intersectionBox[intersects].push_back(Vertex((real)minX, (real)minY, (real)maxZ));
         intersects++;
     }
     if (b.minY < maxY && b.maxY > maxY) { //maxY is intersect
         intersectionBox[intersects].push_back(Vertex((real)minX, (real)maxY, (real)minZ));
         intersectionBox[intersects].push_back(Vertex((real)maxX, (real)maxY, (real)minZ));
         intersectionBox[intersects].push_back(Vertex((real)minX, (real)maxY, (real)maxZ));
         intersects++;
     }
     if (b.minZ < minZ && b.maxZ > minZ) { //minZ is intersect
         intersectionBox[intersects].push_back(Vertex((real)minX, (real)minY, (real)minZ));
         intersectionBox[intersects].push_back(Vertex((real)maxX, (real)minY, (real)minZ));
         intersectionBox[intersects].push_back(Vertex((real)minX, (real)maxY, (real)minZ));
         intersects++;
     }
     if (b.minZ < maxZ && b.maxZ > maxZ) { //maxZ is intersect
         intersectionBox[intersects].push_back(Vertex((real)minX, (real)minY, (real)maxZ));
         intersectionBox[intersects].push_back(Vertex((real)maxX, (real)minY, (real)maxZ));
         intersectionBox[intersects].push_back(Vertex((real)minX, (real)maxY, (real)maxZ));
         intersects++;
     }

     return intersectionBox;
 }

 bool BoundingBox::intersect(const BoundingBox &box) const
 {
     struct Vertex v[8];
     box.getPoints(v);

     for (int i = 0; i < 8; i++) {
         if (isInside(v[i]))
             return true;
     }
     return false;
 }

 void BoundingBox::getPoints(Vertex v[8]) const
 {
     v[0] = Vertex(minX, minY, minZ);
     v[1] = Vertex(maxX, minY, minZ);
     v[2] = Vertex(minX, maxY, minZ);
     v[3] = Vertex(maxX, maxY, minZ);

     v[4] = Vertex(minX, minY, maxZ);
     v[5] = Vertex(maxX, minY, maxZ);
     v[6] = Vertex(minX, maxY, maxZ);
     v[7] = Vertex(maxX, maxY, maxZ);
 }


 void BoundingBox::print() const
 {
     printf("min/max - x: %2.4f/ %2.4f, y: %2.4f, %2.4f, z: %2.4f, %2.4f \n", minX, maxX, minY, maxY, minZ, maxZ);
 }


 bool BoundingBox::operator==(const BoundingBox &box) const
 {
     return vf::Math::equal(minX, box.minX)
         && vf::Math::equal(maxX, box.maxX)
         && vf::Math::equal(minY, box.minY)
         && vf::Math::equal(maxY, box.maxY)
         && vf::Math::equal(minZ, box.minZ)
         && vf::Math::equal(maxZ, box.maxZ);
 }

 void BoundingBox::extend(real delta)
 {
     this->minX -= delta;
     this->minY -= delta;
     this->minZ -= delta;

     this->maxX += delta;
     this->maxY += delta;
     this->maxZ += delta;
 }



//! \}
