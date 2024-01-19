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
//! \addtogroup gpu_utilities utilities
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "TransformatorImp.h"

#include <memory>

#include "geometries/BoundingBox/BoundingBox.h"
#include "geometries/Triangle/Triangle.h"
#include "geometries/TriangularMesh/TriangularMesh.h"
#include "geometries/Vertex/Vertex.h"
#include "geometries/Arrow/Arrow.h"

TransformatorImp::TransformatorImp() 
{
    this->translater = std::shared_ptr<Vertex>(new Vertex());
    this->delta = 1.0;
    this->translater->x = 0;
    this->translater->y = 0;
    this->translater->z = 0;
}

TransformatorImp::TransformatorImp(real delta, const Vertex& translater) : delta(delta), translater(std::make_shared<Vertex>(translater))
{
    this->verifyDelta(delta);
}

TransformatorImp::TransformatorImp(real delta, real dx, real dy, real dz) : TransformatorImp(delta, Vertex(dx,dy,dz))
{

}

TransformatorImp::TransformatorImp(const TransformatorImp& trafo)
{
    this->delta = trafo.delta;
    this->translater = std::shared_ptr<Vertex>(new Vertex(*trafo.translater.get()));
}

TransformatorImp::~TransformatorImp()
{

}

void TransformatorImp::transformWorldToGrid(TriangularMesh &geom) const
{
    for (int i = 0; i < geom.size; i++)
        transformWorldToGrid(geom.triangleVec[i]);
}

void TransformatorImp::transformWorldToGrid(Triangle &value) const
{
    transformWorldToGrid(value.v1);
    transformWorldToGrid(value.v2);
    transformWorldToGrid(value.v3);
}

void TransformatorImp::transformGridToWorld(std::shared_ptr<Arrow> arrow) const
{
    transformGridToWorld(*arrow->getStart());
    transformGridToWorld(*arrow->getEnd());
}

void TransformatorImp::transformWorldToGrid(Vertex &v) const
{
    translateWorldToView(v);
    scaleWorldToView(v);
}


void TransformatorImp::translateWorldToView(Vertex& v) const
{
    v = *translater.get() + v;
}

void TransformatorImp::scaleWorldToView(Vertex& v) const
{
    v = v * (1.0f / this->delta);
}


void TransformatorImp::transformGridToWorld(Triangle & t) const
{
    transformGridToWorld(t.v1);
    transformGridToWorld(t.v2);
    transformGridToWorld(t.v3);
}

void TransformatorImp::transformGridToWorld(Vertex &value) const
{
    scaleGridToWorld(value);
    translateGridToWorld(value);
}

void TransformatorImp::scaleGridToWorld(Vertex & value) const
{
    value = value * this->delta;
}


void TransformatorImp::translateGridToWorld(Vertex & value) const
{
    value = value - *translater.get();
}


void TransformatorImp::transformGridToWorld(BoundingBox &box) const
{
    //scale
    box.minX = (box.minX * this->delta);
    box.minY = (box.minY * this->delta);
    box.minZ = (box.minZ * this->delta);

    box.maxX = (box.maxX * this->delta);
    box.maxY = (box.maxY * this->delta);
    box.maxZ = (box.maxZ * this->delta);

    //translate
    box.minX = (box.minX - this->translater->x);
    box.minY = (box.minY - this->translater->y);
    box.minZ = (box.minZ - this->translater->z);

    box.maxX = (box.maxX - this->translater->x);
    box.maxY = (box.maxY - this->translater->y);
    box.maxZ = (box.maxZ - this->translater->z);
}

void TransformatorImp::transformWorldToGrid(BoundingBox &box) const
{
    //translate
    box.minX += this->translater->x;
    box.minY += this->translater->y;
    box.minZ += this->translater->z;

    box.maxX += this->translater->x;
    box.maxY += this->translater->y;
    box.maxZ += this->translater->z;

    //scale
    box.minX *= (1.0f / this->delta);
    box.minY *= (1.0f / this->delta);
    box.minZ *= (1.0f / this->delta);

    box.maxX *= (1.0f / this->delta);
    box.maxY *= (1.0f / this->delta);
    box.maxZ *= (1.0f / this->delta);
}


void TransformatorImp::verifyDelta(real delta) const
{
    if (delta <= 0.0)
        throw invalidDelta();
}

bool TransformatorImp::operator==(const TransformatorImp& trafo) const
{
    return (this->delta == trafo.delta
        && this->translater->x == trafo.translater->x
        && this->translater->y == trafo.translater->y
        && this->translater->z == trafo.translater->z);
}

//! \}
