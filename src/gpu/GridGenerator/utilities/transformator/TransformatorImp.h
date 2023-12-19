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
#ifndef TransformatorImp_h
#define TransformatorImp_h

#include <exception>
#include <sstream>

#include "global.h"

#include "utilities/transformator/Transformator.h"
#include "utilities/transformator/ArrowTransformator.h"

class BoundingBox;
struct Triangle;
class TriangularMesh;
struct Vertex;

class invalidDelta : public std::exception
{
    public:
    invalidDelta() : error_message ("Delta cant be < Null. To enable no changes change delta to 1.0.")
    {}

    const char* what() const noexcept
    {
        return error_message.c_str();
    }

    private:
    std::string error_message;
};

class TransformatorImp
    : public Transformator, public ArrowTransformator
{
public:
    TransformatorImp();
    TransformatorImp(const TransformatorImp& trafo);
    TransformatorImp(real delta, const Vertex& translater);
    TransformatorImp(real delta, real dx, real dy, real dz);
    virtual ~TransformatorImp();
    
    void transformWorldToGrid(Triangle &value) const override;
    void transformWorldToGrid(TriangularMesh &geom) const override;
    void transformWorldToGrid(Vertex &value) const override;

    void transformGridToWorld(Triangle &t) const override;
    void transformGridToWorld(Vertex &value) const override;

    void transformGridToWorld(BoundingBox &box) const override;
    void transformWorldToGrid(BoundingBox &box) const override;

    bool operator==(const TransformatorImp& trafo) const;

    virtual void transformGridToWorld(std::shared_ptr<Arrow> arrow) const override;

private:
    real delta;
    std::shared_ptr<Vertex> translater;

    void scaleWorldToView(Vertex & v) const;
    void translateWorldToView(Vertex & v) const;

    void translateGridToWorld(Vertex & value) const;
    void scaleGridToWorld(Vertex & value) const;

    void verifyDelta(real delta) const;
};


#endif

//! \}
