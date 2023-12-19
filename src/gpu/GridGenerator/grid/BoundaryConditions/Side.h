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
//! \addtogroup gpu_grid grid
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef SIDE_H
#define SIDE_H

#include <cstddef>
#include <string>
#include <vector>
#include <map>
#include <array>

#include "gpu/GridGenerator/global.h"
#include "lbm/constants/D3Q27.h"

#define X_INDEX 0
#define Y_INDEX 1
#define Z_INDEX 2

#define POSITIVE_DIR 1
#define NEGATIVE_DIR -1

class Grid;

namespace gg
{
class BoundaryCondition;
}

class Side;

enum class SideType
{
    MX, PX, MY, PY, MZ, PZ, GEOMETRY
};

class Side
{
public:
    virtual ~Side() = default;
    virtual void addIndices(const std::vector<SPtr<Grid>> &grid, uint level,
                            SPtr<gg::BoundaryCondition> boundaryCondition) = 0;

    virtual int getCoordinate() const = 0;
    virtual int getDirection() const = 0;

    virtual SideType whoAmI() const = 0;

    std::array<real, 3> getNormal() const;

protected:
    void addIndices(SPtr<Grid> grid, SPtr<gg::BoundaryCondition> boundaryCondition, std::string coord, real constant,
                           real startInner, real endInner, real startOuter, real endOuter);

    static void setPressureNeighborIndices(SPtr<gg::BoundaryCondition> boundaryCondition, SPtr<Grid> grid, const uint index);

    static void setStressSamplingIndices(SPtr<gg::BoundaryCondition> boundaryCondition, SPtr<Grid> grid, const uint index);

    void setQs(SPtr<Grid> grid, SPtr<gg::BoundaryCondition> boundaryCondition, uint index);

    virtual void correctNeighborForPeriodicBoundaries(const Grid *grid, std::array<real, 3>& coords, std::array<real, 3>& neighbors) const;

    virtual bool isAlignedWithMyNormal(const Grid *grid, int dir) const;
    bool isAlignedWithNormal(const Grid *grid, int dir, const std::array<real, 3>& normal) const;

private:
    static uint getIndex(SPtr<Grid> grid, std::string coord, real constant, real v1, real v2);
    void resetDiagonalsInCaseOfOtherBC(Grid *grid, std::vector<real>& qNode, int dir, const std::array<real, 3> &coordinates) const;
    std::array<real, 3> getNeighborCoordinates(Grid *grid, const std::array<real, 3> &coordinates,
                                               size_t direction) const;
    bool neighborNormalToSideIsAStopper(Grid *grid, const std::array<real, 3> &coordinates, SideType side) const;

protected:
    const std::map<SideType, const std::array<real, 3>> normals = {
        { SideType::MX, { NEGATIVE_DIR, 0.0, 0.0 } }, { SideType::PX, { POSITIVE_DIR, 0.0, 0.0 } },
        { SideType::MY, { 0.0, NEGATIVE_DIR, 0.0 } }, { SideType::PY, { 0.0, POSITIVE_DIR, 0.0 } },
        { SideType::MZ, { 0.0, 0.0, NEGATIVE_DIR } }, { SideType::PZ, { 0.0, 0.0, POSITIVE_DIR } }
    };
    const std::map<SideType, size_t> sideToD3Q27 = {
        { SideType::MX, vf::lbm::dir::dM00 }, { SideType::PX, vf::lbm::dir::dP00 },
        { SideType::MY, vf::lbm::dir::d0M0 }, { SideType::PY, vf::lbm::dir::d0P0 },
        { SideType::MZ, vf::lbm::dir::d00M }, { SideType::PZ, vf::lbm::dir::d00P }
    };
};

class Geometry : public Side
{
public:
    void addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<gg::BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return X_INDEX;
    }

    int getDirection() const override
    {
        return NEGATIVE_DIR;
    }

    SideType whoAmI() const override
    {
        return SideType::GEOMETRY;
    }
};

class MX : public Side
{
public:
    void addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<gg::BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return X_INDEX;
    }

    int getDirection() const override
    {
        return NEGATIVE_DIR;
    }

    SideType whoAmI() const override
    {
        return SideType::MX;
    }
};

class PX : public Side
{
public:
    void addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<gg::BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return X_INDEX;
    }

    int getDirection() const override
    {
        return POSITIVE_DIR;
    }

    SideType whoAmI() const override
    {
        return SideType::PX;
    }
};


class MY : public Side
{
public:
    void addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<gg::BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return Y_INDEX;
    }

    int getDirection() const override
    {
        return NEGATIVE_DIR;
    }

    SideType whoAmI() const override
    {
        return SideType::MY;
    }
};

class PY : public Side
{
public:
    void addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<gg::BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return Y_INDEX;
    }

    int getDirection() const override
    {
        return POSITIVE_DIR;
    }

    SideType whoAmI() const override
    {
        return SideType::PY;
    }
};


class MZ : public Side
{
public:
    void addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<gg::BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return Z_INDEX;
    }

    int getDirection() const override
    {
        return NEGATIVE_DIR;
    }

    SideType whoAmI() const override
    {
        return SideType::MZ;
    }
};

class PZ : public Side
{
public:
    void addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<gg::BoundaryCondition> boundaryCondition) override;

    int getCoordinate() const override
    {
        return Z_INDEX;
    }

    int getDirection() const override
    {
        return POSITIVE_DIR;
    }

    SideType whoAmI() const override
    {
        return SideType::PZ;
    }
};


class SideFactory
{
public:
    static SPtr<Side> make(SideType sideType)
    {
        switch (sideType)
        {
        case SideType::MX:
            return SPtr<Side>(new MX());
        case SideType::PX:
            return SPtr<Side>(new PX());
        case SideType::MY:
            return SPtr<Side>(new MY());
        case SideType::PY:
            return SPtr<Side>(new PY());
        case SideType::MZ:
            return SPtr<Side>(new MZ());
        case SideType::PZ:
            return SPtr<Side>(new PZ());
        case SideType::GEOMETRY:
            return SPtr<Side>(new Geometry());
        default:
            throw std::runtime_error("SideFactory::make() - SideType not valid.");
        }
    }
};

#endif

//! \}
