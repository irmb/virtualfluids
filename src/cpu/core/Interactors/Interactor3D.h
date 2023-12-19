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
//! \addtogroup cpu_Interactors Interactors
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef INTERACTOR3D_H
#define INTERACTOR3D_H

#include <PointerDefinitions.h>
#include <vector>

#include "UbSystem.h"
#include "UbTuple.h"
#include "lbm/constants/D3Q27.h"

class Block3D;
class Grid3D;
class GbObject3D;

//! A base class for grid generation.
class Interactor3D : public enableSharedFromThis<Interactor3D>
{
public:
    enum Accuracy { SIMPLE, EDGES, FACES, POINTS };
    Interactor3D();
    Interactor3D(SPtr<Grid3D> grid, int type = Interactor3D::SOLID);
    Interactor3D(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type);
    //! constructor
    //! \param a set accuracy for arePointsInObject() and arePointsNotInObject()
    Interactor3D(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type, Interactor3D::Accuracy a);

    virtual ~Interactor3D();
    virtual void initInteractor(const real &timestep = 0) = 0;
    virtual void updateInteractor(const real &timestep = 0) = 0;

    void setSolidBlock(SPtr<Block3D> block);
    void setBCBlock(SPtr<Block3D> block);

    virtual UbTupleDouble3 getForces();

    void setSolid() { UbSystem::setBit(this->type, SOLID); }
    void setMoveable() { UbSystem::setBit(this->type, MOVEABLE); }

    bool isSolid() { return UbSystem::bitCheck(this->type, SOLID); }
    bool isInverseSolid() { return UbSystem::bitCheck(this->type, INVERSESOLID); }
    bool isTimeDependent() { return UbSystem::bitCheck(this->type, TIMEDEPENDENT); }
    bool isMoveable() { return UbSystem::bitCheck(this->type, MOVEABLE); }

    SPtr<Grid3D> getGrid3D() const { return grid.lock(); }
    void setGrid3D(SPtr<Grid3D> grid) { this->grid = grid; }
    virtual SPtr<GbObject3D> getGbObject3D() const { return geoObject3D; }
    virtual bool setDifferencesToGbObject3D(const SPtr<Block3D>) = 0;
    virtual std::vector<SPtr<Block3D>> &getBcBlocks() { return this->bcBlocks; }
    virtual void removeBcBlocks() { this->bcBlocks.clear(); }
    virtual std::vector<SPtr<Block3D>> &getSolidBlockSet() { return this->solidBlocks; }
    virtual void removeSolidBlocks() { this->solidBlocks.clear(); }

    void setID(int id);
    int getID();

    void setActive();
    void setInactive();
    bool isActive();

protected:
    void setTimeDependent() { UbSystem::setBit(this->type, TIMEDEPENDENT); }
    void unsetTimeDependent() { UbSystem::unsetBit(this->type, TIMEDEPENDENT); }

    //! detect that points are inside object
    //! \param min/max coordinates of bounding box
    //! \param delta is delta x
    bool arePointsInsideGeoObject(real minX1, real minX2, real minX3, real maxX1, real maxX2, real maxX3,
                                  real delta);

    //! detect that points aren't inside object
    //! \param min/max coordinates of bounding box
    //! \param delta is delta x
    bool arePointsOutsideGeoObject(real minX1, real minX2, real minX3, real maxX1, real maxX2, real maxX3,
                                   real delta);

    //! detect that points are cutting object
    //! \param min/max coordinates of bounding box
    //! \param delta is delta x
    bool arePointsCuttingGeoObject(real minX1, real minX2, real minX3, real maxX1, real maxX2, real maxX3,
                                   real delta);

    bool isBlockOutsideGeoObject(real minX1, real minX2, real minX3, real maxX1, real maxX2, real maxX3,
                                 real delta);
    bool isBlockInsideGeoObject(real minX1, real minX2, real minX3, real maxX1, real maxX2, real maxX3,
                                real delta);
    bool isBlockCuttingGeoObject(real minX1, real minX2, real minX3, real maxX1, real maxX2, real maxX3,
                                 real delta);

    void updateBlocks();

    SPtr<GbObject3D> geoObject3D;
    WPtr<Grid3D> grid;
    int type;

    std::vector<SPtr<Block3D>> bcBlocks;
    std::vector<SPtr<Block3D>> solidBlocks;
    int accuracy;

    bool active;
    int id;

public:
    static const int SOLID;              //= (1<<0); //1
    static const int INVERSESOLID;       //= (1<<1); //2
    static const int TIMEDEPENDENT;      //= (1<<2); //4   //zeitlich
    static const int FLUID;              //= (1<<3); //8
    static const int MOVEABLE;           //= (1<<4); //16  // geometrisch
    static const int CHANGENOTNECESSARY; //= (1<<5); //32

private:
};

#endif

//! \}
