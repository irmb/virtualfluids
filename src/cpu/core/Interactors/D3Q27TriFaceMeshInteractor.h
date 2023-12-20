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
//! \author Sören Freudiger
//! \author Sebastian Geller
//! \author Ehsan Kian Far
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef D3Q19AMRTRIFACEMESHINTERACTOR_H
#define D3Q19AMRTRIFACEMESHINTERACTOR_H

#include <PointerDefinitions.h>
#include <map>
#include <string>
#include <vector>

#include "CbArray3D.h"
#include "D3Q27Interactor.h"

class GbObject3D;
class Grid3D;
class BC;
class GbTriFaceMesh3D;
class Block3D;

class D3Q27TriFaceMeshInteractor : public D3Q27Interactor
{
public:
    D3Q27TriFaceMeshInteractor();
    D3Q27TriFaceMeshInteractor(SPtr<Grid3D> grid, std::string name = "D3Q27TriFaceMeshInteractor");
    D3Q27TriFaceMeshInteractor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type);
    D3Q27TriFaceMeshInteractor(SPtr<GbTriFaceMesh3D> triFaceMesh, SPtr<Grid3D> grid, SPtr<BC> BC, int type);
    D3Q27TriFaceMeshInteractor(SPtr<GbTriFaceMesh3D> triFaceMesh, SPtr<Grid3D> grid, SPtr<BC> BC, int type, Interactor3D::Accuracy a);

    ~D3Q27TriFaceMeshInteractor() override;

    void initInteractor(const real &timeStep = 0) override;
    void updateInteractor(const real &timestep = 0) override;
    void updateMovedGeometry(const real &timeStep = 0);
    void setQs(const real &timeStep);
    void refineBlockGridToLevel(int level, real startDistance, real stopDistance);

    bool setDifferencesToGbObject3D(const SPtr<Block3D> block) override;

    void setRegardPointInObjectTest(bool opt) { this->regardPIOTest = opt; }

    ObObject *clone() { throw UbException(UB_EXARGS, "not implemented"); }

    void clearBcNodeIndicesAndQsMap() { this->bcNodeIndicesAndQsMap.clear(); }

    virtual std::string toString();

protected:
    bool useHalfSpace{ true };
    bool regardPIOTest{ true };

    void reinitWithStoredQs(const real &timeStep);
    //   bool reinitWithStoredQsFlag;
    std::map<SPtr<Block3D>, std::map<UbTupleInt3, std::vector<float>>> bcNodeIndicesAndQsMap; 
    //!!! it may be that in this interactor
    // a BC was set to an rpos, but the same node in
    // changed to another type (e.g. Solid) in another
    // became --> there is no longer any BC in the place!
    enum SolidCheckMethod { ScanLine, PointInObject };

    enum FLAGS { BC_FLAG, UNDEF_FLAG, FLUID_FLAG, SOLID_FLAG, OLDSOLID_FLAG };
    void recursiveGridFill(CbArray3D<FLAGS> &flagfield, const short &xs, const short &ys, const short &zs,
                           const FLAGS &type);
    void iterativeGridFill(CbArray3D<FLAGS> &flagfield, const short &xs, const short &ys, const short &zs,
                           const FLAGS &type);
};

#endif

//! \}
