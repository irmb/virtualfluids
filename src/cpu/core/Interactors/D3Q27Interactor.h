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
//! \author Soeren Freudiger
//! \author Sebastian Geller
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef D3Q27INTERACTOR_H
#define D3Q27INTERACTOR_H

#include <PointerDefinitions.h>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "D3Q27System.h"
#include "GbPoint3D.h"
#include "Interactor3D.h"
#include "UbException.h"
#include "UbTuple.h"

class BC;
class Block3D;
class Grid3D;
class GbObject3D;

using BcNodeIndicesMap    = std::map<SPtr<Block3D>, std::set<std::vector<int>>>;
using SolidNodeIndicesMap = std::map<SPtr<Block3D>, std::set<UbTupleInt3>>;

//! \brief A specialized class for grid generation.
//! \details Support standard geometric primitives.
class D3Q27Interactor : public Interactor3D
{
public:
    D3Q27Interactor();
    D3Q27Interactor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type);
    D3Q27Interactor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, SPtr<BC> BC, int type);
    D3Q27Interactor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, SPtr<BC> BC, int type,
                    Interactor3D::Accuracy a);

    ~D3Q27Interactor() override;

    void setRelevantForForces(const bool &value) { this->relevantForForces = value; }
    bool isRelevantForForces() { return this->relevantForForces; }

    virtual void addBC(const SPtr<BC> BC) { BCs.push_back(BC); }
    void deleteBC() { BCs.clear(); }

    void initInteractor(const real &timeStep = 0) override;
    void updateInteractor(const real &timestep = 0) override;

    void setReinitWithStoredQs(bool reinitWithStoredQsFlag) { this->reinitWithStoredQsFlag = reinitWithStoredQsFlag; }

    void removeSolidBlocks() override
    {
        Interactor3D::removeSolidBlocks();
        solidNodeIndicesMap.clear();
    }
    void removeBcBlocks() override
    {
        Interactor3D::removeBcBlocks();
        bcNodeIndicesMap.clear();
    }

    bool setDifferencesToGbObject3D(const SPtr<Block3D> block) override;

    ObObject *clone() { throw UbException(UB_EXARGS, "not implemented"); }

    void writeValidationAVSFile(std::string filename);
    virtual std::vector<std::pair<GbPoint3D, GbPoint3D>> getQsLineSet();

    void addQsLineSet(std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt2> &lines);

    const BcNodeIndicesMap &getBcNodeIndicesMap() const { return bcNodeIndicesMap; }

protected:
    bool relevantForForces;
    bool reinitWithStoredQsFlag;

    std::vector<SPtr<BC>> BCs;

    SolidNodeIndicesMap solidNodeIndicesMap;
    BcNodeIndicesMap bcNodeIndicesMap;

    void initRayVectors();
    real rayX1[D3Q27System::FENDDIR + 1];
    real rayX2[D3Q27System::FENDDIR + 1];
    real rayX3[D3Q27System::FENDDIR + 1];
};

#endif

//! \}
