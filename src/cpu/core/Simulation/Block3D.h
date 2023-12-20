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
//! \addtogroup cpu_Simulation Simulation
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef BLOCK3D_H
#define BLOCK3D_H

#include <PointerDefinitions.h>
#include <map>
#include <string>
#include <vector>
#include "lbm/constants/D3Q27.h"

class Block3DConnector;
class LBMKernel;
class ILBMKernel;

//! A class implements a block structure
class Block3D
{
public:
    Block3D() = default;
    Block3D(int x1, int x2, int x3, int level);
    virtual ~Block3D();
    bool operator==(const Block3D &src) const;
    bool operator!=(const Block3D &src) const;

    int getX1() const;
    int getX2() const;
    int getX3() const;

    void setActive(bool active);
    bool isActive() const;
    bool isNotActive() const;

    void setKernel(SPtr<LBMKernel> kernel);
    SPtr<ILBMKernel> getKernel() const;
    void deleteKernel();

    void setBundle(int bundle);
    int getBundle() const;

    void setRank(int rank);
    int getRank() const;

    void setLocalRank(int rank);
    int getLocalRank() const;

    int getGlobalID() const;
    void setGlobalID(int id);

    int getLocalID() const;
    void setLocalID(int id);

    int getPart() const;
    void setPart(int part);

    int getLevel() const;
    void setLevel(int level);

    // Connector-Section
    void setConnector(SPtr<Block3DConnector> connector);
    SPtr<Block3DConnector> getConnector(int dir) const;
    bool hasConnectors();
    void deleteConnectors();
    void pushBackSameLevelConnectors(std::vector<SPtr<Block3DConnector>> &localSameLevelConnectors,
                                     std::vector<SPtr<Block3DConnector>> &remoteSameLevelConnectors);
    void pushBackLocalSameLevelConnectors(std::vector<SPtr<Block3DConnector>> &localSameLevelConnectors);
    void pushBackRemoteSameLevelConnectors(std::vector<SPtr<Block3DConnector>> &remoteSameLevelConnectors);
    void pushBackLocalInterpolationConnectorsCF(std::vector<SPtr<Block3DConnector>> &localInterpolationConnectors);
    void pushBackRemoteInterpolationConnectorsCF(std::vector<SPtr<Block3DConnector>> &remoteInterpolationConnectors);
    void pushBackLocalInterpolationConnectorsFC(std::vector<SPtr<Block3DConnector>> &localInterpolationConnectors);
    void pushBackRemoteInterpolationConnectorsFC(std::vector<SPtr<Block3DConnector>> &remoteInterpolationConnectors);
    void pushBackLocalSameLevelConnectors(std::vector<SPtr<Block3DConnector>> &localSameLevelConnectors,
                                          const int &dir);
    void pushBackRemoteSameLevelConnectors(std::vector<SPtr<Block3DConnector>> &remoteSameLevelConnectors,
                                           const int &dir);
    int getNumberOfLocalConnectors();
    int getNumberOfRemoteConnectors();
    int getNumberOfLocalConnectorsForSurfaces();
    int getNumberOfRemoteConnectorsForSurfaces();

    void setWeight(int rank, int weight);
    int getWeight(int rank);
    void addWeightForAll(int weight);
    void addWeight(int rank, int weight);
    void clearWeight();
    int getWeightSize();

    // interpolation
    bool hasInterpolationFlag();
    bool hasInterpolationFlag(int dir);
    void deleteInterpolationFlag();

    int getCollectionOfInterpolationFlagCF();
    void setCollectionOfInterpolationFlagCF(int flags);

    void setInterpolationFlagCF(int dir);
    bool hasInterpolationFlagCF(int dir);
    bool hasInterpolationFlagCF();

    int getCollectionOfInterpolationFlagFC();
    void setCollectionOfInterpolationFlagFC(int flags);

    void setInterpolationFlagFC(int dir);
    bool hasInterpolationFlagFC(int dir);
    bool hasInterpolationFlagFC();

    real getWorkLoad();

    std::string toString();

    static int getMaxGlobalID() { return counter; }
    static void setMaxGlobalID(int /*c*/) { counter = 0; } // FIXME ???

private:
    int x1{ 0 };
    int x2{ 0 };
    int x3{ 0 };

    bool active{ true };

    int interpolationFlagCF{ 0 };
    int interpolationFlagFC{ 0 };

    SPtr<LBMKernel> kernel;
    std::vector<SPtr<Block3DConnector>> connectors;
    std::map<int, int> weight;

    int bundle{ -1 };
    int rank{ -1 };
    int lrank{ -1 };
    int globalID{ -1 };
    int localID{ -1 };
    int part{ -1 };
    int level{ -1 };
    static int counter;
};

#endif // BLOCK3D_H

//! \}
