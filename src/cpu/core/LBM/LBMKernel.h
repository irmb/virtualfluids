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
//! \addtogroup cpu_LBM LBM
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef LBMKERNEL_H
#define LBMKERNEL_H

#include "LBMSystem.h"
#include <PointerDefinitions.h>
#include <array>
#include <limits>
#include <muParser.h>

class BCSet;
class DataSet3D;
class Block3D;

//! \brief A base class provides basic functionality for LBM kernel
class LBMKernel
{
public:
    LBMKernel();
    virtual ~LBMKernel() = default;

    virtual SPtr<LBMKernel> clone() = 0;

    virtual void calculate(int step) = 0;
    virtual real getCalculationTime() = 0;

    virtual void setBCSet(SPtr<BCSet> bcp);
    virtual SPtr<BCSet> getBCSet() const;

    virtual void setCollisionFactor(real collFactor);
    virtual real getCollisionFactor() const;

    virtual void setGhostLayerWidth(int witdh);
    virtual int getGhostLayerWidth() const;

    virtual void setDataSet(SPtr<DataSet3D> dataSet);
    virtual SPtr<DataSet3D> getDataSet() const;

    void setForcingX1(real forcingX1);
    void setForcingX2(real forcingX2);
    void setForcingX3(real forcingX3);

    void setForcingX1(const mu::Parser &parser);
    void setForcingX2(const mu::Parser &parser);
    void setForcingX3(const mu::Parser &parser);

    void setForcingX1(const std::string &muParserString);
    void setForcingX2(const std::string &muParserString);
    void setForcingX3(const std::string &muParserString);

    void setIndex(int x1, int x2, int x3);

    real getDeltaT() const;
    void setDeltaT(real dt);

    bool getCompressible() const;
    void setCompressible(bool val);

    bool getWithForcing() const;
    void setWithForcing(bool val);

    void setBlock(SPtr<Block3D> block);
    SPtr<Block3D> getBlock() const;

    bool isInsideOfDomain(const int &x1, const int &x2, const int &x3) const;

    void swapDistributions();

    void setNX(std::array<int, 3> nx);
    std::array<int, 3> getNX();

protected:
    SPtr<DataSet3D> dataSet;
    SPtr<BCSet> bcSet;
    real collFactor;
    int ghostLayerWidth{ 1 };
    bool compressible{ false };

    // forcing
    bool withForcing{ false };
    mu::Parser muForcingX1;
    mu::Parser muForcingX2;
    mu::Parser muForcingX3;
    int ix1, ix2, ix3;
    real deltaT{ 1.0 };

    WPtr<Block3D> block;

    std::array<int, 3> nx;

private:
    void checkFunction(mu::Parser fct);
};

#endif

//! \}
