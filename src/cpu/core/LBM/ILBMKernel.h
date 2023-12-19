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

#ifndef I_LBMKERNEL_H
#define I_LBMKERNEL_H

#include <PointerDefinitions.h>

#include "LBMSystem.h"

class BCSet;
class DataSet3D;

//! Abstract class provides interface for LBM kernel
class ILBMKernel
{
public:
    virtual ~ILBMKernel() = default;

    virtual void calculate(int step)    = 0;
    virtual real getCalculationTime() = 0;
    virtual void swapDistributions()    = 0;

    virtual bool getCompressible() const                                             = 0;
    virtual SPtr<BCSet> getBCSet() const                                 = 0;
    virtual void setBCSet(SPtr<BCSet> BCSet)                       = 0;
    virtual SPtr<DataSet3D> getDataSet() const                                       = 0;
    virtual real getCollisionFactor() const                                        = 0;
    virtual void setCollisionFactor(real collFactor)                               = 0;
    virtual bool isInsideOfDomain(const int &x1, const int &x2, const int &x3) const = 0;
    virtual int getGhostLayerWidth() const                                           = 0;
    virtual real getDeltaT() const                                                = 0;
    virtual bool getWithForcing() const                                              = 0;
};

#endif

//! \}
