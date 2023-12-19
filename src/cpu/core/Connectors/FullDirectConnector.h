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
//! \addtogroup cpu_Connectors Connectors
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef FullDirectConnector_H
#define FullDirectConnector_H

#include "Block3D.h"
#include "LocalBlock3DConnector.h"


//! \brief   Exchange data between blocks.
//! \details Connector send and receive full distributions between two blocks in shared memory.
class FullDirectConnector : public LocalBlock3DConnector
{
public:
    FullDirectConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir);
    virtual void init() override;
    void sendVectors() override;

protected:
    virtual void updatePointers() = 0;
    virtual void exchangeData();
    virtual void exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To) = 0;

    int maxX1;
    int maxX2;
    int maxX3;

private:
};

//////////////////////////////////////////////////////////////////////////

#endif

//! \}
