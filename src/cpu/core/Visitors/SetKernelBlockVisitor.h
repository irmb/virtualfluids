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
//! \addtogroup cpu_Visitors Visitors
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef SetKernelBlockVisitor_h
#define SetKernelBlockVisitor_h

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"
#include "LBMSystem.h"

class Grid3D;
class Block3D;
class LBMKernel;

//! \brief A class generates new LBMKernel and associates it with a block.
class SetKernelBlockVisitor : public Block3DVisitor
{
public:
    enum Action { NewKernel, ChangeKernel, ChangeKernelWithData };

    SetKernelBlockVisitor(SPtr<LBMKernel> kernel, real nue, SetKernelBlockVisitor::Action action = SetKernelBlockVisitor::NewKernel);

    SetKernelBlockVisitor(SPtr<LBMKernel> kernel, real nue, int numberOfProcesses,
                          SetKernelBlockVisitor::Action action = SetKernelBlockVisitor::NewKernel);

    ~SetKernelBlockVisitor() override = default;

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
    SPtr<LBMKernel> kernel;
    real nue;
    Action action;

    int numberOfProcesses{ 1 };

    real getRequiredPhysicalMemory(const SPtr<Grid3D> &grid) const;

    void throwExceptionIfNotEnoughMemory(const SPtr<Grid3D> &grid);
};

#endif

//! \}
