#include "MemoryUtil.h"
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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file SetKernelBlockVisitor.cpp
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================

#include "SetKernelBlockVisitor.h"

#include "BCSet.h"
#include "Block3D.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include "LBMKernel.h"
#include "LBMSystem.h"
#include <utility>

//////////////////////////////////////////////////////////////////////////
SetKernelBlockVisitor::SetKernelBlockVisitor(SPtr<LBMKernel> kernel, real nue, SetKernelBlockVisitor::Action action)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), kernel(std::move(kernel)), nue(nue), action(action)
{
}

SetKernelBlockVisitor::SetKernelBlockVisitor(SPtr<LBMKernel> kernel, real nue, int numberOfProcesses,
                                             SetKernelBlockVisitor::Action action)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), kernel(std::move(kernel)), nue(nue), action(action), numberOfProcesses(numberOfProcesses)
{
}

//////////////////////////////////////////////////////////////////////////
void SetKernelBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    throwExceptionIfNotEnoughMemory(grid);

    if (kernel && (block->getRank() == grid->getRank())) {
        real collFactor = LBMSystem::calcCollisionFactor(nue, block->getLevel());
        kernel->setCollisionFactor(collFactor);
        kernel->setIndex(block->getX1(), block->getX2(), block->getX3());
        kernel->setDeltaT(LBMSystem::getDeltaT(block->getLevel()));
        kernel->setBlock(block);
        UbTupleInt3 blockNX = grid->getBlockNX();
        kernel->setNX(std::array<int, 3>{ { val<1>(blockNX), val<2>(blockNX), val<3>(blockNX) } });
        SPtr<LBMKernel> newKernel = kernel->clone();

        switch (action) {
            case SetKernelBlockVisitor::NewKernel:
                block->setKernel(newKernel);
                break;
            case SetKernelBlockVisitor::ChangeKernel: {
                SPtr<DataSet3D> dataSet = block->getKernel()->getDataSet();
                if (!dataSet) {
                    UB_THROW(UbException(
                        UB_EXARGS, "It is not possible to change a DataSet in kernel! Old DataSet is not exist!"));
                }

                newKernel->setDataSet(dataSet);

                SPtr<BCSet> bcProc = block->getKernel()->getBCSet();
                if (!bcProc) {
                    UB_THROW(UbException(
                        UB_EXARGS,
                        "It is not possible to change a BCSet in kernel! Old BCSet is not exist!"));
                }
                newKernel->setBCSet(bcProc);
                block->setKernel(newKernel);
            } break;

            case SetKernelBlockVisitor::ChangeKernelWithData: {
                SPtr<BCSet> bcProc = block->getKernel()->getBCSet();
                if (!bcProc) {
                    UB_THROW(UbException(
                        UB_EXARGS,
                        "It is not possible to change a BCSet in kernel! Old BCSet is not exist!"));
                }
                newKernel->setBCSet(bcProc);
                block->setKernel(newKernel);
            } break;
        }
    }
}

void SetKernelBlockVisitor::throwExceptionIfNotEnoughMemory(const SPtr<Grid3D> &grid)
{
    auto availableMemory = Utilities::getTotalPhysMem();
    auto requiredMemory  = getRequiredPhysicalMemory(grid);
    if (requiredMemory > availableMemory)
        throw UbException(UB_EXARGS, "SetKernelBlockVisitor: Not enough memory!!!");
}

real SetKernelBlockVisitor::getRequiredPhysicalMemory(const SPtr<Grid3D> &grid) const
{
    using namespace vf::basics::constant;

    unsigned long long numberOfNodesPerBlockWithGhostLayer;
    auto numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
    auto blockNx        = grid->getBlockNX();
    int ghostLayer      = grid->getGhostLayerWidth() * 2 + 1;

    numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (val<1>(blockNx) + ghostLayer) *
                                          (val<2>(blockNx) + ghostLayer) * (val<3>(blockNx) + ghostLayer);

    auto needMemAll =
        real(numberOfNodesPerBlockWithGhostLayer * (c27o1 * sizeof(real) + sizeof(int) + sizeof(float) * c4o1));

    return needMemAll / real(numberOfProcesses);
}
