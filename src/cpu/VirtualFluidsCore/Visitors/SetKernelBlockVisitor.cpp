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

#include "BCProcessor.h"
#include "Block3D.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "Grid3DSystem.h"
#include "LBMKernel.h"
#include "LBMSystem.h"
#include <utility>

//////////////////////////////////////////////////////////////////////////
SetKernelBlockVisitor::SetKernelBlockVisitor(SPtr<LBMKernel> kernel, LBMReal nue, double availMem, double needMem,
                                             SetKernelBlockVisitor::Action action)
    : Block3DVisitor(0, Grid3DSystem::MAXLEVEL), kernel(std::move(kernel)), nue(nue), action(action), dataSetFlag(true)
{
    if (needMem > availMem) {
        throw UbException(UB_EXARGS, "SetKernelBlockVisitor: Not enough memory!!!");
    }
}

SetKernelBlockVisitor::SetKernelBlockVisitor(SPtr<LBMKernel> kernel, LBMReal nue, int numberOfProcesses,
                                             SetKernelBlockVisitor::Action action)
    : Block3DVisitor(0, Grid3DSystem::MAXLEVEL), kernel(std::move(kernel)), nue(nue), action(action), dataSetFlag(true),
      numberOfProcesses(numberOfProcesses)
{
}

//////////////////////////////////////////////////////////////////////////
void SetKernelBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    throwExceptionIfNotEnoughMemory(grid);

    if (kernel && (block->getRank() == grid->getRank())) {
        LBMReal collFactor = LBMSystem::calcCollisionFactor(nue, block->getLevel());
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

                SPtr<BCProcessor> bcProc = block->getKernel()->getBCProcessor();
                if (!bcProc) {
                    UB_THROW(UbException(
                        UB_EXARGS,
                        "It is not possible to change a BCProcessor in kernel! Old BCProcessor is not exist!"));
                }
                newKernel->setBCProcessor(bcProc);
                block->setKernel(newKernel);
            } break;

            case SetKernelBlockVisitor::ChangeKernelWithData: {
                SPtr<BCProcessor> bcProc = block->getKernel()->getBCProcessor();
                if (!bcProc) {
                    UB_THROW(UbException(
                        UB_EXARGS,
                        "It is not possible to change a BCProcessor in kernel! Old BCProcessor is not exist!"));
                }
                newKernel->setBCProcessor(bcProc);
                block->setKernel(newKernel);
            } break;
        }
    }
}

void SetKernelBlockVisitor::setNoDataSetFlag(bool flag) { dataSetFlag = flag; }

void SetKernelBlockVisitor::throwExceptionIfNotEnoughMemory(const SPtr<Grid3D> &grid)
{
    auto availableMemory = Utilities::getTotalPhysMem();
    auto requiredMemory  = getRequiredPhysicalMemory(grid);
    if (requiredMemory > availableMemory)
        throw UbException(UB_EXARGS, "SetKernelBlockVisitor: Not enough memory!!!");
}

double SetKernelBlockVisitor::getRequiredPhysicalMemory(const SPtr<Grid3D> &grid) const
{
    unsigned long long numberOfNodesPerBlockWithGhostLayer;
    auto numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
    auto blockNx        = grid->getBlockNX();
    int ghostLayer      = grid->getGhostLayerWidth() * 2 + 1;

    numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (val<1>(blockNx) + ghostLayer) *
                                          (val<2>(blockNx) + ghostLayer) * (val<3>(blockNx) + ghostLayer);

    auto needMemAll =
        double(numberOfNodesPerBlockWithGhostLayer * (27 * sizeof(double) + sizeof(int) + sizeof(float) * 4));

    return needMemAll / double(numberOfProcesses);
}
