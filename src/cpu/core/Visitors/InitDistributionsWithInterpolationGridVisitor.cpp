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
//! \file InitDistributionsWithInterpolationGridVisitor.cpp
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================
#include "InitDistributionsWithInterpolationGridVisitor.h"

#include "mpi.h"

#include "BCSet.h"
#include "Block3D.h"
#include "EsoSplit.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include "Interpolator.h"
#include "LBMKernel.h"
#include <CbArray2D.h>
#include <basics/utilities/UbFileInputASCII.h>

using namespace std;

InitDistributionsWithInterpolationGridVisitor::InitDistributionsWithInterpolationGridVisitor(
    SPtr<Grid3D> oldGrid, InterpolationProcessorPtr iProcessor, real nu)
    : oldGrid(oldGrid), iProcessor(iProcessor), nu(nu)
{
}
//////////////////////////////////////////////////////////////////////////
InitDistributionsWithInterpolationGridVisitor::~InitDistributionsWithInterpolationGridVisitor() = default;
//////////////////////////////////////////////////////////////////////////
void InitDistributionsWithInterpolationGridVisitor::visit(SPtr<Grid3D> grid)
{
    newGrid          = grid;
    int minInitLevel = newGrid->getCoarsestInitializedLevel();
    int maxInitLevel = newGrid->getFinestInitializedLevel();
    int newGridRank  = newGrid->getRank();

    for (int l = minInitLevel; l <= maxInitLevel; l++) {
        //      int n = 0;
        vector<SPtr<Block3D>> blockVector;
        newGrid->getBlocks(l, blockVector);
        vector<SPtr<Block3D>> tBlockID;

        for (SPtr<Block3D> newBlock : blockVector) {
            if (!newBlock)
                UB_THROW(UbException(UB_EXARGS, "block is not exist"));

            int newBlockRank = newBlock->getRank();
            //         int newBlockLevel = newBlock->getLevel();

            SPtr<Block3D> oldBlock =
                oldGrid->getBlock(newBlock->getX1(), newBlock->getX2(), newBlock->getX3(), newBlock->getLevel());
            if (oldBlock) {
                int oldBlockRank = oldBlock->getRank();
                if (oldBlockRank == newBlockRank && oldBlock->isActive() && newBlockRank == newGridRank &&
                    newBlock->isActive()) {
                    copyLocalBlock(oldBlock, newBlock);
                } else {
                    copyRemoteBlock(oldBlock, newBlock);
                }
            } else {
                int newlevel    = newBlock->getLevel();
                Vector3D coords = newGrid->getNodeCoordinates(newBlock, 1, 1, 1);

                UbTupleInt3 oldGridBlockIndexes =
                    oldGrid->getBlockIndexes(coords[0], coords[1], coords[2], newlevel - 1);
                oldBlock = oldGrid->getBlock(val<1>(oldGridBlockIndexes), val<2>(oldGridBlockIndexes),
                                             val<3>(oldGridBlockIndexes), newlevel - 1);

                if (oldBlock) {
                    int oldBlockRank = oldBlock->getRank();
                    //               int oldBlockLevel = oldBlock->getLevel();

                    if (oldBlockRank == newBlockRank && oldBlock->isActive() && newBlockRank == newGridRank &&
                        newBlock->isActive()) {
                        interpolateLocalBlockCoarseToFine(oldBlock, newBlock);
                    } else {
                        interpolateRemoteBlockCoarseToFine(oldBlock, newBlock);
                    }
                } else {
                    oldGridBlockIndexes = oldGrid->getBlockIndexes(coords[0], coords[1], coords[2], newlevel + 1);
                    oldBlock            = oldGrid->getBlock(val<1>(oldGridBlockIndexes), val<2>(oldGridBlockIndexes),
                                                 val<3>(oldGridBlockIndexes), newlevel + 1);
                    if (oldBlock) {
                        int oldBlockRank = oldBlock->getRank();
                        //                  int oldBlockLevel = oldBlock->getLevel();

                        if (oldBlockRank == newBlockRank && oldBlock->isActive() && newBlockRank == newGridRank &&
                            newBlock->isActive()) {
                            interpolateLocalBlockFineToCoarse(oldBlock, newBlock);
                        } else {
                            interpolateRemoteBlockFineToCoarse(oldBlock, newBlock);
                        }
                    }
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsWithInterpolationGridVisitor::copyLocalBlock(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock)
{
    SPtr<ILBMKernel> oldKernel = oldBlock->getKernel();
    if (!oldKernel)
        throw UbException(UB_EXARGS, "The LBM kernel isn't exist in block: " + oldBlock->toString());
    SPtr<EsoTwist3D> oldDistributions = dynamicPointerCast<EsoTwist3D>(oldKernel->getDataSet()->getFdistributions());

    SPtr<ILBMKernel> kernel = newBlock->getKernel();
    if (!kernel)
        throw UbException(UB_EXARGS, "The LBM kernel isn't exist in new block: " + newBlock->toString());
    kernel->getDataSet()->setFdistributions(oldDistributions);
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsWithInterpolationGridVisitor::copyRemoteBlock(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock)
{
    int newGridRank  = newGrid->getRank();
    int oldBlockRank = oldBlock->getRank();
    int newBlockRank = newBlock->getRank();

    if (oldBlockRank == newGridRank && oldBlock->isActive()) {
        SPtr<ILBMKernel> oldKernel = oldBlock->getKernel();
        if (!oldKernel)
            throw UbException(UB_EXARGS, "The LBM kernel isn't exist in block: " + oldBlock->toString());
        SPtr<EsoTwist3D> oldDistributions =
            dynamicPointerCast<EsoTwist3D>(oldKernel->getDataSet()->getFdistributions());

        CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getLocalDistributions();
        CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getNonLocalDistributions();
        CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getZeroDistributions();

        MPI_Send(localDistributions->getStartAdressOfSortedArray(0, 0, 0, 0),
                 (int)localDistributions->getDataVector().size(), MPI_DOUBLE, newBlockRank, 0, MPI_COMM_WORLD);
        MPI_Send(nonLocalDistributions->getStartAdressOfSortedArray(0, 0, 0, 0),
                 (int)nonLocalDistributions->getDataVector().size(), MPI_DOUBLE, newBlockRank, 0, MPI_COMM_WORLD);
        MPI_Send(zeroDistributions->getStartAdressOfSortedArray(0, 0, 0),
                 (int)zeroDistributions->getDataVector().size(), MPI_DOUBLE, newBlockRank, 0, MPI_COMM_WORLD);
    } else if (newBlockRank == newGridRank && newBlock->isActive()) {
        SPtr<ILBMKernel> newKernel = newBlock->getKernel();
        if (!newKernel)
            throw UbException(UB_EXARGS, "The LBM kernel isn't exist in new block: " + newBlock->toString() +
                                             UbSystem::toString(newGridRank));

        SPtr<EsoTwist3D> newDistributions =
            dynamicPointerCast<EsoTwist3D>(newKernel->getDataSet()->getFdistributions());

        CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributions =
            dynamicPointerCast<EsoSplit>(newDistributions)->getLocalDistributions();
        CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions =
            dynamicPointerCast<EsoSplit>(newDistributions)->getNonLocalDistributions();
        CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributions =
            dynamicPointerCast<EsoSplit>(newDistributions)->getZeroDistributions();

        MPI_Recv(localDistributions->getStartAdressOfSortedArray(0, 0, 0, 0),
                 (int)localDistributions->getDataVector().size(), MPI_DOUBLE, oldBlockRank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Recv(nonLocalDistributions->getStartAdressOfSortedArray(0, 0, 0, 0),
                 (int)nonLocalDistributions->getDataVector().size(), MPI_DOUBLE, oldBlockRank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Recv(zeroDistributions->getStartAdressOfSortedArray(0, 0, 0),
                 (int)zeroDistributions->getDataVector().size(), MPI_DOUBLE, oldBlockRank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsWithInterpolationGridVisitor::interpolateLocalBlockCoarseToFine(SPtr<Block3D> oldBlock,
                                                                                      SPtr<Block3D> newBlock)
{
    using namespace vf::basics::constant;
    
    D3Q27ICell icellC;
    D3Q27ICell icellF;
    real xoff, yoff, zoff;

    real omegaC = LBMSystem::calcCollisionFactor(nu, oldBlock->getLevel());
    real omegaF = LBMSystem::calcCollisionFactor(nu, newBlock->getLevel());

    iProcessor->setOmegas(omegaC, omegaF);

    SPtr<ILBMKernel> oldKernel = oldBlock->getKernel();
    if (!oldKernel)
        throw UbException(UB_EXARGS, "The LBM kernel isn't exist in old block: " + oldBlock->toString());

    SPtr<EsoTwist3D> oldDistributions = dynamicPointerCast<EsoTwist3D>(oldKernel->getDataSet()->getFdistributions());

    SPtr<BCArray3D> bcArrayOldBlock = oldBlock->getKernel()->getBCSet()->getBCArray();

    SPtr<ILBMKernel> newKernel = newBlock->getKernel();
    if (!newKernel)
        throw UbException(UB_EXARGS, "The LBM kernel isn't exist in new block: " + newBlock->toString());

    SPtr<EsoTwist3D> newDistributions = dynamicPointerCast<EsoTwist3D>(newKernel->getDataSet()->getFdistributions());

    int minX1 = 0;
    int minX2 = 0;
    int minX3 = 0;

    int maxX1 = (int)newDistributions->getNX1() - 1;
    int maxX2 = (int)newDistributions->getNX2() - 1;
    int maxX3 = (int)newDistributions->getNX3() - 1;

    int bMaxX1 = (int)newDistributions->getNX1();
    int bMaxX2 = (int)newDistributions->getNX2();
    int bMaxX3 = (int)newDistributions->getNX3();

    for (int ix3 = minX3; ix3 < maxX3; ix3 += 2)
        for (int ix2 = minX2; ix2 < maxX2; ix2 += 2)
            for (int ix1 = minX1; ix1 < maxX1; ix1 += 2) {
                Vector3D coords             = newGrid->getNodeCoordinates(newBlock, ix1, ix2, ix3);
                UbTupleInt3 oldGridIndexMin = oldGrid->getNodeIndexes(oldBlock, coords[0], coords[1], coords[2]);
                int howManySolids           = iProcessor->iCellHowManySolids(bcArrayOldBlock, val<1>(oldGridIndexMin),
                                                                   val<2>(oldGridIndexMin), val<3>(oldGridIndexMin));

                if (howManySolids == 0 || howManySolids == 8) {
                    iProcessor->readICell(oldDistributions, icellC, val<1>(oldGridIndexMin), val<2>(oldGridIndexMin),
                                          val<3>(oldGridIndexMin));
                    iProcessor->interpolateCoarseToFine(icellC, icellF);
                } else {
                    if (iProcessor->findNeighborICell(bcArrayOldBlock, oldDistributions, icellC, bMaxX1, bMaxX2, bMaxX3,
                                                      val<1>(oldGridIndexMin), val<2>(oldGridIndexMin),
                                                      val<3>(oldGridIndexMin), xoff, yoff, zoff)) {
                        // std::string err = "For "+oldBlock->toString()+
                        //   " x1="+UbSystem::toString(val<1>(oldGridIndexMin))+
                        //   ", x2=" + UbSystem::toString(val<2>(oldGridIndexMin))+
                        //   ", x3=" + UbSystem::toString(val<3>(oldGridIndexMin))+
                        //   " interpolation is not implemented for other direction"+
                        //   " by using in: "+(std::string)typeid(*this).name()+
                        //   " or maybe you have a solid on the block boundary";
                        // UB_THROW(UbException(UB_EXARGS, err));
                        iProcessor->interpolateCoarseToFine(icellC, icellF, xoff, yoff, zoff);
                    } else {
                        for (int i = 0; i < 27; i++) {
                            icellF.BSW[i] = c0o1;
                            icellF.BSE[i] = c0o1;
                            icellF.BNW[i] = c0o1;
                            icellF.BNE[i] = c0o1;
                            icellF.TSW[i] = c0o1;
                            icellF.TSE[i] = c0o1;
                            icellF.TNW[i] = c0o1;
                            icellF.TNE[i] = c0o1;
                        }
                        //                     std::string err = "For "+oldBlock->toString()+
                        //   " x1="+UbSystem::toString(val<1>(oldGridIndexMin))+
                        //   ", x2=" + UbSystem::toString(val<2>(oldGridIndexMin))+
                        //   ", x3=" + UbSystem::toString(val<3>(oldGridIndexMin))+
                        //   " interpolation is not implemented for other direction"+
                        //   " by using in: "+(std::string)typeid(*this).name()+
                        //   " or maybe you have a solid on the block boundary";
                        ////UB_THROW(UbException(UB_EXARGS, err));
                        //                     UBLOG(logINFO, err);
                    }
                }

                iProcessor->writeICell(newDistributions, icellF, ix1, ix2, ix3);
                iProcessor->writeICellInv(newDistributions, icellF, ix1, ix2, ix3);
            }
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsWithInterpolationGridVisitor::interpolateRemoteBlockCoarseToFine(SPtr<Block3D> oldBlock,
                                                                                       SPtr<Block3D> newBlock)
{
    using namespace vf::basics::constant;
    
    int newGridRank  = newGrid->getRank();
    int oldBlockRank = oldBlock->getRank();
    int newBlockRank = newBlock->getRank();

    if (oldBlockRank == newGridRank) {
        SPtr<ILBMKernel> oldKernel = oldBlock->getKernel();
        if (!oldKernel)
            throw UbException(UB_EXARGS, "The LBM kernel isn't exist in block: " + oldBlock->toString());
        SPtr<EsoTwist3D> oldDistributions =
            dynamicPointerCast<EsoTwist3D>(oldKernel->getDataSet()->getFdistributions());

        CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getLocalDistributions();
        CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getNonLocalDistributions();
        CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getZeroDistributions();

        MPI_Send(localDistributions->getStartAdressOfSortedArray(0, 0, 0, 0),
                 (int)localDistributions->getDataVector().size(), MPI_DOUBLE, newBlockRank, 0, MPI_COMM_WORLD);
        MPI_Send(nonLocalDistributions->getStartAdressOfSortedArray(0, 0, 0, 0),
                 (int)nonLocalDistributions->getDataVector().size(), MPI_DOUBLE, newBlockRank, 0, MPI_COMM_WORLD);
        MPI_Send(zeroDistributions->getStartAdressOfSortedArray(0, 0, 0),
                 (int)zeroDistributions->getDataVector().size(), MPI_DOUBLE, newBlockRank, 0, MPI_COMM_WORLD);

        SPtr<BCArray3D> bcArrayOldBlock = oldBlock->getKernel()->getBCSet()->getBCArray();
        std::vector<int> &bcDataVector  = bcArrayOldBlock->getBcindexmatrixDataVector();
        MPI_Send(&bcDataVector[0], (int)bcDataVector.size(), MPI_INT, newBlockRank, 0, MPI_COMM_WORLD);
    } else if (newBlockRank == newGridRank && newBlock->isActive()) {
        D3Q27ICell icellC;
        D3Q27ICell icellF;
        real xoff, yoff, zoff;

        real omegaC = LBMSystem::calcCollisionFactor(nu, oldBlock->getLevel());
        real omegaF = LBMSystem::calcCollisionFactor(nu, newBlock->getLevel());

        iProcessor->setOmegas(omegaC, omegaF);

        SPtr<ILBMKernel> newKernel = newBlock->getKernel();
        if (!newKernel)
            throw UbException(UB_EXARGS, "The LBM kernel isn't exist in new block: " + newBlock->toString());

        SPtr<EsoTwist3D> newDistributions =
            dynamicPointerCast<EsoTwist3D>(newKernel->getDataSet()->getFdistributions());

        int minX1 = 0;
        int minX2 = 0;
        int minX3 = 0;

        int maxX1 = (int)newDistributions->getNX1() - 1;
        int maxX2 = (int)newDistributions->getNX2() - 1;
        int maxX3 = (int)newDistributions->getNX3() - 1;

        int bMaxX1 = (int)newDistributions->getNX1();
        int bMaxX2 = (int)newDistributions->getNX2();
        int bMaxX3 = (int)newDistributions->getNX3();

        SPtr<EsoTwist3D> oldDistributions(new EsoSplit(bMaxX1, bMaxX2, bMaxX3, 0));

        CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getLocalDistributions();
        CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getNonLocalDistributions();
        CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getZeroDistributions();

        MPI_Recv(localDistributions->getStartAdressOfSortedArray(0, 0, 0, 0),
                 (int)localDistributions->getDataVector().size(), MPI_DOUBLE, oldBlockRank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Recv(nonLocalDistributions->getStartAdressOfSortedArray(0, 0, 0, 0),
                 (int)nonLocalDistributions->getDataVector().size(), MPI_DOUBLE, oldBlockRank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Recv(zeroDistributions->getStartAdressOfSortedArray(0, 0, 0),
                 (int)zeroDistributions->getDataVector().size(), MPI_DOUBLE, oldBlockRank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        SPtr<BCArray3D> bcArrayOldBlock(new BCArray3D(bMaxX1, bMaxX2, bMaxX3, BCArray3D::FLUID));
        std::vector<int> &bcDataVector = bcArrayOldBlock->getBcindexmatrixDataVector();
        MPI_Recv(&bcDataVector[0], (int)bcDataVector.size(), MPI_INT, oldBlockRank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        for (int ix3 = minX3; ix3 < maxX3; ix3 += 2)
            for (int ix2 = minX2; ix2 < maxX2; ix2 += 2)
                for (int ix1 = minX1; ix1 < maxX1; ix1 += 2) {
                    Vector3D coords             = newGrid->getNodeCoordinates(newBlock, ix1, ix2, ix3);
                    UbTupleInt3 oldGridIndexMin = oldGrid->getNodeIndexes(oldBlock, coords[0], coords[1], coords[2]);

                    int howManySolids = iProcessor->iCellHowManySolids(
                        bcArrayOldBlock, val<1>(oldGridIndexMin), val<2>(oldGridIndexMin), val<3>(oldGridIndexMin));

                    if (howManySolids == 0 || howManySolids == 8) {
                        iProcessor->readICell(oldDistributions, icellC, val<1>(oldGridIndexMin),
                                              val<2>(oldGridIndexMin), val<3>(oldGridIndexMin));
                        iProcessor->interpolateCoarseToFine(icellC, icellF);
                    } else {
                        if (iProcessor->findNeighborICell(bcArrayOldBlock, oldDistributions, icellC, bMaxX1, bMaxX2,
                                                          bMaxX3, val<1>(oldGridIndexMin), val<2>(oldGridIndexMin),
                                                          val<3>(oldGridIndexMin), xoff, yoff, zoff)) {
                            // std::string err = "For "+oldBlock->toString()+
                            //   " x1="+UbSystem::toString(val<1>(oldGridIndexMin))+
                            //   ", x2=" + UbSystem::toString(val<2>(oldGridIndexMin))+
                            //   ", x3=" + UbSystem::toString(val<3>(oldGridIndexMin))+
                            //   " interpolation is not implemented for other direction"+
                            //   " by using in: "+(std::string)typeid(*this).name()+
                            //   " or maybe you have a solid on the block boundary";
                            // UB_THROW(UbException(UB_EXARGS, err));
                            iProcessor->interpolateCoarseToFine(icellC, icellF, xoff, yoff, zoff);
                        } else {
                            for (int i = 0; i < 27; i++) {
                                icellF.BSW[i] = c0o1;
                                icellF.BSE[i] = c0o1;
                                icellF.BNW[i] = c0o1;
                                icellF.BNE[i] = c0o1;
                                icellF.TSW[i] = c0o1;
                                icellF.TSE[i] = c0o1;
                                icellF.TNW[i] = c0o1;
                                icellF.TNE[i] = c0o1;
                            }
                            //                     std::string err = "For "+oldBlock->toString()+
                            //   " x1="+UbSystem::toString(val<1>(oldGridIndexMin))+
                            //   ", x2=" + UbSystem::toString(val<2>(oldGridIndexMin))+
                            //   ", x3=" + UbSystem::toString(val<3>(oldGridIndexMin))+
                            //   " interpolation is not implemented for other direction"+
                            //   " by using in: "+(std::string)typeid(*this).name()+
                            //   " or maybe you have a solid on the block boundary";
                            ////UB_THROW(UbException(UB_EXARGS, err));
                            //                     UBLOG(logINFO, err);
                        }
                    }

                    iProcessor->writeICell(newDistributions, icellF, ix1, ix2, ix3);
                    iProcessor->writeICellInv(newDistributions, icellF, ix1, ix2, ix3);
                }
    }
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsWithInterpolationGridVisitor::interpolateLocalBlockFineToCoarse(SPtr<Block3D> oldBlock,
                                                                                      SPtr<Block3D> newBlock)
{
    using namespace vf::basics::constant;
    
    real icellC[27];
    D3Q27ICell icellF;
    real xoff, yoff, zoff;

    real omegaF = LBMSystem::calcCollisionFactor(nu, oldBlock->getLevel());
    real omegaC = LBMSystem::calcCollisionFactor(nu, newBlock->getLevel());

    iProcessor->setOmegas(omegaC, omegaF);

    SPtr<ILBMKernel> oldKernel = oldBlock->getKernel();
    if (!oldKernel)
        throw UbException(UB_EXARGS, "The LBM kernel isn't exist in old block: " + oldBlock->toString());

    SPtr<EsoTwist3D> oldDistributions = dynamicPointerCast<EsoTwist3D>(oldKernel->getDataSet()->getFdistributions());

    SPtr<BCArray3D> bcArrayOldBlock = oldBlock->getKernel()->getBCSet()->getBCArray();

    SPtr<ILBMKernel> newKernel = newBlock->getKernel();
    if (!newKernel)
        throw UbException(UB_EXARGS, "The LBM kernel isn't exist in new block: " + newBlock->toString());

    SPtr<EsoTwist3D> newDistributions = dynamicPointerCast<EsoTwist3D>(newKernel->getDataSet()->getFdistributions());

    int minX1 = 0;
    int minX2 = 0;
    int minX3 = 0;

    int maxX1 = (int)newDistributions->getNX1() - 1;
    int maxX2 = (int)newDistributions->getNX2() - 1;
    int maxX3 = (int)newDistributions->getNX3() - 1;

    int bMaxX1 = (int)newDistributions->getNX1();
    int bMaxX2 = (int)newDistributions->getNX2();
    int bMaxX3 = (int)newDistributions->getNX3();

    for (int ix3 = minX3; ix3 < maxX3; ix3 += 2)
        for (int ix2 = minX2; ix2 < maxX2; ix2 += 2)
            for (int ix1 = minX1; ix1 < maxX1; ix1 += 2) {
                Vector3D coords             = newGrid->getNodeCoordinates(newBlock, ix1, ix2, ix3);
                UbTupleInt3 oldGridIndexMin = oldGrid->getNodeIndexes(oldBlock, coords[0], coords[1], coords[2]);

                int howManySolids = iProcessor->iCellHowManySolids(bcArrayOldBlock, val<1>(oldGridIndexMin),
                                                                   val<2>(oldGridIndexMin), val<3>(oldGridIndexMin));

                if (howManySolids == 0 || howManySolids == 8) {
                    iProcessor->readICell(oldDistributions, icellF, val<1>(oldGridIndexMin), val<2>(oldGridIndexMin),
                                          val<3>(oldGridIndexMin));
                    iProcessor->interpolateFineToCoarse(icellF, icellC);
                } else {
                    if (iProcessor->findNeighborICell(bcArrayOldBlock, oldDistributions, icellF, bMaxX1, bMaxX2, bMaxX3,
                                                      val<1>(oldGridIndexMin), val<2>(oldGridIndexMin),
                                                      val<3>(oldGridIndexMin), xoff, yoff, zoff)) {
                        // std::string err = "For "+oldBlock->toString()+
                        //   " x1="+UbSystem::toString(val<1>(oldGridIndexMin))+
                        //   ", x2=" + UbSystem::toString(val<2>(oldGridIndexMin))+
                        //   ", x3=" + UbSystem::toString(val<3>(oldGridIndexMin))+
                        //   " interpolation is not implemented for other direction"+
                        //   " by using in: "+(std::string)typeid(*this).name()+
                        //   " or maybe you have a solid on the block boundary";
                        // UB_THROW(UbException(UB_EXARGS, err));
                        iProcessor->interpolateFineToCoarse(icellF, icellC, xoff, yoff, zoff);
                    } else {
                        for (int i = 0; i < 27; i++) {
                            icellF.BSW[i] = c0o1;
                            icellF.BSE[i] = c0o1;
                            icellF.BNW[i] = c0o1;
                            icellF.BNE[i] = c0o1;
                            icellF.TSW[i] = c0o1;
                            icellF.TSE[i] = c0o1;
                            icellF.TNW[i] = c0o1;
                            icellF.TNE[i] = c0o1;
                        }
                        //                     std::string err = "For "+oldBlock->toString()+
                        //   " x1="+UbSystem::toString(val<1>(oldGridIndexMin))+
                        //   ", x2=" + UbSystem::toString(val<2>(oldGridIndexMin))+
                        //   ", x3=" + UbSystem::toString(val<3>(oldGridIndexMin))+
                        //   " interpolation is not implemented for other direction"+
                        //   " by using in: "+(std::string)typeid(*this).name()+
                        //   " or maybe you have a solid on the block boundary";
                        ////UB_THROW(UbException(UB_EXARGS, err));
                        //                     UBLOG(logINFO, err);
                    }
                }

                iProcessor->writeINode(newDistributions, icellC, ix1, ix2, ix3);
                // iProcessor->writeINodeInv(newDistributions, icellC, ix1, ix2, ix3);
            }
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsWithInterpolationGridVisitor::interpolateRemoteBlockFineToCoarse(SPtr<Block3D> oldBlock,
                                                                                       SPtr<Block3D> newBlock)
{
    using namespace vf::basics::constant;
    
    int newGridRank  = newGrid->getRank();
    int oldBlockRank = oldBlock->getRank();
    int newBlockRank = newBlock->getRank();

    if (oldBlockRank == newGridRank) {
        SPtr<ILBMKernel> oldKernel = oldBlock->getKernel();
        if (!oldKernel)
            throw UbException(UB_EXARGS, "The LBM kernel isn't exist in block: " + oldBlock->toString());
        SPtr<EsoTwist3D> oldDistributions =
            dynamicPointerCast<EsoTwist3D>(oldKernel->getDataSet()->getFdistributions());

        CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getLocalDistributions();
        CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getNonLocalDistributions();
        CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getZeroDistributions();

        MPI_Send(localDistributions->getStartAdressOfSortedArray(0, 0, 0, 0),
                 (int)localDistributions->getDataVector().size(), MPI_DOUBLE, newBlockRank, 0, MPI_COMM_WORLD);
        MPI_Send(nonLocalDistributions->getStartAdressOfSortedArray(0, 0, 0, 0),
                 (int)nonLocalDistributions->getDataVector().size(), MPI_DOUBLE, newBlockRank, 0, MPI_COMM_WORLD);
        MPI_Send(zeroDistributions->getStartAdressOfSortedArray(0, 0, 0),
                 (int)zeroDistributions->getDataVector().size(), MPI_DOUBLE, newBlockRank, 0, MPI_COMM_WORLD);

        SPtr<BCArray3D> bcArrayOldBlock = oldBlock->getKernel()->getBCSet()->getBCArray();
        std::vector<int> &bcDataVector  = bcArrayOldBlock->getBcindexmatrixDataVector();
        MPI_Send(&bcDataVector[0], (int)bcDataVector.size(), MPI_INT, newBlockRank, 0, MPI_COMM_WORLD);
    } else if (newBlockRank == newGridRank && newBlock->isActive()) {
        real icellC[27];
        D3Q27ICell icellF;
        real xoff, yoff, zoff;

        real omegaF = LBMSystem::calcCollisionFactor(nu, oldBlock->getLevel());
        real omegaC = LBMSystem::calcCollisionFactor(nu, newBlock->getLevel());

        iProcessor->setOmegas(omegaC, omegaF);

        SPtr<ILBMKernel> newKernel = newBlock->getKernel();
        if (!newKernel)
            throw UbException(UB_EXARGS, "The LBM kernel isn't exist in new block: " + newBlock->toString());

        SPtr<EsoTwist3D> newDistributions =
            dynamicPointerCast<EsoTwist3D>(newKernel->getDataSet()->getFdistributions());

        int minX1 = 0;
        int minX2 = 0;
        int minX3 = 0;

        int maxX1 = (int)newDistributions->getNX1() - 1;
        int maxX2 = (int)newDistributions->getNX2() - 1;
        int maxX3 = (int)newDistributions->getNX3() - 1;

        int bMaxX1 = (int)newDistributions->getNX1();
        int bMaxX2 = (int)newDistributions->getNX2();
        int bMaxX3 = (int)newDistributions->getNX3();

        SPtr<EsoTwist3D> oldDistributions(new EsoSplit(bMaxX1, bMaxX2, bMaxX3, 0));

        CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getLocalDistributions();
        CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getNonLocalDistributions();
        CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributions =
            dynamicPointerCast<EsoSplit>(oldDistributions)->getZeroDistributions();

        MPI_Recv(localDistributions->getStartAdressOfSortedArray(0, 0, 0, 0),
                 (int)localDistributions->getDataVector().size(), MPI_DOUBLE, oldBlockRank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Recv(nonLocalDistributions->getStartAdressOfSortedArray(0, 0, 0, 0),
                 (int)nonLocalDistributions->getDataVector().size(), MPI_DOUBLE, oldBlockRank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Recv(zeroDistributions->getStartAdressOfSortedArray(0, 0, 0),
                 (int)zeroDistributions->getDataVector().size(), MPI_DOUBLE, oldBlockRank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        SPtr<BCArray3D> bcArrayOldBlock(new BCArray3D(bMaxX1, bMaxX2, bMaxX3, BCArray3D::FLUID));
        std::vector<int> &bcDataVector = bcArrayOldBlock->getBcindexmatrixDataVector();
        MPI_Recv(&bcDataVector[0], (int)bcDataVector.size(), MPI_INT, oldBlockRank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        for (int ix3 = minX3; ix3 < maxX3; ix3 += 2)
            for (int ix2 = minX2; ix2 < maxX2; ix2 += 2)
                for (int ix1 = minX1; ix1 < maxX1; ix1 += 2) {
                    Vector3D coords             = newGrid->getNodeCoordinates(newBlock, ix1, ix2, ix3);
                    UbTupleInt3 oldGridIndexMin = oldGrid->getNodeIndexes(oldBlock, coords[0], coords[1], coords[2]);

                    int howManySolids = iProcessor->iCellHowManySolids(
                        bcArrayOldBlock, val<1>(oldGridIndexMin), val<2>(oldGridIndexMin), val<3>(oldGridIndexMin));

                    if (howManySolids == 0 || howManySolids == 8) {
                        iProcessor->readICell(oldDistributions, icellF, val<1>(oldGridIndexMin),
                                              val<2>(oldGridIndexMin), val<3>(oldGridIndexMin));
                        iProcessor->interpolateFineToCoarse(icellF, icellC);
                    } else {
                        if (iProcessor->findNeighborICell(bcArrayOldBlock, oldDistributions, icellF, bMaxX1, bMaxX2,
                                                          bMaxX3, val<1>(oldGridIndexMin), val<2>(oldGridIndexMin),
                                                          val<3>(oldGridIndexMin), xoff, yoff, zoff)) {
                            // std::string err = "For "+oldBlock->toString()+
                            //   " x1="+UbSystem::toString(val<1>(oldGridIndexMin))+
                            //   ", x2=" + UbSystem::toString(val<2>(oldGridIndexMin))+
                            //   ", x3=" + UbSystem::toString(val<3>(oldGridIndexMin))+
                            //   " interpolation is not implemented for other direction"+
                            //   " by using in: "+(std::string)typeid(*this).name()+
                            //   " or maybe you have a solid on the block boundary";
                            // UB_THROW(UbException(UB_EXARGS, err));
                            iProcessor->interpolateFineToCoarse(icellF, icellC, xoff, yoff, zoff);
                        } else {
                            for (int i = 0; i < 27; i++) {
                                icellF.BSW[i] = c0o1;
                                icellF.BSE[i] = c0o1;
                                icellF.BNW[i] = c0o1;
                                icellF.BNE[i] = c0o1;
                                icellF.TSW[i] = c0o1;
                                icellF.TSE[i] = c0o1;
                                icellF.TNW[i] = c0o1;
                                icellF.TNE[i] = c0o1;
                            }
                            //                     std::string err = "For "+oldBlock->toString()+
                            //   " x1="+UbSystem::toString(val<1>(oldGridIndexMin))+
                            //   ", x2=" + UbSystem::toString(val<2>(oldGridIndexMin))+
                            //   ", x3=" + UbSystem::toString(val<3>(oldGridIndexMin))+
                            //   " interpolation is not implemented for other direction"+
                            //   " by using in: "+(std::string)typeid(*this).name()+
                            //   " or maybe you have a solid on the block boundary";
                            ////UB_THROW(UbException(UB_EXARGS, err));
                            //                     UBLOG(logINFO, err);
                        }
                    }

                    iProcessor->writeINode(newDistributions, icellC, ix1, ix2, ix3);
                    // iProcessor->writeINodeInv(newDistributions, icellC, ix1, ix2, ix3);
                }
    }
}

//////////////////////////////////////////////////////////////////////////
