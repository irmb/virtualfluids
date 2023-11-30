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
//! \file SetInterpolationConnectorsBlockVisitor.cpp
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================

#include "SetInterpolationConnectorsBlockVisitor.h"
#include "CoarseToFineVectorConnector.h"
#include "FineToCoarseVectorConnector.h"
#include "TwoDistributionsFullDirectConnector.h"
#include "TwoDistributionsFullVectorConnector.h"
#include "D3Q27System.h"

#include <parallel/Communicator.h>
#include <parallel/transmitter/TbTransmitterLocal.h>

#include "Interpolator.h"

SetInterpolationConnectorsBlockVisitor::SetInterpolationConnectorsBlockVisitor(std::shared_ptr<vf::parallel::Communicator> comm, real nue, SPtr<Interpolator> iProcessor) :
Block3DVisitor(0, D3Q27System::MAXLEVEL), 
    comm(comm),
    nue(nue),
    iProcessor(iProcessor)
{
}
//////////////////////////////////////////////////////////////////////////
SetInterpolationConnectorsBlockVisitor::~SetInterpolationConnectorsBlockVisitor(void)
{
}
//////////////////////////////////////////////////////////////////////////
void SetInterpolationConnectorsBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    if(!block) return;

    UbTupleInt3 blockNX = grid->getBlockNX();
    if (val<1>(blockNX) % 2 == 0 && val<2>(blockNX) % 2 == 0 && val<3>(blockNX) % 2 == 0)
        UB_THROW(UbException(UB_EXARGS, "Odd number of nodes: The number of nodes in each dimension (x,y,z) has to be even."));

    UBLOG(logDEBUG5, "SetInterpolationConnectorsBlockVisitor::visit() - start");
    UBLOG(logDEBUG5, block->toString());

    gridRank = comm->getProcessID();
    grid->setRank(gridRank);

    if(grid->getFinestInitializedLevel() > grid->getCoarsestInitializedLevel())
        setInterpolationConnectors(grid, block);

    UBLOG(logDEBUG5, "SetInterpolationConnectorsBlockVisitor::visit() - end");
}
//////////////////////////////////////////////////////////////////////////
void SetInterpolationConnectorsBlockVisitor::setInterpolationConnectors(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    using namespace vf::lbm::dir;

   UBLOG(logDEBUG5, "SetInterpolationConnectorsBlockVisitor::setInterpolationConnectors() - start");

    //search for all blocks with different ranks
    if (block->hasInterpolationFlagCF() && block->isActive())
    {
        int fbx1 = block->getX1() << 1;
        int fbx2 = block->getX2() << 1;
        int fbx3 = block->getX3() << 1;
        int level = block->getLevel() + 1;

        if( block->hasInterpolationFlagCF(dP00))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1,fbx2,fbx3,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);
            SPtr<Block3D> fblockNW = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);
            SPtr<Block3D> fblockNE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dP00);
        }
        if( block->hasInterpolationFlagCF(dM00))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2,fbx3,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1,fbx2+1,fbx3,level);
            SPtr<Block3D> fblockNW = grid->getBlock(fbx1,fbx2,fbx3+1,level);
            SPtr<Block3D> fblockNE = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dM00);
        }
        if( block->hasInterpolationFlagCF(d0P0))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2+1,fbx3,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);
            SPtr<Block3D> fblockNW = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);
            SPtr<Block3D> fblockNE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, d0P0);
        }
        if( block->hasInterpolationFlagCF(d0M0))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2,fbx3,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3,level);
            SPtr<Block3D> fblockNW = grid->getBlock(fbx1,fbx2,fbx3+1,level);
            SPtr<Block3D> fblockNE = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, d0M0);
        }
        if( block->hasInterpolationFlagCF(d00P))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2,fbx3+1,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);
            SPtr<Block3D> fblockNW = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);
            SPtr<Block3D> fblockNE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, d00P);
        }
        if( block->hasInterpolationFlagCF(d00M))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2,fbx3,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3,level);
            SPtr<Block3D> fblockNW = grid->getBlock(fbx1,fbx2+1,fbx3,level);
            SPtr<Block3D> fblockNE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, d00M);
        }

        //////NE-NW-SE-SW
        if( block->hasInterpolationFlagCF(dPP0)&&!block->hasInterpolationFlagCF(d0P0) && !block->hasInterpolationFlagCF(dP00))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1,fbx2+1,fbx3+0,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);
            SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
            SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dPP0);
        }
        if( block->hasInterpolationFlagCF(dMM0)&& !block->hasInterpolationFlagCF(dM00) && !block->hasInterpolationFlagCF(d0M0))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2,fbx3,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1,fbx2,fbx3+1,level);
            SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);
            SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dMM0);
        }
        if( block->hasInterpolationFlagCF(dPM0)&& !block->hasInterpolationFlagCF(dP00) && !block->hasInterpolationFlagCF(d0M0))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1,fbx2,fbx3+0,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);
            SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
            SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dPM0);
        }
        if( block->hasInterpolationFlagCF(dMP0)&& !block->hasInterpolationFlagCF(d0P0) && !block->hasInterpolationFlagCF(dM00))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2+1,fbx3,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);
            SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
            SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dMP0);
        }

        /////////TE-BW-BE-TW 1-0
        if( block->hasInterpolationFlagCF(dP0P)&& !block->hasInterpolationFlagCF(dP00) && !block->hasInterpolationFlagCF(d00P))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1,fbx2+0,fbx3+1,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);
            SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2+0, fbx3+1, level);
            SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dP0P);
        }
        if( block->hasInterpolationFlagCF(dM0M)&& !block->hasInterpolationFlagCF(dM00) && !block->hasInterpolationFlagCF(d00M))
        {

            SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2+0,fbx3,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1,fbx2+1,fbx3,level);
            SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2+0, fbx3, level);
            SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2+1, fbx3, level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dM0M);
        }
        if( block->hasInterpolationFlagCF(dP0M)&& !block->hasInterpolationFlagCF(dP00) && !block->hasInterpolationFlagCF(d00M))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1,fbx2+0,fbx3,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);
            SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2+0, fbx3, level);
            SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3, level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dP0M);
        }
        if( block->hasInterpolationFlagCF(dM0P)&& !block->hasInterpolationFlagCF(dM00) && !block->hasInterpolationFlagCF(d00P))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2+0,fbx3+1,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);
            SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2+0, fbx3+1, level);
            SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dM0P);
        }

        //////TN-BS-BN-TS
        if( block->hasInterpolationFlagCF(d0PP)&& !block->hasInterpolationFlagCF(d0P0) && !block->hasInterpolationFlagCF(d00P))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1+0,fbx2+1,fbx3+1,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);
            SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+0, fbx2+1, fbx3+1, level);
            SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, d0PP);
        }
        if( block->hasInterpolationFlagCF(d0MM)&& !block->hasInterpolationFlagCF(d0M0) && !block->hasInterpolationFlagCF(d00M))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1+0,fbx2,fbx3,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3,level);
            SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+0, fbx2, fbx3, level);
            SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2, fbx3, level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, d0MM);
        }
        if( block->hasInterpolationFlagCF(d0PM)&& !block->hasInterpolationFlagCF(d0P0) && !block->hasInterpolationFlagCF(d00M))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1+0,fbx2+1,fbx3,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);
            SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+0, fbx2+1, fbx3, level);
            SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3, level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, d0PM);
        }
        if( block->hasInterpolationFlagCF(d0MP)&& !block->hasInterpolationFlagCF(d0M0) && !block->hasInterpolationFlagCF(d00P))
        {
            SPtr<Block3D> fblockSW = grid->getBlock(fbx1+0,fbx2,fbx3+1,level);
            SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);
            SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+0, fbx2, fbx3+1, level);
            SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);

            setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, d0MP);
        }




      //////corners
      if (block->hasInterpolationFlagCF(dPPP)&&!block->hasInterpolationFlagCF(dP0P)&&!block->hasInterpolationFlagCF(d0PP)&&!block->hasInterpolationFlagCF(dPP0)&&!block->hasInterpolationFlagCF(d00P)&&!block->hasInterpolationFlagCF(d0P0) && !block->hasInterpolationFlagCF(dP00))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dPPP);
      }
      if (block->hasInterpolationFlagCF(dMMP)&&!block->hasInterpolationFlagCF(dM0P)&&!block->hasInterpolationFlagCF(d0MP)&& !block->hasInterpolationFlagCF(dMM0)&& !block->hasInterpolationFlagCF(d00P)&& !block->hasInterpolationFlagCF(dM00) && !block->hasInterpolationFlagCF(d0M0))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1, fbx2, fbx3, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dMMP);
      }
      if (block->hasInterpolationFlagCF(dPMP)&&!block->hasInterpolationFlagCF(dP0P)&&!block->hasInterpolationFlagCF(d0MP)&& !block->hasInterpolationFlagCF(dPM0)&& !block->hasInterpolationFlagCF(d00P)&& !block->hasInterpolationFlagCF(dP00) && !block->hasInterpolationFlagCF(d0M0))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1+1, fbx2, fbx3+0, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dPMP);
      }
      if (block->hasInterpolationFlagCF(dMPP)&&!block->hasInterpolationFlagCF(dM0P)&&!block->hasInterpolationFlagCF(d0PP)&& !block->hasInterpolationFlagCF(dMP0)&& !block->hasInterpolationFlagCF(d00P)&& !block->hasInterpolationFlagCF(d0P0) && !block->hasInterpolationFlagCF(dM00))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dMPP);
      }
      if (block->hasInterpolationFlagCF(dPPM)&&!block->hasInterpolationFlagCF(dP0M)&&!block->hasInterpolationFlagCF(d0PM)&& !block->hasInterpolationFlagCF(dPP0)&&!block->hasInterpolationFlagCF(d00M)&&!block->hasInterpolationFlagCF(d0P0) && !block->hasInterpolationFlagCF(dP00))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dPPM);
      }
      if (block->hasInterpolationFlagCF(dMMM)&& !block->hasInterpolationFlagCF(d0MM)&& !block->hasInterpolationFlagCF(dM0M)&& !block->hasInterpolationFlagCF(dMM0)&& !block->hasInterpolationFlagCF(d00M)&& !block->hasInterpolationFlagCF(dM00) && !block->hasInterpolationFlagCF(d0M0))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1, fbx2, fbx3+0, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1, fbx2, fbx3, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dMMM);
      }
      if (block->hasInterpolationFlagCF(dPMM)&& !block->hasInterpolationFlagCF(d0MM)&& !block->hasInterpolationFlagCF(dP0M)&& !block->hasInterpolationFlagCF(dPM0)&& !block->hasInterpolationFlagCF(d00M)&& !block->hasInterpolationFlagCF(dP00) && !block->hasInterpolationFlagCF(d0M0))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1, fbx2, fbx3, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1+1, fbx2, fbx3+0, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dPMM);
      }
      if (block->hasInterpolationFlagCF(dMPM)&& !block->hasInterpolationFlagCF(d0PM)&& !block->hasInterpolationFlagCF(dM0M)&& !block->hasInterpolationFlagCF(dMP0)&& !block->hasInterpolationFlagCF(d00M)&& !block->hasInterpolationFlagCF(d0P0) && !block->hasInterpolationFlagCF(dM00))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, dMPM);
      }

    }
   UBLOG(logDEBUG5, "SetInterpolationConnectorsBlockVisitor::setInterpolationConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void SetInterpolationConnectorsBlockVisitor::setInterpolationConnectors(SPtr<Block3D> fBlockSW, SPtr<Block3D> fBlockSE, SPtr<Block3D> fBlockNW, SPtr<Block3D> fBlockNE, SPtr<Block3D> cBlock, int dir)
{
   UBLOG(logDEBUG5, "SetInterpolationConnectorsBlockVisitor::setInterpolationConnectors(...) - start");
    int fBlockSWRank = -999, fBlockSERank = -999, fBlockNWRank = -999, fBlockNERank = -999;
    if(fBlockSW) fBlockSWRank = fBlockSW->getRank();
    if(fBlockNW) fBlockNWRank = fBlockNW->getRank();
    if(fBlockSE) fBlockSERank = fBlockSE->getRank();
    if(fBlockNE) fBlockNERank = fBlockNE->getRank();
    int cBlockRank   = cBlock->getRank();

    real omegaF {0.0};
    if(fBlockSW) omegaF =LBMSystem::calcCollisionFactor(nue, fBlockSW->getLevel());
    if(fBlockNW) omegaF =LBMSystem::calcCollisionFactor(nue, fBlockNW->getLevel());
    if(fBlockSE) omegaF =LBMSystem::calcCollisionFactor(nue, fBlockSE->getLevel());
    if(fBlockNE) omegaF =LBMSystem::calcCollisionFactor(nue, fBlockNE->getLevel());
    real omegaC = LBMSystem::calcCollisionFactor(nue, cBlock->getLevel());
    iProcessor->setOmegas(omegaC, omegaF);

    InterpolationProcessorPtr cIProcessor(iProcessor->clone());
    InterpolationProcessorPtr fIProcessorSW(iProcessor->clone());
    InterpolationProcessorPtr fIProcessorSE(iProcessor->clone());
    InterpolationProcessorPtr fIProcessorNW(iProcessor->clone());
    InterpolationProcessorPtr fIProcessorNE(iProcessor->clone());

    CreateTransmittersHelper::TransmitterPtr senderCFevenEvenSW, receiverCFevenEvenSW, 
        senderCFevenOddNW,  receiverCFevenOddNW, 
        senderCFoddEvenSE,  receiverCFoddEvenSE, 
        senderCFoddOddNE,   receiverCFoddOddNE,
        senderFCevenEvenSW, receiverFCevenEvenSW, 
        senderFCevenOddNW,  receiverFCevenOddNW, 
        senderFCoddEvenSE,  receiverFCoddEvenSE, 
        senderFCoddOddNE,   receiverFCoddOddNE;

    if(fBlockSW) createTransmitters(cBlock, fBlockSW, dir, CreateTransmittersHelper::SW, senderCFevenEvenSW, receiverCFevenEvenSW, senderFCevenEvenSW, receiverFCevenEvenSW);
    if(fBlockNW) createTransmitters(cBlock, fBlockNW, dir, CreateTransmittersHelper::NW, senderCFevenOddNW, receiverCFevenOddNW, senderFCevenOddNW, receiverFCevenOddNW);
    if(fBlockSE) createTransmitters(cBlock, fBlockSE, dir, CreateTransmittersHelper::SE, senderCFoddEvenSE, receiverCFoddEvenSE, senderFCoddEvenSE, receiverFCoddEvenSE);
    if(fBlockNE) createTransmitters(cBlock, fBlockNE, dir, CreateTransmittersHelper::NE, senderCFoddOddNE, receiverCFoddOddNE, senderFCoddOddNE, receiverFCoddOddNE);

    if(cBlockRank == gridRank)
    {
      SPtr<Block3DConnector> connector(new CoarseToFineVectorConnector< TbTransmitter< CbVector< real > > >(cBlock,
            senderCFevenEvenSW, receiverCFevenEvenSW, senderCFevenOddNW,  receiverCFevenOddNW, 
            senderCFoddEvenSE,  receiverCFoddEvenSE,  senderCFoddOddNE,   receiverCFoddOddNE, 
            dir, cIProcessor) );
        cBlock->setConnector(connector);
    }
    if(fBlockSW && fBlockSWRank == gridRank)
    {
        SPtr<Block3DConnector> connector( new FineToCoarseVectorConnector< TbTransmitter< CbVector< real > > >(fBlockSW, 
            senderFCevenEvenSW, receiverFCevenEvenSW, dir, fIProcessorSW, EvenEvenSW) );
        fBlockSW->setConnector(connector);
    }
    if(fBlockNW && fBlockNWRank == gridRank)
    {
        SPtr<Block3DConnector> connector( new FineToCoarseVectorConnector< TbTransmitter< CbVector< real > > >(fBlockNW, 
            senderFCevenOddNW, receiverFCevenOddNW, dir, fIProcessorNW, EvenOddNW) );
        fBlockNW->setConnector(connector);
    }
    if(fBlockSE && fBlockSERank == gridRank)
    {
        SPtr<Block3DConnector> connector( new FineToCoarseVectorConnector< TbTransmitter< CbVector< real > > >(fBlockSE, 
            senderFCoddEvenSE, receiverFCoddEvenSE, dir, fIProcessorSE, OddEvenSE) );
        fBlockSE->setConnector(connector);
    }
    if(fBlockNE && fBlockNERank == gridRank)
    {
        SPtr<Block3DConnector> connector( new FineToCoarseVectorConnector< TbTransmitter< CbVector< real > > >(fBlockNE, 
            senderFCoddOddNE, receiverFCoddOddNE, dir, fIProcessorNE, OddOddNE) );
        fBlockNE->setConnector(connector);
    }
   UBLOG(logDEBUG5, "SetInterpolationConnectorsBlockVisitor::setInterpolationConnectors(...) - end");
}
//////////////////////////////////////////////////////////////////////////
void SetInterpolationConnectorsBlockVisitor::createTransmitters(SPtr<Block3D> cBlock, SPtr<Block3D> fBlock, int dir, 
                                                        CreateTransmittersHelper::IBlock ib, 
                                                                      CreateTransmittersHelper::TransmitterPtr& senderCF, 
                                                                      CreateTransmittersHelper::TransmitterPtr& receiverCF, 
                                                                      CreateTransmittersHelper::TransmitterPtr& senderFC, 
                                                                      CreateTransmittersHelper::TransmitterPtr& receiverFC)
{
   UBLOG(logDEBUG5, "SetInterpolationConnectorsBlockVisitor::createTransmitters(...) - start");
    CreateTransmittersHelper helper;
    int fBlockRank = fBlock->getRank();
    int cBlockRank = cBlock->getRank();
    if(fBlockRank == cBlockRank && fBlockRank == gridRank)
    {
        senderCF = receiverFC = CreateTransmittersHelper::TransmitterPtr( new TbLocalTransmitter< CbVector< real > >());
        senderFC = receiverCF = CreateTransmittersHelper::TransmitterPtr( new TbLocalTransmitter< CbVector< real > >());
    }
    else if(cBlockRank == gridRank)
    {
        helper.createTransmitters(cBlock, fBlock, dir, ib, senderCF, receiverCF, comm, CreateTransmittersHelper::MPI);
    }
    else if(fBlockRank == gridRank)
    {
        helper.createTransmitters(fBlock, cBlock, dir, ib, senderFC, receiverFC, comm, CreateTransmittersHelper::MPI);
    }
   UBLOG(logDEBUG5, "SetInterpolationConnectorsBlockVisitor::createTransmitters(...) - end");
}

