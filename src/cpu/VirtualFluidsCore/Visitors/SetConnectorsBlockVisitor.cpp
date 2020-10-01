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
//! \file SetConnectorsBlockVisitor.cpp
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================

#include "SetConnectorsBlockVisitor.h"
#include "D3Q27ETFullDirectConnector.h"
#include "D3Q27ETFullVectorConnector.h"
#include "D3Q27ETCFOffVectorConnector.h"
#include "D3Q27ETFCOffVectorConnector.h"
#include "Grid3DSystem.h"
#include <basics/transmitter/TbTransmitterLocal.h>

#include "Communicator.h"
#include "InterpolationProcessor.h"

SetConnectorsBlockVisitor::SetConnectorsBlockVisitor(SPtr<Communicator> comm, bool fullConnector, int dirs, 
															   LBMReal nue, InterpolationProcessorPtr iProcessor) :
Block3DVisitor(0, Grid3DSystem::MAXLEVEL), 
	comm(comm),
	fullConnector(fullConnector),
	dirs(dirs),
	nue(nue),
	iProcessor(iProcessor)
{
}
//////////////////////////////////////////////////////////////////////////
SetConnectorsBlockVisitor::~SetConnectorsBlockVisitor(void)
{
}
//////////////////////////////////////////////////////////////////////////
void SetConnectorsBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
	if(!block) return;

	UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::visit() - start");
   UBLOG(logDEBUG5, block->toString());

	gridRank = comm->getProcessID();
	grid->setRank(gridRank);

	setSameLevelConnectors(grid, block);

	if(grid->getFinestInitializedLevel() > grid->getCoarsestInitializedLevel())
		setInterpolationConnectors(grid, block);

	UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::visit() - end");
}
//////////////////////////////////////////////////////////////////////////
void SetConnectorsBlockVisitor::setSameLevelConnectors(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::setSameLevelConnectors() - start");
	int blockRank = block->getRank();
	if (gridRank == blockRank && block->isActive())
	{
		block->clearWeight();
		std::vector<SPtr<Block3D>> neighbors; 
		int ix1 = block->getX1();
		int ix2 = block->getX2();
		int ix3 = block->getX3();
		int level = block->getLevel();
		//grid->getAllNeighbors(ix1, ix2, ix3, level, level, neighbors);

      //if (block->getGlobalID()==2512)
      //{
      //   int test = 0;
      //}

		for( int dir = 0; dir < dirs; dir++)
		{
			SPtr<Block3D> neighBlock = grid->getNeighborBlock(dir, ix1, ix2, ix3, level);

			if(neighBlock)
			{
				int neighBlockRank = neighBlock->getRank();
				if(blockRank == neighBlockRank && neighBlock->isActive())
				{
					SPtr<Block3DConnector> connector;
               connector = SPtr<Block3DConnector>(new D3Q27ETFullDirectConnector( block, neighBlock, dir));
					block->setConnector(connector);
				}
				else if(blockRank != neighBlockRank && neighBlock->isActive())
				{
					setRemoteConnectors(block, neighBlock, dir, fullConnector);  

					if(dir >=0 && dir<=5)
					{
						int weight = block->getWeight(neighBlockRank);
						weight++;
						block->setWeight(neighBlockRank, weight);
					}
				}
			}
		}
      
      //if (block->getGlobalID()==2794)
      //{
      //   UBLOG(logINFO, block->toString());
      //}
		
      int weight = block->getNumberOfLocalConnectorsForSurfaces();
		weight = 6 - weight;
		block->addWeightForAll(weight);
	}
   UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::setSameLevelConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void SetConnectorsBlockVisitor::setRemoteConnectors(SPtr<Block3D> sblock, SPtr<Block3D> tblock, int dir, bool fullConnector)
{
   UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::setRemoteConnectors() - start");
	CreateTransmittersHelper helper;
	CreateTransmittersHelper::TransmitterPtr sender, receiver;
	helper.createTransmitters(sblock, tblock, dir, CreateTransmittersHelper::NONE, sender, receiver, comm, CreateTransmittersHelper::MPI);


	SPtr<Block3DConnector> connector;
	connector = SPtr<Block3DConnector>(new D3Q27ETFullVectorConnector(sblock, sender, receiver, dir));
	sblock->setConnector(connector);
   UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::setRemoteConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void SetConnectorsBlockVisitor::setInterpolationConnectors(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::setInterpolationConnectors() - start");
	int blockRank = block->getRank();
	if (block->getGlobalID()==394)
	{
		int test=0;
	}

	//search for all blocks with different ranks
	if (block->hasInterpolationFlagCF() && block->isActive())
	{
		int fbx1 = block->getX1() << 1;
		int fbx2 = block->getX2() << 1;
		int fbx3 = block->getX3() << 1;
		int level = block->getLevel() + 1;

		if( block->hasInterpolationFlagCF(D3Q27System::E))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1,fbx2,fbx3,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);
			SPtr<Block3D> fblockNW = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);
			SPtr<Block3D> fblockNE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::E);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::W))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2,fbx3,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1,fbx2+1,fbx3,level);
			SPtr<Block3D> fblockNW = grid->getBlock(fbx1,fbx2,fbx3+1,level);
			SPtr<Block3D> fblockNE = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::W);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::N))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2+1,fbx3,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);
			SPtr<Block3D> fblockNW = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);
			SPtr<Block3D> fblockNE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::N);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::S))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2,fbx3,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3,level);
			SPtr<Block3D> fblockNW = grid->getBlock(fbx1,fbx2,fbx3+1,level);
			SPtr<Block3D> fblockNE = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::S);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::T))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2,fbx3+1,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);
			SPtr<Block3D> fblockNW = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);
			SPtr<Block3D> fblockNE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::T);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::B))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2,fbx3,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3,level);
			SPtr<Block3D> fblockNW = grid->getBlock(fbx1,fbx2+1,fbx3,level);
			SPtr<Block3D> fblockNE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::B);
		}

		//////NE-NW-SE-SW
		if( block->hasInterpolationFlagCF(D3Q27System::NE)&&!block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::E))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1,fbx2+1,fbx3+0,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::NE);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::SW)&& !block->hasInterpolationFlagCF(D3Q27System::W) && !block->hasInterpolationFlagCF(D3Q27System::S))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2,fbx3,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1,fbx2,fbx3+1,level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::SW);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::SE)&& !block->hasInterpolationFlagCF(D3Q27System::E) && !block->hasInterpolationFlagCF(D3Q27System::S))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1,fbx2,fbx3+0,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::SE);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::NW)&& !block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::W))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2+1,fbx3,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::NW);
		}

		/////////TE-BW-BE-TW 1-0
		if( block->hasInterpolationFlagCF(D3Q27System::TE)&& !block->hasInterpolationFlagCF(D3Q27System::E) && !block->hasInterpolationFlagCF(D3Q27System::T))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1,fbx2+0,fbx3+1,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2+0, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::TE);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::BW)&& !block->hasInterpolationFlagCF(D3Q27System::W) && !block->hasInterpolationFlagCF(D3Q27System::B))
		{

			SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2+0,fbx3,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1,fbx2+1,fbx3,level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2+0, fbx3, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2+1, fbx3, level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::BW);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::BE)&& !block->hasInterpolationFlagCF(D3Q27System::E) && !block->hasInterpolationFlagCF(D3Q27System::B))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1,fbx2+0,fbx3,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2+0, fbx3, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3, level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::BE);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::TW)&& !block->hasInterpolationFlagCF(D3Q27System::W) && !block->hasInterpolationFlagCF(D3Q27System::T))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1,fbx2+0,fbx3+1,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2+0, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::TW);
		}

		//////TN-BS-BN-TS
		if( block->hasInterpolationFlagCF(D3Q27System::TN)&& !block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::T))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1+0,fbx2+1,fbx3+1,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+0, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::TN);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::BS)&& !block->hasInterpolationFlagCF(D3Q27System::S) && !block->hasInterpolationFlagCF(D3Q27System::B))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1+0,fbx2,fbx3,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3,level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+0, fbx2, fbx3, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2, fbx3, level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::BS);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::BN)&& !block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::B))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1+0,fbx2+1,fbx3,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+0, fbx2+1, fbx3, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3, level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::BN);
		}
		if( block->hasInterpolationFlagCF(D3Q27System::TS)&& !block->hasInterpolationFlagCF(D3Q27System::S) && !block->hasInterpolationFlagCF(D3Q27System::T))
		{
			SPtr<Block3D> fblockSW = grid->getBlock(fbx1+0,fbx2,fbx3+1,level);
			SPtr<Block3D> fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+0, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);

			setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::TS);
		}




      //////corners
      if (block->hasInterpolationFlagCF(D3Q27System::TNE)&&!block->hasInterpolationFlagCF(D3Q27System::TE)&&!block->hasInterpolationFlagCF(D3Q27System::TN)&&!block->hasInterpolationFlagCF(D3Q27System::NE)&&!block->hasInterpolationFlagCF(D3Q27System::T)&&!block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::E))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::TNE);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::TSW)&&!block->hasInterpolationFlagCF(D3Q27System::TW)&&!block->hasInterpolationFlagCF(D3Q27System::TS)&& !block->hasInterpolationFlagCF(D3Q27System::SW)&& !block->hasInterpolationFlagCF(D3Q27System::T)&& !block->hasInterpolationFlagCF(D3Q27System::W) && !block->hasInterpolationFlagCF(D3Q27System::S))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1, fbx2, fbx3, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::TSW);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::TSE)&&!block->hasInterpolationFlagCF(D3Q27System::TE)&&!block->hasInterpolationFlagCF(D3Q27System::TS)&& !block->hasInterpolationFlagCF(D3Q27System::SE)&& !block->hasInterpolationFlagCF(D3Q27System::T)&& !block->hasInterpolationFlagCF(D3Q27System::E) && !block->hasInterpolationFlagCF(D3Q27System::S))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1+1, fbx2, fbx3+0, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::TSE);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::TNW)&&!block->hasInterpolationFlagCF(D3Q27System::TW)&&!block->hasInterpolationFlagCF(D3Q27System::TN)&& !block->hasInterpolationFlagCF(D3Q27System::NW)&& !block->hasInterpolationFlagCF(D3Q27System::T)&& !block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::W))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::TNW);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::BNE)&&!block->hasInterpolationFlagCF(D3Q27System::BE)&&!block->hasInterpolationFlagCF(D3Q27System::BN)&& !block->hasInterpolationFlagCF(D3Q27System::NE)&&!block->hasInterpolationFlagCF(D3Q27System::B)&&!block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::E))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::BNE);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::BSW)&& !block->hasInterpolationFlagCF(D3Q27System::BS)&& !block->hasInterpolationFlagCF(D3Q27System::BW)&& !block->hasInterpolationFlagCF(D3Q27System::SW)&& !block->hasInterpolationFlagCF(D3Q27System::B)&& !block->hasInterpolationFlagCF(D3Q27System::W) && !block->hasInterpolationFlagCF(D3Q27System::S))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1, fbx2, fbx3+0, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1, fbx2, fbx3, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::BSW);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::BSE)&& !block->hasInterpolationFlagCF(D3Q27System::BS)&& !block->hasInterpolationFlagCF(D3Q27System::BE)&& !block->hasInterpolationFlagCF(D3Q27System::SE)&& !block->hasInterpolationFlagCF(D3Q27System::B)&& !block->hasInterpolationFlagCF(D3Q27System::E) && !block->hasInterpolationFlagCF(D3Q27System::S))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1+1, fbx2, fbx3, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1+1, fbx2, fbx3+0, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::BSE);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::BNW)&& !block->hasInterpolationFlagCF(D3Q27System::BN)&& !block->hasInterpolationFlagCF(D3Q27System::BW)&& !block->hasInterpolationFlagCF(D3Q27System::NW)&& !block->hasInterpolationFlagCF(D3Q27System::B)&& !block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::W))
      {
         SPtr<Block3D> fblockSW = grid->getBlock(fbx1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblockSE;// = grid->getBlock(fbx1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblockNW;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblockNE;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::BNW);
      }

	}
   UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::setInterpolationConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void SetConnectorsBlockVisitor::setInterpolationConnectors(SPtr<Block3D> fBlockSW, SPtr<Block3D> fBlockSE, SPtr<Block3D> fBlockNW, SPtr<Block3D> fBlockNE, SPtr<Block3D> cBlock, int dir)
{
   UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::setInterpolationConnectors(...) - start");
	int fBlockSWRank = -999, fBlockSERank = -999, fBlockNWRank = -999, fBlockNERank = -999;
	if(fBlockSW) fBlockSWRank = fBlockSW->getRank();
	if(fBlockNW) fBlockNWRank = fBlockNW->getRank();
	if(fBlockSE) fBlockSERank = fBlockSE->getRank();
	if(fBlockNE) fBlockNERank = fBlockNE->getRank();
	int cBlockRank   = cBlock->getRank();

	LBMReal omegaF;
	if(fBlockSW) omegaF =LBMSystem::calcCollisionFactor(nue, fBlockSW->getLevel());
	if(fBlockNW) omegaF =LBMSystem::calcCollisionFactor(nue, fBlockNW->getLevel());
	if(fBlockSE) omegaF =LBMSystem::calcCollisionFactor(nue, fBlockSE->getLevel());
	if(fBlockNE) omegaF =LBMSystem::calcCollisionFactor(nue, fBlockNE->getLevel());
	LBMReal omegaC = LBMSystem::calcCollisionFactor(nue, cBlock->getLevel());
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
      SPtr<Block3DConnector> connector(new D3Q27ETCFOffVectorConnector< TbTransmitter< CbVector< LBMReal > > >(cBlock,
			senderCFevenEvenSW, receiverCFevenEvenSW, senderCFevenOddNW,  receiverCFevenOddNW, 
			senderCFoddEvenSE,  receiverCFoddEvenSE,  senderCFoddOddNE,   receiverCFoddOddNE, 
			dir, cIProcessor) );
		cBlock->setConnector(connector);
	}
	if(fBlockSW && fBlockSWRank == gridRank)
	{
		SPtr<Block3DConnector> connector( new D3Q27ETFCOffVectorConnector< TbTransmitter< CbVector< LBMReal > > >(fBlockSW, 
			senderFCevenEvenSW, receiverFCevenEvenSW, dir, fIProcessorSW, EvenEvenSW) );
		fBlockSW->setConnector(connector);
	}
	if(fBlockNW && fBlockNWRank == gridRank)
	{
		SPtr<Block3DConnector> connector( new D3Q27ETFCOffVectorConnector< TbTransmitter< CbVector< LBMReal > > >(fBlockNW, 
			senderFCevenOddNW, receiverFCevenOddNW, dir, fIProcessorNW, EvenOddNW) );
		fBlockNW->setConnector(connector);
	}
	if(fBlockSE && fBlockSERank == gridRank)
	{
		SPtr<Block3DConnector> connector( new D3Q27ETFCOffVectorConnector< TbTransmitter< CbVector< LBMReal > > >(fBlockSE, 
			senderFCoddEvenSE, receiverFCoddEvenSE, dir, fIProcessorSE, OddEvenSE) );
		fBlockSE->setConnector(connector);
	}
	if(fBlockNE && fBlockNERank == gridRank)
	{
		SPtr<Block3DConnector> connector( new D3Q27ETFCOffVectorConnector< TbTransmitter< CbVector< LBMReal > > >(fBlockNE, 
			senderFCoddOddNE, receiverFCoddOddNE, dir, fIProcessorNE, OddOddNE) );
		fBlockNE->setConnector(connector);
	}
   UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::setInterpolationConnectors(...) - end");
}
//////////////////////////////////////////////////////////////////////////
void SetConnectorsBlockVisitor::createTransmitters(SPtr<Block3D> cBlock, SPtr<Block3D> fBlock, int dir, 
                                                        CreateTransmittersHelper::IBlock ib, 
														              CreateTransmittersHelper::TransmitterPtr& senderCF, 
														              CreateTransmittersHelper::TransmitterPtr& receiverCF, 
														              CreateTransmittersHelper::TransmitterPtr& senderFC, 
														              CreateTransmittersHelper::TransmitterPtr& receiverFC)
{
   UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::createTransmitters(...) - start");
	CreateTransmittersHelper helper;
	bool MPIpool = true;
	bool orthogonal = false;
	int fBlockRank = fBlock->getRank();
	int cBlockRank = cBlock->getRank();
	if(fBlockRank == cBlockRank && fBlockRank == gridRank)
	{
		senderCF = receiverFC = CreateTransmittersHelper::TransmitterPtr( new TbLocalTransmitter< CbVector< LBMReal > >());
		senderFC = receiverCF = CreateTransmittersHelper::TransmitterPtr( new TbLocalTransmitter< CbVector< LBMReal > >());
	}
	else if(cBlockRank == gridRank)
	{
		helper.createTransmitters(cBlock, fBlock, dir, ib, senderCF, receiverCF, comm, CreateTransmittersHelper::MPI);
	}
	else if(fBlockRank == gridRank)
	{
		helper.createTransmitters(fBlock, cBlock, dir, ib, senderFC, receiverFC, comm, CreateTransmittersHelper::MPI);
	}
   UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::createTransmitters(...) - end");
}

