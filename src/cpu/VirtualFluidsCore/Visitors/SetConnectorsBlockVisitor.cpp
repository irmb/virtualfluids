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
#include "Grid3DSystem.h"
#include "Communicator.h"
#include "Grid3D.h"

SetConnectorsBlockVisitor::SetConnectorsBlockVisitor(SPtr<Communicator> comm, bool fullConnector, int dirs, LBMReal nu) : Block3DVisitor(0, Grid3DSystem::MAXLEVEL),	comm(comm), fullConnector(fullConnector), dirs(dirs), nu(nu)
{
}
//////////////////////////////////////////////////////////////////////////
SetConnectorsBlockVisitor::~SetConnectorsBlockVisitor(void)
= default;
//////////////////////////////////////////////////////////////////////////
void SetConnectorsBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
	if(!block) return;

	UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::visit() - start");
   UBLOG(logDEBUG5, block->toString());

	gridRank = comm->getProcessID();
	grid->setRank(gridRank);

	setSameLevelConnectors(grid, block);

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
				//else if(blockRank != neighBlockRank && neighBlock->isActive())
				//{
				//	setRemoteConnectors(block, neighBlock, dir, fullConnector);  

				//	if(dir >=0 && dir<=5)
				//	{
				//		int weight = block->getWeight(neighBlockRank);
				//		weight++;
				//		block->setWeight(neighBlockRank, weight);
				//	}
				//}
			}
		}
      
  //    int weight = block->getNumberOfLocalConnectorsForSurfaces();
		//weight = 6 - weight;
		//block->addWeightForAll(weight);
	}
   UBLOG(logDEBUG5, "D3Q27SetConnectorsBlockVisitor::setSameLevelConnectors() - end");
}

