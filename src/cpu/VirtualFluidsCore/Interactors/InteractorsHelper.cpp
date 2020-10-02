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
//! \file InteractorsHelper.cpp
//! \ingroup Interactor
//! \author Konstantin Kutscher
//=======================================================================================

#include "InteractorsHelper.h"

#include <Grid3DVisitor.h>
#include <Grid3D.h>
#include <Interactor3D.h>
#include "Block3D.h"
#include "Communicator.h"
#include "SetSolidBlocksBlockVisitor.h"
#include "SetBcBlocksBlockVisitor.h"


InteractorsHelper::InteractorsHelper(SPtr<Grid3D> grid, SPtr<Grid3DVisitor> visitor, bool deleteBlocks) :
                                     grid(grid), visitor(visitor), deleteBlocks(deleteBlocks)
{

}
//////////////////////////////////////////////////////////////////////////
InteractorsHelper::~InteractorsHelper()
{

}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::addInteractor( SPtr<Interactor3D> interactor )
{
   interactors.push_back(interactor);
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::setBC()
{
    for(SPtr<Interactor3D> i : interactors)
        i->initInteractor();
}

void InteractorsHelper::sendDomainDecompositionVisitor() const
{
    grid->accept( visitor );
}

//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::selectBlocks()
{
   sendDomainDecompositionVisitor();
   deleteSolidBlocks();

   sendDomainDecompositionVisitor();
   setBcBlocks();
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::deleteSolidBlocks()
{
    for(SPtr<Interactor3D> interactor : interactors)
    {
        SetSolidBlocksBlockVisitor v(interactor);
        grid->accept(v);
        if (deleteBlocks)
        {
           std::vector<SPtr<Block3D>>& sb = interactor->getSolidBlockSet();
           solidBlocks.insert(solidBlocks.end(), sb.begin(), sb.end());
           interactor->removeSolidBlocks();
        }
    }

   if (deleteBlocks) updateGrid();
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::setBcBlocks()
{
    for(const auto& interactor : interactors)
    {
       SetBcBlocksBlockVisitor v(interactor);
       grid->accept(v);
    }
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::updateGrid()
{
    std::vector<int> ids;

    for(const SPtr<Block3D> block : solidBlocks)
        ids.push_back(block->getGlobalID());

    std::vector<int> rids;
    Communicator::getInstance()->allGather(ids, rids);
    grid->deleteBlocks(rids);
}

