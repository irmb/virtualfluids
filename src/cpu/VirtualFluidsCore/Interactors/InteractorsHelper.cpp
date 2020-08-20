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
#include "SetSolidBlocksBlockVisitor.h"
#include "SetBcBlocksBlockVisitor.h"


InteractorsHelper::InteractorsHelper(SPtr<Grid3D> grid) :grid(grid)
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
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::selectBlocks()
{
   setBcBlocks();
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::setBcBlocks()
{
    for(const SPtr<Interactor3D> interactor : interactors)
    {
       SetBcBlocksBlockVisitor v(interactor);
       grid->accept(v);
    }
}

