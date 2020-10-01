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
//! \file InteractorsHelper.h
//! \ingroup Interactor
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef InteractorsHelper_h 
#define InteractorsHelper_h

#include <vector>
#include <PointerDefinitions.h>


class Interactor3D;
class Block3D;
class Grid3D;
class Grid3DVisitor;

//! A helper class for grid generation.
class InteractorsHelper
{
public:
   InteractorsHelper(SPtr<Grid3D> grid, SPtr<Grid3DVisitor> visitor, bool deleteBlocks=true);
   ~InteractorsHelper();

   void addInteractor(SPtr<Interactor3D> interactor);
   void selectBlocks();
   void setBC();
   void sendDomainDecompositionVisitor() const;

protected:
   void deleteSolidBlocks();
   void setBcBlocks();

private:
   void updateGrid();

   std::vector<SPtr<Interactor3D> > interactors;
   SPtr<Grid3D> grid;
   std::vector<SPtr<Block3D> > solidBlocks;
   SPtr<Grid3DVisitor> visitor;
   bool deleteBlocks;
};

#endif
