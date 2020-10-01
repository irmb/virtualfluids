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
//! \file Block3DVisitor.h
//! \ingroup Visitors
//! \author Konstantin Kutscher, Soeren Freudiger, Sebastian Geller
//=======================================================================================

#ifndef Block3DVisitor_h
#define Block3DVisitor_h

#include <PointerDefinitions.h>

class Block3D;
class Grid3D;

//! Abstract class provides interface for visitor design pettern
class Block3DVisitor
{
public:
   Block3DVisitor() : startLevel(-1), stopLevel(-1)
   {
   }

   Block3DVisitor(int startLevel, int stopLevel) : startLevel(startLevel), stopLevel(stopLevel)
   {
   }

	virtual ~Block3DVisitor()
   {
   }
	
   virtual void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) = 0;
   
   int  getStartLevel() const; 
   int  getStopLevel() const;
   void setStartLevel(int level);
   void setStopLevel(int level);

private:
   int  startLevel;
   int  stopLevel;
};
//////////////////////////////////////////////////////////////////////////
inline int  Block3DVisitor::getStartLevel() const
{ 
   return this->startLevel;  
}
//////////////////////////////////////////////////////////////////////////
inline int  Block3DVisitor::getStopLevel() const
{ 
   return this->stopLevel;   
}
//////////////////////////////////////////////////////////////////////////
inline void Block3DVisitor::setStartLevel(int level)
{ 
   this->startLevel = level; 
}
//////////////////////////////////////////////////////////////////////////
inline void Block3DVisitor::setStopLevel(int level) 
{ 
   this->stopLevel = level;  
}

#endif 
