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
//! \file NoSlipBCAdapter.cpp
//! \ingroup BoundarConditions
//! \author Sören Freudiger
//=======================================================================================

#ifndef NoSlipBCAdapter_H
#define NoSlipBCAdapter_H

#include "BCAdapter.h"

//! A class provides an interface for no-slip boundary condition in grid generator
class NoSlipBCAdapter : public BCAdapter
{
public:
   NoSlipBCAdapter()
    : BCAdapter()
   {
   }
   NoSlipBCAdapter(const short& secondaryBcOption)
      : BCAdapter(secondaryBcOption)
   {
   }

   void init(const D3Q27Interactor* const& interactor, const double& time=0) {}
   void update(const D3Q27Interactor* const& interactor, const double& time=0) {}

   void adaptBCForDirection( const D3Q27Interactor& interactor, SPtr<BoundaryConditions> bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& q, const int& fdirection, const double& time=0 )
   {
      bc->setNoSlipBoundaryFlag(D3Q27System::INVDIR[fdirection],secondaryBcOption);
      bc->setQ((float)q,fdirection);
   }
   void adaptBC( const D3Q27Interactor& interactor, SPtr<BoundaryConditions> bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& time=0 ) 
   {
      bc->setBcAlgorithmType(algorithmType);
   }

private:

};
#endif //NoSlipBCAdapter_H
