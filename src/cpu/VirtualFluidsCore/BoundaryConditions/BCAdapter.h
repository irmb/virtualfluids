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
//! \file BCAdapter.h
//! \ingroup BoundarConditions
//! \author Sören Freudiger
//=======================================================================================
#ifndef BCAdapter_H
#define BCAdapter_H

#include <PointerDefinitions.h>

#include "BoundaryConditions.h"
#include "BCAlgorithm.h"

class D3Q27Interactor;

//! \brief Abstract class of baundary conditions adapter
//! \details  BCAdapter supports the definition of boundary conditions in grid generation
class BCAdapter
{
public:
   BCAdapter() = default;

   //! \param secondaryBcOption additional option of boundary conditions
   BCAdapter(const short& secondaryBcOption) 
      :  secondaryBcOption(secondaryBcOption)
   {
   }
   virtual ~BCAdapter() = default;

   //methods
   bool isTimeDependent() { return((this->type & TIMEDEPENDENT) ==  TIMEDEPENDENT); }

   virtual short getSecondaryBcOption() { return this->secondaryBcOption; }
   virtual void  setSecondaryBcOption(const short& val) { this->secondaryBcOption=val; }

   virtual void init(const D3Q27Interactor* const& interactor, const double& time=0) = 0;
   virtual void update(const D3Q27Interactor* const& interactor, const double& time=0) = 0;

   virtual void adaptBC( const D3Q27Interactor& interactor, SPtr<BoundaryConditions> bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& time=0 ) = 0;
   virtual void adaptBCForDirection( const D3Q27Interactor& interactor, SPtr<BoundaryConditions> bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& q, const int& fdirection, const double& time=0 ) = 0;

   void setBcAlgorithm(SPtr<BCAlgorithm> alg) {algorithmType = alg->getType(); algorithm = alg;}
   SPtr<BCAlgorithm> getAlgorithm() {return algorithm;} 
   char getBcAlgorithmType() {return algorithmType;}

protected:
   short secondaryBcOption{0};

   char  type{0};

   SPtr<BCAlgorithm> algorithm;
   char algorithmType{-1};

   static const char   TIMEDEPENDENT = 1<<0;//'1';
   static const char   TIMEPERIODIC  = 1<<1;//'2';
};


#endif //D3Q27BOUNDARYCONDITIONADAPTER_H
