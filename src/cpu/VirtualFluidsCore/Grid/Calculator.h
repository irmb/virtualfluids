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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \file Calculator.h
//! \ingroup Grid
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <PointerDefinitions.h>
#include <vector>

class Grid3D;
class UbScheduler;
class Block3D;
class Block3DConnector;
class CoProcessor;

//! \class Calculator 
//! \brief A base class for main calculation loop  

class Calculator 
{
public:
   Calculator(SPtr<Grid3D> grid, SPtr<UbScheduler> additionalGhostLayerUpdateScheduler, int numberOfTimeSteps);
   virtual ~Calculator();
   //! control of coProcessors
   void addCoProcessor(SPtr<CoProcessor> coProcessor);
   void coProcess(double step);

   virtual void calculate()=0;
protected:
   virtual void initLocalConnectors();
   virtual void initRemoteConnectors();
   void initConnectors(std::vector<SPtr<Block3DConnector> >& connectors);
   void deleteBlocks();
   void deleteConnectors();
   void deleteConnectors(std::vector< std::vector< SPtr<Block3DConnector> > >& conns);

   int minLevel, maxLevel;
   int startTimeStep;
   int numberOfTimeSteps;
   std::vector< std::vector< SPtr<Block3DConnector> > > localConns;
   std::vector< std::vector< SPtr<Block3DConnector> > > remoteConns;

   bool refinement;
   SPtr<Grid3D> grid;
   SPtr<UbScheduler> additionalGhostLayerUpdateScheduler;
   std::vector< std::vector<SPtr<Block3D> > > blocks;

   //localInterConns and remoteInterConns save interpolation connectors 
   //every element save CF connectors for current level and FC connectors for next level
   //e.g. 
   //localInterConns[0] = CF(0), FC(1)
   //localInterConns[1] = CF(1), FC(2)
   //localInterConns[2] 
   std::vector< std::vector< SPtr<Block3DConnector> > > localInterConns;
   std::vector< std::vector< SPtr<Block3DConnector> > > remoteInterConns;

   std::vector< SPtr<CoProcessor> > coProcessors;
};

#endif
