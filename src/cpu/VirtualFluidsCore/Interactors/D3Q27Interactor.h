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
//! \file D3Q27Interactor.h
//! \ingroup Interactor
//! \author Soeren Freudiger
//! \author Sebastian Geller
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef D3Q27INTERACTOR_H
#define D3Q27INTERACTOR_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <PointerDefinitions.h>

#include "UbException.h"
#include "UbTuple.h"
#include "GbPoint3D.h"
#include "Interactor3D.h"
#include "D3Q27System.h"

class BCAdapter;
class Block3D;
class Grid3D;
class GbObject3D;

typedef std::map<SPtr<Block3D>, std::set< std::vector<int> > > BcNodeIndicesMap;
typedef std::map<SPtr<Block3D>, std::set< UbTupleInt3 > > SolidNodeIndicesMap;

//! \brief A specialized class for grid generation.
//! \details Support standard geometric primitives.
class D3Q27Interactor : public Interactor3D 
{
public:
   D3Q27Interactor();
   D3Q27Interactor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type);
   D3Q27Interactor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, SPtr<BCAdapter> bcAdapter,  int type);
   D3Q27Interactor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, SPtr<BCAdapter> bcAdapter,  int type, Interactor3D::Accuracy a);

   virtual ~D3Q27Interactor();

   void setRelevantForForces(const bool& value) {  this->relevantForForces = value; }
   bool isRelevantForForces() { return this->relevantForForces; }

   virtual void addBCAdapter(const SPtr<BCAdapter> bcAdapter) { bcAdapters.push_back(bcAdapter); }
   void deleteBCAdapter() { bcAdapters.clear(); }

 
   virtual void initInteractor(const double& timeStep=0);
   void updateInteractor(const double& timestep=0); 

   void setReinitWithStoredQs(bool reinitWithStoredQsFlag) { this->reinitWithStoredQsFlag = reinitWithStoredQsFlag; }
   
   void removeSolidBlocks() { Interactor3D::removeSolidBlocks(); solidNodeIndicesMap.clear(); }
   void removeBcBlocks() { Interactor3D::removeBcBlocks(); bcNodeIndicesMap.clear(); }

   bool setDifferencesToGbObject3D(const SPtr<Block3D> block);

   ObObject* clone() { throw UbException(UB_EXARGS,"not implemented");	}

   void writeValidationAVSFile(std::string filename);  
   virtual std::vector< std::pair<GbPoint3D,GbPoint3D> >  getQsLineSet();

   void addQsLineSet(std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt2 >& lines);

   const BcNodeIndicesMap& getBcNodeIndicesMap() const { return bcNodeIndicesMap; }

protected:
   bool relevantForForces;
   bool reinitWithStoredQsFlag;

   std::vector<SPtr<BCAdapter> > bcAdapters;

   SolidNodeIndicesMap solidNodeIndicesMap;
   BcNodeIndicesMap bcNodeIndicesMap;
   
   void   initRayVectors();
   double rayX1[D3Q27System::FENDDIR+1];
   double rayX2[D3Q27System::FENDDIR+1];
   double rayX3[D3Q27System::FENDDIR+1];

};


#endif
