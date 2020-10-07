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
//! \file BCProcessor.h
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================

#include "BCProcessor.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "DataSet3D.h"
#include "ILBMKernel.h"
#include "BCArray3D.h"
#include "BCAlgorithm.h"

BCProcessor::BCProcessor()
= default;
//////////////////////////////////////////////////////////////////////////
BCProcessor::BCProcessor(SPtr<ILBMKernel> kernel)
{
   SPtr<DistributionArray3D> distributions = std::dynamic_pointer_cast<EsoTwist3D>(kernel->getDataSet()->getFdistributions());
   bcArray = std::make_shared<BCArray3D>(distributions->getNX1(), distributions->getNX2(), distributions->getNX3(), BCArray3D::FLUID);
}
//////////////////////////////////////////////////////////////////////////
BCProcessor::~BCProcessor()
= default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCProcessor> BCProcessor::clone(SPtr<ILBMKernel> kernel)
{
   SPtr<BCProcessor> bcProcessor(new BCProcessor(kernel));
   return bcProcessor;
}
//////////////////////////////////////////////////////////////////////////
SPtr<BCArray3D> BCProcessor::getBCArray()
{ 
   return bcArray; 
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::setBCArray(SPtr<BCArray3D> bcarray)
{
   bcArray = bcarray;
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::addBC(SPtr<BCAlgorithm> bc)
{
   if (bc->isPreCollision())
   {
      preBC.push_back(bc);
   }
   else
   {
      postBC.push_back(bc);
   }
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::applyPreCollisionBC()
{
   for(SPtr<BCAlgorithm> bc : preBC)
      bc->applyBC();
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::applyPostCollisionBC()
{
    for (SPtr<BCAlgorithm> bc : postBC)
        bc->applyBC();
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::clearBC()
{
   preBC.clear();
   postBC.clear();
}

