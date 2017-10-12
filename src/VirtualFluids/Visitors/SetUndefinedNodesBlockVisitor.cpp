#include "SetUndefinedNodesBlockVisitor.h"
#include "BCArray3D.h"
#include "BoundaryConditions.h"
#include "BCProcessor.h"
#include "Grid3DSystem.h"
#include "D3Q27System.h"

#include <boost/pointer_cast.hpp>

SetUndefinedNodesBlockVisitor::SetUndefinedNodesBlockVisitor() : 
                                    Block3DVisitor(0, Grid3DSystem::MAXLEVEL) 
{

}
//////////////////////////////////////////////////////////////////////////
void SetUndefinedNodesBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   if(!block->hasInterpolationFlag()) return;

   LBMKernelPtr kernel = block->getKernel();

   if(!kernel && (block->getRank() != grid->getRank())) return;

   //width of ghost layer 
   //int gl = kernel->getGhostLayerWidth();
   int gl = 0;
   
   BCArray3D& bcMatrix = kernel->getBCProcessor()->getBCArray();

   int minX1 = gl;
   int minX2 = gl;
   int minX3 = gl;

   int maxX1 = static_cast<int>(bcMatrix.getNX1())-1-gl;
   int maxX2 = static_cast<int>(bcMatrix.getNX2())-1-gl;
   int maxX3 = static_cast<int>(bcMatrix.getNX3())-1-gl;

   //int offset = 2;
   int offset = 3;

   if(block->hasInterpolationFlag(D3Q27System::E))
   {
      int startix1 = maxX1;
      int endix1   = maxX1;
      if(block->hasInterpolationFlagCF()) startix1 = startix1-offset;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::W))
   {
      int startix1 = minX1;
      int endix1   = minX1;
      if(block->hasInterpolationFlagCF()) endix1 = endix1+offset;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::N))
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = maxX2;
      int endix2   = maxX2;
      if(block->hasInterpolationFlagCF()) startix2 = startix2-offset;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::S))
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = minX2;
      if(block->hasInterpolationFlagCF()) endix2 = endix2+offset;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::T))
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = maxX3;
      int endix3   = maxX3;
      if(block->hasInterpolationFlagCF()) startix3 = startix3-offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::B))
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = minX3;
      int endix3   = minX3;
      if(block->hasInterpolationFlagCF()) endix3 = endix3+offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::NE))
   {
      int startix1 = maxX1;
      int endix1   = maxX1;
      if(block->hasInterpolationFlagCF()) startix1 = startix1-offset;
      int startix2 = maxX2;
      int endix2   = maxX2;
      if(block->hasInterpolationFlagCF()) startix2 = startix2-offset;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::SW))
   {
      int startix1 = minX1;
      int endix1   = minX1;
      if(block->hasInterpolationFlagCF()) endix1 = endix1+offset;
      int startix2 = minX2;
      int endix2   = minX2;
      if(block->hasInterpolationFlagCF()) endix2 = endix2+offset;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::SE))
   {
      int startix1 = maxX1;
      int endix1   = maxX1;
      if(block->hasInterpolationFlagCF()) startix1 = startix1-offset;
      int startix2 = minX2;
      int endix2   = minX2;
      if(block->hasInterpolationFlagCF()) endix2 = endix2+offset;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::NW))  
   {
      int startix1 = minX1;
      int endix1   = minX1;
      if(block->hasInterpolationFlagCF()) endix1 = endix1+offset;
      int startix2 = maxX2;
      int endix2   = maxX2;
      if(block->hasInterpolationFlagCF()) startix2 = startix2-offset;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::TE))  
   {
      int startix1 = maxX1;
      int endix1   = maxX1;
      if(block->hasInterpolationFlagCF()) startix1 = startix1-offset;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = maxX3;
      int endix3   = maxX3;
      if(block->hasInterpolationFlagCF()) startix3 = startix3-offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::BW))  
   {
      int startix1 = minX1;
      int endix1   = minX1;
      if(block->hasInterpolationFlagCF()) endix1 = endix1+offset;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = minX3;
      int endix3   = minX3;
      if(block->hasInterpolationFlagCF()) endix3 = endix3+offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::BE)) 
   {
      int startix1 = maxX1;
      int endix1   = maxX1;
      if(block->hasInterpolationFlagCF()) startix1 = startix1-offset;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = minX3;
      int endix3   = minX3;
      if(block->hasInterpolationFlagCF()) endix3 = endix3+offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::TW))  
   {
      int startix1 = minX1;
      int endix1   = minX1;
      if(block->hasInterpolationFlagCF()) endix1 = endix1+offset;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = maxX3;
      int endix3   = maxX3;
      if(block->hasInterpolationFlagCF()) startix3 = startix3-offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::TN)) 
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = maxX2;
      int endix2   = maxX2;
      if(block->hasInterpolationFlagCF()) startix2 = startix2-offset;
      int startix3 = maxX3;
      int endix3   = maxX3;
      if(block->hasInterpolationFlagCF()) startix3 = startix3-offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::BS))  
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = minX2;
      if(block->hasInterpolationFlagCF()) endix2 = endix2+offset;
      int startix3 = minX3;
      int endix3   = minX3;
      if(block->hasInterpolationFlagCF()) endix3 = endix3+offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::BN))  
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = maxX2;
      int endix2   = maxX2;
      if(block->hasInterpolationFlagCF()) startix2 = startix2-offset;
      int startix3 = minX3;
      int endix3   = minX3;
      if(block->hasInterpolationFlagCF()) endix3 = endix3+offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::TS))  
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = minX2;
      if(block->hasInterpolationFlagCF()) endix2 = endix2+offset;
      int startix3 = maxX3;
      int endix3   = maxX3;
      if(block->hasInterpolationFlagCF()) startix3 = startix3-offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::TNE)) 
   {
      int startix1 = maxX1;
      int endix1   = maxX1;
      if(block->hasInterpolationFlagCF()) startix1 = startix1-offset;
      int startix2 = maxX2;
      int endix2   = maxX2;
      if(block->hasInterpolationFlagCF()) startix2 = startix2-offset;
      int startix3 = maxX3;
      int endix3   = maxX3;
      if(block->hasInterpolationFlagCF()) startix3 = startix3-offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::TNW)) 
   {
      int startix1 = minX1;
      int endix1   = minX1;
      if(block->hasInterpolationFlagCF()) endix1 = endix1+offset;
      int startix2 = maxX2;
      int endix2   = maxX2;
      if(block->hasInterpolationFlagCF()) startix2 = startix2-offset;
      int startix3 = maxX3;
      int endix3   = maxX3;
      if(block->hasInterpolationFlagCF()) startix3 = startix3-offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::TSE)) 
   {
      int startix1 = maxX1;
      int endix1   = maxX1;
      if(block->hasInterpolationFlagCF()) startix1 = startix1-offset;
      int startix2 = minX2;
      int endix2   = minX2;
      if(block->hasInterpolationFlagCF()) endix2 = endix2+offset;
      int startix3 = maxX3;
      int endix3   = maxX3;
      if(block->hasInterpolationFlagCF()) startix3 = startix3-offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::TSW)) 
   {
      int startix1 = minX1;
      int endix1   = minX1;
      if(block->hasInterpolationFlagCF()) endix1 = endix1+offset;
      int startix2 = minX2;
      int endix2   = minX2;
      if(block->hasInterpolationFlagCF()) endix2 = endix2+offset;
      int startix3 = maxX3;
      int endix3   = maxX3;
      if(block->hasInterpolationFlagCF()) startix3 = startix3-offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::BNE)) 
   {
      int startix1 = maxX1;
      int endix1   = maxX1;
      if(block->hasInterpolationFlagCF()) startix1 = startix1-offset;
      int startix2 = maxX2;
      int endix2   = maxX2;
      if(block->hasInterpolationFlagCF()) startix2 = startix2-offset;
      int startix3 = minX3;
      int endix3   = minX3;
      if(block->hasInterpolationFlagCF()) endix3 = endix3+offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::BNW))
   {
      int startix1 = minX1;
      int endix1   = minX1;
      if(block->hasInterpolationFlagCF()) endix1 = endix1+offset;
      int startix2 = maxX2;
      int endix2   = maxX2;
      if(block->hasInterpolationFlagCF()) startix2 = startix2-offset;
      int startix3 = minX3;
      int endix3   = minX3;
      if(block->hasInterpolationFlagCF()) endix3 = endix3+offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::BSE)) 
   {
      int startix1 = maxX1;
      int endix1   = maxX1;
      if(block->hasInterpolationFlagCF()) startix1 = startix1-offset;
      int startix2 = minX2;
      int endix2   = minX2;
      if(block->hasInterpolationFlagCF()) endix2 = endix2+offset;
      int startix3 = minX3;
      int endix3   = minX3;
      if(block->hasInterpolationFlagCF()) endix3 = endix3+offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlag(D3Q27System::BSW)) 
   {
      int startix1 = minX1;
      int endix1   = minX1;
      if(block->hasInterpolationFlagCF()) endix1 = endix1+offset;
      int startix2 = minX2;
      int endix2   = minX2;
      if(block->hasInterpolationFlagCF()) endix2 = endix2+offset;
      int startix3 = minX3;
      int endix3   = minX3;
      if(block->hasInterpolationFlagCF()) endix3 = endix3+offset;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }



	//////////////////////////////////////////////////////////////////////////
   int offset2 = 1;
   int ll = 0;

   minX1 = ll;
   minX2 = ll;
   minX3 = ll;

   maxX1 = static_cast<int>(bcMatrix.getNX1())-1-ll;
   maxX2 = static_cast<int>(bcMatrix.getNX2())-1-ll;
   maxX3 = static_cast<int>(bcMatrix.getNX3())-1-ll;

   if(block->hasInterpolationFlagFC(D3Q27System::E))
   {
      int startix1 = maxX1-offset2;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::W))
   {
      int startix1 = minX1;
      int endix1   = minX1+offset2;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::N))
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = maxX2-offset2;
      int endix2   = maxX2;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::S))
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = minX2+offset2;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::T))
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = maxX3-offset2;
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::B))
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = minX3;
      int endix3   = minX3+offset2;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::NE))
   {
      int startix1 = maxX1-offset2;
      int endix1   = maxX1;
      int startix2 = maxX2-offset2;
      int endix2   = maxX2;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::SW))
   {
      int startix1 = minX1;
      int endix1   = minX1+offset2;
      int startix2 = minX2;
      int endix2   = minX2+offset2;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::SE))
   {
      int startix1 = maxX1-offset2;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = minX2+offset2;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::NW))  
   {
      int startix1 = minX1;
      int endix1   = minX1+offset2;
      int startix2 = maxX2-offset2;
      int endix2   = maxX2;
      int startix3 = minX3; 
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::TE))  
   {
      int startix1 = maxX1-offset2;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = maxX3-offset2;
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::BW))  
   {
      int startix1 = minX1;
      int endix1   = minX1+offset2;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = minX3;
      int endix3   = minX3+offset2;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::BE)) 
   {
      int startix1 = maxX1-offset2;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = minX3;
      int endix3   = minX3+offset2;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::TW))  
   {
      int startix1 = minX1;
      int endix1   = minX1+offset2;
      int startix2 = minX2;
      int endix2   = maxX2;
      int startix3 = maxX3-offset2;
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::TN)) 
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = maxX2-offset2;
      int endix2   = maxX2;
      int startix3 = maxX3-offset2;
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::BS))  
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = minX2+offset2;
      int startix3 = minX3;
      int endix3   = minX3+offset2;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::BN))  
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = maxX2-offset2;
      int endix2   = maxX2;
      int startix3 = minX3;
      int endix3   = minX3+offset2;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::TS))  
   {
      int startix1 = minX1;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = minX2+offset2;
      int startix3 = maxX3-offset2;
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::TNE)) 
   {
      int startix1 = maxX1-offset2;
      int endix1   = maxX1;
      int startix2 = maxX2-offset2;
      int endix2   = maxX2;
      int startix3 = maxX3-offset2;
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::TNW)) 
   {
      int startix1 = minX1;
      int endix1   = minX1+offset2;
      int startix2 = maxX2-offset2;
      int endix2   = maxX2;
      int startix3 = maxX3-offset2;
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::TSE)) 
   {
      int startix1 = maxX1-offset2;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = minX2+offset2;
      int startix3 = maxX3-offset2;
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::TSW)) 
   {
      int startix1 = minX1;
      int endix1   = minX1+offset2;
      int startix2 = minX2;
      int endix2   = minX2+offset2;
      int startix3 = maxX3-offset2;
      int endix3   = maxX3;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::BNE)) 
   {
      int startix1 = maxX1-offset2;
      int endix1   = maxX1;
      int startix2 = maxX2-offset2;
      int endix2   = maxX2;
      int startix3 = minX3;
      int endix3   = minX3+offset2;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::BNW))
   {
      int startix1 = minX1;
      int endix1   = minX1+offset2;
      int startix2 = maxX2-offset2;
      int endix2   = maxX2;
      int startix3 = minX3;
      int endix3   = minX3+offset2;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::BSE)) 
   {
      int startix1 = maxX1-offset2;
      int endix1   = maxX1;
      int startix2 = minX2;
      int endix2   = minX2+offset2;
      int startix3 = minX3;
      int endix3   = minX3+offset2;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   if(block->hasInterpolationFlagFC(D3Q27System::BSW)) 
   {
      int startix1 = minX1;
      int endix1   = minX1+offset2;
      int startix2 = minX2;
      int endix2   = minX2+offset2;
      int startix3 = minX3;
      int endix3   = minX3+offset2;
      this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
   }
   
   //invert scaleCF blocks
   if(block->hasInterpolationFlagCF())
   {
      if(block->hasInterpolationFlagFC()) 
      {
         for (int i = D3Q27System::E; i <= D3Q27System::BSW; i++)
         {
             UBLOG(logINFO, "FC in dir="<<i<<" "<<block->hasInterpolationFlagFC(i));
         }
         for (int i = D3Q27System::E; i<=D3Q27System::BSW; i++)
         {
            UBLOG(logINFO, "CF in dir="<<i<<" "<<block->hasInterpolationFlagCF(i));
         }
         throw UbException(UB_EXARGS, "block "+block->toString()+" has CF and FC");
      }

      minX1 = gl;
      minX2 = gl;
      minX3 = gl;

      maxX1 = static_cast<int>(bcMatrix.getNX1())-1-gl;
      maxX2 = static_cast<int>(bcMatrix.getNX2())-1-gl;
      maxX3 = static_cast<int>(bcMatrix.getNX3())-1-gl;

      for (int ix3=minX3; ix3<=maxX3; ix3++)
         for (int ix2=minX2; ix2<=maxX2; ix2++)
            for (int ix1=minX1; ix1<=maxX1; ix1++)
            {
               if(bcMatrix.isUndefined(ix1, ix2, ix3)) bcMatrix.setFluid(ix1, ix2, ix3);
               else                                    bcMatrix.setUndefined(ix1, ix2, ix3);
            }
            return;
   }

}
//////////////////////////////////////////////////////////////////////////
void SetUndefinedNodesBlockVisitor::setNodesUndefined( int startix1, int endix1, int startix2, int endix2, int startix3, int endix3, BCArray3D& bcMatrix )
{
   for (int ix3=startix3; ix3<=endix3; ix3++)
      for (int ix2=startix2; ix2<=endix2; ix2++)
         for (int ix1=startix1; ix1<=endix1; ix1++)
         {
            bcMatrix.setUndefined(ix1, ix2, ix3);
         }
}
