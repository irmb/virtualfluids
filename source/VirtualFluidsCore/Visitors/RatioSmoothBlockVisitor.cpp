#include "RatioSmoothBlockVisitor.h"
#include "Grid3DSystem.h"
#include "Block3D.h"
#include "Grid3D.h"

RatioSmoothBlockVisitor::RatioSmoothBlockVisitor(int levelDepth, bool includeNotActiveBlocks)
   :   Block3DVisitor(Grid3DSystem::MAXLEVEL, 0)
   , maxLevelRatio(1)
   , expandBlocks(true)
   , levelDepth(levelDepth)
   , includeNotActiveBlocks(includeNotActiveBlocks)
{
}
//////////////////////////////////////////////////////////////////////////
void RatioSmoothBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   int ix1, ix2, ix3, level;
   ix1 = block->getX1();
   ix2 = block->getX2();
   ix3 = block->getX3();
   level = block->getLevel();

   if( block->isActive()  || includeNotActiveBlocks )
   {
      if(this->expandBlocks)
      {
         if(this->lookForExpand(grid,ix1, ix2, ix3, level))
         {
            grid->expandBlock(ix1, ix2, ix3, level);
         }
      }
      else
      {
         if(this->lookForCollapse(grid,ix1, ix2, ix3, level))
         {
            grid->collapseBlock(ix1, ix2, ix3, level, levelDepth);
         }
      }
   }
}

//////////////////////////////////////////////////////////////////////////
bool RatioSmoothBlockVisitor::lookForExpand(Grid3DPtr grid, const int& ix1, const int& ix2, const int& ix3, const int& level)
{
   std::vector<Block3DPtr> neighbors;
   grid->getAllNeighbors(ix1, ix2, ix3, level, this->levelDepth, neighbors);
   int nix1, nix2,nix3, nlev;
   for(size_t i=0; i<neighbors.size(); i++)
   {
      if(   ( neighbors[i]->isActive() || includeNotActiveBlocks )
         && neighbors[i]->getLevel() > level) 
      {
         nix1 = (neighbors)[i]->getX1();
         nix2 = (neighbors)[i]->getX2();
         nix3 = (neighbors)[i]->getX3();
         nlev = (neighbors)[i]->getLevel();

         std::vector<Block3DPtr> neighbors1;
         grid->getAllNeighbors(nix1, nix2, nix3, nlev, nlev+1, neighbors1);
         for(size_t j=0; j<neighbors1.size(); j++)
         {
            if(   ( neighbors1[j]->isActive() || includeNotActiveBlocks )
               && neighbors1[j]->getLevel() > level+this->maxLevelRatio) 
            {
               return true;
            }
         }
      }
   }
   return false;
}
//////////////////////////////////////////////////////////////////////////
bool RatioSmoothBlockVisitor::lookForCollapse(Grid3DPtr grid, const int& ix1, const int& ix2, const int& ix3, const int& level)
{
   std::vector<Block3DPtr> neighbors;
   grid->getAllNeighbors(ix1, ix2,ix3, level, this->levelDepth, neighbors);
   for(size_t i=0; i<neighbors.size(); i++)
   {     
      if(    ( neighbors[i]->isActive() || includeNotActiveBlocks )
         &&  neighbors[i]->getLevel() < level-this->maxLevelRatio) 
      {
         throw UbException(UB_EXARGS," not implemented till now");
         return true;
      }
   }

   return false;
}
//////////////////////////////////////////////////////////////////////////
void RatioSmoothBlockVisitor::setExpandByAdaptation(bool expandBlocks)
{
   if(this->expandBlocks != expandBlocks)
   {
      this->expandBlocks = expandBlocks;

      int l1 = Block3DVisitor::getStartLevel();
      int l2 = Block3DVisitor::getStopLevel();

      if(expandBlocks) { if(l1 < l2) { Block3DVisitor::setStartLevel(l2); Block3DVisitor::setStopLevel(l1); } }
      else             { if(l2 < l1) { Block3DVisitor::setStartLevel(l2); Block3DVisitor::setStopLevel(l1); } }
   }
}
//////////////////////////////////////////////////////////////////////////
void RatioSmoothBlockVisitor::setLevelRatio(int ratio)
{
   if(ratio < 1) throw UbException(UB_EXARGS,"illegal ratio specified");
   this->maxLevelRatio = ratio;
}
//////////////////////////////////////////////////////////////////////////
std::string RatioSmoothBlockVisitor::getSpecificDescription()
{
   std::string str("Ratio:");
   return str;
}
//////////////////////////////////////////////////////////////////////////
int RatioSmoothBlockVisitor::getStartLevel()
{
   int l1 = Block3DVisitor::getStartLevel();
   int l2 = Block3DVisitor::getStopLevel();

   if(this->expandBlocks) { if(l2 < l1) return(l1); else return(l2); }
   else                   { if(l2 < l1) return(l2); else return(l1); }
}
//////////////////////////////////////////////////////////////////////////
int RatioSmoothBlockVisitor::getStopLevel()
{
   int l1 = Block3DVisitor::getStartLevel();
   int l2 = Block3DVisitor::getStopLevel();

   if(this->expandBlocks) { if(l2 < l1) return(l2); else return(l1); }
   else                   { if(l2 < l1) return(l1); else return(l2); }
}
//////////////////////////////////////////////////////////////////////////
void RatioSmoothBlockVisitor::setStartLevel(int level)
{
   if(this->expandBlocks) { if(level >= Block3DVisitor::getStopLevel()) Block3DVisitor::setStartLevel(level); }
   else                   { if(level <= Block3DVisitor::getStopLevel()) Block3DVisitor::setStartLevel(level); }
}
//////////////////////////////////////////////////////////////////////////
void RatioSmoothBlockVisitor::setStopLevel(int level)
{
   if(this->expandBlocks) { if(level <= Block3DVisitor::getStartLevel()) Block3DVisitor::setStopLevel(level); }
   else                   { if(level >= Block3DVisitor::getStartLevel()) Block3DVisitor::setStopLevel(level); }
}

